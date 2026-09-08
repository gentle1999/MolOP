"""
Author: TMJ
Date: 2024-06-17 20:42:47
LastEditors: TMJ
LastEditTime: 2026-05-05 19:33:30
Description: 请填写简介
"""

import json
import sys
import threading
from collections.abc import Iterable, Mapping, Sequence
from contextlib import ExitStack
from dataclasses import asdict
from itertools import islice
from typing import Any, ClassVar, Literal, Optional, cast

import numpy as np
from molgr import ReconstructionBatchRequest, ReconstructionBatchResult, iter_xyz_to_rdmol_batch
from molgr.config import CONFIG as MOLGR_CONFIG
from molgr.interface import xyz_to_rdmol
from openbabel import pybel
from pint._typing import UnitLike
from pydantic import Field, PrivateAttr, computed_field, model_validator
from rdkit import Chem
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from molop.config import molopconfig, moloplogger
from molop.io.base_models._format_transform import (
    FrameFormatTransformMixin,  # pyright: ignore[reportAttributeAccessIssue]
)
from molop.io.base_models.summary import SummaryDict, summary_column
from molop.structure.editing import replace_substituent
from molop.structure.FormatConverter import rdmol_to_omol
from molop.structure.GeometryTransformation import get_geometry_info, standard_orient
from molop.structure.substitution_models import StereoPolicy
from molop.structure.topology import (
    build_mol_from_atoms_and_bonds,
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
    get_total_charge,
    get_total_multiplicity,
    reset_atom_index,
)
from molop.structure.utils import canonical_smiles
from molop.unit import atom_ureg
from molop.utils.progressbar import (
    NativeReconstructionConcurrencyError,
    is_loky_worker,
    native_reconstruction_guard,
)
from molop.utils.types import PintArrayNx3, RdMol

from .Bases import BaseDataClassWithUnit
from .DataClasses import InternalCoords
from .source import canonical_json_sha256


pt = Chem.GetPeriodicTable()
_ReconstructionOptions = tuple[
    Literal["cpp", "python"],
    Literal["raise", "return_suspicious"],
    bool,
    bool,
]


def _reconstruction_diagnostics_dict(status: object | None) -> dict[str, Any] | None:
    """Return MolGR diagnostics as JSON-compatible model evidence."""

    if status is None:
        return None
    as_dict = getattr(status, "as_dict", None)
    if callable(as_dict):
        serialized = as_dict()
        if isinstance(serialized, Mapping):
            return {str(key): value for key, value in serialized.items()}
    return {"message": str(status)}


class _PickleSafeRLock:
    """A per-molecule lock that is recreated when a model crosses loky."""

    __slots__ = ("_lock",)

    def __init__(self) -> None:
        self._lock = threading.RLock()

    def __enter__(self) -> "_PickleSafeRLock":
        self._lock.acquire()
        return self

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        self._lock.release()

    def __getstate__(self) -> dict[str, object]:
        return {}

    def __setstate__(self, _state: dict[str, object]) -> None:
        self._lock = threading.RLock()


EXCLUDE_FIELDS_IF_NO_BOND = {
    "smiles",
    "canonical_smiles",
    "bonds",
    "formal_charges",
    "formal_num_radicals",
}


class Molecule(FrameFormatTransformMixin, BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {"coords": atom_ureg.angstrom}
    atoms: list[int] = Field(default_factory=list, description="atom numbers", title="Atom numbers")
    coords: PintArrayNx3 = Field(
        default=np.zeros((0, 3)) * atom_ureg.angstrom,
        description="Atom coordinates, unit is `angstrom`",
        title="Atom coordinates",
    )
    charge: int = Field(default=0, description="Molecule total charge")
    multiplicity: int = Field(default=1, description="Molecule total multiplicity")

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def smiles(self) -> str:
        return self.to_SMILES()

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def canonical_smiles(self) -> str:
        return self.to_canonical_SMILES()

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def topology_v3000_molblock(self) -> str | None:
        """Return a map-free V3000 graph carrier for a trusted topology."""

        if (
            self.rdmol is None
            or self.source_to_topology_atom_permutation is None
            or self.topology_reconstruction_status in {"failed", "suspicious_fallback"}
        ):
            return None
        topology = self.rdmol_no_conformer
        if topology.HasProp("_Name"):
            topology.ClearProp("_Name")
        for atom in topology.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToMolBlock(
            topology,
            includeStereo=True,
            confId=-1,
            kekulize=True,
            forceV3000=True,
        )

    bonds: list[tuple[int, int, int, int]] = Field(
        default_factory=list,
        description="Bond information, each bond is represented by a tuple of three integers, "
        "where the first two integers represent the indices of the two atoms involved in the bond, "
        "the third integer represents the bond order (follow the RDKit convention)"
        "the fourth integer represents the stereo configuration (follow the RDKit convention)",
    )
    formal_charges: list[int] = Field(
        default_factory=list,
        description="Formal charges of each atom",
    )
    formal_num_radicals: list[int] = Field(
        default_factory=list,
        description="Number of radical electrons of each atom",
    )
    source_to_topology_atom_permutation: list[int] | None = Field(
        default=None,
        description=(
            "Mapping from each source atom index to its topology atom index; identity is "
            "recorded only after source-order verification"
        ),
        exclude_if=lambda value: value is None,
    )
    topology_reconstruction_backend: Literal["cpp", "python"] | None = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    topology_reconstruction_failure_policy: Literal["raise", "return_suspicious"] | None = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    topology_make_dative_bonds: bool | None = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    topology_make_stereochemistry: bool | None = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    topology_reconstruction_diagnostics: dict[str, Any] | None = Field(
        default=None,
        description="Structured MolGR diagnostics for a failed reconstruction",
        exclude_if=lambda value: value is None,
    )
    topology_reconstruction_config_sha256: str | None = Field(
        default=None,
        pattern=r"^[0-9a-f]{64}$",
        exclude_if=lambda value: value is None,
    )
    topology_reconstruction_status: (
        Literal["provided", "succeeded", "failed", "suspicious_fallback"] | None
    ) = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    _rdmol: RdMol | None = PrivateAttr(default=None)
    _smiles_cache: str | None = PrivateAttr(default=None)
    _canonical_smiles_cache: str | None = PrivateAttr(default=None)
    _topology_lock: _PickleSafeRLock = PrivateAttr(default_factory=_PickleSafeRLock)

    @model_validator(mode="after")
    def validate_source_to_topology_atom_permutation(self) -> "Molecule":
        permutation = self.source_to_topology_atom_permutation
        if permutation is not None and sorted(permutation) != list(range(len(self.atoms))):
            raise ValueError(
                "source_to_topology_atom_permutation must be a permutation of atom indices"
            )
        return self

    @property
    def atom_symbols(self) -> list[str]:
        """
        Get the atom symbols.

        Returns:
            List[str]: A list of atom symbols.
        """
        return [Chem.Atom(atom).GetSymbol() for atom in self.atoms]

    @property
    def total_electrons(self) -> int:
        """
        Get the total electrons.

        Returns:
            int: The total electrons.
        """
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self.atoms) - self.charge

    @property
    def elements(self) -> list[str]:
        """
        Get the elements set.

        Returns:
            List[str]: A list of element symbols dropped duplicates.
        """
        return list(set(self.atom_symbols))

    @property
    def formula(self) -> str:
        """
        Get the formula.

        Returns:
            str: The formula.
        """
        if self.rdmol:
            return CalcMolFormula(self.rdmol)
        return "".join(
            [
                f"{pt.GetElementSymbol(i)}{self.atoms.count(i)}"
                for i in range(1, 119)
                if self.atoms.count(i) != 0
            ]
        )

    def _resolve_reconstruction_options(
        self,
        *,
        backend: Literal["cpp", "python"] | None = None,
        reconstruction_failure_policy: Literal["raise", "return_suspicious"] | None = None,
        make_dative_bonds: bool | None = None,
        make_stereochemistry: bool | None = None,
    ) -> _ReconstructionOptions:
        resolved_backend = (
            backend
            or self.topology_reconstruction_backend
            or molopconfig.graph_reconstruction_backend
        )
        resolved_failure_policy = (
            reconstruction_failure_policy
            or self.topology_reconstruction_failure_policy
            or molopconfig.reconstruction_failure_policy
        )
        resolved_dative_bonds = (
            make_dative_bonds if make_dative_bonds is not None else self.topology_make_dative_bonds
        )
        if resolved_dative_bonds is None:
            resolved_dative_bonds = molopconfig.make_dative_bonds
        resolved_stereochemistry = (
            make_stereochemistry
            if make_stereochemistry is not None
            else self.topology_make_stereochemistry
        )
        if resolved_stereochemistry is None:
            resolved_stereochemistry = molopconfig.make_stereochemistry
        return (
            resolved_backend,
            resolved_failure_policy,
            resolved_dative_bonds,
            resolved_stereochemistry,
        )

    def _record_topology_reconstruction_provenance(
        self,
        options: _ReconstructionOptions | None = None,
    ) -> None:
        backend, failure_policy, make_dative_bonds, make_stereochemistry = (
            options if options is not None else self._resolve_reconstruction_options()
        )
        molopconfig.apply_molgr_reconstruction_policy(failure_policy)
        config = {
            "backend": backend,
            "reconstruction_failure_policy": failure_policy,
            "make_dative_bonds": make_dative_bonds,
            "make_stereochemistry": make_stereochemistry,
            "molgr": asdict(MOLGR_CONFIG),
        }
        self.topology_reconstruction_backend = backend
        self.topology_reconstruction_failure_policy = failure_policy
        self.topology_make_dative_bonds = make_dative_bonds
        self.topology_make_stereochemistry = make_stereochemistry
        self.topology_reconstruction_config_sha256 = canonical_json_sha256(config)

    @staticmethod
    def _rdmol_reconstruction_diagnostics(rdmol: RdMol) -> dict[str, Any] | None:
        if not rdmol.HasProp("_MolGRReconstructionDiagnostics"):
            return None
        raw_diagnostics = rdmol.GetProp("_MolGRReconstructionDiagnostics")
        try:
            diagnostics = json.loads(raw_diagnostics)
        except (TypeError, ValueError):
            diagnostics = {"message": raw_diagnostics}
        return diagnostics if isinstance(diagnostics, dict) else None

    def _apply_reconstruction_result(
        self,
        reconstructed: RdMol | None,
        *,
        status: object | None = None,
        error: BaseException | None = None,
    ) -> None:
        if status is not None:
            self.topology_reconstruction_diagnostics = _reconstruction_diagnostics_dict(status)
        if reconstructed is None:
            self.topology_reconstruction_status = "failed"
            if error is not None:
                error_diagnostics = getattr(error, "diagnostics", None)
                if error_diagnostics is not None:
                    self.topology_reconstruction_diagnostics = _reconstruction_diagnostics_dict(
                        error_diagnostics
                    )
                moloplogger.error(f"{error}")
            elif status is not None:
                moloplogger.error(getattr(status, "message", str(status)))
            self.__get_topology()
            return

        self._rdmol = reconstructed
        diagnostics = self._rdmol_reconstruction_diagnostics(reconstructed)
        if diagnostics is not None:
            self.topology_reconstruction_diagnostics = diagnostics
        self.topology_reconstruction_status = (
            "suspicious_fallback"
            if reconstructed.HasProp("_MolGRReconstructionStatus")
            and reconstructed.GetProp("_MolGRReconstructionStatus") == "suspicious_fallback"
            else "succeeded"
        )
        self.__get_topology()

    def _materialize_topology_from_fields(
        self,
        status: Literal["provided", "succeeded"],
    ) -> RdMol:
        if not self.formal_charges or not self.formal_num_radicals:
            raise ValueError("If bonds given, formal charges and spins must be provided.")
        rdmol = build_mol_from_atoms_and_bonds(
            self.atoms,
            self.bonds,
            self.formal_charges,
            self.formal_num_radicals,
            coords=self.coords.m,
        )
        if rdmol is None:
            raise ValueError("Building the provided topology returned no molecule.")
        rdmol = cast(RdMol, rdmol)
        self._rdmol = rdmol
        self.topology_reconstruction_status = status
        self.__get_topology()
        return rdmol

    def _record_source_to_topology_atom_permutation(self) -> None:
        rdmol = self._rdmol
        self.source_to_topology_atom_permutation = None
        if rdmol is None or self.topology_reconstruction_status == "suspicious_fallback":
            return
        topology_atoms = [atom.GetAtomicNum() for atom in rdmol.GetAtoms()]
        if topology_atoms != self.atoms or rdmol.GetNumConformers() == 0:
            return
        source_coords = np.asarray(
            self.coords.to(atom_ureg.angstrom).magnitude,
            dtype=float,
        )
        topology_coords = np.asarray(rdmol.GetConformer().GetPositions(), dtype=float)
        # The Python MolGR backend round-trips coordinates through five-decimal XYZ text.
        if source_coords.shape != topology_coords.shape or not np.allclose(
            source_coords,
            topology_coords,
            rtol=0.0,
            atol=1e-5,
        ):
            return
        self.source_to_topology_atom_permutation = list(range(len(self.atoms)))

    def __get_topology(self) -> None:
        rdmol = self._rdmol
        if rdmol is None:
            return
        self._record_source_to_topology_atom_permutation()
        if self.source_to_topology_atom_permutation is None:
            return
        self.bonds = get_bond_pairs(rdmol)
        self.formal_charges = get_formal_charges(rdmol)
        self.formal_num_radicals = get_formal_num_radicals(rdmol)

    def _materialize_unitless_dump_with_unit_keys(self) -> None:
        """Build lazy topology before Pydantic snapshots regular fields."""

        super()._materialize_unitless_dump_with_unit_keys()
        _ = self.rdmol

    @property
    def omol(self) -> pybel.Molecule:
        """
        Get the openbabel molecule object.

        Returns:
            pybel.Molecule: The openbabel molecule object.
        """
        return rdmol_to_omol(self.rdmol)

    @property
    def rdmol(self) -> Chem.rdchem.Mol | None:
        """
        Get the rdkit molecule object.

        If reconstruction failed, return None.

        Returns:
            Union[Chem.rdchem.Mol,None]:
                The rdkit molecule object. If reconstruction failed, return None.
        """
        with self._topology_lock:
            if self._rdmol is not None:
                return self._rdmol
            if not self.atoms:
                return None
            if self.bonds:
                return self._materialize_topology_from_fields("provided")
            if self.topology_reconstruction_status in {"failed", "suspicious_fallback"}:
                # A suspicious MolGR graph is raw evidence, not a trusted
                # topology. Never rebuild a replacement from stale fields.
                return None
            status = self.topology_reconstruction_status
            if status == "provided" or status == "succeeded":
                try:
                    # A prewarmed graph may cross a process boundary without
                    # its private RDKit cache. Rebuild that cache locally and
                    # do not re-enter MolGR in the worker.
                    return self._materialize_topology_from_fields(status)
                except Exception as error:
                    self._apply_reconstruction_result(None, error=error)
                    return None
            return self._reconstruct_single_topology()

    def _reconstruct_single_topology(self) -> RdMol | None:
        options = self._resolve_reconstruction_options()
        attempt_started = False
        try:
            with native_reconstruction_guard():
                # A failed process-boundary check leaves status unset so a
                # later access can retry. Once the guard is acquired, any
                # interruption becomes a terminal failed reconstruction.
                attempt_started = True
                self._record_topology_reconstruction_provenance(options)
                reconstructed = xyz_to_rdmol(
                    self.to_XYZ(),
                    self.charge,
                    self.multiplicity,
                    backend=options[0],
                    make_dative_bonds=options[2],
                    make_stereochemistry=options[3],
                    config=MOLGR_CONFIG,
                )
            if reconstructed is None:
                raise ValueError("MolGR topology reconstruction returned no molecule")
            self._apply_reconstruction_result(reconstructed)
        except NativeReconstructionConcurrencyError:
            raise
        except Exception as error:
            self._apply_reconstruction_result(None, error=error)
        finally:
            if (
                attempt_started
                and self._rdmol is None
                and self.topology_reconstruction_status is None
            ):
                self.topology_reconstruction_status = "failed"
                self.__get_topology()
        return self._rdmol

    def _apply_batch_reconstruction_result(
        self,
        result: ReconstructionBatchResult,
    ) -> None:
        """Apply one MolGR batch result to this molecule's lazy topology state."""
        self._apply_reconstruction_result(
            result.result,
            status=getattr(result, "status", None),
        )

    @property
    def rdmol_no_conformer(self) -> Chem.rdchem.Mol:
        """
        Get the rdkit molecule object without conformer.

        Returns:
            Chem.rdchem.Mol: The rdkit molecule object without conformer.
        """
        if self.rdmol is None:
            raise ValueError("No RDKit molecule found.")
        rdmol = Chem.RWMol(self.rdmol)
        rdmol.RemoveAllConformers()
        return rdmol.GetMol()

    def to_SMILES(self) -> str:
        """
        Get the SMILES with explicit hydrogens.

        Returns:
            str: The SMILES.
        """
        if self._smiles_cache is not None:
            return self._smiles_cache
        if self.rdmol is None:
            return ""
        if smi := Chem.MolToSmiles(self.rdmol):
            self._smiles_cache = smi
            return smi
        moloplogger.error("SMILES building failed.")
        return ""

    def to_canonical_SMILES(self) -> str:
        """
        Get the SMILES with standardization.

        Returns:
            str: The standard SMILES.
        """
        if self._canonical_smiles_cache is not None:
            return self._canonical_smiles_cache
        smi = self.to_SMILES()
        if not smi:
            return ""
        try:
            canonical = canonical_smiles(smi)
        except Exception as e:
            # RDKit's CanonSmiles reparses the string and may pass None to
            # MolToSmiles when the reconstructed graph is not sanitizable.
            moloplogger.error(f"Canonical SMILES building failed: {e}")
            self._canonical_smiles_cache = ""
            return ""
        self._canonical_smiles_cache = canonical
        return canonical

    def to_InChI(self) -> str:
        """
        Get the InChI.

        Returns:
            str: The InChI.
        """
        if self.rdmol is None:
            return ""
        if inchi := Chem.MolToInchi(self.rdmol):
            return inchi  # type: ignore
        moloplogger.error("InChI building failed.")
        return ""

    def __hash__(self) -> int:
        return hash(str(self))

    def __len__(self) -> int:
        return len(self.atoms)

    def geometry_analysis(self, atom_idxs: Sequence[int], one_start=False) -> float:
        """
        Get the geometry infos among the atoms

        Parameters:
            atom_idxs (Sequence[int]):
                A Sequence of index of the atoms, starts from 0
            one_start (bool):
                If true, consider atom index starts from 1, so let index value subtracts 1 for all the atoms

        Returns:
            float:
                - If the length of atom_idxs is 2, the bond length with unit Angstrom between the two atoms will be returned.
                - If the length of atom_idxs is 3, the angle with unit degree between  the three atoms will be returned.
                - If the length of atom_idxs is 4, the dihedral angle with unit degree between the four atoms will be returned.
        """
        if one_start:
            atom_idxs = [atom_idx - 1 for atom_idx in atom_idxs]
        return get_geometry_info(self.rdmol, atom_idxs)

    def compare_rmsd(
        self,
        other: "Molecule",
        *,
        ignore_H: bool = False,
        prbId: int = -1,
        refId: int = -1,
        map: Sequence[tuple[int, int]] | None = None,
        maxMatches: int = 1000000,
        symmetrizeConjugatedTerminalGroups: bool = True,
        weights: Sequence[float] = [],
        numThreads: int = 1,
    ) -> float:
        """
        Calculate the RMSD between two molecules.

        Parameters:
            other (BaseMolFrame):
                The other molecule to compare.
            ignore_H (bool):
                If True, ignore the H atoms in the calculation. Default is False.
            prbId (int):
                The probe(self) molecule id. Default is -1.
            refId (int):
                The reference(other) molecule id. Default is -1.
            map (Optional[Sequence[tuple[int, int]]]):
                a list of lists of (probeAtomId,refAtomId) tuples with the atom-atom mappings of the two molecules.
                If not provided, these will be generated using a substructure search.
            maxMatches (int):
                if map isn't specified, this will be the max number of matches found in a SubstructMatch()
            symmetrizeConjugatedTerminalGroups (bool):
                if set, conjugated terminal functional groups (like nitro or carboxylate) will be considered symmetrically
            weights (Sequence[float]):
                weights for mapping
            numThreads (int):
                number of threads to use
        Returns:
            float: The RMSD value.
        """
        assert isinstance(other, Molecule), "The other object is not a BaseMolFrame object."

        return GetBestRMS(
            Chem.RemoveHs(self.rdmol) if ignore_H else self.rdmol,
            Chem.RemoveHs(other.rdmol) if ignore_H else other.rdmol,
            prbId=prbId,
            refId=refId,
            map=map,
            maxMatches=maxMatches,
            symmetrizeConjugatedTerminalGroups=symmetrizeConjugatedTerminalGroups,
            weights=weights,
            numThreads=numThreads,
        )

    @staticmethod
    def from_rdmol(rdmol: RdMol) -> "Molecule":
        return Molecule.model_validate(
            {
                "atoms": [atom.GetAtomicNum() for atom in rdmol.GetAtoms()],
                "coords": rdmol.GetConformer().GetPositions() * atom_ureg.angstrom,
                "charge": get_total_charge(rdmol),
                "multiplicity": get_total_multiplicity(rdmol),
                "bonds": get_bond_pairs(rdmol),
                "formal_charges": get_formal_charges(rdmol),
                "formal_num_radicals": get_formal_num_radicals(rdmol),
            }
        )

    @staticmethod
    def from_coords(
        atom_symbols: Sequence[str | int],
        coords: np.ndarray,
        charge: int = 0,
        multiplicity: int = 1,
        bonds: list[tuple[int, int, int, int]] | None = None,
        formal_charges: list[int] | None = None,
        formal_num_radicals: list[int] | None = None,
    ) -> "Molecule":
        return Molecule.model_validate(
            {
                "atoms": [Chem.Atom(atom).GetAtomicNum() for atom in atom_symbols],
                "coords": coords * atom_ureg.angstrom,
                "charge": charge,
                "multiplicity": multiplicity,
            }
            | ({"bonds": bonds} if bonds else {})
            | ({"formal_charges": formal_charges} if formal_charges else {})
            | ({"formal_num_radicals": formal_num_radicals} if formal_num_radicals else {})
        )

    def replace_substituent(
        self,
        query: str | RdMol,
        replacement: str | RdMol,
        bind_idx: int | None = None,
        replace_all=False,
        attempt_num=10,
        crowding_threshold=0.75,
        angle_split=10,
        randomSeed=114514,
        start_idx: int | None = None,
        end_idx: int | None = None,
        *,
        stereo_policy: StereoPolicy = "preserve",
    ) -> "Molecule":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi,
        and the new substituent is defined by the replacement_smi.

        Parameters:
            query (str | RdMol):
                The SMARTS or Mol object to query the substituent in the original molecule.
            replacement (str | RdMol):
                The SMILES or Mol object of the new substituent. It must
                contain one RDKit ``[*]`` attachment marker.
            bind_idx (int):
                The atom index on the retained scaffold to which the substituent is attached.
                Use it to disambiguate multiple pendant matches.
            replace_all (bool):
                If True, replace all the substituent queried in the original molecule.
            attempt_num (int):
                Max attempt times to replace the substituent. Each time a new substituent conformation
                will be used for substitution.
            crowding_threshold (float):
                The threshold of crowding. If the new substituent is too crowded
                `d(a-b) < threshold * (R(a)+R(b))`, the substitution will be rejected.
            angle_split (int):
                Decide how many equal parts of 360° you want to divide. The larger the number the finer
                the rotation will be attempted but the slower the calculation will be.
            randomSeed (int):
                The random seed.
            start_idx (int):
                If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
                key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
            end_idx (int):
                If both `start_idx` and `end_idx` are specified, simply ignore the `query`, break the
                key between `start_idx` and `end_idx` and replace the base group where `end_idx` is located
            stereo_policy (str):
                Double-bond stereo policy: "preserve", "E", "Z", or "any".

        Returns:
            Molecule: The new molecule with the substituent replaced.
        """
        new_mol = replace_substituent(
            self.rdmol,
            query=query,
            replacement=replacement,
            bind_idx=bind_idx,
            replace_all=replace_all,
            attempt_num=attempt_num,
            crowding_threshold=crowding_threshold,
            angle_split=angle_split,
            randomSeed=randomSeed,
            start_idx=start_idx,
            end_idx=end_idx,
            stereo_policy=stereo_policy,
        )
        return self.from_rdmol(new_mol)

    def reset_atom_index(
        self,
        mapping_smarts: str | None = None,
        *,
        mapping_indice: Sequence[int] | None = None,
    ) -> Optional["Molecule"]:
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts (str):
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according
                to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.
            mapping_indice (Sequence[int]):
                The indices of the atoms to be renumbered. The relative order of the remaining atoms
                remains unchanged.
                e.g. atoms = [0, 1, 2, 3, 4, 5]; mapping = [3, 5, 4] means the first atom in the new
                molecule is mapped to the third atom in the original molecule,
                the second atom is mapped to the first atom, and the third atom is mapped to the second
                atom. Result is [3, 5, 4, 0, 1, 2].
        Returns:
            Molecule: The new Molecule with the atom index reset.
        """
        if self.rdmol is None:
            raise ValueError("No RDKit molecule found.")
        if mapping_smarts is not None:
            smarts = Chem.MolFromSmarts(mapping_smarts)
            if not self.rdmol.HasSubstructMatch(smarts):
                raise ValueError(
                    f"Failed to match {self.to_canonical_SMILES()} with {mapping_smarts}"
                )
            mappings: list[Sequence[int]] = self.rdmol.GetSubstructMatches(smarts)
            if len(mappings) > 1:
                moloplogger.warning(
                    f"Multiple matches found in {self.to_canonical_SMILES()} with {mapping_smarts}, using the first one"
                )
            mapping = mappings[0]
            rdmol = reset_atom_index(self.rdmol, mapping)
            return self.from_rdmol(rdmol)
        elif mapping_indice is not None:
            assert max(mapping_indice) < self.rdmol.GetNumAtoms(), "Invalid mapping index"
            rdmol = reset_atom_index(self.rdmol, mapping_indice)
            return self.from_rdmol(rdmol)
        return None

    def standard_orient(
        self,
        anchor_list: Sequence[int],
    ) -> "Molecule":
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and
        `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches
                the origin.
            - `rotate_anchor_to_axis`: Rotate the specified second atom along the axis passing through
                the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_plane`: Rotate along the axis passing through the origin so that the
                specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list (Sequence[int]):
                A Sequence of indices of the atoms to be translated to origin, rotated to X axis,
                    and rotated again to XY face:

            - If length is 1, execute `translate_anchor`
            - If length is 2, execute `translate_anchor` and `rotate_anchor_to_axis`
            - If length is 3, execute `translate_anchor`, `rotate_anchor_to_axis` and `rotate_anchor_to_plane`
            - If the length of the input `anchor_list` is greater than 3, subsequent atomic numbers are not considered.
        Returns:
            BaseMolFrameParser: The new parser.
        """
        if self.rdmol is None:
            raise ValueError("No RDKit molecule found.")
        mol = Chem.RWMol(self.rdmol)
        standard_orient(mol, anchor_list)
        return self.from_rdmol(mol)

    def to_summary_dict(self, brief: bool = True, **kwargs) -> SummaryDict:
        return {
            summary_column("General", "Charge"): self.charge,
            summary_column("General", "Multiplicity"): self.multiplicity,
            summary_column("General", "CanonicalSMILES"): self.to_canonical_SMILES(),
            summary_column("General", "NumAtoms"): len(self.atoms),
        }

    def _render(self, **kwargs) -> str:
        """
        Render the Molecule as a string.

        Returns:
            str: The rendered Molecule.
        """
        return self.to_canonical_SMILES()

    def to_internal_coords(self) -> InternalCoords:
        """
        Convert the cartesian coordinates to internal coordinates.

        Returns:
            Molecule: The new Molecule with internal coordinates.
        """
        return InternalCoords.from_cartesian_coords(self.atom_symbols, self.coords)

    def to_XYZ(self) -> str:
        return f"{len(self.atoms)}\n\n" + "\n".join(
            [
                f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                for atom, (x, y, z) in zip(self.atom_symbols, self.coords.m, strict=True)
            ]
        )


_DEFAULT_TOPOLOGY_BATCH_SIZE = 256


def _resolve_batch_worker_count(
    max_workers: int | None,
    backend: Literal["cpp", "python"],
) -> int | None:
    if backend == "python":
        return 1
    # MolGR owns the native pool, so its automatic worker count must use the
    # stricter native-runtime budget rather than the general joblib limit.
    worker_count = (
        molopconfig.effective_molgr_max_jobs
        if max_workers is None
        else molopconfig.set_molgr_n_jobs(max_workers)
    )
    configured_limit = MOLGR_CONFIG.cpp_backend.max_threads
    if configured_limit is not None:
        worker_count = (
            configured_limit if worker_count is None else min(worker_count, configured_limit)
        )
    # MolGR deliberately defaults to one Open Babel thread on Windows.  Keep
    # that safety boundary even when MolOP receives an explicit n_jobs value.
    if sys.platform == "win32" and configured_limit == 1:
        worker_count = 1
    return worker_count


def reconstruct_topologies_batch(
    molecules: Iterable[Molecule],
    *,
    backend: Literal["cpp", "python"] | None = None,
    reconstruction_failure_policy: Literal["raise", "return_suspicious"] | None = None,
    make_dative_bonds: bool | None = None,
    make_stereochemistry: bool | None = None,
    max_workers: int | None = None,
    queue_size: int = 16,
    ordered: bool = False,
    raise_on_error: bool = False,
    batch_size: int = _DEFAULT_TOPOLOGY_BATCH_SIZE,
    retain_results: bool = True,
) -> list[ReconstructionBatchResult]:
    """Reconstruct coordinate-only molecules through MolGR's native batch API.

    The native C++ backend owns its worker pool.  This helper deliberately
    gathers all eligible molecules in the calling process instead of invoking
    :mod:`joblib` once per molecule, which would create nested native pools.
    Results are returned in the input molecule order even when the native
    iterator is configured for unordered completion.
    """

    if max_workers is not None and max_workers < 1:
        raise ValueError("max_workers must be >= 1 when provided")
    if queue_size < 1:
        raise ValueError("queue_size must be >= 1")
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    if is_loky_worker():
        raise RuntimeError(
            "MolGR topology reconstruction is forbidden in a loky worker; "
            "prewarm the topology in the parent process first."
        )

    completed: list[tuple[int, ReconstructionBatchResult]] = []
    molecule_iter = iter(molecules)
    input_index = 0
    while chunk := list(islice(molecule_iter, batch_size)):
        indexed_chunk = list(enumerate(chunk, start=input_index))
        input_index += len(chunk)
        unique_molecules = {
            id(molecule): molecule
            for _, molecule in indexed_chunk
            if isinstance(molecule, Molecule)
        }
        grouped: dict[_ReconstructionOptions, list[tuple[int, Molecule]]] = {}

        # Lock acquisition is deterministic across concurrent batch calls.
        # The locks remain held while native work applies its results, so a
        # concurrent lazy reader waits for the completed cache instead of
        # observing a transient ``None``.
        with ExitStack() as molecule_locks:
            for molecule in sorted(unique_molecules.values(), key=id):
                molecule_locks.enter_context(molecule._topology_lock)

            for index, molecule in indexed_chunk:
                if not isinstance(molecule, Molecule):
                    continue
                # Existing graphs and previously attempted lazy reconstructions
                # remain untouched; callers can construct a fresh model to retry.
                if molecule._rdmol is not None or molecule.bonds or not molecule.atoms:
                    continue
                if molecule.topology_reconstruction_status is not None:
                    continue

                options = molecule._resolve_reconstruction_options(
                    backend=backend,
                    reconstruction_failure_policy=reconstruction_failure_policy,
                    make_dative_bonds=make_dative_bonds,
                    make_stereochemistry=make_stereochemistry,
                )
                grouped.setdefault(
                    options,
                    [],
                ).append((index, molecule))

            if not grouped:
                continue
            for (
                group_backend,
                group_failure_policy,
                group_dative_bonds,
                group_stereochemistry,
            ), entries in grouped.items():
                worker_count = _resolve_batch_worker_count(max_workers, group_backend)
                received_indices: set[int] = set()
                molecules_by_index = dict(entries)
                attempt_started = False
                try:
                    with native_reconstruction_guard():
                        # Do not mutate lazy reconstruction state until the
                        # process-wide native boundary is acquired.  A
                        # temporary contention error must leave all models
                        # retryable instead of sealing them as failed.
                        attempt_started = True
                        for _, molecule in entries:
                            molecule._record_topology_reconstruction_provenance(
                                (
                                    group_backend,
                                    group_failure_policy,
                                    group_dative_bonds,
                                    group_stereochemistry,
                                )
                            )
                        requests = [
                            ReconstructionBatchRequest(
                                molecule.to_XYZ(),
                                total_charge=molecule.charge,
                                spin_multiplicity=molecule.multiplicity,
                            )
                            for _, molecule in entries
                        ]
                        request_indices = {
                            id(request): item_index
                            for request, (item_index, _) in zip(requests, entries, strict=True)
                        }
                        for result in iter_xyz_to_rdmol_batch(
                            requests,
                            backend=group_backend,
                            max_workers=worker_count,
                            queue_size=queue_size,
                            ordered=ordered,
                            make_dative_bonds=group_dative_bonds,
                            make_stereochemistry=group_stereochemistry,
                            config=MOLGR_CONFIG,
                            raise_on_error=raise_on_error,
                        ):
                            item_index = request_indices.get(id(result.input))
                            if item_index is None:
                                raise RuntimeError(
                                    "MolGR batch returned a result for an unknown request"
                                )
                            if item_index in received_indices:
                                raise RuntimeError(
                                    "MolGR batch returned a duplicate result for one request"
                                )
                            molecule = molecules_by_index[item_index]
                            molecule._apply_batch_reconstruction_result(result)
                            received_indices.add(item_index)
                            if retain_results:
                                completed.append((item_index, result))
                except BaseException:
                    if attempt_started:
                        for item_index, molecule in entries:
                            if item_index not in received_indices:
                                molecule.topology_reconstruction_status = "failed"
                    raise

                # A malformed or interrupted iterator must not leave an
                # unprocessed coordinate-only model able to call native MolGR
                # later from a worker.
                for item_index, molecule in entries:
                    if item_index not in received_indices:
                        molecule.topology_reconstruction_status = "failed"

    if not retain_results:
        return []
    completed.sort(key=lambda item: item[0])
    return [result for _, result in completed]
