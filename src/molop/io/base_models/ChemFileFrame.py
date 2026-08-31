"""
Author: TMJ
Date: 2025-07-28 18:43:45
LastEditors: TMJ
LastEditTime: 2026-08-30 23:26:00
Description: 请填写简介
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from io import StringIO
from pathlib import Path
from typing import Any, ClassVar, Generic, Literal, Protocol, TypeVar, cast

import numpy as np
import numpy.typing as npt
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field, model_validator
from rdkit import Chem
from scipy.sparse import coo_matrix
from typing_extensions import Self

from molop.config import molopconfig, moloplogger
from molop.io.base_models.DataClasses import (
    NMR,
    BondOrders,
    ChargeSpinPopulations,
    ElectronicStates,
    Energies,
    ExcitedStateRequest,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    MolecularOrbitals,
    MultireferenceRequest,
    MultireferenceResult,
    Polarizability,
    QMModelChemistry,
    QMResourceRequest,
    QMTaskRequest,
    SinglePointProperties,
    Status,
    ThermalInformations,
    TotalSpin,
    Vibration,
    Vibrations,
)
from molop.io.base_models.Molecule import Molecule, reconstruct_topologies_batch
from molop.io.base_models.source import (
    ParseCompleteness,
    ParseDiagnostic,
    ParsePresence,
    SourceSpan,
)
from molop.io.base_models.summary import SummaryDict, summary_column, summary_item
from molop.structure.StructureTransformation import check_crowding
from molop.unit import atom_ureg
from molop.utils.progressbar import NativeReconstructionConcurrencyError
from molop.utils.types import OMol, PintArrayN, PintArrayNx3, PintSquareMatrix, RdMol
from molop.visualization.animation import AnimationFormat, render_molecule_animation


class _HasCoords(Protocol):
    frame_id: int
    atoms: list[int]

    @property
    def atom_symbols(self) -> list[str]: ...

    coords: PintArrayNx3
    charge: int
    multiplicity: int
    _default_units: dict[str, UnitLike]

    @property
    def rdmol(self) -> RdMol | None: ...

    @property
    def omol(self) -> OMol | None: ...

    def to_SMILES(self) -> str: ...

    def to_canonical_SMILES(self) -> str: ...


class _HasKeywords(_HasCoords):
    keywords: str
    basis_set: str
    functional: str
    method: str


class _HasVibrations(_HasKeywords):
    vibrations: Vibrations | None


def _legacy_method_from_model_chemistry(model: QMModelChemistry) -> str:
    if model.method_family:
        if model.method:
            if model.method.lower() == model.method_family.lower():
                return model.method_family
            if model.method_family.upper() == "DFT":
                return "DFT"
            return model.method
        return model.method_family
    if model.functional:
        return "DFT"
    if model.method:
        return model.method
    return ""


def _canonical_smiles_from_rdmol(rdmol: RdMol | None) -> str:
    if rdmol is None:
        return ""
    smiles = Chem.MolToSmiles(rdmol)
    if not smiles:
        return ""
    try:
        return Chem.CanonSmiles(smiles)
    except Exception:
        return smiles


def _topology_frequency_key(rdmol: RdMol) -> str | bytes:
    """Return a conformer-independent key for TS endpoint voting."""

    try:
        smiles = Chem.MolToSmiles(rdmol, canonical=True, isomericSmiles=False)
        if smiles:
            return smiles
    except Exception:
        pass
    topology = Chem.Mol(rdmol)
    topology.RemoveAllConformers()
    return topology.ToBinary()


def _most_frequent_topology(candidates: Sequence[RdMol], *, side: str) -> RdMol:
    """Select a side's topology mode and retain its largest-amplitude conformer."""

    if not candidates:
        raise ValueError(f"Failed to reconstruct any {side}-space TS endpoint candidates")

    # Candidates arrive in ascending amplitude order. Updating the representative
    # on every hit preserves the largest-amplitude conformer for the winning graph.
    grouped: dict[str | bytes, tuple[int, int, RdMol]] = {}
    for index, rdmol in enumerate(candidates):
        key = _topology_frequency_key(rdmol)
        count = grouped[key][0] + 1 if key in grouped else 1
        grouped[key] = (count, index, Chem.Mol(rdmol))

    _, _, representative = max(grouped.values(), key=lambda item: (item[0], item[1]))
    return representative


def _bond_change_pairs(first: RdMol, second: RdMol) -> set[tuple[int, int]]:
    """Return atom pairs whose bond presence or type differs between two graphs."""

    if first.GetNumAtoms() != second.GetNumAtoms():
        return set()

    def bond_types(molecule: RdMol) -> dict[tuple[int, int], Chem.BondType]:
        result: dict[tuple[int, int], Chem.BondType] = {}
        for bond in molecule.GetBonds():
            start_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            if start_atom_idx > end_atom_idx:
                start_atom_idx, end_atom_idx = end_atom_idx, start_atom_idx
            result[(start_atom_idx, end_atom_idx)] = bond.GetBondType()
        return result

    first_bonds = bond_types(first)
    second_bonds = bond_types(second)
    return {
        atom_pair
        for atom_pair in first_bonds.keys() | second_bonds.keys()
        if first_bonds.get(atom_pair) != second_bonds.get(atom_pair)
    }


def _bond_change_atom_indices(first: RdMol, second: RdMol) -> set[int]:
    """Return atoms incident to bonds that differ between two endpoint graphs."""

    return {
        atom_index for atom_pair in _bond_change_pairs(first, second) for atom_index in atom_pair
    }


def _method_family_allows_functional(method_family: str | None) -> bool:
    if method_family is None:
        return False
    normalized = method_family.upper().replace("_", "-").replace(" ", "-")
    return normalized in {"DFT", "DOUBLE-HYBRID", "DOUBLE-HYBRID-DFT"}


def _legacy_functional_can_backfill(
    model: QMModelChemistry,
    legacy_method: str,
) -> bool:
    if _method_family_allows_functional(model.method_family):
        return True
    if model.method_family is not None:
        return False
    return _method_family_allows_functional(legacy_method)


def _backfill_common_qm_containers_from_legacy(target: Any) -> None:
    model = target.model_chemistry
    legacy_method = getattr(target, "method", "")
    legacy_functional = getattr(target, "functional", "")
    legacy_basis_set = getattr(target, "basis_set", "")
    legacy_keywords = getattr(target, "keywords", "")

    if legacy_method:
        if model.method is None:
            model.method = legacy_method
        if model.method_family is None and legacy_method.upper() in {
            "DFT",
            "HF",
            "MP2",
            "CCSD",
        }:
            model.method_family = legacy_method.upper()
    if (
        legacy_functional
        and model.functional is None
        and _legacy_functional_can_backfill(model, legacy_method)
    ):
        model.functional = legacy_functional
    if legacy_basis_set and model.basis_set is None:
        model.basis_set = legacy_basis_set
    if legacy_keywords and not model.raw_keywords:
        model.raw_keywords = legacy_keywords

    auxiliary_basis_set = getattr(target, "auxiliary_basis_set", "")
    if auxiliary_basis_set and model.auxiliary_basis_set is None:
        model.auxiliary_basis_set = auxiliary_basis_set

    dispersion_correction = getattr(target, "dispersion_correction", "")
    if dispersion_correction and model.dispersion_correction is None:
        model.dispersion_correction = dispersion_correction

    resource = target.resource_request
    if getattr(target, "request_num_cpu", None) is not None and resource.num_cpu is None:
        resource.num_cpu = target.request_num_cpu
    if getattr(target, "request_memory", None) is not None and resource.memory is None:
        resource.memory = target.request_memory
    resources_raw = getattr(target, "resources_raw", "")
    if resources_raw and not resource.raw:
        resource.raw = resources_raw


def _project_common_qm_fields(target: Any) -> None:
    model = target.model_chemistry
    if model.raw_keywords:
        target.keywords = model.raw_keywords
    target.method = _legacy_method_from_model_chemistry(model)
    target.functional = model.functional or ""
    target.basis_set = model.basis_set or ""

    if hasattr(target, "auxiliary_basis_set"):
        target.auxiliary_basis_set = model.auxiliary_basis_set or ""
    if hasattr(target, "dispersion_correction"):
        target.dispersion_correction = model.dispersion_correction or ""

    resource = target.resource_request
    target.request_num_cpu = resource.num_cpu
    target.request_memory = resource.memory
    if resource.raw:
        target.resources_raw = resource.raw


def _process_bond_helper(
    rwmol: Chem.RWMol,
    start_atom_idx: int,
    end_atom_idx: int,
    other_rdmol: RdMol,
) -> None:
    """
    Helper function to process bonds and set bond types to zero if necessary.

    Parameters:
        rwmol (Chem.RWMol):
            The RDKit molecule to modify.
        start_atom_idx (int):
            The index of the start atom in the bond.
        end_atom_idx (int):
            The index of the end atom in the bond.
        other_rdmol (RdMol):
            The other RDKit molecule to compare against.
    """
    bond_1 = rwmol.GetBondBetweenAtoms(start_atom_idx, end_atom_idx)
    bond_2 = other_rdmol.GetBondBetweenAtoms(start_atom_idx, end_atom_idx)
    if bond_1 is None:
        rwmol.AddBond(start_atom_idx, end_atom_idx, Chem.BondType.ZERO)
    elif bond_2 is None or bond_1.GetBondType() != bond_2.GetBondType():
        bond_1.SetBondType(Chem.BondType.ZERO)


ChemFileFrame = TypeVar("ChemFileFrame", bound="BaseChemFileFrame")


class BaseChemFileFrame(Molecule, Generic[ChemFileFrame]):
    frame_id: int = Field(default=0, description="Frame ID")
    frame_content: str = Field(default="", repr=False, exclude=True)
    source_span: SourceSpan | None = Field(
        default=None,
        description="Optional half-open source byte, character, and line offsets",
        exclude_if=lambda value: value is None,
    )
    source_block_sha256: str | None = Field(
        default=None,
        pattern=r"^[0-9a-f]{64}$",
        description="Optional SHA-256 digest of the source block",
        exclude_if=lambda value: value is None,
    )
    segment_index: int | None = Field(
        default=None,
        ge=0,
        exclude_if=lambda value: value is None,
    )
    segment_frame_index: int | None = Field(
        default=None,
        ge=0,
        exclude_if=lambda value: value is None,
    )
    file_frame_index: int | None = Field(
        default=None,
        ge=0,
        description="Stable frame ordinal in the complete located source artifact",
        exclude_if=lambda value: value is None,
    )
    parse_presence: dict[str, ParsePresence] = Field(
        default_factory=dict,
        description="Presence state for optional scientific fields assessed during parsing",
        exclude_if=lambda value: len(value) == 0,
    )
    parse_diagnostics: list[ParseDiagnostic] = Field(
        default_factory=list,
        description="Structured frame-scoped parser diagnostics",
        exclude_if=lambda value: len(value) == 0,
    )
    parse_completeness: ParseCompleteness = Field(
        default=ParseCompleteness.NOT_ASSESSED,
        description="Completeness of the requested scientific parsing work for this frame",
        exclude_if=lambda value: value is ParseCompleteness.NOT_ASSESSED,
    )
    _frame_type: str = PrivateAttr(default="")
    _next_frame: ChemFileFrame | None = PrivateAttr(default=None)
    _prev_frame: ChemFileFrame | None = PrivateAttr(default=None)

    @property
    def next(self) -> ChemFileFrame | None:
        return self._next_frame

    @property
    def prev(self) -> ChemFileFrame | None:
        return self._prev_frame

    @property
    def is_error(self) -> bool | None:
        """
        Abstrcact method to check if the current frame is an error.
        The details are implemented in the derived classes.
        """

    @property
    def is_normal(self) -> bool | None: ...

    @property
    def is_TS(self) -> bool | None:
        """
        Check if the molecule is a TS. Can not check if the molecule
        is a TS without frequency information. Thus this function returns False.
        """

    @property
    def is_optimized(self) -> bool | None:
        """
        Check if the molecule is optimized.
        """

    def to_summary_dict(self, brief: bool = True, **kwargs) -> SummaryDict:
        return {
            **super().to_summary_dict(brief=brief, **kwargs),
            summary_column("General", "FrameID"): self.frame_id,
        }

    def log_with_file_info(self, content: str, level: str = "info"):
        if file_name := getattr(self, "filename", None):
            getattr(moloplogger, level)(f"{file_name} - Frame {self.frame_id}: {content}")

    def release_frame_content(self) -> None:
        self.frame_content = ""


class BaseCoordsFrame(BaseChemFileFrame[ChemFileFrame]):
    coordinate_source: str | None = Field(default=None, exclude_if=lambda value: value is None)
    coordinate_provenance: str | None = Field(
        default=None,
        exclude_if=lambda value: value is None,
    )
    coordinate_decimal_places: int | None = Field(
        default=None,
        ge=0,
        le=18,
        exclude_if=lambda value: value is None,
    )


class BaseQMInputFrame(BaseCoordsFrame[ChemFileFrame]):
    """Frame type for *QM input* files that contain coordinates + input metadata.

    This sits conceptually between `BaseCoordsFrame` and `BaseCalcFrame`:
    - Has coordinates (like coords formats)
    - Carries lightweight method/keyword/resource metadata (like QM calculations)
    - Does NOT imply that QM output properties (energies/forces/vibrations/...) exist
    """

    # QM software
    qm_software: str = Field(
        default="",
        description="QM software used for this input (e.g., gaussian/orca)",
    )
    qm_software_version: str = Field(
        default="",
        description="QM software version (if known)",
    )

    # QM input parameters
    keywords: str = Field(
        default="",
        description="Input keywords / route section (raw or normalized)",
    )
    method: str = Field(
        default="",
        description="QM method (best-effort, may be empty for raw-only inputs)",
    )
    basis_set: str = Field(
        default="",
        description="Basis set (best-effort, may be empty for raw-only inputs)",
    )
    functional: str = Field(
        default="",
        description="Functional (best-effort, may be empty for raw-only inputs)",
    )
    model_chemistry: QMModelChemistry = Field(
        default_factory=QMModelChemistry,
        description="Structured model chemistry information",
    )
    task_requests: list[QMTaskRequest] = Field(
        default_factory=list,
        description="Structured calculation task requests",
    )
    excited_state_requests: list[ExcitedStateRequest] = Field(
        default_factory=list,
        description="Structured excited-state task requests",
    )
    multireference_requests: list[MultireferenceRequest] = Field(
        default_factory=list,
        description="Structured multi-reference task requests",
    )

    # Resources
    resources_raw: str = Field(
        default="",
        description="Raw resource directives from the input (preservation-only)",
    )
    request_num_cpu: int | None = Field(
        default=None,
        description="Number of CPUs used for the QM calculation",
    )
    request_memory: PlainQuantity | None = Field(
        default=None,
        description="Memory used for the QM calculation",
    )
    resource_request: QMResourceRequest = Field(
        default_factory=QMResourceRequest,
        description="Structured resource request",
    )

    def backfill_common_qm_containers_from_legacy(self) -> None:
        """Fill structured QM containers from legacy flat fields when they are empty."""
        _backfill_common_qm_containers_from_legacy(self)

    def project_common_qm_fields(self) -> None:
        """Project authoritative structured QM containers onto legacy flat fields."""
        _project_common_qm_fields(self)

    def refresh_common_qm_containers(self) -> None:
        self.backfill_common_qm_containers_from_legacy()


class BaseCalcFrame(BaseQMInputFrame[ChemFileFrame]):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "coords": atom_ureg.angstrom,
        "forces": atom_ureg.Unit("hartree / bohr"),
        "rotation_constants": atom_ureg.Unit("gigahertz"),
        "running_time": atom_ureg.Unit("second"),
        "temperature": atom_ureg.Unit("K"),
        "electron_temperature": atom_ureg.Unit("K"),
    }
    # Note: QM input metadata (keywords/method/basis_set/functional/resources_raw)
    # lives on BaseQMInputFrame.
    # solvation
    solvent: ImplicitSolvation | None = Field(
        default=None,
        description="Solvent used in the QM calculation",
    )
    # physical settings
    temperature: PlainQuantity | None = Field(
        default=None,
        description="Temperature used in the QM calculation, unit is `K`",
    )
    pressure: PlainQuantity | None = Field(
        default=None,
        description="Pressure used in the QM calculation, unit is `atm`",
    )
    # QM properties
    forces: PintArrayNx3 | None = Field(
        default=None,
        description="Forces of each atom, unit is `hartree/bohr`.\n"
        "In Gaussian, the extracted forces data are all calculated using the "
        "input coordinates as a reference.",
    )
    hessian: PintSquareMatrix | None = Field(
        default=None,
        description="Hessian matrix of the QM calculation, unit is `hartree/bohr^2`.\n"
        "In Gaussian, the extracted hessian data are all calculated using the "
        "input coordinates as a reference.",
    )
    forces_axis_order: tuple[Literal["atom"], Literal["cartesian"]] | None = Field(
        default=None,
        description="Axis order of forces, with Cartesian components ordered x, y, z",
        exclude_if=lambda value: value is None,
    )
    forces_atom_order: Literal["source"] | None = Field(
        default=None,
        description="Atom ordering used by forces",
        exclude_if=lambda value: value is None,
    )
    forces_orientation: Literal["input", "standard", "source", "unknown"] | None = Field(
        default=None,
        description="Cartesian orientation used by forces",
        exclude_if=lambda value: value is None,
    )
    hessian_axis_order: tuple[Literal["atom_cartesian"], Literal["atom_cartesian"]] | None = Field(
        default=None,
        description="Axis order of the flattened atom-major Cartesian Hessian",
        exclude_if=lambda value: value is None,
    )
    hessian_atom_order: Literal["source"] | None = Field(
        default=None,
        description="Atom ordering used by the Hessian",
        exclude_if=lambda value: value is None,
    )
    hessian_orientation: Literal["input", "standard", "source", "unknown"] | None = Field(
        default=None,
        description="Cartesian orientation used by the Hessian",
        exclude_if=lambda value: value is None,
    )
    rotation_constants: PintArrayN | None = Field(
        default=np.array([]) * atom_ureg.gigahertz,
        description="Rotational constants, unit is `gigahertz`",
    )
    energies: Energies | None = Field(
        default=None,
        description="Energies",
    )
    thermal_informations: ThermalInformations | None = Field(
        default=None,
        description="Thermal Energies",
    )
    molecular_orbitals: MolecularOrbitals | None = Field(
        default=None,
        description="Molecular Orbitals",
    )
    vibrations: Vibrations | None = Field(default=None, description="vibrations")
    charge_spin_populations: ChargeSpinPopulations | None = Field(
        default=None, description="Charge and spin populations"
    )
    polarizability: Polarizability | None = Field(
        default=None,
        description="Polarizability of the molecule.\n"
        "In Gaussian, the extracted polarization-related data are all calculated using the "
        "input coordinates as a reference.",
    )
    nmr: NMR | None = Field(default=None, description="NMR shielding and coupling properties")
    bond_orders: BondOrders | None = Field(default=None, description="Bond orders of the molecule")
    total_spin: TotalSpin | None = Field(default=None, description="Total spin of the molecule")
    single_point_properties: SinglePointProperties | None = Field(
        default=None,
        description="Single point properties of the molecule",
    )
    electronic_states: ElectronicStates | None = Field(
        default=None,
        description="Electronic-state resolved properties",
    )
    multireference_result: MultireferenceResult | None = Field(
        default=None,
        description="Multi-reference calculation results",
    )
    status: Status | None = Field(default=None, description="Status of the frame")
    geometry_optimization_status: GeometryOptimizationStatus | None = Field(
        default=None,
        description="Geometry optimization status",
    )
    running_time: PlainQuantity | None = Field(
        default=None,
        description="Running time of the QM calculation, unit is `second`",
    )
    frame_role: str | None = Field(default=None, exclude_if=lambda value: value is None)
    force_source_field: str | None = Field(default=None, exclude_if=lambda value: value is None)
    force_transformation: str | None = Field(default=None, exclude_if=lambda value: value is None)

    @model_validator(mode="after")
    def validate_cartesian_array_conventions(self) -> Self:
        atom_count = len(self.atoms)
        force_metadata = (
            self.forces_axis_order,
            self.forces_atom_order,
            self.forces_orientation,
        )
        if self.forces is None:
            if any(value is not None for value in force_metadata):
                raise ValueError("forces convention metadata require forces")
        else:
            if tuple(self.forces.shape) != (atom_count, 3):
                raise ValueError("forces must have shape (N, 3) in source atom order")
            self.forces_axis_order = self.forces_axis_order or ("atom", "cartesian")
            self.forces_atom_order = self.forces_atom_order or "source"
            self.forces_orientation = self.forces_orientation or "unknown"

        hessian_metadata = (
            self.hessian_axis_order,
            self.hessian_atom_order,
            self.hessian_orientation,
        )
        if self.hessian is None:
            if any(value is not None for value in hessian_metadata):
                raise ValueError("Hessian convention metadata require hessian")
        else:
            cartesian_size = atom_count * 3
            if tuple(self.hessian.shape) != (cartesian_size, cartesian_size):
                raise ValueError("hessian must have shape (3N, 3N) in source atom order")
            self.hessian_axis_order = self.hessian_axis_order or (
                "atom_cartesian",
                "atom_cartesian",
            )
            self.hessian_atom_order = self.hessian_atom_order or "source"
            self.hessian_orientation = self.hessian_orientation or "unknown"

        if self.nmr is not None:
            for shielding in self.nmr.shielding_tensors:
                if shielding.atom_index >= atom_count:
                    raise ValueError("NMR shielding atom index is outside the frame atom range")
                expected_symbol = self.atom_symbols[shielding.atom_index]
                if shielding.atom_symbol != expected_symbol:
                    raise ValueError(
                        "NMR shielding atom symbol does not match the frame source atom order"
                    )
            for atom_index in self.nmr.coupling_atom_indices:
                if atom_index < 0 or atom_index >= atom_count:
                    raise ValueError("NMR coupling atom index is outside the frame atom range")
        return self

    def qm_embedded_rdmol(
        self,
        embed_populations: bool = True,
        embed_bond_orders: bool = True,
        embed_nmr: bool = True,
    ) -> RdMol | None:
        """
        Store atom- and bond-resolved QM properties in an RDKit molecule.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        If `embed_populations` is True, this function will use all population properties in the `charge_spin_populations` field
        to generate the population embedded rdkit molecule object.
        If `embed_bond_orders` is True, this function will use all bond order properties in the `bond_orders` field
        to generate the bond order embedded rdkit molecule object.
        If `embed_nmr` is True, per-atom shielding scalars, principal values, and tensor components
        are embedded in ppm. Spin-spin coupling matrices remain frame-level atom-pair data.

        Parameters:
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.
            embed_nmr (bool): If True, embed per-atom NMR shielding properties. Defaults to True.

        Returns:
            Optional[RdMol]: The population embedded rdkit molecule object.
        """
        if self.rdmol is None:
            return None
        rwmol = Chem.RWMol(self.rdmol)
        if embed_populations and self.charge_spin_populations is not None:
            for population, series in self.charge_spin_populations.population_items():
                for atom_idx, pop in enumerate(series.values):
                    rwmol.GetAtomWithIdx(atom_idx).SetDoubleProp(
                        f"{population}_by_{self.qm_software}".upper(), pop
                    )
                Chem.CreateAtomDoublePropertyList(
                    rwmol, f"{population}_by_{self.qm_software}".upper()
                )
        if embed_bond_orders and self.bond_orders is not None:
            bond_orders: dict[str, npt.NDArray[np.floating]] = self.bond_orders.model_dump(
                exclude_defaults=True
            )
            for bond_order, bond_order_matrix in bond_orders.items():
                coo = coo_matrix(bond_order_matrix)
                for i, j, val in zip(coo.row, coo.col, coo.data, strict=True):
                    if rwmol.GetBondBetweenAtoms(i, j) is None:
                        rwmol.AddBond(i, j, Chem.BondType.UNSPECIFIED)
                    rwmol.GetBondBetweenAtoms(i, j).SetDoubleProp(
                        f"{bond_order}_by_{self.qm_software}".upper(), val
                    )
                Chem.CreateBondDoublePropertyList(
                    rwmol, f"{bond_order}_by_{self.qm_software}".upper()
                )
        if embed_nmr and self.nmr is not None and self.nmr.shielding_tensors:
            software_label = self.qm_software or "UNKNOWN"
            property_suffix = f"BY_{software_label}".upper()
            numeric_property_names: set[str] = set()
            string_property_names: set[str] = set()
            tensor_labels = (
                ("XX", "XY", "XZ"),
                ("YX", "YY", "YZ"),
                ("ZX", "ZY", "ZZ"),
            )

            if self.nmr.gauge:
                rwmol.SetProp(f"NMR_GAUGE_{property_suffix}", self.nmr.gauge)

            for shielding in self.nmr.shielding_tensors:
                atom = rwmol.GetAtomWithIdx(shielding.atom_index)
                numeric_properties: dict[str, float] = {}
                if shielding.isotropic is not None:
                    numeric_properties["NMR_SHIELDING_ISOTROPIC_PPM"] = float(
                        shielding.isotropic.m_as("ppm")
                    )
                if shielding.anisotropy is not None:
                    numeric_properties["NMR_SHIELDING_ANISOTROPY_PPM"] = float(
                        shielding.anisotropy.m_as("ppm")
                    )
                if shielding.principal_values is not None:
                    for principal_index, value in enumerate(
                        shielding.principal_values.m_as("ppm"), start=1
                    ):
                        numeric_properties[
                            f"NMR_SHIELDING_PRINCIPAL_VALUE_{principal_index}_PPM"
                        ] = float(value)
                tensor = np.asarray(shielding.shielding_tensor.m_as("ppm"), dtype=float)
                for row_index, row_labels in enumerate(tensor_labels):
                    for column_index, component_label in enumerate(row_labels):
                        numeric_properties[f"NMR_SHIELDING_TENSOR_{component_label}_PPM"] = float(
                            tensor[row_index, column_index]
                        )

                for property_name, value in numeric_properties.items():
                    full_name = f"{property_name}_{property_suffix}"
                    atom.SetDoubleProp(full_name, value)
                    numeric_property_names.add(full_name)

                string_properties = {
                    "NMR_SHIELDING_ANISOTROPY_CONVENTION": shielding.anisotropy_convention,
                    "NMR_SHIELDING_ORIENTATION": shielding.orientation,
                }
                for property_name, value in string_properties.items():
                    if value is None:
                        continue
                    full_name = f"{property_name}_{property_suffix}"
                    atom.SetProp(full_name, value)
                    string_property_names.add(full_name)

            for property_name in sorted(numeric_property_names):
                Chem.CreateAtomDoublePropertyList(rwmol, property_name)
            for property_name in sorted(string_property_names):
                Chem.CreateAtomStringPropertyList(rwmol, property_name)
        return rwmol.GetMol()

    def to_population_embedded_SDF_block(
        self,
        embed_populations: bool = True,
        embed_bond_orders: bool = True,
        embed_nmr: bool = True,
    ) -> str:
        """
        Write the SDF block with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        If `embed_populations` is True, this function will use all population properties in the `charge_spin_populations` field
        to generate the population embedded rdkit molecule object.
        If `embed_bond_orders` is True, this function will use all bond order properties in the `bond_orders` field
        to generate the bond order embedded rdkit molecule object.
        If `embed_nmr` is True, per-atom NMR shielding properties are included.

        Parameters:
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.
            embed_nmr (bool): If True, embed per-atom NMR shielding properties. Defaults to True.

        Returns:
            str: The SDF block with population embedded properties.
        """
        sio = StringIO()
        with Chem.SDWriter(sio) as w:
            w.write(self.qm_embedded_rdmol(embed_populations, embed_bond_orders, embed_nmr))
        return sio.getvalue()

    def to_population_embedded_SDF_file(
        self,
        filepath: os.PathLike | str,
        embed_populations: bool = True,
        embed_bond_orders: bool = True,
        embed_nmr: bool = True,
    ):
        """
        Write the SDF block to a file with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        Parameters:
            filepath (os.PathLike| str): The path to the output file.
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.
            embed_nmr (bool): If True, embed per-atom NMR shielding properties. Defaults to True.
        """
        with open(filepath, "w") as f:
            f.write(
                self.to_population_embedded_SDF_block(
                    embed_populations, embed_bond_orders, embed_nmr
                )
            )

    def vibrate(
        self,
        vibration_id: int | None = None,
        vibration: Vibration | None = None,
        *,
        ratio: float = 1.75,
        steps: int = 7,
    ) -> list[Molecule]:
        """
        Generate a list of base block parsers for vibration calculations.

        Parameters:
            vibration_id (Union[int, None]):
                The index of the vibration to be calculated. If not specified, the first vibration will be used.
            vibration (Union[Vibration, None]):
                The Vibration object to be calculated. If not specified, the first vibration will be used.
            ratio (float):
                The ratio to force the geometry to vibrate.
            steps (int):
                The number of steps to generate.

        Returns:
            List[BaseMolFrameParser]: A list of base block parsers for vibration calculations.
        """

        if vibration is None:
            if vibration_id is None:
                vibration_id = 0
            if self.vibrations is None:
                raise ValueError("No vibrations found in this frame")
            if vibration_id < 0 or vibration_id >= len(self.vibrations):
                raise IndexError(f"Invalid vibration id {vibration_id}")
            vibration = self.vibrations[vibration_id]
        assert vibration.vibration_mode.m.shape == self.coords.m.shape, "Invalid vibration mode"

        temp_moleculues = []  # Initialize a list of base block parsers

        # Iterate over a list of ratios
        for r in np.linspace(-ratio, ratio, num=steps, endpoint=True):
            # Calculate extreme coordinates based on current ratio
            extreme_coords = cast(np.ndarray, self.coords.m - vibration.vibration_mode.m * r)

            # Convert extreme coordinates to rdkit molecule object
            rdmol = Chem.MolFromXYZBlock(
                f"{len(self.atoms)}\n"
                + f"charge {self.charge} multiplicity {self.multiplicity}\n"
                + "\n".join(
                    [
                        f"{Chem.Atom(atom).GetSymbol():10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                        for atom, x, y, z in zip(
                            self.atoms,
                            *zip(*extreme_coords, strict=True),
                            strict=True,
                        )
                    ]
                )
            )
            # Rebuild using the rdkit molecule object
            if rdmol is None:
                continue
            if not check_crowding(rdmol):
                continue
            molecule = Molecule.from_coords(
                atom_symbols=self.atom_symbols,
                coords=extreme_coords,
                charge=self.charge,
                multiplicity=self.multiplicity,
            )
            # Check if the molecule satisfies crowding conditions and append it to the list
            temp_moleculues.append(molecule)
        return temp_moleculues

    def draw_vibration_animation(
        self,
        vibration_id: int | None = None,
        vibration: Vibration | None = None,
        *,
        ratio: float = 1.75,
        steps: int = 7,
        image_format: AnimationFormat = "gif",
        file_path: os.PathLike[str] | str | None = None,
        duration: int | Sequence[int] = 200,
        loop: int = 0,
        legends: Sequence[str | None] | None = None,
        **kwargs: Any,
    ) -> Any:
        """Render structures displaced along one normal mode as an animation."""

        candidates = self.vibrate(
            vibration_id=vibration_id,
            vibration=vibration,
            ratio=ratio,
            steps=steps,
        )
        # Displaced vibration geometries are coordinate-only molecules.  An
        # optional native batch avoids repeated lazy reconstruction when the
        # caller has enabled parent-process topology prewarming.
        if molopconfig.prewarm_topologies:
            reconstruct_topologies_batch(candidates, retain_results=False)
        candidate_legends = legends
        if candidate_legends is None:
            candidate_legends = self._vibration_animation_legends(
                candidates,
                vibration_id=vibration_id,
                vibration=vibration,
            )
        return render_molecule_animation(
            candidates,
            image_format=image_format,
            file_path=file_path,
            duration=duration,
            loop=loop,
            legends=candidate_legends,
            **kwargs,
        )

    def _vibration_animation_legends(
        self,
        candidates: Sequence[Molecule],
        *,
        vibration_id: int | None,
        vibration: Vibration | None,
    ) -> list[str]:
        selected_vibration = vibration
        selected_id = vibration_id
        if selected_vibration is None:
            selected_id = 0 if selected_id is None else selected_id
            assert self.vibrations is not None
            selected_vibration = self.vibrations[selected_id]

        mode_label = "Mode" if selected_id is None else f"Mode {selected_id}"
        frequency_label = "frequency unavailable"
        if selected_vibration.frequency is not None:
            try:
                frequency = float(selected_vibration.frequency.m_as("cm^-1"))
                frequency_label = f"frequency = {frequency:.2f} cm^-1"
            except (AttributeError, TypeError, ValueError):
                pass
        return [
            f"{mode_label} | {frequency_label} | geometry {index}/{len(candidates)}"
            for index in range(1, len(candidates) + 1)
        ]

    def get_QRC(self, ratio: float = 1.75, vibration_id: int = 0) -> list[Molecule]:
        assert self.vibrations is not None and self.vibrations[vibration_id].is_imaginary, (
            "Must be an imaginary vibration"
        )
        res = []
        for r in (-ratio, ratio):
            temp_moleculues = self.vibrate(vibration_id=vibration_id, ratio=r, steps=1)
            res.extend(temp_moleculues)
        return res

    def ts_vibration(self, ratio: float = 1.75, steps: int = 7) -> list[Molecule]:
        """
        Generate a list of base block parsers for transition state vibration calculations.

        Parameters:
            ratio (float):
                The ratio to force the geometry to vibrate.
            steps (int):
                The number of steps to generate.

        Returns:
            List[BaseMolFrameParser]: A list of base block parsers for transition state vibration calculations.
        """
        return self.vibrate(vibration_id=self._ts_vibration_id(), ratio=ratio, steps=steps)

    def draw_ts_vibration_animation(
        self,
        *,
        ratio: float = 1.75,
        steps: int = 7,
        image_format: AnimationFormat = "gif",
        file_path: os.PathLike[str] | str | None = None,
        duration: int | Sequence[int] = 200,
        loop: int = 0,
        legends: Sequence[str | None] | None = None,
        **kwargs: Any,
    ) -> Any:
        """Render the unique imaginary mode of a transition-state frame."""

        return self.draw_vibration_animation(
            vibration_id=self._ts_vibration_id(),
            ratio=ratio,
            steps=steps,
            image_format=image_format,
            file_path=file_path,
            duration=duration,
            loop=loop,
            legends=legends,
            **kwargs,
        )

    def _ts_vibration_id(self) -> int:
        if not self.is_TS:
            raise ValueError("Must be a TS frame")
        assert self.vibrations is not None
        imaginary_idxs = self.vibrations.imaginary_idxs
        if len(imaginary_idxs) != 1:
            raise ValueError("A TS frame must have exactly one imaginary vibration")
        return imaginary_idxs[0]

    def possible_pre_post_ts(
        self,
        show_3D: bool = False,
        *,
        min_ratio: float = 0.2,
        max_ratio: float = 1.8,
        steps: int = 9,
    ) -> tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]:
        """Infer possible pre- and post-TS molecules from stable side topologies.

        Parameters:
            show_3D (bool):
                Whether to show 3D coordinates. Defaults to False.
            min_ratio (float):
                The smallest displacement amplitude to sample. Defaults to 0.2.
            max_ratio (float):
                The largest displacement amplitude to sample. Defaults to 1.8.
            steps (int):
                The number of amplitudes sampled on each side. Defaults to 9.

        Returns:
            Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]: A tuple containing the possible pre- and post-transition state molecules.
        """
        if not np.isfinite(min_ratio) or min_ratio <= 0:
            raise ValueError("min_ratio must be a finite value greater than 0")
        if not np.isfinite(max_ratio) or max_ratio < min_ratio:
            raise ValueError("max_ratio must be finite and greater than or equal to min_ratio")
        if steps < 1:
            raise ValueError("steps must be >= 1")

        amplitudes = np.linspace(min_ratio, max_ratio, num=steps, endpoint=True)
        selected_sides: list[RdMol] = []
        for side_name, direction in (("negative", -1.0), ("positive", 1.0)):
            side_candidates: list[RdMol] = []
            for amplitude in amplitudes:
                # With one vibration step, ``ratio`` selects the signed extreme.
                displaced = self.ts_vibration(ratio=direction * float(amplitude), steps=1)
                for molecule in displaced:
                    if (rdmol := molecule.rdmol) is not None:
                        side_candidates.append(rdmol)
                        break
            selected_sides.append(_most_frequent_topology(side_candidates, side=side_name))

        # The side with more disconnected fragments is treated as the precursor.
        # Equal fragment counts retain the deterministic negative/positive ordering.
        selected_sides.sort(key=lambda mol: len(Chem.GetMolFrags(mol)), reverse=True)
        reactant_rdmol, product_rdmol = selected_sides
        if not show_3D:
            reactant_rdmol.RemoveAllConformers()
            product_rdmol.RemoveAllConformers()
        return reactant_rdmol, product_rdmol

    def _additional_pre_post_ts_side(
        self,
        *,
        side: str,
        direction: float,
        endpoint: RdMol,
        fixed_atom_indices: set[int],
        amplitudes: npt.NDArray[np.floating[Any]],
    ) -> RdMol | None:
        """Resample one mode direction while holding reaction-center atoms fixed."""

        if endpoint.GetNumConformers() == 0:
            raise ValueError(
                "Additional TS endpoint sampling requires endpoint conformers; "
                "call possible_pre_post_ts(show_3D=True) first."
            )

        endpoint_coords = np.asarray(endpoint.GetConformer().GetPositions(), dtype=float)
        ts_coords = np.asarray(self.coords.m, dtype=float)
        vibrations = self.vibrations
        if vibrations is None:
            raise ValueError("Additional TS endpoint sampling requires vibration data")
        vibration = vibrations[self._ts_vibration_id()]
        mode_coords = np.asarray(vibration.vibration_mode.m, dtype=float)
        if endpoint_coords.shape != ts_coords.shape or mode_coords.shape != ts_coords.shape:
            raise ValueError("TS endpoint and vibration coordinates must have matching shapes")

        fixed_indices = sorted(fixed_atom_indices)
        side_candidates: list[RdMol] = []
        for amplitude in amplitudes:
            displaced_coords = np.array(
                ts_coords + direction * mode_coords * float(amplitude),
                dtype=float,
                copy=True,
            )
            displaced_coords[fixed_indices, :] = endpoint_coords[fixed_indices, :]
            coordinate_rdmol = Chem.MolFromXYZBlock(
                f"{len(self.atoms)}\n"
                + f"charge {self.charge} multiplicity {self.multiplicity}\n"
                + "\n".join(
                    [
                        f"{Chem.Atom(atom).GetSymbol():10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                        for atom, (x, y, z) in zip(
                            self.atoms,
                            displaced_coords,
                            strict=True,
                        )
                    ]
                )
            )
            if coordinate_rdmol is None or not check_crowding(coordinate_rdmol):
                continue

            molecule = Molecule.from_coords(
                atom_symbols=self.atom_symbols,
                coords=displaced_coords,
                charge=self.charge,
                multiplicity=self.multiplicity,
            )
            if (rdmol := molecule.rdmol) is not None:
                side_candidates.append(rdmol)

        if not side_candidates:
            return None
        return _most_frequent_topology(side_candidates, side=f"{side}-additional")

    def additional_pre_post_ts(
        self,
        pre_rdmol: RdMol,
        post_rdmol: RdMol,
        *,
        min_ratio: float = 0.2,
        max_ratio: float = 1.8,
        steps: int = 9,
    ) -> tuple[RdMol, RdMol]:
        """Generate an additional-sampling TS endpoint representation.

        Atoms incident to bonds that differ between ``pre_rdmol`` and ``post_rdmol``
        retain their endpoint coordinates. The remaining atoms are resampled along
        the imaginary mode, and the resulting side-wise topology modes are returned.
        The standard :meth:`possible_pre_post_ts` result is never replaced or mutated.
        If no bond-change atoms are present, the input endpoints are returned directly.

        The input endpoints must retain 3D conformers, for example by calling
        ``possible_pre_post_ts(show_3D=True)``.
        """

        if not np.isfinite(min_ratio) or min_ratio <= 0:
            raise ValueError("min_ratio must be a finite value greater than 0")
        if not np.isfinite(max_ratio) or max_ratio < min_ratio:
            raise ValueError("max_ratio must be finite and greater than or equal to min_ratio")
        if steps < 1:
            raise ValueError("steps must be >= 1")
        if pre_rdmol.GetNumAtoms() != post_rdmol.GetNumAtoms() or pre_rdmol.GetNumAtoms() != len(
            self.atoms
        ):
            raise ValueError("TS endpoints must contain the same atoms as the source frame")

        changed_atom_indices = _bond_change_atom_indices(pre_rdmol, post_rdmol)
        if not changed_atom_indices:
            return pre_rdmol, post_rdmol

        ts_coords = np.asarray(self.coords.m, dtype=float)
        vibrations = self.vibrations
        if vibrations is None:
            raise ValueError("Additional TS endpoint sampling requires vibration data")
        vibration = vibrations[self._ts_vibration_id()]
        mode_coords = np.asarray(vibration.vibration_mode.m, dtype=float)
        if mode_coords.shape != ts_coords.shape:
            raise ValueError("TS vibration mode and source coordinates must have matching shapes")

        endpoint_directions: list[float] = []
        for endpoint in (pre_rdmol, post_rdmol):
            if endpoint.GetNumConformers() == 0:
                raise ValueError(
                    "Additional TS endpoint sampling requires endpoint conformers; "
                    "call possible_pre_post_ts(show_3D=True) first."
                )
            endpoint_coords = np.asarray(endpoint.GetConformer().GetPositions(), dtype=float)
            if endpoint_coords.shape != ts_coords.shape:
                raise ValueError("TS endpoint and source coordinates must have matching shapes")
            projection = float(np.sum((endpoint_coords - ts_coords) * mode_coords))
            if np.isclose(projection, 0.0):
                raise ValueError("Could not determine the vibration direction of a TS endpoint")
            endpoint_directions.append(1.0 if projection > 0.0 else -1.0)

        amplitudes = np.linspace(min_ratio, max_ratio, num=steps, endpoint=True)
        additional_endpoints: list[RdMol] = []
        for side, direction, endpoint in zip(
            ("pre", "post"),
            endpoint_directions,
            (pre_rdmol, post_rdmol),
            strict=True,
        ):
            additional = self._additional_pre_post_ts_side(
                side=side,
                direction=direction,
                endpoint=endpoint,
                fixed_atom_indices=changed_atom_indices,
                amplitudes=amplitudes,
            )
            additional_endpoints.append(additional if additional is not None else endpoint)

        additional_endpoints.sort(
            key=lambda molecule: len(Chem.GetMolFrags(molecule)),
            reverse=True,
        )
        return additional_endpoints[0], additional_endpoints[1]

    def save_pre_post_ts(
        self,
        output_dir: os.PathLike[str] | str,
        *,
        prefix: str | None = None,
        format: Literal["xyz", "sdf"] = "xyz",
        min_ratio: float = 0.2,
        max_ratio: float = 1.8,
        steps: int = 9,
    ) -> tuple[Path, Path]:
        """Save inferred pre- and post-TS endpoint candidates.

        The endpoints are geometric candidates from :meth:`possible_pre_post_ts`,
        not separately optimized reactant and product structures.

        Parameters:
            output_dir: Directory for the exported endpoint files.
            prefix: Filename prefix. Defaults to the source filename stem and
                frame ID when the frame has file metadata, otherwise the frame ID.
            format: Endpoint format, either ``"xyz"`` or ``"sdf"``. SDF
                preserves the reconstructed molecular graph and 3D conformer.
            min_ratio: Smallest displacement amplitude to sample.
            max_ratio: Largest displacement amplitude to sample.
            steps: Number of amplitudes sampled on each side.
        """

        normalized_format = format.lower()
        if normalized_format not in {"xyz", "sdf"}:
            raise ValueError(f"Unsupported endpoint format: {format!r}. Use 'xyz' or 'sdf'.")

        pre_rdmol, post_rdmol = self.possible_pre_post_ts(
            show_3D=True,
            min_ratio=min_ratio,
            max_ratio=max_ratio,
            steps=steps,
        )
        source_prefix = getattr(self, "pure_filename", None)
        if not source_prefix:
            filename = getattr(self, "filename", None)
            if filename:
                source_prefix = Path(str(filename)).stem
        if not source_prefix:
            file_path = getattr(self, "file_path", None)
            if file_path:
                source_prefix = Path(str(file_path)).stem
        name_prefix = prefix or (
            f"{source_prefix}_frame_{self.frame_id:03d}"
            if source_prefix
            else f"ts_frame_{self.frame_id:03d}"
        )
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)
        pre_path = destination / f"{name_prefix}_pre.{normalized_format}"
        post_path = destination / f"{name_prefix}_post.{normalized_format}"

        if normalized_format == "xyz":
            pre_path.write_text(Molecule.from_rdmol(pre_rdmol).to_XYZ(), encoding="utf-8")
            post_path.write_text(Molecule.from_rdmol(post_rdmol).to_XYZ(), encoding="utf-8")
        else:
            for endpoint_name, rdmol, path in (
                ("pre-TS endpoint candidate", pre_rdmol, pre_path),
                ("post-TS endpoint candidate", post_rdmol, post_path),
            ):
                endpoint_rdmol = Chem.Mol(rdmol)
                endpoint_rdmol.SetProp("_Name", endpoint_name)
                writer = Chem.SDWriter(str(path))
                try:
                    writer.write(endpoint_rdmol)
                finally:
                    writer.close()
        return pre_path, post_path

    def to_diff_rdmol(
        self,
        *,
        min_ratio: float = 0.2,
        max_ratio: float = 1.8,
        steps: int = 9,
    ) -> RdMol | None:
        """
        Generate a rdkit molecule object for the transition state with bond-breaking.

        Parameters:
            min_ratio (float):
                The smallest displacement amplitude to sample. Defaults to 0.2.
            max_ratio (float):
                The largest displacement amplitude to sample. Defaults to 1.8.
            steps (int):
                The number of amplitudes sampled on each side. Defaults to 9.

        Returns:
            Optional[RdMol]: The rdkit molecule object for the transition state with bond-breaking.
        """
        try:
            assert self.is_TS, "Must be a TS frame"

            reactant_rdmol, product_rdmol = self.possible_pre_post_ts(
                show_3D=True,
                min_ratio=min_ratio,
                max_ratio=max_ratio,
                steps=steps,
            )
            assert not (
                reactant_rdmol.HasSubstructMatch(product_rdmol)
                or product_rdmol.HasSubstructMatch(reactant_rdmol)
            ), (
                "The inferred reactant and product rdmol objects are consistent, thus it is not a bond-breaking transition state."
            )

            rwmol = Chem.RWMol(reactant_rdmol)
            for start_atom_idx, end_atom_idx in sorted(
                _bond_change_pairs(reactant_rdmol, product_rdmol)
            ):
                _process_bond_helper(
                    rwmol,
                    start_atom_idx,
                    end_atom_idx,
                    product_rdmol,
                )
            for atom_idx in range(rwmol.GetNumAtoms()):
                atom = rwmol.GetAtomWithIdx(atom_idx)
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)
            return rwmol.GetMol()

        except AssertionError as e:
            moloplogger.error(f"Assertion failed: {e}")
            return None
        except Exception as e:
            moloplogger.error(f"Unexpected error occurred: {e}")
            return None

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_error(self) -> bool | None:
        """
        Abstrcact method to check if the current frame is an error. The details are implemented in the derived classes.
        """
        if self.energies is None:
            return True
        if self.energies.total_energy is None:
            return True
        if self.status is None:
            return None
        if self.status.scf_converged is False or self.status.normal_terminated is False:
            return True
        if self.status.scf_converged is True or self.status.normal_terminated is True:
            return False
        return None

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_normal(self) -> bool | None:
        is_error = self.is_error
        return None if is_error is None else not is_error

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_TS(self) -> bool | None:
        """
        Abstract method to check if the current frame is a transition state. The details are implemented in the derived classes.
        """
        if self.is_error:
            return False
        if self.vibrations is None:
            return False
        if len(self.vibrations.frequencies) == 0:
            return False
        return len(self.vibrations.imaginary_idxs) == 1

    @computed_field()  # type: ignore[prop-decorator]
    @property
    def is_optimized(self) -> bool | None:
        """
        Check if the molecule is optimized.
        """
        if self.geometry_optimization_status is None:
            return False
        if not self.geometry_optimization_status.geometry_optimized:
            return False
        return self.vibrations is None or self.vibrations.num_imaginary <= 1

    def to_summary_dict(self, brief: bool = True, **kwargs) -> SummaryDict:
        brief_dict = super().to_summary_dict(brief=brief, **kwargs) | {
            summary_column("Calc Parameter", "Software"): self.qm_software,
            summary_column("Calc Parameter", "Version"): self.qm_software_version,
            summary_column("Calc Parameter", "Method"): self.method,
            summary_column("Calc Parameter", "BasisSet"): self.basis_set,
            summary_column("Calc Parameter", "Functional"): self.functional,
            summary_column("Calc Parameter", "Keywords"): self.keywords,
            summary_column("Environment", "SolventModel"): (
                self.solvent.solvent_model if self.solvent else None
            ),
            summary_column("Environment", "Solvent"): (
                self.solvent.solvent if self.solvent else None
            ),
            summary_column("Status", "IsError"): self.is_error,
            summary_column("Status", "IsNormal"): self.is_normal,
            summary_column("Status", "IsTS"): self.is_TS,
            summary_column("Status", "IsOptimized"): self.is_optimized,
        }
        if self.is_TS:
            try:
                pre, post = self.possible_pre_post_ts()
                pre_smiles = _canonical_smiles_from_rdmol(pre)
                post_smiles = _canonical_smiles_from_rdmol(post)
            except NativeReconstructionConcurrencyError:
                raise
            except Exception as e:
                moloplogger.error(f"Error in possible_pre_post_ts: {e}")
                pre_smiles, post_smiles = "", ""
            brief_dict |= {
                summary_column("General", "PreCanonicalSMILES"): pre_smiles,
                summary_column("General", "PostCanonicalSMILES"): post_smiles,
            }
        if self.temperature:
            item = summary_item("Environment", "Temperature", self.temperature)
            if item is not None:
                column, value = item
                brief_dict[column] = value
        if self.pressure:
            item = summary_item("Environment", "Pressure", self.pressure)
            if item is not None:
                column, value = item
                brief_dict[column] = value

        if not brief:
            brief_dict |= self.energies.to_summary_dict() if self.energies else {}
            brief_dict |= (
                self.thermal_informations.to_summary_dict() if self.thermal_informations else {}
            )
            brief_dict |= (
                self.geometry_optimization_status.to_summary_dict()
                if self.geometry_optimization_status
                else {}
            )
            brief_dict |= self.vibrations.to_summary_dict() if self.vibrations else {}

        return brief_dict

    @model_validator(mode="after")
    def physical_check(self) -> Self:
        num_atoms = len(self.atoms)
        if self.vibrations and len(self.vibrations) not in (
            num_atoms * 3 - 6,
            num_atoms * 3 - 5,
            num_atoms * 3 - 3,
        ):
            raise ValueError(
                f"Invalid vibrationa mode count: {len(self.vibrations)} for {num_atoms} atoms"
            )
        return self


calc_frame = TypeVar("calc_frame", bound="BaseCalcFrame")
coords_frame = TypeVar("coords_frame", bound="BaseCoordsFrame")
