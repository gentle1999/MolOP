"""
Author: TMJ
Date: 2024-06-17 20:42:47
LastEditors: TMJ
LastEditTime: 2024-06-17 20:47:40
Description: 请填写简介
"""

import os
from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union, Any

import numpy as np
from openbabel import pybel
from pint.facets.numpy.quantity import NumpyQuantity
from pydantic import Field, PrivateAttr, computed_field
from rdkit import Chem
from rdkit.Chem.rdMolAlign import GetBestRMS
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from molop.descriptor.spms import SPMSCalculator
from molop.logger.logger import moloplogger
from molop.structure.FormatConverter import rdmol_to_omol
from molop.structure.GeometryTransformation import get_geometry_info, standard_orient
from molop.structure.GraphReconstruction import xyz_to_rdmol
from molop.structure.StructureTransformation import (
    attempt_replacement,
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
from molop.utils.types import RdMol

from .Bases import BaseDataClassWithUnit

pt = Chem.GetPeriodicTable()


class Molecule(BaseDataClassWithUnit):

    atoms: Sequence[int] = Field(
        default_factory=tuple, description="atom numbers", title="Atom numbers"
    )
    coords: NumpyQuantity = Field(
        default=np.array([[]]) * atom_ureg.angstrom,
        description="Atom coordinates, unit is `angstrom`",
        title="Atom coordinates",
    )
    charge: int = Field(default=0, description="Molecule total charge")
    multiplicity: int = Field(default=1, description="Molecule total multiplicity")

    def _add_default_units(self) -> None:
        self._default_units.update({"coords": atom_ureg.angstrom})

    @computed_field
    @property
    def smiles(self) -> str:
        return self.to_SMILES()

    @computed_field()
    @property
    def canonical_smiles(self) -> str:
        return self.to_canonical_SMILES()

    bonds: List[Tuple[int, int, int, int]] = Field(
        default=[],
        description="Bond information, each bond is represented by a tuple of three integers, "
        "where the first two integers represent the indices of the two atoms involved in the bond, "
        "the third integer represents the bond order (follow the RDKit convention)"
        "the fourth integer represents the stereo configuration (follow the RDKit convention)",
    )
    formal_charges: List[int] = Field(
        default=[],
        description="Formal charges of each atom",
    )
    formal_num_radicals: List[int] = Field(
        default=[],
        description="Number of radical electrons of each atom",
    )
    _rdmol: Optional[RdMol] = PrivateAttr(default=None)

    @property
    def atom_symbols(self) -> List[str]:
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
        return sum(Chem.Atom(atom).GetAtomicNum() for atom in self.atoms) + self.charge

    @property
    def elements(self) -> List[str]:
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
                f"{pt.GetElementSymbol(i)}{self.atoms.count(pt.GetElementSymbol(i))}"
                for i in range(1, 119)
                if self.atoms.count(pt.GetElementSymbol(i)) != 0
            ]
        )

    def __get_topology(self):
        if self._rdmol is None:
            return
        self.bonds = get_bond_pairs(self.rdmol)
        self.formal_charges = get_formal_charges(self.rdmol)
        self.formal_num_radicals = get_formal_num_radicals(self.rdmol)

    @property
    def omol(self) -> pybel.Molecule:
        """
        Get the openbabel molecule object.

        Returns:
            pybel.Molecule: The openbabel molecule object.
        """
        return rdmol_to_omol(self.rdmol)

    @property
    def rdmol(self) -> Union[Chem.rdchem.Mol, None]:
        """
        Get the rdkit molecule object.

        If reconstruction failed, return None.

        Returns:
            Union[Chem.rdchem.Mol,None]:
                The rdkit molecule object. If reconstruction failed, return None.
        """
        if self._rdmol is None:
            if len(self.bonds) == 0:
                try:
                    self._rdmol = xyz_to_rdmol(
                        self.to_XYZ_block(), self.charge, self.multiplicity - 1
                    )
                except Exception as e:
                    moloplogger.error(f"{e}")
                finally:
                    self.__get_topology()
            else:
                assert (
                    self.formal_charges and self.formal_num_radicals
                ), "If bonds given, formal charges and spins must be provided."
                self._rdmol = build_mol_from_atoms_and_bonds(
                    self.atoms,
                    self.bonds,
                    self.formal_charges,
                    self.formal_num_radicals,
                    coords=self.coords.m,
                )
        return self._rdmol

    @property
    def rdmol_no_conformer(self) -> Chem.rdchem.Mol:
        """
        Get the rdkit molecule object without conformer.

        Returns:
            Chem.rdchem.Mol: The rdkit molecule object without conformer.
        """
        rdmol = Chem.RWMol(self.rdmol)
        rdmol.RemoveAllConformers()
        return rdmol.GetMol()

    def to_XYZ_block(self) -> str:
        """
        Get the XYZ block.

        Returns:
            str: The XYZ block.
        """
        if self.coords is None:
            raise ValueError("No coordinates found!")
        return (
            f"{len(self.atoms)}\n"
            + f"charge {self.charge} multiplicity {self.multiplicity}\n"
            + "\n".join(
                [
                    f"{atom:10s}{x:10.5f}{y:10.5f}{z:10.5f}"
                    for atom, (x, y, z) in zip(
                        self.atom_symbols, self.coords.m, strict=True
                    )
                ]
            )
        )

    def to_SDF_block(self, engine: Literal["rdkit", "openbabel"] = "rdkit") -> str:
        """
        Get the SDF block.

        Parameters:
            engine (Literal["rdkit", "openbabel"]):
                The engine to generate the SDF block. default is "rdkit".

        Returns:
            str: The SDF block.
        """
        if self.rdmol is None:
            moloplogger.warning("SDF building failed. No RDKit molecule found.")
            return ""
        if engine == "rdkit":
            return Chem.MolToMolBlock(self.rdmol)
        elif engine == "openbabel":
            return self.omol.write("sdf")  # type: ignore
        else:
            raise ValueError(f"Unsupported engine: {engine}")

    def to_SDF_file(self, filepath: os.PathLike, **kwargs):
        """
        Write the SDF block to a file.

        Parameters:
            filepath (os.PathLike): The path to the file.
            **kwargs: The keyword arguments for the `to_SDF_block` method.
        """
        with open(filepath, "w") as f:
            f.write(self.to_SDF_block(**kwargs))

    def to_CML_block(self, engine: Literal["rdkit", "openbabel"] = "rdkit") -> str:
        """
        Get the CML block.

        Returns:
            str: The CML block.
        """
        if self.rdmol is None:
            moloplogger.warning("CML building failed. No RDKit molecule found.")
            return ""
        if engine == "rdkit":
            return Chem.MolToMrvBlock(self.rdmol)
        elif engine == "openbabel":
            return self.omol.write("cml")  # type: ignore
        else:
            raise ValueError(f"Unsupported engine: {engine}")

    def to_SMILES(self) -> str:
        """
        Get the SMILES with explicit hydrogens.

        Returns:
            str: The SMILES.
        """
        if self.rdmol is None:
            return ""
        if smi := Chem.MolToSmiles(self.rdmol):
            return smi
        moloplogger.error("SMILES building failed.")
        return ""

    def to_canonical_SMILES(self) -> str:
        """
        Get the SMILES with standardization.

        Returns:
            str: The standard SMILES.
        """
        return canonical_smiles(self.to_SMILES())

    def to_InChI(self) -> str:
        """
        Get the InChI.

        Returns:
            str: The InChI.
        """
        return Chem.MolToInchi(self.rdmol)  # type: ignore

    def calc_spms_descriptor(
        self,
        anchor_list: Union[Sequence[int], None] = None,
        sphere_radius: Union[float, None] = None,
        atom_radius: Literal["vdw", "covalent"] = "vdw",
        latitudinal_resolution: int = 40,
        longitude_resolution: int = 40,
        precision: int = 8,
        *,
        custom_first_anchors: Union[Sequence[int], None] = None,
        custom_second_anchors: Union[Sequence[int], None] = None,
        custom_third_anchors: Union[Sequence[int], None] = None,
    ) -> SPMSCalculator:
        """
        Re-implementation of SPMS descriptor [A Molecular Stereostructure Descriptor based on
        Spherical Projection](https://www.thieme-connect.de/products/ejournals/abstract/10.1055/s-0040-1705977).
        GitHub repository: https://github.com/licheng-xu-echo/SPMS.git

        Parameters:
            anchor_list (Sequence[int]):
                List of anchor atom ids to keep the invariance of 3D geometry. Default is None, which means using Centroid to
                origin and the atom nearest to the centroid to Z axis, and the atom farthest from the centroid to ZY face. If
                user provides a list of anchor atom ids, use the center of those atoms as the anchor, and do the same rotation -
                the atom nearest to the centroid to Z axis, and the atom farthest from the centroid to ZY face.
            sphere_radius (Union[float, None]):
                Sphere radius. Default is None to use the largest possible radius for each input molecule.
                If you want to use a fixed radius, set this parameter to a float value.
            atom_radius (Literal["vdw", "covalent"]):
                Atom radius type. Default is 'vdw', which means using Van der Waals radii as atom radii. If you want to use
                covalent radii, set this parameter to 'covalent'.
            latitudinal_resolution (int):
                Number of splits on the latitudinal axis. Default is 40.
            longitude_resolution (int):
                Number of splits on the longitude axis. Default is 40.
            precision (int):
                Precision of the SPMS descriptor. Default is 8.
            custom_first_anchors (Union[Sequence[int], None]):
                List of atom ids to use as the first anchor. Default is None.
            custom_second_anchors (Union[Sequence[int], None]):
                List of atom ids to use as the second anchor. Default is None.
            custom_third_anchors (Union[Sequence[int], None]):
                List of atom ids to use as the third anchor. Default is None.

        Returns:
            SPMSCalculator: The SPMS calculator object.
        """
        return SPMSCalculator(
            rdmol=self.rdmol,
            anchor_list=anchor_list,
            sphere_radius=sphere_radius,
            atom_radius=atom_radius,
            latitudinal_resolution=latitudinal_resolution,
            longitude_resolution=longitude_resolution,
            precision=precision,
            custom_first_anchors=custom_first_anchors,
            custom_second_anchors=custom_second_anchors,
            custom_third_anchors=custom_third_anchors,
        )

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
        map: Optional[Sequence[tuple[int, int]]] = None,
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
        assert isinstance(
            other, Molecule
        ), "The other object is not a BaseMolFrame object."

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
        return Molecule(
            atoms=[atom.GetAtomicNum() for atom in rdmol.GetAtoms()],
            coords=rdmol.GetConformer().GetPositions() * atom_ureg.angstrom,
            charge=get_total_charge(rdmol),
            multiplicity=get_total_multiplicity(rdmol),
            bonds=get_bond_pairs(rdmol),
            formal_charges=get_formal_charges(rdmol),
            formal_num_radicals=get_formal_num_radicals(rdmol),
        )

    def replace_substituent(
        self,
        query: Union[str, RdMol],
        replacement: Union[str, RdMol],
        bind_idx: Optional[int] = None,
        replace_all=False,
        attempt_num=10,
        crowding_threshold=0.75,
        angle_split=10,
        randomSeed=114514,
        start_idx: Optional[int] = None,
        end_idx: Optional[int] = None,
        *,
        replacement_relative_idx: int = 0,
        replacement_absolute_idx: Optional[int] = None,
        prefer_ZE: str = "Z",
    ) -> "Molecule":
        """
        Replace the substituent with the given SMARTS. The substituent is defined by the query_smi,
        and the new substituent is defined by the replacement_smi.

        Parameters:
            query (str | RdMol):
                The SMARTS or Mol object to query the substituent in the original molecule.
            replacement (str | RdMol):
                The SMARTS or Mol object of new substituent.
            bind_idx (int):
                The index of the atom to bind the new substituent. The default is None, which means
                to replace the first legal atom in original molecule.
                If specified, try to replace the legal substruct where the atom in it. User should
                meke sure the atom is legal.
                Detail example in (Repalce Substituent)[Repalce Substituent]
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
            replacement_relative_idx (int):
                The relative index of the radical atom in the replacement molecule to be
                transformed to the first atom.
            replacement_absolute_idx (Union[int, None]):
                Priority is higher than replacement_relative_idx.
                The absolute index of the radical atom in the replacement molecule to be
                transformed to the first atom.
                If None, the function will try to find the first atom in the replacement
                molecule that is a radical atom.
            prefer_ZE (str):
                The preferred stereochemistry of the bond to be replaced.
                If "Z", the function will try to replace the bond with Z stereochemistry.
                If "E", the function will try to replace the bond with E stereochemistry.
                only works for bond type DOUBLE.

        Returns:
            Molecule: The new molecule with the substituent replaced.
        """
        new_mol = attempt_replacement(
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
            replacement_relative_idx=replacement_relative_idx,
            replacement_absolute_idx=replacement_absolute_idx,
            prefer_ZE=prefer_ZE,
        )
        return self.from_rdmol(new_mol)

    def reset_atom_index(
        self,
        mapping_smarts: Union[str, None] = None,
        *,
        mapping_indice: Union[Sequence[int], None] = None,
    ) -> Optional["Molecule"]:
        """
        Reset the atom index of the molecule according to the mapping SMARTS.

        Parameters:
            mapping_smarts (str):
                The SMARTS to query the molecule substructure.
                The queried atoms will be renumbered and placed at the beginning of all atoms according to the order of the atoms in SMARTS. The relative order of the remaining atoms remains unchanged.
            mapping_indice (Sequence[int]):
                The indices of the atoms to be renumbered. The relative order of the remaining atoms remains unchanged.
                e.g. atoms = [0, 1, 2, 3, 4, 5]; mapping = [3, 5, 4] means the first atom in the new molecule is mapped to the third atom in the original molecule,
                the second atom is mapped to the first atom, and the third atom is mapped to the second atom. Result is [3, 5, 4, 0, 1, 2].
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
            mappings: List[Sequence[int]] = self.rdmol.GetSubstructMatches(smarts)
            if len(mappings) > 1:
                moloplogger.warning(
                    f"Multiple matches found in {self.to_canonical_SMILES()} with {mapping_smarts}, using the first one"
                )
            mapping = mappings[0]
            rdmol = reset_atom_index(self.rdmol, mapping)
            return self.from_rdmol(rdmol)
        elif mapping_indice is not None:
            assert max(mapping_indice) < len(
                self.rdmol.GetAtoms()
            ), "Invalid mapping index"
            rdmol = reset_atom_index(self.rdmol, mapping_indice)
            return self.from_rdmol(rdmol)
        return None

    def standard_orient(
        self,
        anchor_list: Sequence[int],
    ) -> "Molecule":
        """
        Depending on the input `idx_list`, `translate_anchor`, `rotate_anchor_to_axis`, and `rotate_anchor_to_plane` are executed in order to obtain the normalized oriented molecule.

        Sub-functions:
            - `translate_anchor`: Translate the entire molecule so that the specified atom reaches the origin.
            - `rotate_anchor_to_axis`: Rotate the specified second atom along the axis passing through the origin so that it reaches the positive half-axis of the X-axis.
            - `rotate_anchor_to_plane`: Rotate along the axis passing through the origin so that the specified third atom reaches quadrant 1 or 2 of the XY plane.

        Parameters:
            anchor_list (Sequence[int]):
                A Sequence of indices of the atoms to be translated to origin, rotated to X axis, and rotated again to XY face:

                - If length is 1, execute `translate_anchor`
                - If length is 2, execute `translate_anchor` and `rotate_anchor_to_axis`
                - If length is 3, execute `translate_anchor`, `rotate_anchor_to_axis` and `rotate_anchor_to_plane`
                - If the length of the input `anchor_list` is greater than 3, subsequent atomic numbers are not considered.
        Returns:
            BaseMolFrameParser: The new parser.
        """
        mol = Chem.RWMol(self.rdmol)
        standard_orient(mol, anchor_list)
        return self.from_rdmol(mol)

    def to_summary_dict(self) -> Dict[str, Any]:
        return {
            "SMILES": self.to_SMILES(),
            "CanonicalSMILES": self.to_canonical_SMILES(),
            "InChI": self.to_InChI(),
            "Charge": self.charge,
            "Multiplicity": self.multiplicity,
            "NumAtoms": len(self.atoms),
        }
        
