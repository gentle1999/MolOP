"""
Author: TMJ
Date: 2025-07-28 18:43:45
LastEditors: TMJ
LastEditTime: 2025-12-11 00:10:26
Description: 请填写简介
"""

import os
from io import StringIO
from typing import Any, Dict, Generic, List, Optional, Tuple, TypeVar, Union

import numpy as np
import numpy.typing as npt
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field
from rdkit import Chem
from scipy.sparse import coo_matrix

from molop.config import moloplogger
from molop.io.base_models.DataClasses import (
    BondOrders,
    ChargeSpinPopulations,
    Energies,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    MolecularOrbitals,
    Polarizability,
    SinglePointProperties,
    Status,
    ThermalInformations,
    TotalSpin,
    Vibration,
    Vibrations,
)
from molop.io.base_models.Molecule import Molecule
from molop.structure.GraphReconstruction import xyz_to_rdmol
from molop.structure.StructureTransformation import check_crowding
from molop.unit import atom_ureg
from molop.utils.types import RdMol

ChemFileFrame = TypeVar("ChemFileFrame", bound="BaseChemFileFrame")


class BaseChemFileFrame(Molecule, Generic[ChemFileFrame]):
    frame_id: int = Field(default=0, description="Frame ID")
    frame_content: str = Field(default="", repr=False, exclude=True)
    _frame_type: str = PrivateAttr(default="")
    _next_frame: Optional[ChemFileFrame] = PrivateAttr(default=None)
    _prev_frame: Optional[ChemFileFrame] = PrivateAttr(default=None)

    @property
    def next(self) -> Optional[ChemFileFrame]:
        return self._next_frame

    @property
    def prev(self) -> Optional[ChemFileFrame]:
        return self._prev_frame

    @property
    def is_error(self) -> Optional[bool]:
        """
        Abstrcact method to check if the current frame is an error.
        The details are implemented in the derived classes.
        """

    @property
    def is_normal(self) -> Optional[bool]: ...

    @property
    def is_TS(self) -> Optional[bool]:
        """
        Check if the molecule is a TS. Can not check if the molecule
        is a TS without frequency information. Thus this function returns False.
        """

    @property
    def is_optimized(self) -> Optional[bool]:
        """
        Check if the molecule is optimized.
        """

    def to_summary_dict(self, **kwargs) -> Dict[tuple[str, str], Any]:
        return {**super().to_summary_dict(), ("General", "FrameID"): self.frame_id}

    def log_with_file_info(self, content: str, level: str = "info"):
        if hasattr(self, "filename"):
            getattr(moloplogger, level)(f"{self.filename}: {content}")  # type: ignore
        else:
            getattr(moloplogger, level)(content)  # type: ignore


class BaseCoordsFrame(BaseChemFileFrame[ChemFileFrame]): ...


class BaseCalcFrame(BaseChemFileFrame[ChemFileFrame]):
    # QM software
    qm_software: str = Field(
        default="",
        description="QM software used to perform the calculation",
    )
    qm_software_version: str = Field(
        default="",
        description="QM software version used to perform the calculation",
    )
    # QM parameters
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )
    method: str = Field(
        default="",
        description="QM method used to perform the calculation. "
        "e.g. DFT or SEMI-EMPIRICAL or HF et. al.",
    )
    basis_set: str = Field(
        default="",
        description="Basis set used in the QM calculation, only for DFT calculations",
    )
    functional: str = Field(
        default="",
        description="Functional used in the QM calculation, only for DFT calculations",
    )
    # solvation
    solvent: Optional[ImplicitSolvation] = Field(
        default=None,
        description="Solvent used in the QM calculation",
    )
    # physical settings
    temperature: Optional[PlainQuantity] = Field(
        default=None,
        description="Temperature used in the QM calculation, unit is `K`",
    )
    pressure: Optional[PlainQuantity] = Field(
        default=None,
        description="Pressure used in the QM calculation, unit is `atm`",
    )
    # QM properties
    forces: Optional[NumpyQuantity] = Field(
        default=None,
        description="Forces of each atom, unit is `hartree/bohr`.\n"
        "In Gaussian, the extracted forces data are all calculated using the "
        "input coordinates as a reference.",
    )
    hessian: Optional[NumpyQuantity] = Field(
        default=None,
        description="Hessian matrix of the QM calculation, unit is `hartree/bohr^2`.\n"
        "In Gaussian, the extracted hessian data are all calculated using the "
        "input coordinates as a reference.",
    )
    rotation_constants: Optional[NumpyQuantity] = Field(
        default=np.array([[]]) * atom_ureg.gigahertz,
        description="Rotational constants, unit is `gigahertz`",
    )
    energies: Optional[Energies] = Field(
        default=None,
        description="Energies",
    )
    thermal_informations: Optional[ThermalInformations] = Field(
        default=None,
        description="Thermal Energies",
    )
    molecular_orbitals: Optional[MolecularOrbitals] = Field(
        default=None,
        description="Molecular Orbitals",
    )
    vibrations: Optional[Vibrations] = Field(default=None, description="vibrations")
    charge_spin_populations: Optional[ChargeSpinPopulations] = Field(
        default=None, description="Charge and spin populations"
    )
    polarizability: Optional[Polarizability] = Field(
        default=None,
        description="Polarizability of the molecule.\n"
        "In Gaussian, the extracted polarization-related data are all calculated using the "
        "input coordinates as a reference.",
    )
    bond_orders: Optional[BondOrders] = Field(
        default=None, description="Bond orders of the molecule"
    )
    total_spin: Optional[TotalSpin] = Field(
        default=None, description="Total spin of the molecule"
    )
    single_point_properties: Optional[SinglePointProperties] = Field(
        default=None,
        description="Single point properties of the molecule",
    )
    status: Optional[Status] = Field(default=None, description="Status of the frame")
    geometry_optimization_status: Optional[GeometryOptimizationStatus] = Field(
        default=None,
        description="Geometry optimization status",
    )
    running_time: Optional[PlainQuantity] = Field(
        default=None,
        description="Running time of the QM calculation, unit is `second`",
    )

    def qm_embedded_rdmol(
        self, embed_populations: bool = True, embed_bond_orders: bool = True
    ) -> Optional[RdMol]:
        """
        Store the population embedded rdkit molecule object.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        If `embed_populations` is True, this function will use all population properties in the `charge_spin_populations` field
        to generate the population embedded rdkit molecule object.
        If `embed_bond_orders` is True, this function will use all bond order properties in the `bond_orders` field
        to generate the bond order embedded rdkit molecule object.

        Parameters:
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.

        Returns:
            Optional[RdMol]: The population embedded rdkit molecule object.
        """
        if self.rdmol is None:
            return None
        rwmol = Chem.RWMol(self.rdmol)
        if embed_populations and self.charge_spin_populations is not None:
            populations: dict[str, List[float]] = (
                self.charge_spin_populations.model_dump(exclude_defaults=True)
            )
            for population, pop_list in populations.items():
                for atom_idx, pop in enumerate(pop_list):
                    rwmol.GetAtomWithIdx(atom_idx).SetDoubleProp(
                        f"{population}_by_{self.qm_software}".upper(), pop
                    )
                Chem.CreateAtomDoublePropertyList(
                    rwmol, f"{population}_by_{self.qm_software}".upper()
                )
        if embed_bond_orders and self.bond_orders is not None:
            bond_orders: dict[str, npt.NDArray[np.floating]] = (
                self.bond_orders.model_dump(exclude_defaults=True)
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
        return rwmol.GetMol()

    def to_population_embedded_SDF_block(
        self, embed_populations: bool = True, embed_bond_orders: bool = True
    ) -> str:
        """
        Write the SDF block with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        If `embed_populations` is True, this function will use all population properties in the `charge_spin_populations` field
        to generate the population embedded rdkit molecule object.
        If `embed_bond_orders` is True, this function will use all bond order properties in the `bond_orders` field
        to generate the bond order embedded rdkit molecule object.

        Parameters:
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.

        Returns:
            str: The SDF block with population embedded properties.
        """
        sio = StringIO()
        with Chem.SDWriter(sio) as w:
            w.write(self.qm_embedded_rdmol(embed_populations, embed_bond_orders))
        return sio.getvalue()

    def to_population_embedded_SDF_file(
        self,
        filepath: os.PathLike | str,
        embed_populations: bool = True,
        embed_bond_orders: bool = True,
    ):
        """
        Write the SDF block to a file with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        Parameters:
            filepath (os.PathLike| str): The path to the output file.
            embed_populations (bool): If True, embed the population properties. Defaults to True.
            embed_bond_orders (bool): If True, embed the bond order properties. Defaults to True.
        """
        with open(filepath, "w") as f:
            f.write(
                self.to_population_embedded_SDF_block(
                    embed_populations, embed_bond_orders
                )
            )

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update(
            {
                "coords": atom_ureg.angstrom,
                "forces": atom_ureg.Unit("hartree / bohr"),
                "rotation_constants": atom_ureg.Unit("gigahertz"),
                "running_time": atom_ureg.Unit("second"),
                "temperature": atom_ureg.Unit("K"),
                "electron_temperature": atom_ureg.Unit("K"),
            }
        )

    def vibrate(
        self,
        vibration_id: Union[int, None] = None,
        vibration: Union[Vibration, None] = None,
        *,
        ratio: float = 1.75,
        steps: int = 15,
        ignore_dative=True,
    ) -> List[Molecule]:
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
            ignore_dative (bool):
                Whether to ignore dative bonds.

        Returns:
            List[BaseMolFrameParser]: A list of base block parsers for vibration calculations.
        """

        if vibration is None:
            if vibration_id is None:
                vibration_id = 0
            if self.vibrations is None:
                raise ValueError("No vibrations found in this frame")
            assert len(self.vibrations) > vibration_id, (
                f"Invalid vibration id {vibration_id}"
            )
            vibration = self.vibrations[vibration_id]
        assert vibration.vibration_mode.m.shape == self.coords.m.shape, (
            "Invalid vibration mode"
        )

        temp_moleculues = []  # Initialize a list of base block parsers

        # Iterate over a list of ratios
        for r in np.linspace(-ratio, ratio, num=steps, endpoint=True):
            # Calculate extreme coordinates based on current ratio
            extreme_coords = self.coords.m - vibration.vibration_mode.m * r

            try:
                # Convert extreme coordinates to rdkit molecule object
                rdmol = xyz_to_rdmol(
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
                    ),
                    total_charge=self.charge,
                    total_radical_electrons=self.multiplicity - 1,
                    make_dative=not ignore_dative,
                )
                # Rebuild using the rdkit molecule object
                if rdmol is None:
                    continue
                if not check_crowding(rdmol):
                    continue
                Chem.SanitizeMol(rdmol)
                molecule = Molecule.from_rdmol(rdmol)
            except Exception as e:
                moloplogger.warning(
                    f"Failed to rebuild for {self.to_SMILES()} with error {e}"
                )
                continue
            # Check if the molecule satisfies crowding conditions and append it to the list
            temp_moleculues.append(molecule)
        return temp_moleculues

    def ts_vibration(
        self, ratio: float = 1.75, steps: int = 15, ignore_dative=True
    ) -> List[Molecule]:
        """
        Generate a list of base block parsers for transition state vibration calculations.

        Parameters:
            ratio (float):
                The ratio to force the geometry to vibrate.
            steps (int):
                The number of steps to generate.
            ignore_dative (bool):
                Whether to ignore dative bonds.

        Returns:
            List[BaseMolFrameParser]: A list of base block parsers for transition state vibration calculations.
        """
        assert self.is_TS, "Must be a TS frame"

        return self.vibrate(
            vibration_id=0,
            ratio=ratio,
            steps=steps,
            ignore_dative=ignore_dative,
        )

    def possible_pre_post_ts(
        self,
        show_3D: bool = False,
        *,
        ratio: float = 1.75,
        steps: int = 15,
        ignore_dative=True,
    ) -> Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]:
        """
        This method returns the possible pre- and post-transition state molecules.

        Parameters:
            show_3D (bool):
                Whether to show 3D coordinates. Defaults to False.
            ratio (float):
                The ratio to force the geometry to vibrate. Defaults to 1.75.
            steps (int):
                The number of steps to generate. Defaults to 15.
            ignore_dative (bool):
                Whether to ignore dative bonds. Defaults to True.

        Returns:
            Tuple[Chem.rdchem.Mol, Chem.rdchem.Mol]: A tuple containing the possible pre- and post-transition state molecules.
        """
        temp_moleculues = self.ts_vibration(
            ratio=ratio, steps=steps, ignore_dative=ignore_dative
        )
        assert len(temp_moleculues) > 0, "Failed to generate TS vibrations"
        assert temp_moleculues[0].rdmol and temp_moleculues[-1].rdmol, (
            "Failed to generate TS vibrations"
        )
        reactant_rdmol, product_rdmol = (
            temp_moleculues[0].rdmol,
            temp_moleculues[-1].rdmol,
        )
        if not show_3D:
            reactant_rdmol.RemoveAllConformers()
            product_rdmol.RemoveAllConformers()
        return reactant_rdmol, product_rdmol

    def to_diff_rdmol(self, *, ratio: float = 1.75, steps: int = 15) -> Optional[RdMol]:
        """
        Generate a rdkit molecule object for the transition state with bond-breaking.

        Parameters:
            ratio (float):
                The ratio to force the geometry to vibrate. Defaults to 1.75.
            steps (int):
                The number of steps to generate. Defaults to 15.

        Returns:
            Optional[RdMol]: The rdkit molecule object for the transition state with bond-breaking.
        """
        try:
            assert self.is_TS, "Must be a TS frame"

            reactant_rdmol, product_rdmol = self.possible_pre_post_ts(
                show_3D=True, ratio=ratio, steps=steps
            )
            assert not (
                reactant_rdmol.HasSubstructMatch(product_rdmol)
                or product_rdmol.HasSubstructMatch(reactant_rdmol)
            ), (
                "The inferred reactant and product rdmol objects are consistent, thus it is not a bond-breaking transition state."
            )

            rwmol = Chem.RWMol(reactant_rdmol)

            for bond_idx in range(reactant_rdmol.GetNumBonds()):
                bond = reactant_rdmol.GetBondWithIdx(bond_idx)
                start_atom_idx, end_atom_idx = (
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                )
                self._process_bond(
                    rwmol,
                    start_atom_idx,
                    end_atom_idx,
                    product_rdmol,
                    bond.GetBondType(),
                )

            for bond_idx in range(product_rdmol.GetNumBonds()):
                bond = product_rdmol.GetBondWithIdx(bond_idx)
                start_atom_idx, end_atom_idx = (
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                )
                self._process_bond(
                    rwmol,
                    start_atom_idx,
                    end_atom_idx,
                    reactant_rdmol,
                    bond.GetBondType(),
                )

            return rwmol.GetMol()

        except AssertionError as e:
            moloplogger.error(f"Assertion failed: {e}")
            return None
        except Exception as e:
            moloplogger.error(f"Unexpected error occurred: {e}")
            return None

    def _process_bond(
        self,
        rwmol: Chem.RWMol,
        start_atom_idx: int,
        end_atom_idx: int,
        other_rdmol: RdMol,
        bond_type: Chem.BondType,
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
            bond_type (Chem.BondType):
                The type of the bond in the reactant or product molecule.
        """
        if rwmol.GetBondBetweenAtoms(start_atom_idx, end_atom_idx) is None:
            rwmol.AddBond(start_atom_idx, end_atom_idx, Chem.BondType.ZERO)
        elif (
            other_rdmol.GetBondBetweenAtoms(start_atom_idx, end_atom_idx) is None
            or bond_type
            != other_rdmol.GetBondBetweenAtoms(
                start_atom_idx, end_atom_idx
            ).GetBondType()
        ):
            rwmol.GetBondBetweenAtoms(start_atom_idx, end_atom_idx).SetBondType(
                Chem.BondType.ZERO
            )

    @computed_field
    @property
    def is_error(self) -> Optional[bool]:
        """
        Abstrcact method to check if the current frame is an error. The details are implemented in the derived classes.
        """
        if self.status is None:
            return True
        if self.energies is None:
            return True
        if self.energies.total_energy is None:
            return True
        return not self.status.normal_terminated

    @computed_field()
    @property
    def is_normal(self) -> bool:
        return not self.is_error

    @computed_field
    @property
    def is_TS(self) -> Optional[bool]:
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

    @computed_field
    @property
    def is_optimized(self) -> Optional[bool]:
        """
        Check if the molecule is optimized.
        """
        if self.geometry_optimization_status is None:
            return False
        return self.geometry_optimization_status.geometry_optimized

    def to_summary_dict(
        self, brief: bool = True, **kwargs
    ) -> Dict[Tuple[str, str], Any]:
        try:
            brief_dict = super().to_summary_dict(**kwargs) | {
                ("Calc Parameter", "Software"): self.qm_software,
                ("Calc Parameter", "Version"): self.qm_software_version,
                ("Calc Parameter", "Method"): self.method,
                ("Calc Parameter", "BasisSet"): self.basis_set,
                ("Calc Parameter", "Functional"): self.functional,
                ("Calc Parameter", "Keywords"): self.keywords,
                ("Environment", "SolventModel"): self.solvent.solvent_model
                if self.solvent
                else None,
                ("Environment", "Solvent"): self.solvent.solvent
                if self.solvent
                else None,
                ("Status", "IsError"): self.is_error,
                ("Status", "IsNormal"): self.is_normal,
                ("Status", "IsTS"): self.is_TS,
                ("Status", "IsOptimized"): self.is_optimized,
            }
            if self.temperature:
                brief_dict = brief_dict | {
                    (
                        "Environment",
                        f"Temperature ({self.temperature.units})",
                    ): self.temperature
                }
            if self.pressure:
                brief_dict = brief_dict | {
                    ("Environment", f"Pressure ({self.pressure.units})"): self.pressure
                }

            if not brief:
                brief_dict |= self.energies.to_summary_dict() if self.energies else {}
                brief_dict |= (
                    self.thermal_informations.to_summary_dict()
                    if self.thermal_informations
                    else {}
                )
                brief_dict |= (
                    self.geometry_optimization_status.to_summary_dict()
                    if self.geometry_optimization_status
                    else {}
                )
                brief_dict |= (
                    self.vibrations.to_summary_dict() if self.vibrations else {}
                )

            return brief_dict

        except Exception as e:
            moloplogger.error(f"Error in to_summary_dict: {e}")
            return {}


calc_frame = TypeVar("calc_frame", bound="BaseCalcFrame")
coords_frame = TypeVar("coords_frame", bound="BaseCoordsFrame")
