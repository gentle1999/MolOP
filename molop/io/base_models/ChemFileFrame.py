"""
Author: TMJ
Date: 2025-07-28 18:43:45
LastEditors: TMJ
LastEditTime: 2025-07-29 14:17:47
Description: 请填写简介
"""

import os
from copy import deepcopy
from io import StringIO
from typing import Generic, List, Optional, Tuple, TypeVar, Union

import numpy as np
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr, computed_field
from rdkit import Chem

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
    basis: str = Field(
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
        default=None, description="Forces of each atom, unit is `hartree/bohr`"
    )
    hessian: Optional[NumpyQuantity] = Field(
        default=None,
        description="Hessian matrix of the QM calculation, unit is `hartree/bohr^2`",
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
        default=None, description="Polarizability of the molecule"
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

    @property
    def population_embedded_rdmol(self) -> Optional[RdMol]:
        """
        Store the population embedded rdkit molecule object.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        This function will use all population properties in the `charge_spin_populations` field to generate the population embedded rdkit molecule object.
        """
        if self.rdmol is None:
            return None
        if self.charge_spin_populations is None:
            return self.rdmol
        rdmol_copy = deepcopy(self.rdmol)
        populations: dict[str, List[float]] = self.charge_spin_populations.model_dump(
            exclude_defaults=True
        )
        for population in populations:
            for atom_idx, pop in enumerate(populations[population]):
                rdmol_copy.GetAtomWithIdx(atom_idx).SetDoubleProp(
                    f"{population}_by_{self.qm_software}".upper(), pop
                )
            Chem.CreateAtomDoublePropertyList(
                rdmol_copy, f"{population}_by_{self.qm_software}".upper()
            )
        return rdmol_copy

    def to_population_embedded_SDF_block(self) -> str:
        """
        Write the SDF block with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html
        """
        sio = StringIO()
        with Chem.SDWriter(sio) as w:
            w.write(self.population_embedded_rdmol)
        return sio.getvalue()

    def to_population_embedded_SDF_file(self, filepath: os.PathLike):
        """
        Write the SDF block to a file with population embedded properties.

        Follow the guide in https://greglandrum.github.io/rdkit-blog/posts/2025-07-24-writing-partial-charges-to-sd-files.html

        Parameters:
            filepath (os.PathLike): The path to the output file.
        """
        with open(filepath, "w") as f:
            f.write(self.to_population_embedded_SDF_block())

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
            assert (
                len(self.vibrations) > vibration_id
            ), f"Invalid vibration id {vibration_id}"
            vibration = self.vibrations[vibration_id]
        assert (
            vibration.vibration_mode.m.shape == self.coords.m.shape
        ), "Invalid vibration mode"

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
        moloplogger.info(f"Generated {len(temp_moleculues)} TS vibration calculations")
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
        if not self.is_TS:
            raise RuntimeError("This is not a TS")

        return self.vibrate(
            vibration_id=0,
            ratio=ratio,
            steps=steps,
            ignore_dative=ignore_dative,
        )

    def possible_pre_post_ts(
        self,
        show_3D=False,
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
        assert (
            temp_moleculues[0].rdmol and temp_moleculues[-1].rdmol
        ), "Failed to generate TS vibrations"
        if not show_3D:
            temp_moleculues[0].rdmol.RemoveAllConformers()
            temp_moleculues[-1].rdmol.RemoveAllConformers()
        return temp_moleculues[0].rdmol, temp_moleculues[-1].rdmol

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
        if not self.is_optimized:
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


calc_frame = TypeVar("calc_frame", bound="BaseCalcFrame")
coords_frame = TypeVar("coords_frame", bound="BaseCoordsFrame")
