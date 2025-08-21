"""
Author: TMJ
Date: 2025-07-28 18:44:12
LastEditors: TMJ
LastEditTime: 2025-07-28 22:39:18
Description: 请填写简介
"""

from typing import Generic, List, Optional, Sequence, TypeVar, Union, overload

import pandas as pd
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import ChemFileFrame, calc_frame, coords_frame
from molop.io.base_models.DataClasses import (
    GeometryOptimizationStatus,
    ImplicitSolvation,
    Status,
)
from molop.unit import atom_ureg


class BaseChemFile(BaseDataClassWithUnit, Generic[ChemFileFrame]):
    # inside frames parsed from the file
    _frames_: List[ChemFileFrame] = PrivateAttr(default_factory=list)
    _index_: int = PrivateAttr(default=0)

    file_content: str = Field(
        default="", description="File content.", repr=False, exclude=True
    )
    charge: int = Field(default=0, description="charge")
    multiplicity: int = Field(default=1, description="multiplicity")

    @overload
    def __getitem__(self, frameID: int) -> ChemFileFrame: ...
    @overload
    def __getitem__(self, frameID: slice) -> List[ChemFileFrame]: ...
    @overload
    def __getitem__(self, frameID: Sequence[int]) -> List[ChemFileFrame]: ...
    def __getitem__(
        self, frameID: Union[int, slice, Sequence[int]]
    ) -> Union[ChemFileFrame, List[ChemFileFrame]]:
        try:
            if isinstance(frameID, int):
                return self._frames_[frameID]
            if isinstance(frameID, slice):
                return self._frames_[frameID]
            if isinstance(frameID, Sequence):
                return [self[i] for i in frameID]
        except IndexError as e:
            raise IndexError(
                "Invalid index type. Only `int`, `slice`, and `Sequence[int]` are allowed."
            ) from e

    def __iter__(self):
        self._index_ = 0
        return self

    def __next__(
        self,
    ) -> ChemFileFrame:
        if self._index_ >= len(self):
            raise StopIteration
        else:
            self._index_ += 1
            return self._frames_[self._index_ - 1]

    def __len__(self) -> int:
        return len(self._frames_)

    def append(
        self,
        frame: ChemFileFrame,
    ):
        """
        Append a frame to the list of frames.

        Parameters:
            frame (chem_file_frame): The frame to be appended.
        """
        frame.frame_id = len(self._frames_)
        if frame.frame_id > 0:
            self._frames_[frame.frame_id - 1]._next_frame = frame  # type: ignore
            frame._prev_frame = self._frames_[frame.frame_id - 1]  # type: ignore
        self._frames_.append(frame)

    @property
    def frames(
        self,
    ) -> List[ChemFileFrame]:
        """
        Get the list of parsed frames.

        Returns:
            List[chem_file_frame]: The list of parsed frames.
        """
        return self._frames_

    def _add_default_units(self) -> None: ...


ChemFile = TypeVar("ChemFile", bound="BaseChemFile")


class CoordsFileMixin(BaseChemFile[coords_frame]): ...


class CalcFileMixin(BaseChemFile[calc_frame]):
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
        description="QM method used to perform the calculation. e.g. DFT or GFN2-xTB",
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
    running_time: Optional[PlainQuantity] = Field(
        default=None,
        description="Running time of the QM calculation, unit is `second`",
    )

    status: Status = Field(default=Status(), description="Status of the last frame")
    geometry_optimization_status: Optional[GeometryOptimizationStatus] = Field(
        default=None,
        description="Geometry optimization status of the last frame",
    )

    def _add_default_units(self) -> None:
        super()._add_default_units()
        self._default_units.update(
            {
                "temperature": atom_ureg.K,
                "electron_temperature": atom_ureg.K,
                "running_time": atom_ureg.second,
            }
        )

    @property
    def sort_by_optimization(self) -> Sequence[calc_frame]:
        """
        Sort the frames by the optimization status. The closer the frame is to the optimized state, the higher the priority.

        Returns:
            List[QMMolFrameType]: A list of frames sorted by the optimization status.
        """
        frames = self.frames
        for frame in frames:
            if frame.geometry_optimization_status is None:
                raise ValueError(f"Frame {frame.frame_id} has no optimization status.")
        return sorted(frames, key=lambda frame: frame.geometry_optimization_status)  # type: ignore

    @property
    def closest_optimized_frame(self) -> calc_frame:
        """
        Get the frame with the closest optimization status to the optimized state.

        Returns:
            QMMolFrameType: The frame with the closest optimization status to the optimized state.
        """
        return self.sort_by_optimization[0]

    @property
    def is_error(self) -> bool:
        """
        Abstrcact method to check if the current frame is an error. The details are implemented in the derived classes.
        """
        if self[-1].energies is None:
            return True
        if self[-1].energies.total_energy is None:
            return True
        return not self.status.normal_terminated

    def draw_energy_curve(self):
        try:
            import seaborn as sns  # type: ignore # lazy import  # noqa: I001
        except ImportError as e:
            raise ImportError(
                "Seaborn is required for drawing energy curve. Please install it first by `pip install seaborn`."
            ) from e
        energies = {
            frame.frame_id: frame.energies.total_energy.m
            for frame in self.frames
            if frame.energies is not None and frame.energies.total_energy is not None
        }
        temp_df = pd.DataFrame(
            {"frame_id": list(energies.keys()), "total_energy": list(energies.values())}
        )
        return sns.lineplot(x="frame_id", y="total_energy", data=temp_df)
