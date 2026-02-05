"""
Author: TMJ
Date: 2025-07-28 18:44:12
LastEditors: TMJ
LastEditTime: 2026-02-04 10:51:49
Description: 请填写简介
"""

import os
from collections.abc import Iterator, Sequence
from sys import getsizeof
from typing import Any, ClassVar, Generic, Literal, TypeVar, cast, overload

import pandas as pd
from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field, PrivateAttr

from molop.config import moloplogger
from molop.io import codec_registry
from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseCalcFrame, BaseChemFileFrame, BaseCoordsFrame
from molop.io.base_models.DataClasses import (
    GeometryOptimizationStatus,
    ImplicitSolvation,
    Status,
)
from molop.unit import atom_ureg


FrameT = TypeVar("FrameT", bound=BaseChemFileFrame)
CoordsFrameT = TypeVar("CoordsFrameT", bound=BaseCoordsFrame)
CalcFrameT = TypeVar("CalcFrameT", bound=BaseCalcFrame)


class BaseChemFile(BaseDataClassWithUnit, Sequence[FrameT], Generic[FrameT]):
    # inside frames parsed from the file
    _frames_: list[FrameT] = PrivateAttr(default_factory=list)
    _index_: int = PrivateAttr(default=0)

    file_content: str = Field(default="", description="File content.", repr=False, exclude=True)
    charge: int = Field(default=0, description="charge")
    multiplicity: int = Field(default=1, description="multiplicity")

    @property
    def file_size(self) -> int:
        return getsizeof(self.file_content)

    def _format_file_size(self) -> str:
        """
        Format file size with adaptive units.

        Returns:
            str: Formatted file size with appropriate unit (B, KB, MB, GB, etc.)
        """
        size = float(self.file_size)
        for unit in ["B", "KB", "MB", "GB", "TB"]:
            if size < 1024.0:
                return f"{size:.2f} {unit}"
            size /= 1024.0
        return f"{size:.2f} PB"

    def __repr__(self) -> str:
        return f"frames={len(self)}, {super().__repr__()}"

    @overload
    def __getitem__(self, frameID: int) -> FrameT: ...
    @overload
    def __getitem__(self, frameID: slice) -> list[FrameT]: ...
    @overload
    def __getitem__(self, frameID: Sequence[int]) -> list[FrameT]: ...
    def __getitem__(self, frameID: int | slice | Sequence[int]) -> FrameT | list[FrameT]:
        try:
            if isinstance(frameID, int):
                return self._frames_[frameID]
            if isinstance(frameID, slice):
                return self._frames_[frameID]
            if isinstance(frameID, Sequence):
                return [self._frames_[i] for i in frameID]
        except IndexError as e:
            raise e

    def __iter__(self) -> Iterator[FrameT]:  # type: ignore[override]
        return iter(self._frames_)

    def __next__(
        self,
    ) -> FrameT:
        if self._index_ >= len(self):
            raise StopIteration
        else:
            self._index_ += 1
            return self._frames_[self._index_ - 1]

    def __len__(self) -> int:
        return len(self._frames_)

    def append(
        self,
        frame: FrameT,
    ):
        """
        Append a frame to the list of frames.

        Parameters:
            frame (chem_file_frame): The frame to be appended.
        """
        frame.frame_id = len(self._frames_)
        if frame.frame_id > 0:
            self._frames_[frame.frame_id - 1]._next_frame = frame
            frame._prev_frame = self._frames_[frame.frame_id - 1]
        self._frames_.append(frame)

    @property
    def frames(
        self,
    ) -> list[FrameT]:
        """
        Get the list of parsed frames.

        Returns:
            List[chem_file_frame]: The list of parsed frames.
        """
        return self._frames_

    def format_transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        frameID: Sequence[int] | int | Literal["all"] | slice = -1,
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str | list[str]:
        """
        Transform the selected frames to the specified format.

        Parameters:
            format (str): The format to be transformed to.
            frameID (int | slice | Sequence[int] | Literal["all"]): The frame(s) to be transformed.
            file_path (os.PathLike| str | None): If given, the file path to save the transformed block(s); otherwise,
                no file will be saved.
            **kwargs: Additional keyword arguments for the transformation.
        Returns:
            (str | List[str]): The transformed block(s).
        """
        assert format in ("xyz", "sdf", "cml", "gjf", "smi"), (
            "Only 'xyz', 'sdf', 'gjf', 'smi' supported"
        )
        assert file_path is None or not os.path.isdir(file_path), (
            "file_path should be a file path or None"
        )
        if isinstance(frameID, int):
            frame_ids = [frameID if frameID >= 0 else len(self.frames) + frameID]
        elif frameID == "all":
            frame_ids = list(range(len(self.frames)))
        elif isinstance(frameID, slice):
            frame_ids = list(range(len(self.frames)))[frameID]
        elif isinstance(frameID, Sequence):
            frame_ids = []
            for i in frameID:
                assert isinstance(i, int), "frameID should be a sequence of integers"
                frame_ids.append(i)
        else:
            raise ValueError("frameID should be an integer, a sequence of integers, or 'all'")

        graph_policy = kwargs.pop("graph_policy", "prefer")
        rendered = cast(
            str | list[str],
            codec_registry.write(
                format,
                self,
                frameID=frame_ids,
                embed_in_one_file=embed_in_one_file,
                graph_policy=graph_policy,
                **kwargs,
            ),
        )
        if file_path:
            base = os.path.basename(file_path).split(".")[0]
            if isinstance(rendered, str):
                with open(base + f".{format}", "w") as f:
                    f.write(rendered)
            elif isinstance(rendered, list):
                for idx, frame_content in zip(frame_ids, rendered, strict=True):
                    with open(base + f"{idx:03d}.{format}", "w") as f:
                        f.write(frame_content)
        return rendered

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("General", "NumberOfFrames"): len(self),
            ("General", "FileSize"): self._format_file_size(),
        }

    def to_summary_df(self, **kwargs) -> pd.DataFrame:
        df = pd.concat(
            [frame.to_summary_series(**kwargs) for frame in self._frames_],
            axis=1,
        ).T
        top_level_order = df.columns.get_level_values(0).unique()
        return df.loc[:, top_level_order]

    def release_file_content(self) -> None:
        self.file_content = ""

    def log_with_file_info(self, content: str, level: str = "info"):
        filename = getattr(self, "filename", None)
        if filename:
            getattr(moloplogger, level)(f"{filename}: {content}")
        else:
            getattr(moloplogger, level)(content)


ChemFile = TypeVar("ChemFile", bound="BaseChemFile")


class BaseCoordsFile(BaseChemFile[CoordsFrameT], Generic[CoordsFrameT]): ...


class BaseCalcFile(BaseChemFile[CalcFrameT], Generic[CalcFrameT]):
    default_units: ClassVar[dict[str, UnitLike]] = {
        "temperature": atom_ureg.K,
        "electron_temperature": atom_ureg.K,
        "running_time": atom_ureg.second,
    }
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
    basis_set: str = Field(
        default="",
        description="Basis set used in the QM calculation, only for DFT calculations",
    )
    functional: str = Field(
        default="",
        description="Functional used in the QM calculation, only for DFT calculations",
    )
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
    running_time: PlainQuantity | None = Field(
        default=None,
        description="Running time of the QM calculation, unit is `second`",
    )

    status: Status = Field(default=Status(), description="Status of the last frame")
    geometry_optimization_status: GeometryOptimizationStatus | None = Field(
        default=None,
        description="Geometry optimization status of the last frame",
    )

    @property
    def sort_by_optimization(self) -> Sequence[BaseCalcFrame]:
        """
        Sort the frames by the optimization status. The closer the frame is to the optimized state, the higher the priority.

        Returns:
            List[QMMolFrameType]: A list of frames sorted by the optimization status.
        """
        frames = self.frames
        frames_with_opt = [
            frame for frame in frames if frame.geometry_optimization_status is not None
        ]
        return sorted(
            frames_with_opt,
            key=lambda frame: cast(
                GeometryOptimizationStatus,
                frame.geometry_optimization_status,
            ),
        )

    @property
    def closest_optimized_frame(self) -> BaseCalcFrame:
        """
        Get the frame with the closest optimization status to the optimized state.
        Optimized to O(N) complexity using min().
        """
        frames_with_opt = [
            frame for frame in self.frames if frame.geometry_optimization_status is not None
        ]

        if not frames_with_opt:
            raise ValueError("No frames with optimization status found.")
        return min(
            frames_with_opt,
            key=lambda frame: cast(
                GeometryOptimizationStatus,
                frame.geometry_optimization_status,
            ),
        )

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

    def draw_energy_curve(self, unit: str = "hartree"):
        """
        Draw the energy curve of the QM calculation.

        Parameters:
            unit (str, optional): The unit of the energy. Defaults to "hartree".

        Returns:
            matplotlib.axes.Axes: The axes of the energy curve plot.
        """
        try:
            import seaborn as sns
        except ImportError as e:
            raise ImportError(
                "Seaborn is required for drawing energy curve. Please install it first by `pip install seaborn`."
            ) from e
        energies = {
            frame.frame_id: frame.energies.total_energy.to(unit).m
            for frame in self.frames
            if frame.energies is not None and frame.energies.total_energy is not None
        }
        if len(energies) == 0:
            raise ValueError("No valid energy data found.")
        temp_df = pd.DataFrame(
            {
                "frame_id": list(energies.keys()),
                f"total_energy ({unit})": list(energies.values()),
            }
        )
        return sns.lineplot(x="frame_id", y=f"total_energy ({unit})", data=temp_df)

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        brief_dict = super().to_summary_dict(**kwargs) | {
            ("Calc Parameter", "Software"): self.qm_software,
            ("Calc Parameter", "Version"): self.qm_software_version,
            ("Calc Parameter", "Method"): self.method,
            ("Calc Parameter", "BasisSet"): self.basis_set,
            ("Calc Parameter", "Functional"): self.functional,
            ("Calc Parameter", "Keywords"): self.keywords,
            ("Environment", "SolventModel"): self.solvent.solvent_model if self.solvent else None,
            ("Status", "IsError"): self.is_error,
        }
        if self.temperature:
            brief_dict = brief_dict | {
                (
                    "Environment",
                    f"Temperature ({self.temperature.units})",
                ): self.temperature.m
            }
        if self.pressure:
            brief_dict = brief_dict | {
                ("Environment", f"Pressure ({self.pressure.units})"): self.pressure.m
            }
        if self.running_time:
            brief_dict = brief_dict | {
                (
                    "Status",
                    f"RuningTime ({self.running_time.units})",
                ): self.running_time.m
            }
        return brief_dict
