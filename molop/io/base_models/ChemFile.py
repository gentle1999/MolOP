"""
Author: TMJ
Date: 2025-07-28 18:44:12
LastEditors: TMJ
LastEditTime: 2025-11-27 16:42:41
Description: 请填写简介
"""

import os
from sys import getsizeof
from typing import (
    Any,
    Dict,
    Generic,
    List,
    Literal,
    Optional,
    Sequence,
    TypeVar,
    Union,
    overload,
)

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

    @property
    def file_size(self) -> int:
        return getsizeof(self.file_content)

    def _format_file_size(self) -> str:
        """
        Format file size with adaptive units.

        Returns:
            str: Formatted file size with appropriate unit (B, KB, MB, GB, etc.)
        """
        size = self.file_size
        for unit in ["B", "KB", "MB", "GB", "TB"]:
            if size < 1024.0:
                return f"{size:.2f} {unit}"
            size /= 1024.0
        return f"{size:.2f} PB"

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
            raise e

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
            self._frames_[frame.frame_id - 1]._next_frame = frame
            frame._prev_frame = self._frames_[frame.frame_id - 1]
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

    @overload
    def format_transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        frameID: int,
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str: ...
    @overload
    def format_transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        frameID: Sequence[int],
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> List[str]: ...
    def format_transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        frameID: Sequence[int] | int = -1,
        file_path: os.PathLike | str | None = None,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> str | List[str]:
        """
        Transform the selected frames to the specified format.

        Parameters:
            format (str): The format to be transformed to.
            frameID (int | slice | Sequence[int]): The frame(s) to be transformed.
            file_path (os.PathLike| str | None): If given, the file path to save the transformed block(s); otherwise,
                no file will be saved.
            **kwargs: Additional keyword arguments for the transformation.
        Returns:
            (str | List[str]): The transformed block(s).
        """
        assert format in (
            "xyz",
            "sdf",
            "cml",
            "gjf",
            "smi",
        ), "Only 'xyz', 'sdf', 'cml', 'gjf', 'smi' supported"
        assert file_path is None or not os.path.isdir(file_path), (
            "file_path should be a file path or None"
        )
        if isinstance(frameID, int):
            assert (
                file_path is None or os.path.splitext(file_path)[1] == f".{format}"
            ), "file_path should have the same extension as format"
            keywords = kwargs.copy()
            if format == "smi":
                block = self[frameID].to_canonical_SMILES(**kwargs)
            elif format == "gjf":
                block = self[frameID].to_GJF_block(
                    title_card=keywords.pop(
                        "title_card",
                        os.path.splitext(file_path)[0] if file_path else "title",
                    ),
                    **keywords,
                )
            else:
                block: str = getattr(self[frameID], f"to_{format.upper()}_block")(
                    **keywords
                )
            if file_path is not None:
                with open(file_path, "w") as f:
                    f.write(block)
            return block
        elif isinstance(frameID, Sequence):
            if embed_in_one_file:
                assert (
                    file_path is None or os.path.splitext(file_path)[1] == f".{format}"
                ), "file_path should have the same extension as format"
                sep = {
                    "xyz": "\n",
                    "sdf": "$$$$\n",
                    "cml": "\n",
                    "gjf": "\n",
                    "smi": "\n",
                }
                blocks = sep[format].join(
                    [
                        self.format_transform(format, idx, file_path=None, **kwargs)
                        for idx in frameID
                    ]
                )
                if file_path is not None:
                    with open(file_path, "w") as f:
                        f.write(blocks)
                return blocks
            return [
                self.format_transform(
                    format,
                    idx,
                    file_path=(
                        f"{os.path.splitext(file_path)[0]}-{idx:03d}.{format}"
                        if file_path is not None
                        else None
                    ),  # type: ignore
                    **kwargs,
                )
                for idx in frameID
            ]

    def _add_default_units(self) -> None: ...

    def to_summary_dict(self, **kwargs) -> Dict[tuple[str, str], Any]:
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
        return df[top_level_order]

    def release_file_content(self) -> None:
        self.file_content = ""


ChemFile = TypeVar("ChemFile", bound="BaseChemFile")


class BaseCoordsFile(BaseChemFile[coords_frame]): ...


class BaseCalcFile(BaseChemFile[calc_frame]):
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
        frames_with_opt = [
            frame for frame in frames if frame.geometry_optimization_status is not None
        ]
        return sorted(
            frames_with_opt,
            key=lambda frame: frame.geometry_optimization_status,  # type: ignore
        )

    @property
    def closest_optimized_frame(self) -> calc_frame:
        """
        Get the frame with the closest optimization status to the optimized state.
        Optimized to O(N) complexity using min().
        """
        frames_with_opt = [
            frame
            for frame in self.frames
            if frame.geometry_optimization_status is not None
        ]

        if not frames_with_opt:
            raise ValueError("No frames with optimization status found.")
        return min(
            frames_with_opt,
            key=lambda frame: frame.geometry_optimization_status,  # type: ignore
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
            import seaborn as sns  # type: ignore # lazy import  # noqa: I001
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

    def to_summary_dict(self, **kwargs) -> Dict[tuple[str, str], Any]:
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
                ("Status", f"RuningTime ({self.running_time.units})"): self.running_time.m
            }
        return brief_dict
