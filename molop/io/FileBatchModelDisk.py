import os
from collections import OrderedDict
from collections.abc import MutableMapping
from typing import Iterable, List, Literal, Sequence, Union, overload

from tqdm import tqdm

from molop.config import molopconfig, moloplogger
from molop.io.coords_models import GJFFileDisk, SDFFileDisk, XYZFileDisk
from molop.io.QM_models import G16LogFileDisk

FileDiskType = G16LogFileDisk | GJFFileDisk | XYZFileDisk | SDFFileDisk
QMFileDiskType = G16LogFileDisk


class FileBatchModelDisk(MutableMapping):
    """
    A class to store a batch of files and their corresponding models.
    """

    __dickfiles: OrderedDict[str, FileDiskType]

    def __init__(self, diskfiles: Iterable[FileDiskType] | None = None):
        self.__dickfiles: OrderedDict[str, FileDiskType] = OrderedDict()
        if diskfiles is not None:
            self.add_diskfiles(diskfiles)

    def add_diskfiles(self, diskfiles: Iterable[FileDiskType]) -> None:
        for diskfile in sorted(diskfiles):
            if not isinstance(diskfile, FileDiskType):
                raise TypeError(
                    f"file_parsers should be a list in which each element is a subclass of {FileDiskType}, got {type(diskfile)}"
                )
            if diskfile.file_path not in self.__dickfiles:
                self[diskfile.file_path] = diskfile
            else:
                moloplogger.warning(
                    f"File {diskfile.file_path} already exists in the batch, skipped"
                )

    @classmethod
    def new_batch(cls, parsers: Iterable[FileDiskType]):
        new_batch = cls()
        new_batch.add_diskfiles(parsers)
        return new_batch

    def __contains__(self, key: Union[str, FileDiskType]) -> bool:
        if isinstance(key, str):
            return key in self.__dickfiles.keys()
        else:
            return key in self.__dickfiles.values()

    @overload
    def __getitem__(self, key: int) -> FileDiskType: ...
    @overload
    def __getitem__(self, key: str) -> FileDiskType: ...
    @overload
    def __getitem__(self, key: slice) -> "FileBatchModelDisk": ...
    @overload
    def __getitem__(self, key: Sequence[int | str]) -> "FileBatchModelDisk": ...

    def __getitem__(
        self, key: int | str | slice | Sequence[int | str]
    ) -> Union[FileDiskType, "FileBatchModelDisk"]:
        if isinstance(key, int):
            return list(self.__dickfiles.values())[key]
        if isinstance(key, str):
            return self.__dickfiles[key]
        if isinstance(key, slice):
            return self.new_batch(
                [self.__dickfiles[k] for k in list(self.__dickfiles.keys())[key]]
            )
        if isinstance(key, Sequence):
            return self.new_batch([self[k] for k in key])
        raise TypeError(f"Invalid key type: {type(key)}")

    def __setitem__(self, key: str, value: FileDiskType) -> None:
        if not isinstance(key, str):
            raise TypeError(f"Key should be a string, got {type(key)}")
        if not isinstance(value, FileDiskType):
            raise TypeError(
                f"file_parsers should be a list in which each element is a subclass of {FileDiskType}, got {type(value)}"
            )
        if key in self.__dickfiles:
            moloplogger.warning(
                f"File {key} already exists in the batch, replaced with the new one"
            )
        self.__dickfiles[key] = value

    def __delitem__(self, key: str | int) -> None:
        if isinstance(key, int):
            key = list(self.__dickfiles.keys())[key]
        del self.__dickfiles[key]

    def __len__(self) -> int:
        return len(self.__dickfiles)

    def __iter__(self):
        self.__index = 0
        return self

    def __next__(self) -> FileDiskType:
        if self.__index >= len(self.__dickfiles):
            raise StopIteration
        self.__index += 1
        return self[self.__index - 1]

    @property
    def file_paths(self) -> List[str]:
        return [diskfile.file_path for diskfile in self]

    @property
    def file_names(self) -> List[str]:
        return [diskfile.filename for diskfile in self]

    def filter_state(
        self, state: Literal["ts", "error", "opt", "normal"], negate: bool = False
    ) -> "FileBatchModelDisk":
        """
        Filter the files based on their state.

        Parameters:
            state ("ts" | "error" | "opt" | "normal"): The state to filter.
            negate (bool): Whether to negate the filter.
        Returns:
            FileBatchModelDisk: A new batch of files with the filtered files.
        """

        def judge_func(diskfile: FileDiskType) -> bool:
            if not isinstance(diskfile, QMFileDiskType):
                return False
            if state == "ts":
                return diskfile[-1].is_TS != negate
            elif state == "error":
                return diskfile[-1].is_error != negate
            elif state == "opt":
                return (
                    diskfile[-1].is_normal != negate
                    and diskfile.closest_optimized_frame.is_optimized != negate
                )
            elif state == "normal":
                return diskfile[-1].is_normal != negate
            else:
                raise ValueError(f"Invalid state: {state}")

        return self.new_batch(filter(judge_func, self))

    def filter_value(
        self,
        target: Literal["charge", "multiplicity", "format"],
        value: Union[str, float, int],
        negate: bool = False,
    ):
        """
        Filter the files based on their properties.

        Parameters:
            target ("charge" | "multiplicity" | "format"): The property to filter.
            value (str | float | int): The value to filter.
            negate (bool): Whether to negate the filter.
        Returns:
            FileBatchModelDisk: A new batch of files with the filtered files.
        """

        def judge_func(diskfile: FileDiskType) -> bool:
            if target == "charge":
                return (diskfile[-1].charge == value) != negate
            elif target == "multiplicity":
                return (diskfile[-1].multiplicity == value) != negate
            elif target == "format":
                return (diskfile.file_format == value) != negate

        return self.new_batch(filter(judge_func, self))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"

    def format_transform(
        self,
        format: Literal["xyz", "sdf", "cml", "gjf", "smi"],
        output_dir: str = os.getcwd(),
        frameID: int | Literal["all"] = -1,
        embed_in_one_file: bool = True,
        **kwargs,
    ) -> dict[str, str | List[str]]:
        """
        Transform the file format of the batch of files.

        Parameters:
            format ("xyz" | "sdf" | "cml" | "gjf" | "smi"): The target format to transform.
            output_dir (str): The directory to store the transformed files.
            frameID (int | "all"): The frame ID to transform. If "all", all frames will be transformed.
            embed_in_one_file (bool): Whether to embed the transformed data in one file.
            **kwargs: Additional keyword arguments to pass to the format_transform method of the file parsers.
        Returns:
            dict[str, str | List[str]]: A dictionary containing the file paths and the transformed data.
        """
        assert os.path.isdir(output_dir), f"{output_dir} is not a directory"

        def transform_func(diskfile: FileDiskType) -> str | List[str]:
            file_path = os.path.join(
                output_dir, f"{os.path.splitext(diskfile.filename)[0]}.{format}"
            )
            frame_id = frameID if isinstance(frameID, int) else range(len(diskfile))
            try:
                return diskfile.format_transform(
                    format,
                    frameID=frame_id,  # type: ignore
                    embed_in_one_file=embed_in_one_file,
                    file_path=file_path,  # type: ignore
                )
            except ValueError as e:
                moloplogger.warning(
                    f"Format transform failed for {diskfile.filename}: {e}"
                )
                return "" if embed_in_one_file else []

        desc = "MolOP parsing with single process"
        if molopconfig.show_progress_bar:
            return {
                diskfile.file_path: transform_func(diskfile)
                for diskfile in tqdm(self, desc=desc)
            }
        else:
            return {diskfile.file_path: transform_func(diskfile) for diskfile in self}
