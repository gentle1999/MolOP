from collections import OrderedDict
from collections.abc import MutableMapping
from typing import List, Literal, Sequence, Union, overload

from molop.io.coords_models import GJFFileDisk, SDFFileDisk, XYZFileDisk
from molop.io.QM_models import G16LogFileDisk
from molop.logger.logger import moloplogger

FileDiskType = G16LogFileDisk | GJFFileDisk | XYZFileDisk | SDFFileDisk
QMFileDiskType = G16LogFileDisk


class FileBatchModelDisk(MutableMapping):
    """
    A class to store a batch of files and their corresponding models.
    """

    __dickfiles: OrderedDict[str, FileDiskType]

    def __init__(self):
        self.__dickfiles: OrderedDict[str, FileDiskType] = OrderedDict()

    def add_diskfiles(self, diskfiles: Sequence[FileDiskType]) -> None:
        for diskfile in diskfiles:
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
    def new_batch(cls, parsers: Sequence[FileDiskType]):
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

    def filter(
        self, condition: Literal["ts", "error", "opt", "normal"]
    ) -> "FileBatchModelDisk":
        if condition == "ts":
            return self.new_batch(
                [
                    diskfile
                    for diskfile in self
                    if isinstance(diskfile, QMFileDiskType) and diskfile[-1].is_TS
                ]
            )
        elif condition == "error":
            return self.new_batch(
                [
                    diskfile
                    for diskfile in self
                    if isinstance(diskfile, QMFileDiskType) and diskfile[-1].is_error
                ]
            )
        elif condition == "opt":
            return self.new_batch(
                [
                    diskfile
                    for diskfile in self
                    if isinstance(diskfile, QMFileDiskType)
                    and diskfile[-1].is_normal
                    and diskfile.closest_optimized_frame.is_optimized
                ]
            )
        elif condition == "normal":
            return self.new_batch(
                [
                    diskfile
                    for diskfile in self
                    if isinstance(diskfile, QMFileDiskType) and diskfile[-1].is_normal
                ]
            )
        else:
            raise ValueError(f"Invalid condition: {condition}")

    def filter_by_charge(self, charge: int) -> "FileBatchModelDisk":
        return self.new_batch(
            [diskfile for diskfile in self if diskfile[-1].charge == charge]
        )

    def filter_by_multiplicity(self, multiplicity: int) -> "FileBatchModelDisk":
        return self.new_batch(
            [diskfile for diskfile in self if diskfile[-1].multiplicity == multiplicity]
        )

    def filter_by_format(self, format: str) -> "FileBatchModelDisk":
        return self.new_batch(
            [diskfile for diskfile in self if diskfile.file_format == format]
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"
