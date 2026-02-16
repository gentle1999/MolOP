from __future__ import annotations

import operator
import os
import random
from collections import OrderedDict
from collections.abc import Callable, Iterable, MutableMapping, Sequence
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeAlias, TypeVar, cast, overload

import pandas as pd

from molop.config import molopconfig, moloplogger
from molop.io._batch_format_transform import BatchFormatTransformMixin
from molop.utils.progressbar import parallel_map


if TYPE_CHECKING:
    from molop.io._typing_catalog import FileDiskObj
else:
    FileDiskObj: TypeAlias = Any


def _looks_like_disk_file(obj: object) -> bool:
    """Compromise runtime validation.

    We avoid importing concrete file model classes here to prevent import graph
    coupling. Instead we check for a minimal, file-like surface.
    """

    if not hasattr(obj, "file_path"):
        return False
    file_path = getattr(obj, "file_path", None)
    if not isinstance(file_path, str) or not file_path:
        return False
    return all(
        hasattr(obj, name)
        for name in (
            "filename",
            "file_format",
            "__len__",
            "__getitem__",
            "format_transform",
            "to_summary_series",
            "release_file_content",
        )
    )


R = TypeVar("R")
TFile = TypeVar("TFile")
TFileDisk = TypeVar("TFileDisk", bound=FileDiskObj)


class FileBatchModelDisk(BatchFormatTransformMixin, MutableMapping, Generic[TFileDisk]):
    """
    A class to store a batch of files and their corresponding models.
    """

    def __init__(self, diskfiles: Iterable[TFileDisk] | None = None):
        self.__diskfiles: OrderedDict[str, TFileDisk] = OrderedDict()
        if diskfiles is not None:
            self.add_diskfiles(diskfiles)

    def parallel_execute(
        self, func: Callable[[TFileDisk], R], desc: str = "", n_jobs: int = 1
    ) -> list[R]:
        """
        Internal helper to execute a function over the batch in parallel or serial.
        Refactored to support Adaptive UI (Rich/Ipywidgets) and Type Hints.
        """
        return parallel_map(
            func,
            cast(Iterable[TFileDisk], self),
            n_jobs=molopconfig.set_n_jobs(n_jobs),
            desc=desc,
            total=len(self),
            disable=not molopconfig.show_progress_bar,
            maxtasks_per_child=50,
            max_nbytes=molopconfig.parallel_max_size,
        )

    def add_diskfiles(self, diskfiles: Iterable[TFileDisk]) -> None:
        """
        Add disk files to the batch.

        Parameters:
            diskfiles (Iterable[FileDiscObj]): The disk files to add.
        """
        for diskfile in sorted(diskfiles, key=lambda d: d.file_path):
            if not _looks_like_disk_file(diskfile):
                raise TypeError(f"diskfiles must be MolOP file-like objects, got {type(diskfile)}")
            if diskfile.file_path not in self.__diskfiles:
                self[diskfile.file_path] = diskfile
            else:
                moloplogger.warning(
                    f"File {diskfile.file_path} already exists in the batch, skipped"
                )

    @classmethod
    def new_batch(cls, diskfiles: Iterable[TFileDisk]):
        new_batch: FileBatchModelDisk[TFileDisk] = cls()
        new_batch.add_diskfiles(diskfiles)
        return new_batch

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self)})"

    def __contains__(self, key: object) -> bool:
        if isinstance(key, str):
            return key in self.__diskfiles
        else:
            return key in self.__diskfiles.values()

    @overload
    def __getitem__(self, key: int) -> TFileDisk: ...
    @overload
    def __getitem__(self, key: str) -> TFileDisk: ...
    @overload
    def __getitem__(self, key: slice) -> FileBatchModelDisk[TFileDisk]: ...
    @overload
    def __getitem__(self, key: Sequence[int]) -> FileBatchModelDisk[TFileDisk]: ...

    def __getitem__(
        self,
        key: int | str | slice | Sequence[int | str],
    ) -> TFileDisk | FileBatchModelDisk[TFileDisk]:
        if isinstance(key, int):
            return list(self.__diskfiles.values())[key]
        if isinstance(key, str):
            return self.__diskfiles[key]
        if isinstance(key, slice):
            return self.new_batch(list(self.__diskfiles.values())[key])
        if isinstance(key, Sequence):
            return self.new_batch([self[k] for k in key])
        raise TypeError(f"Invalid key type: {type(key)}")

    @overload
    def require(self, key: int, expected_cls: type[TFile]) -> TFile: ...

    @overload
    def require(self, key: str, expected_cls: type[TFile]) -> TFile: ...

    def require(self, key: int | str, expected_cls: type[TFile]) -> TFile:
        """Return a file with precise type hints, validating at runtime.

        This avoids hardcoding a global list of concrete file classes while still
        allowing callers to opt into format-specific APIs.
        """

        value = self[key]
        if not isinstance(value, expected_cls):
            raise TypeError(
                f"Expected {expected_cls.__name__} for key {key!r}, got {type(value).__name__}"
            )
        return value

    def __setitem__(self, key: str, value: FileDiskObj) -> None:
        if not isinstance(key, str):
            raise TypeError(f"Key should be a string, got {type(key)}")
        if not _looks_like_disk_file(value):
            raise TypeError(f"Value must be a MolOP file-like object, got {type(value)}")
        if key in self.__diskfiles:
            moloplogger.warning(
                f"File {key} already exists in the batch, replaced with the new one"
            )
        self.__diskfiles[key] = cast(TFileDisk, value)

    def __delitem__(self, key: str | int) -> None:
        if isinstance(key, int):
            key = list(self.__diskfiles.keys())[key]
        del self.__diskfiles[key]

    def __len__(self) -> int:
        return len(self.__diskfiles)

    def __iter__(self) -> FileBatchModelDisk[TFileDisk]:  # type: ignore[override]
        self.__index = 0
        return self

    def __next__(self) -> TFileDisk:
        if self.__index >= len(self.__diskfiles):
            raise StopIteration
        self.__index += 1
        return self[self.__index - 1]

    def keys(self):
        return self.__diskfiles.keys()

    def values(self):
        return self.__diskfiles.values()

    def items(self):
        return self.__diskfiles.items()

    def __add__(self, other: FileBatchModelDisk) -> FileBatchModelDisk:
        """
        Operator +: Merge two batches (Union).

        Supports passing another FileBatchModelDisk instance.

        Parameters:
            other (FileBatchModelDisk): The other batch to merge with.

        Returns:
            FileBatchModelDisk: A new batch containing all unique files from both batches.
        """
        if not isinstance(other, FileBatchModelDisk):
            return NotImplemented
        new_batch = self.new_batch(self.__diskfiles.values())
        new_batch.add_diskfiles(other.__diskfiles.values())
        return new_batch

    def __sub__(self, other: FileBatchModelDisk | Iterable[str]) -> FileBatchModelDisk:
        """
        Operator -: Remove files from the batch (Difference).

        Supports passing another batch OR a simple list of file paths.
        This avoids the IO overhead of parsing files just to exclude them.

        Parameters:
            other (Union[FileBatchModelDisk, Iterable[str]]): The other batch or list of file paths to exclude.

        Returns:
            FileBatchModelDisk: A new batch containing files not present in the other batch.
        """
        other_paths: set[str]

        if isinstance(other, FileBatchModelDisk):
            other_paths = set(other.file_paths)
        elif isinstance(other, str):
            other_paths = {other}
        elif isinstance(other, Iterable):
            other_paths = set(other)
        else:
            return NotImplemented
        filtered_files = [f for f in self.__diskfiles.values() if f.file_path not in other_paths]
        return self.new_batch(filtered_files)

    def __and__(self, other: FileBatchModelDisk | Iterable[str]) -> FileBatchModelDisk:
        """
        Operator &: Keep only files present in both (Intersection).

        Also supports passing a list of file paths.

        Parameters:
            other (Union[FileBatchModelDisk, Iterable[str]]): The other batch or list of file paths to intersect with.

        Returns:
            FileBatchModelDisk: A new batch containing files present in both batches.
        """
        other_paths: set[str]

        if isinstance(other, FileBatchModelDisk):
            other_paths = set(other.file_paths)
        elif isinstance(other, str):
            other_paths = {other}
        elif isinstance(other, Iterable):
            other_paths = set(other)
        else:
            return NotImplemented
        common_files = [f for f in self.__diskfiles.values() if f.file_path in other_paths]
        return self.new_batch(common_files)

    @property
    def file_paths(self) -> list[str]:
        """
        Get the file paths of the batch.

        Returns:
            List[str]: The file paths.
        """
        return [diskfile.file_path for diskfile in self]

    @property
    def file_names(self) -> list[str]:
        """
        Get the file names of the batch.

        Returns:
            List[str]: The file names.
        """
        return [diskfile.filename for diskfile in self]

    def filter_state(
        self,
        state: Literal["ts", "error", "opt", "normal"],
        negate: bool = False,
        n_jobs: int = 1,
    ) -> FileBatchModelDisk:
        """
        Filter the files based on their state using the new execution engine.

        Parameters:
            state (Literal["ts", "error", "opt", "normal"]): The state to filter.
            negate (bool): Whether to negate the filter.
            n_jobs (int): Number of parallel jobs to filter.

        Returns:
            FileBatchModelDisk: The filtered batch.
        """

        def judge_func(diskfile: TFileDisk) -> bool:
            last_frame = diskfile[-1]
            if state == "ts":
                return last_frame.is_TS is not None and last_frame.is_TS != negate
            elif state == "error":
                return last_frame.is_error is not None and last_frame.is_error != negate
            elif state == "opt":
                if (
                    any(
                        getattr(frame, "geometry_optimization_status", None) is not None
                        for frame in diskfile
                    )
                    and (opt_frame := getattr(diskfile, "closest_optimized_frame", None))
                    is not None
                ):
                    return (
                        last_frame.is_normal is not None
                        and last_frame.is_normal != negate
                        and opt_frame.is_optimized != negate
                    )
                return False
            elif state == "normal":
                return last_frame.is_normal is not None and last_frame.is_normal != negate
            else:
                raise ValueError(f"Invalid state: {state}")

        desc = f"Filtering {state} files {'negated' if negate else ''} with {n_jobs} jobs"
        mask = self.parallel_execute(judge_func, desc, n_jobs)
        filtered_files = [diskfile for diskfile, keep in zip(self, mask, strict=True) if keep]
        return self.new_batch(filtered_files)

    def filter_value(
        self,
        target: Literal["charge", "multiplicity", "format"],
        value: str | float | int,
        compare: Literal["==", "!=", ">", "<", ">=", "<="] = "==",
        n_jobs: int = 1,
    ) -> FileBatchModelDisk[TFileDisk]:
        """
        Filter with updated operator module and parallel execution.

        Parameters:
            target (Literal["charge", "multiplicity", "format"]): The target property to filter.
            value (Union[str, float, int]): The value to compare against.
            compare (Literal["==", "!=", ">", "<", ">=", "<="]): The comparison operator.
            n_jobs (int): Number of parallel jobs to filter.

        Returns:
            FileBatchModelDisk: The filtered batch.
        """
        compare_ops = {
            "==": operator.eq,
            "!=": operator.ne,
            ">": operator.gt,
            "<": operator.lt,
            ">=": operator.ge,
            "<=": operator.le,
        }
        op_func = compare_ops[compare]

        def judge_func(diskfile: TFileDisk) -> bool:
            if target == "charge":
                return op_func(diskfile[-1].charge, value)
            elif target == "multiplicity":
                return op_func(diskfile[-1].multiplicity, value)
            elif target == "format":
                return op_func(diskfile.file_format, value)
            return False

        desc = f"Filtering files with {target} {compare} {value} with {n_jobs} jobs"
        mask = self.parallel_execute(judge_func, desc, n_jobs)
        filtered_files = [diskfile for diskfile, keep in zip(self, mask, strict=True) if keep]
        return self.new_batch(filtered_files)

    def filter_custom(
        self,
        condition: Callable[[TFileDisk], bool],
        n_jobs: int = 1,
    ) -> FileBatchModelDisk[TFileDisk]:
        """
        Filter the files based on a custom callable condition.

        This allows for complex filtering logic that goes beyond simple state or value checks.

        Parameters:
            condition (Callable[[FileDiscObj], bool]): A function that accepts a diskfile object and returns True (keep) or False (discard).
            n_jobs (int): Number of parallel jobs.

        Returns:
            FileBatchModelDisk: A new batch containing only the files that satisfy the condition.
        """
        desc = f"Filtering with custom function with {n_jobs} jobs"
        mask = self.parallel_execute(condition, desc, n_jobs)
        filtered_files = [diskfile for diskfile, keep in zip(self, mask, strict=True) if keep]
        return self.new_batch(filtered_files)

    def groupby(
        self, key_func: Callable[[TFileDisk], str], n_jobs: int = 1
    ) -> dict[str, FileBatchModelDisk[TFileDisk]]:
        """
        Group the files into multiple batches based on a key function.

        Parameters:
            key_func (Callable[[FileDiscObj], str]): A function that generates a group key for each file.
            n_jobs (int): Number of parallel jobs to calculate keys.

        Returns:
            (dict[str, FileBatchModelDisk]): A dictionary where keys are group names and values are new batch objects.
        """
        desc = f"Grouping files with {n_jobs} jobs"
        keys = self.parallel_execute(key_func, desc, n_jobs)
        temp_groups: dict[str, list[TFileDisk]] = {}
        for diskfile, key in zip(self, keys, strict=True):
            if key not in temp_groups:
                temp_groups[key] = []
            temp_groups[key].append(diskfile)
        return {group_key: self.new_batch(files) for group_key, files in temp_groups.items()}

    def filter_by_codec_id(
        self,
        codec_id: str,
        *,
        negate: bool = False,
        n_jobs: int = 1,
        on_missing: Literal["keep", "drop", "error"] = "drop",
    ) -> FileBatchModelDisk[TFileDisk]:
        """Filter batch by detected reader codec id.

        This differs from `filter_value(target="format", ...)` which filters by file suffix.

        Parameters:
            codec_id: Reader format id (e.g. "xyz", "g16log"). Compared against each diskfile's
                `detected_format_id` populated during parsing.
            negate: If True, invert the predicate.
            n_jobs: Parallel jobs.
            on_missing: Behavior when `detected_format_id` is missing/None on an item.
                - "drop": drop the item
                - "keep": keep the item
                - "error": raise ValueError
        """

        wanted = codec_id.strip().lower()
        if not wanted:
            raise ValueError("codec_id must be non-empty")

        def judge_func(diskfile: TFileDisk) -> bool:
            detected = getattr(diskfile, "detected_format_id", None)
            if detected is None:
                if on_missing == "keep":
                    return True
                if on_missing == "drop":
                    return False
                raise ValueError(
                    f"Missing detected_format_id for {getattr(diskfile, 'file_path', '?')}"
                )
            hit = str(detected).strip().lower() == wanted
            return (not hit) if negate else hit

        desc = f"Filtering files with codec_id == {wanted} with {n_jobs} jobs"
        mask = self.parallel_execute(judge_func, desc, n_jobs)
        filtered_files = [diskfile for diskfile, keep in zip(self, mask, strict=True) if keep]
        return self.new_batch(filtered_files)

    def to_summary_df(
        self,
        mode: Literal["file", "frame"] = "frame",
        frameIDs: int | Sequence[int] = -1,
        n_jobs: int = 1,
        **kwargs: Any,
    ) -> pd.DataFrame:
        """
        Generate summary DataFrame using the parallel execution engine.

        Parameters:
            mode (Literal["file", "frame"]): Whether to generate file-level or frame-level summary.
            frameIDs (int | Sequence[int]): Frame IDs to process. Use -1 for all frames.
            n_jobs (int): Number of parallel jobs.
            **kwargs: Additional arguments for to_summary_series.

        Returns:
            pd.DataFrame: A DataFrame containing the summary information.
        """
        if mode not in ["file", "frame"]:
            raise ValueError(f"Invalid mode: {mode}")
        current_frameIDs = [frameIDs] if isinstance(frameIDs, int) else frameIDs

        def process_file_summary(diskfile: FileDiskObj) -> list[pd.Series]:
            if mode == "file":
                return [diskfile.to_summary_series(**kwargs)]
            elif mode == "frame":
                return [
                    diskfile[fid].to_summary_series(**kwargs)
                    for fid in current_frameIDs
                    if fid < len(diskfile) and fid >= -len(diskfile)
                ]
            return []

        desc = f"MolOP processing {mode} summary with {n_jobs} jobs"
        nested_results = self.parallel_execute(process_file_summary, desc, n_jobs)
        series_list: list[pd.Series] = [s for sublist in nested_results for s in sublist]
        if not series_list:
            return pd.DataFrame()
        df = pd.concat(series_list, axis=1).T
        top_level_order = df.columns.get_level_values(0).unique()
        return pd.DataFrame(df[top_level_order])

    def release_file_content(self) -> None:
        """
        Release the content of the files in the batch.
        """
        for diskfile in self:
            diskfile.release_file_content()

    def draw_grid_image(
        self,
        molsPerRow: int = 4,
        subImgSize: tuple[int, int] = (200, 200),
        maxMols: int = 16,
        useSVG: bool = True,
        n_jobs: int = 1,
        **kwargs: Any,
    ):
        """
        Quick preview of the structures of the first N molecules in a Jupyter Notebook.

        Parameters:
            molsPerRow (int): Number of molecules per row in the grid.
            subImgSize (tuple[int, int]): Size of each sub-image in the grid.
            maxMols (int): Maximum number of molecules to visualize.
            useSVG (bool): Whether to use SVG format for visualization.
            n_jobs (int): Number of parallel jobs.
            **kwargs: Additional arguments for Draw.MolsToGridImage.

        Returns:
            (PIL.Image.Image | None): A grid image of the molecules, or None if RDKit is not installed.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
        except ImportError:
            moloplogger.warning("RDKit is not installed, visualization skipped.")
            return None

        mols = []
        legends = []
        for file_path, smi in (
            self[:maxMols].format_transform("smi", output_dir=None, n_jobs=n_jobs).items()
        ):
            smi = smi[-1] if isinstance(smi, list) and smi else smi
            if isinstance(smi, str) and smi:
                mol = Chem.MolFromSmiles(smi.strip())
                if mol:
                    mols.append(mol)
                    legends.append(os.path.basename(file_path))
                else:
                    mols.append(Chem.MolFromSmiles(""))
                    legends.append(os.path.basename(file_path))
        if not mols:
            moloplogger.warning("No valid molecules found for visualization.")
            return None
        return Draw.MolsToGridImage(
            mols,
            molsPerRow=molsPerRow,
            subImgSize=subImgSize,
            legends=legends,
            useSVG=useSVG,
            maxMols=maxMols,
            **kwargs,
        )

    def sample(self, n: int = 10, seed: int | None = None) -> FileBatchModelDisk:
        """
        Randomly sample n files from the batch.

        Parameters:
            n (int): Number of files to sample.
            seed (int): Random seed for reproducibility.

        Returns:
            FileBatchModelDisk: A new batch containing the sampled files.
        """
        if seed is not None:
            random.seed(seed)
        n = min(n, len(self))
        sampled_files = random.sample(list(self.__diskfiles.values()), n)
        return self.new_batch(sampled_files)
