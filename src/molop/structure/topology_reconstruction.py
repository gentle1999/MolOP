"""Coordinate-to-topology reconstruction services.

The service layer owns MolGR calls and their native execution boundary.  Model
classes only provide explicit coordinate data and apply the returned molecule
to their own cache/state, which keeps reconstruction independent from parser
and batch-model implementations.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import AbstractContextManager
from dataclasses import dataclass
from typing import Literal, TypeAlias

from molgr import ReconstructionBatchRequest, ReconstructionBatchResult, iter_xyz_to_rdmol_batch
from molgr.config import MolGRConfig
from molgr.interface import xyz_to_rdmol

from molop.structure.topology import build_mol_from_atoms_and_bonds
from molop.utils.progressbar import (
    NativeReconstructionConcurrencyError,
    is_loky_worker,
    native_reconstruction_guard,
)
from molop.utils.types import RdMol


ReconstructionOptions: TypeAlias = tuple[
    Literal["cpp", "python"],
    Literal["raise", "return_suspicious"],
    bool,
    bool,
]
NativeGuardFactory: TypeAlias = Callable[[], AbstractContextManager[None]]
WorkerProcessChecker: TypeAlias = Callable[[], bool]
SingleReconstructor: TypeAlias = Callable[..., RdMol | None]
BatchIterator: TypeAlias = Callable[..., Iterable[ReconstructionBatchResult]]
BatchStartedHook: TypeAlias = Callable[[Sequence["TopologyReconstructionTask"]], None]
SingleStartedHook: TypeAlias = Callable[[], None]


@dataclass(frozen=True, slots=True)
class TopologyReconstructionTask:
    """Explicit input for one coordinate-to-topology reconstruction."""

    index: int
    xyz_block: str
    total_charge: int
    spin_multiplicity: int
    options: ReconstructionOptions


def materialize_topology_from_fields(
    atoms: Sequence[int | str],
    bonds: Sequence[tuple[int, int, int, int]],
    formal_charges: Sequence[int] | None,
    formal_num_radicals: Sequence[int] | None,
    *,
    coords: Sequence[tuple[float, float, float]] | None = None,
) -> RdMol:
    """Materialize a trusted RDKit graph from serialized topology fields."""

    if not formal_charges or not formal_num_radicals:
        raise ValueError("If bonds given, formal charges and spins must be provided.")
    rdmol = build_mol_from_atoms_and_bonds(
        atoms,
        bonds,
        formal_charges,
        formal_num_radicals,
        coords=coords,
    )
    if rdmol is None:
        raise ValueError("Building the provided topology returned no molecule.")
    return rdmol


def reconstruct_single_topology(
    xyz_block: str,
    total_charge: int,
    spin_multiplicity: int,
    *,
    backend: Literal["cpp", "python"],
    make_dative_bonds: bool,
    make_stereochemistry: bool,
    config: MolGRConfig,
    reconstructor: SingleReconstructor = xyz_to_rdmol,
    native_guard: NativeGuardFactory = native_reconstruction_guard,
    on_reconstruction_started: SingleStartedHook | None = None,
) -> RdMol:
    """Reconstruct one topology inside the process-wide native boundary."""

    with native_guard():
        if on_reconstruction_started is not None:
            on_reconstruction_started()
        reconstructed = reconstructor(
            xyz_block,
            total_charge,
            spin_multiplicity,
            backend=backend,
            make_dative_bonds=make_dative_bonds,
            make_stereochemistry=make_stereochemistry,
            config=config,
        )
    if reconstructed is None:
        raise ValueError("MolGR topology reconstruction returned no molecule")
    return reconstructed


def resolve_reconstruction_worker_count(
    max_workers: int | None,
    backend: Literal["cpp", "python"],
    *,
    default_worker_count: int,
    configured_max_threads: int | None,
    platform_name: str = sys.platform,
) -> int | None:
    """Resolve a MolGR worker count from explicit runtime configuration."""

    if backend == "python":
        return 1

    worker_count = default_worker_count if max_workers is None else max_workers
    if configured_max_threads is not None:
        worker_count = min(worker_count, configured_max_threads)
    if platform_name == "win32" and configured_max_threads == 1:
        worker_count = 1
    return worker_count


def iter_reconstruct_topology_tasks(
    tasks: Iterable[TopologyReconstructionTask],
    *,
    max_workers: int | None,
    queue_size: int,
    ordered: bool,
    raise_on_error: bool,
    config: MolGRConfig,
    default_worker_count: int,
    configured_max_threads: int | None,
    platform_name: str = sys.platform,
    batch_iterator: BatchIterator = iter_xyz_to_rdmol_batch,
    native_guard: NativeGuardFactory = native_reconstruction_guard,
    worker_process_checker: WorkerProcessChecker = is_loky_worker,
    on_batch_started: BatchStartedHook | None = None,
) -> Iterator[tuple[int, ReconstructionBatchResult]]:
    """Run one bounded set of reconstruction tasks grouped by options.

    Results are yielded with their caller-provided task index.  The service
    deliberately does not retain results or know about model objects; callers
    can apply each result immediately and decide whether to preserve it.
    """

    if max_workers is not None and max_workers < 1:
        raise ValueError("max_workers must be >= 1 when provided")
    if queue_size < 1:
        raise ValueError("queue_size must be >= 1")
    if worker_process_checker():
        raise RuntimeError(
            "MolGR topology reconstruction is forbidden in a loky worker; "
            "prewarm the topology in the parent process first."
        )

    task_list = list(tasks)
    grouped: dict[ReconstructionOptions, list[TopologyReconstructionTask]] = defaultdict(list)
    for task in task_list:
        grouped[task.options].append(task)

    for (
        backend,
        _failure_policy,
        make_dative_bonds,
        make_stereochemistry,
    ), entries in grouped.items():
        worker_count = resolve_reconstruction_worker_count(
            max_workers,
            backend,
            default_worker_count=default_worker_count,
            configured_max_threads=configured_max_threads,
            platform_name=platform_name,
        )
        requests = [
            ReconstructionBatchRequest(
                task.xyz_block,
                total_charge=task.total_charge,
                spin_multiplicity=task.spin_multiplicity,
            )
            for task in entries
        ]
        task_by_request_id = {
            id(request): task for request, task in zip(requests, entries, strict=True)
        }
        received_request_ids: set[int] = set()

        # The guard covers both request submission and result consumption.  A
        # model-state hook can therefore freeze provenance only after the
        # native boundary is successfully acquired.
        with native_guard():
            if on_batch_started is not None:
                on_batch_started(entries)
            for result in batch_iterator(
                requests,
                backend=backend,
                max_workers=worker_count,
                queue_size=queue_size,
                ordered=ordered,
                make_dative_bonds=make_dative_bonds,
                make_stereochemistry=make_stereochemistry,
                config=config,
                raise_on_error=raise_on_error,
            ):
                request_id = id(result.input)
                matched_task = task_by_request_id.get(request_id)
                if matched_task is None:
                    raise RuntimeError("MolGR batch returned a result for an unknown request")
                if request_id in received_request_ids:
                    raise RuntimeError("MolGR batch returned a duplicate result for one request")
                received_request_ids.add(request_id)
                yield matched_task.index, result


__all__ = [
    "BatchIterator",
    "NativeReconstructionConcurrencyError",
    "ReconstructionOptions",
    "TopologyReconstructionTask",
    "iter_reconstruct_topology_tasks",
    "is_loky_worker",
    "materialize_topology_from_fields",
    "native_reconstruction_guard",
    "reconstruct_single_topology",
    "resolve_reconstruction_worker_count",
    "xyz_to_rdmol",
]
