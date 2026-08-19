from __future__ import annotations

import importlib
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from molgr import (
    ReconstructionBatchResult,
    ReconstructionDiagnostics,
    ReconstructionFailureCode,
)
from molgr.config import CONFIG as MOLGR_CONFIG
from rdkit import Chem

from molop.config import molopconfig
from molop.io.base_models.Molecule import Molecule, reconstruct_topologies_batch
from molop.io.base_models.source import canonical_json_sha256
from molop.io.codec_types import ParseOptions
from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserMemory
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory
from molop.utils.progressbar import loky_parallel_guard, parallel_map


FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files"
CALCULATION_CASES = (
    (G16LogFileParserMemory, FIXTURE_ROOT / "g16log" / "H2O.log"),
    (
        ORCALogFileParserMemory,
        FIXTURE_ROOT / "orca" / "output_files" / "local" / "H2_sp_orca.out",
    ),
)
WATER_ATOMS = [8, 1, 1]
WATER_COORDS = np.asarray(
    [
        [0.0, 0.0, 0.0],
        [0.9572, 0.0, 0.0],
        [-0.2399872, 0.927297, 0.0],
    ]
)


def _fresh_worker_rdmol(_value: int) -> bool:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    return molecule.rdmol is not None


@pytest.mark.parametrize(("parser_cls", "fixture"), CALCULATION_CASES)
def test_real_calc_frames_capture_trusted_topology_evidence(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    parsed = parser_cls(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]
    provenance = parsed.parser_provenance
    rdmol = frame.rdmol

    assert provenance is not None
    assert rdmol is not None
    assert [atom.GetAtomicNum() for atom in rdmol.GetAtoms()] == frame.atoms
    assert frame.source_to_topology_atom_permutation == list(range(len(frame.atoms)))
    assert frame.topology_reconstruction_status == "succeeded"
    assert (
        frame.topology_reconstruction_backend
        == (provenance.effective_config["molop"]["graph_reconstruction_backend"])
    )
    assert (
        frame.topology_make_dative_bonds
        == (provenance.effective_config["molop"]["make_dative_bonds"])
    )
    assert frame.topology_reconstruction_config_sha256 == canonical_json_sha256(
        {
            "backend": frame.topology_reconstruction_backend,
            "make_dative_bonds": frame.topology_make_dative_bonds,
            "make_stereochemistry": frame.topology_make_stereochemistry,
            "molgr": provenance.effective_config["molgr"],
        }
    )
    assert frame.parse_presence["topology"] == "parsed"
    assert not any(
        diagnostic.code.startswith("MOL.PARSE.TOPOLOGY_") for diagnostic in frame.parse_diagnostics
    )

    payload = frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    assert payload["source_to_topology_atom_permutation"] == list(range(len(frame.atoms)))
    assert payload["topology_reconstruction_backend"] == frame.topology_reconstruction_backend
    assert payload["topology_make_dative_bonds"] == frame.topology_make_dative_bonds
    assert payload["topology_make_stereochemistry"] == frame.topology_make_stereochemistry
    assert payload["topology_reconstruction_config_sha256"] == (
        frame.topology_reconstruction_config_sha256
    )
    assert payload["topology_v3000_molblock"] == frame.topology_v3000_molblock
    assert "V3000" in payload["topology_v3000_molblock"]


@pytest.mark.parametrize("backend", ["cpp", "python"])
def test_molgr_backends_preserve_source_atom_order(
    backend: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", backend)
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    rdmol = molecule.rdmol

    assert rdmol is not None
    assert [atom.GetAtomicNum() for atom in rdmol.GetAtoms()] == WATER_ATOMS
    np.testing.assert_allclose(rdmol.GetConformer().GetPositions(), WATER_COORDS, atol=1e-5)
    assert molecule.source_to_topology_atom_permutation == [0, 1, 2]
    assert molecule.topology_reconstruction_backend == backend
    assert molecule.topology_reconstruction_status == "succeeded"


def test_molecule_reads_global_topology_config_at_reconstruction_time(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    captured: dict[str, Any] = {}

    def capture_config(
        xyz_block: str,
        _charge: int,
        _multiplicity: int,
        *,
        backend: str,
        make_dative_bonds: bool,
        make_stereochemistry: bool,
        config: Any,
    ) -> Chem.Mol:
        captured.update(
            backend=backend,
            make_dative_bonds=make_dative_bonds,
            make_stereochemistry=make_stereochemistry,
            config=config,
        )
        rdmol = Chem.MolFromXYZBlock(xyz_block)
        assert rdmol is not None
        return rdmol

    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", "python")
    monkeypatch.setattr(molopconfig, "make_dative_bonds", False)
    monkeypatch.setattr(molopconfig, "make_stereochemistry", False)
    monkeypatch.setattr(
        MOLGR_CONFIG.resonance,
        "max_depth",
        MOLGR_CONFIG.resonance.max_depth + 1,
    )
    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", capture_config)

    assert molecule.rdmol is not None
    assert captured == {
        "backend": "python",
        "make_dative_bonds": False,
        "make_stereochemistry": False,
        "config": MOLGR_CONFIG,
    }
    assert molecule.topology_reconstruction_backend == "python"
    assert molecule.topology_make_dative_bonds is False
    assert molecule.topology_make_stereochemistry is False


def test_parser_options_freeze_lazy_topology_settings(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    captured: dict[str, Any] = {}

    def capture_config(
        xyz_block: str,
        _charge: int,
        _multiplicity: int,
        *,
        backend: str,
        make_dative_bonds: bool,
        make_stereochemistry: bool,
        config: Any,
    ) -> Chem.Mol:
        captured.update(
            backend=backend,
            make_dative_bonds=make_dative_bonds,
            make_stereochemistry=make_stereochemistry,
            config=config,
        )
        rdmol = Chem.MolFromXYZBlock(xyz_block)
        assert rdmol is not None
        return rdmol

    parsed = XYZFileParserMemory(
        parse_options=ParseOptions(
            graph_reconstruction_backend="python",
            make_dative_bonds=False,
            make_stereochemistry=False,
        )
    ).parse("3\nwater\nO 0.0 0.0 0.0\nH 0.9572 0.0 0.0\nH -0.2399872 0.927297 0.0\n")
    frame = parsed[0]
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", "cpp")
    monkeypatch.setattr(molopconfig, "make_dative_bonds", True)
    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", capture_config)

    assert frame.rdmol is not None
    assert captured == {
        "backend": "python",
        "make_dative_bonds": False,
        "make_stereochemistry": False,
        "config": MOLGR_CONFIG,
    }


def test_batch_reconstruction_uses_native_iterator_and_restores_input_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
    ]
    molecules[2].topology_reconstruction_backend = "python"
    calls: list[dict[str, Any]] = []

    def fake_batch_iterator(requests: Any, **kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append({"count": len(request_list), **kwargs})
        for request in reversed(request_list):
            rdmol = Chem.MolFromXYZBlock(request.xyz_block)
            assert rdmol is not None
            yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", fake_batch_iterator)

    results = reconstruct_topologies_batch(molecules, max_workers=2, ordered=False)

    assert len(results) == 3
    assert len(calls) == 2
    assert [call["backend"] for call in calls] == ["cpp", "python"]
    assert all(call["max_workers"] == 1 for call in calls if call["backend"] == "python")
    assert all(call["make_stereochemistry"] is True for call in calls)
    assert all(molecule.topology_reconstruction_status == "succeeded" for molecule in molecules)
    assert all(molecule.rdmol is not None for molecule in molecules)


def test_batch_reconstruction_skips_existing_and_attempted_topologies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    existing = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    existing._rdmol = Chem.MolFromSmiles("O")
    attempted = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    attempted.topology_reconstruction_status = "failed"
    fresh = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    calls: list[int] = []

    def fake_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append(len(request_list))
        for request in request_list:
            rdmol = Chem.MolFromXYZBlock(request.xyz_block)
            assert rdmol is not None
            yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", fake_batch_iterator)

    results = reconstruct_topologies_batch([existing, attempted, fresh], max_workers=1)

    assert calls == [1]
    assert len(results) == 1
    assert fresh.topology_reconstruction_status == "succeeded"
    assert attempted.topology_reconstruction_status == "failed"


def test_batch_reconstruction_marks_missing_iterator_results_failed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
    ]

    def truncated_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        rdmol = Chem.MolFromXYZBlock(request_list[0].xyz_block)
        assert rdmol is not None
        yield ReconstructionBatchResult(request_list[0], rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", truncated_batch_iterator)

    results = reconstruct_topologies_batch(molecules, max_workers=1)

    assert len(results) == 1
    assert molecules[0].topology_reconstruction_status == "succeeded"
    assert molecules[1].topology_reconstruction_status == "failed"


def test_batch_reconstruction_rejects_duplicate_iterator_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
    ]

    def duplicate_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        rdmol = Chem.MolFromXYZBlock(request_list[0].xyz_block)
        assert rdmol is not None
        yield ReconstructionBatchResult(request_list[0], rdmol)
        yield ReconstructionBatchResult(request_list[0], rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", duplicate_batch_iterator)

    with pytest.raises(RuntimeError, match="duplicate result"):
        reconstruct_topologies_batch(molecules, max_workers=1)

    assert molecules[1].topology_reconstruction_status == "failed"


def test_batch_reconstruction_marks_request_build_failures_failed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    def invalid_xyz(_self: Molecule) -> str:
        raise ValueError("synthetic XYZ serialization failure")

    monkeypatch.setattr(Molecule, "to_XYZ", invalid_xyz)

    with pytest.raises(ValueError, match="synthetic XYZ serialization failure"):
        reconstruct_topologies_batch([molecule], max_workers=1)

    assert molecule.topology_reconstruction_status == "failed"


def test_batch_reconstruction_marks_pending_items_when_iterator_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
    ]

    def failing_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        rdmol = Chem.MolFromXYZBlock(request_list[0].xyz_block)
        assert rdmol is not None
        yield ReconstructionBatchResult(request_list[0], rdmol)
        raise RuntimeError("synthetic iterator failure")

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", failing_batch_iterator)

    with pytest.raises(RuntimeError, match="synthetic iterator failure"):
        reconstruct_topologies_batch(molecules, max_workers=1)

    assert molecules[0].topology_reconstruction_status == "succeeded"
    assert molecules[1].topology_reconstruction_status == "failed"


def test_batch_reconstruction_marks_pending_items_after_iterator_cancelled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
        Molecule.from_coords(WATER_ATOMS, WATER_COORDS),
    ]

    def cancelled_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        rdmol = Chem.MolFromXYZBlock(request_list[0].xyz_block)
        assert rdmol is not None
        yield ReconstructionBatchResult(request_list[0], rdmol)
        raise GeneratorExit

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", cancelled_batch_iterator)

    with pytest.raises(GeneratorExit):
        reconstruct_topologies_batch(molecules, max_workers=1)

    assert molecules[0].topology_reconstruction_status == "succeeded"
    assert molecules[1].topology_reconstruction_status == "failed"


def test_batch_reconstruction_rejects_loky_worker_before_mutating_models(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    monkeypatch.setattr(molecule_module, "is_loky_worker", lambda: True)

    with pytest.raises(RuntimeError, match="forbidden in a loky worker"):
        reconstruct_topologies_batch([molecule], max_workers=1)

    assert molecule.topology_reconstruction_status is None


def test_batch_reconstruction_contention_leaves_models_retryable() -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    with (
        loky_parallel_guard(),
        pytest.raises(RuntimeError, match="MolOP-managed.*task or result generator"),
    ):
        reconstruct_topologies_batch([molecule], max_workers=1)

    assert molecule.topology_reconstruction_status is None
    assert molecule.topology_reconstruction_backend is None
    assert molecule.topology_make_dative_bonds is None


def test_unmanaged_child_conflict_leaves_lazy_molecule_retryable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    progressbar_module = importlib.import_module("molop.utils.progressbar")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    monkeypatch.setattr(
        progressbar_module.multiprocessing,
        "active_children",
        lambda: [SimpleNamespace(name="external-worker")],
    )

    with pytest.raises(RuntimeError, match="unmanaged child processes"):
        _ = molecule.rdmol

    assert molecule.topology_reconstruction_status is None
    assert molecule.topology_reconstruction_backend is None
    assert molecule.topology_make_dative_bonds is None


def test_prewarmed_graph_cache_is_available_in_loky_worker() -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    reconstruct_topologies_batch([molecule], max_workers=1)

    results = parallel_map(
        lambda item: item.rdmol is not None and item.topology_reconstruction_status == "succeeded",
        [molecule],
        n_jobs=2,
        disable=True,
        return_results=True,
    )

    assert results == [True]


def test_lost_suspicious_private_graph_fails_closed_without_inventing_bonds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    def suspicious_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request = list(requests)[0]
        rdmol = Chem.MolFromXYZBlock(request.xyz_block)
        assert rdmol is not None
        rdmol.SetProp("_MolGRReconstructionStatus", "suspicious_fallback")
        yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", suspicious_batch_iterator)
    reconstruct_topologies_batch([molecule], max_workers=1, retain_results=False)

    assert molecule.topology_reconstruction_status == "suspicious_fallback"
    assert molecule.bonds == []
    molecule._rdmol = None

    assert molecule.rdmol is None


def test_batch_failure_diagnostics_survive_loky_serialization(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    diagnostic = ReconstructionDiagnostics(
        code=ReconstructionFailureCode.RDKIT_POSTPROCESS_FAILED,
        stage="rdkit.conversion",
        backend="cpp",
        message="synthetic postprocess failure",
        details={"candidate": 2},
    )

    def failed_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request = list(requests)[0]
        yield ReconstructionBatchResult(request, None, diagnostic)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", failed_batch_iterator)
    reconstruct_topologies_batch([molecule], max_workers=1, retain_results=False)

    serialized = parallel_map(
        lambda item: item.topology_reconstruction_diagnostics,
        [molecule],
        n_jobs=2,
        disable=True,
        return_results=True,
    )

    assert molecule.topology_reconstruction_status == "failed"
    assert molecule.topology_reconstruction_diagnostics == diagnostic.as_dict()
    assert list(serialized or []) == [diagnostic.as_dict()]


def test_unprewarmed_graph_is_reconstructed_in_fresh_loky_worker() -> None:
    results = parallel_map(
        _fresh_worker_rdmol,
        [0],
        n_jobs=2,
        disable=True,
        return_results=True,
    )

    assert results == [True]


def test_batch_reconstruction_is_chunked_and_can_discard_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecules = [Molecule.from_coords(WATER_ATOMS, WATER_COORDS) for _ in range(5)]
    calls: list[int] = []

    def fake_batch_iterator(requests: Any, **_kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append(len(request_list))
        for request in request_list:
            rdmol = Chem.MolFromXYZBlock(request.xyz_block)
            assert rdmol is not None
            yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", fake_batch_iterator)

    assert (
        reconstruct_topologies_batch(
            molecules,
            max_workers=1,
            batch_size=2,
            retain_results=False,
        )
        == []
    )
    assert calls == [2, 2, 1]
    assert all(molecule.rdmol is not None for molecule in molecules)


def test_batch_worker_limit_respects_windows_single_thread_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    calls: list[dict[str, Any]] = []

    def fake_batch_iterator(requests: Any, **kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append(kwargs)
        for request in request_list:
            rdmol = Chem.MolFromXYZBlock(request.xyz_block)
            assert rdmol is not None
            yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", fake_batch_iterator)
    monkeypatch.setattr(molecule_module.sys, "platform", "win32")
    monkeypatch.setattr(MOLGR_CONFIG.cpp_backend, "max_threads", 1)

    reconstruct_topologies_batch([molecule], max_workers=8, retain_results=False)

    assert calls[0]["max_workers"] == 1


def test_lazy_topology_reader_waits_for_concurrent_reconstruction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    entered = threading.Event()
    original_reconstruction = molecule_module.xyz_to_rdmol

    def slow_reconstruction(*args: Any, **kwargs: Any) -> Any:
        entered.set()
        time.sleep(0.1)
        return original_reconstruction(*args, **kwargs)

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", slow_reconstruction)
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(lambda: molecule.rdmol)
        assert entered.wait(2)
        concurrent_value = molecule.rdmol
        first_value = future.result()

    assert concurrent_value is not None
    assert first_value is not None


def test_lazy_topology_reader_fails_fast_inside_thread_parallel_region() -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    def read_graph() -> bool:
        return molecule.rdmol is not None

    with loky_parallel_guard(), ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(read_graph)
        with pytest.raises(RuntimeError, match="MolOP-managed.*task or result generator"):
            future.result(timeout=2)

    assert molecule.topology_reconstruction_status is None


def test_lazy_topology_reader_records_failure_after_base_exception(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    def cancel_reconstruction(*_args: Any, **_kwargs: Any) -> Any:
        raise KeyboardInterrupt

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", cancel_reconstruction)

    with pytest.raises(KeyboardInterrupt):
        _ = molecule.rdmol

    assert molecule.topology_reconstruction_status == "failed"


def test_chem_file_summary_implicitly_uses_native_batch_reconstruction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    parsed = XYZFileParserMemory().parse(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.9572 0.0 0.0\nH -0.2399872 0.927297 0.0\n"
    )
    calls: list[dict[str, Any]] = []

    def fake_batch_iterator(requests: Any, **kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append({"count": len(request_list), **kwargs})
        for request in request_list:
            rdmol = Chem.MolFromXYZBlock(request.xyz_block)
            assert rdmol is not None
            yield ReconstructionBatchResult(request, rdmol)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", fake_batch_iterator)

    summary = parsed.to_summary_df(frame="all")

    assert len(summary) == 1
    assert len(calls) == 1
    assert calls[0]["backend"] == "cpp"
    assert parsed[0].rdmol is not None
    assert parsed[0].topology_reconstruction_status == "succeeded"


def test_topology_v3000_molblock_is_map_free_and_does_not_mutate_cached_rdmol() -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    rdmol = molecule.rdmol
    assert rdmol is not None
    for atom_index, atom in enumerate(rdmol.GetAtoms(), start=1):
        atom.SetAtomMapNum(atom_index)
    cached_maps = [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()]
    cached_conformer_count = rdmol.GetNumConformers()

    first = molecule.topology_v3000_molblock
    second = molecule.topology_v3000_molblock

    assert first is not None
    assert first == second
    assert "V3000" in first
    portable = Chem.MolFromMolBlock(
        first,
        sanitize=False,
        removeHs=False,
        strictParsing=True,
    )
    assert portable is not None
    assert [atom.GetAtomicNum() for atom in portable.GetAtoms()] == WATER_ATOMS
    assert [atom.GetAtomMapNum() for atom in portable.GetAtoms()] == [0, 0, 0]
    assert portable.GetNumBonds() == rdmol.GetNumBonds()
    assert rdmol.GetNumConformers() == cached_conformer_count
    assert [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()] == cached_maps


def test_calc_frame_topology_reconstruction_failure_is_structured(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")

    def fail_reconstruction(*_args: Any, **_kwargs: Any) -> Any:
        raise ValueError("synthetic reconstruction failure")

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", fail_reconstruction)
    fixture = FIXTURE_ROOT / "g16log" / "H2O.log"
    parsed = G16LogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]

    assert frame.rdmol is None
    assert frame.to_canonical_SMILES() == ""
    assert frame.topology_v3000_molblock is None
    assert frame.topology_reconstruction_status == "failed"
    assert frame.topology_reconstruction_backend is not None
    assert frame.topology_reconstruction_config_sha256 is not None
    assert frame.source_to_topology_atom_permutation is None
    assert frame.parse_presence["topology"] == "parse_failed"
    assert frame.parse_completeness == "partial"
    assert parsed.parse_completeness == "partial"
    assert any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_RECONSTRUCTION_FAILED"
        for diagnostic in frame.parse_diagnostics
    )


def test_ambiguous_same_element_reordering_does_not_guess_permutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")

    def reorder_equivalent_atoms(xyz_block: str, *_args: Any, **_kwargs: Any) -> Chem.Mol:
        rdmol = Chem.MolFromXYZBlock(xyz_block)
        assert rdmol is not None
        return Chem.RenumberAtoms(rdmol, [0, 2, 1])

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", reorder_equivalent_atoms)
    fixture = FIXTURE_ROOT / "g16log" / "H2O.log"
    parsed = G16LogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]

    assert frame.rdmol is not None
    assert [atom.GetAtomicNum() for atom in frame.rdmol.GetAtoms()] == frame.atoms
    assert frame.source_to_topology_atom_permutation is None
    assert frame.topology_v3000_molblock is None
    assert frame.bonds == []
    assert frame.formal_charges == []
    assert frame.formal_num_radicals == []
    assert frame.parse_presence["topology"] == "parse_failed"
    assert frame.parse_completeness == "partial"
    assert parsed.parse_completeness == "partial"
    assert any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_ATOM_ORDER_MISMATCH"
        for diagnostic in frame.parse_diagnostics
    )
    assert not any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_RECONSTRUCTION_FAILED"
        for diagnostic in frame.parse_diagnostics
    )
