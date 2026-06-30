import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from molop.io import AutoParser


ORCA_OUTPUT_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "orca" / "output_files"
MANIFEST_PATH = ORCA_OUTPUT_FIXTURE_DIR / "manifest.json"


def _manifest() -> dict[str, Any]:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _fixture_entries() -> list[dict[str, Any]]:
    return list(_manifest()["fixtures"])


def _fixture_path(entry: dict[str, Any]) -> Path:
    return ORCA_OUTPUT_FIXTURE_DIR / str(entry["path"])


def _fixture_text(entry: dict[str, Any]) -> str:
    return _fixture_path(entry).read_text(encoding="utf-8", errors="replace")


def _entries_with_feature(feature: str) -> Iterator[dict[str, Any]]:
    for entry in _fixture_entries():
        if feature in entry["features"]:
            yield entry


def test_orca_output_fixture_manifest_inventory_is_broad() -> None:
    entries = _fixture_entries()
    assert len(entries) == 31
    assert sum(1 for entry in entries if entry["source"] == "local") == 9
    assert sum(1 for entry in entries if entry["source"] == "cclib") == 22

    required_features = {
        "single_point",
        "gradient",
        "numerical_gradient",
        "optimization",
        "frequency",
        "vibrations",
        "ir",
        "raman",
        "polarizability",
        "nmr",
        "spin_spin_coupling",
        "solvation",
        "cpcm",
        "smd",
        "dispersion",
        "mp2",
        "mp3",
        "ccsd",
        "ccsd_t",
        "excited_state",
        "tddft",
        "adc2",
        "eom_ccsd",
        "steom_ccsd",
        "steom_dlpno_ccsd",
        "rocis",
    }
    available_features = {feature for entry in entries for feature in entry["features"]}
    assert required_features <= available_features


def test_orca_output_cclib_license_is_preserved() -> None:
    license_path = ORCA_OUTPUT_FIXTURE_DIR / _manifest()["sources"]["cclib"]["license_file"]
    text = license_path.read_text(encoding="utf-8")
    assert "BSD 3-Clause License" in text
    assert "Copyright (c) 2024, the cclib development team" in text


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["path"])
def test_orca_output_fixture_files_exist_and_are_orca_outputs(entry: dict[str, Any]) -> None:
    path = _fixture_path(entry)
    assert path.is_file()
    assert path.stat().st_size > 1000

    text = _fixture_text(entry)
    assert "Program Version" in text
    assert "ORCA" in text[:12000]


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["path"])
def test_orca_output_fixture_raw_anchors_are_present(entry: dict[str, Any]) -> None:
    text = _fixture_text(entry)
    for anchor in entry["anchors"]:
        assert anchor in text


def test_orca_output_fixture_versions_cover_legacy_and_modern_orca() -> None:
    texts = [_fixture_text(entry) for entry in _fixture_entries()]
    assert any("Program Version 4.1.1" in text for text in texts)
    assert any("Program Version 4.2.1" in text for text in texts)
    assert any("Program Version 5.0.3" in text for text in texts)
    assert any("Program Version 6.0.0" in text for text in texts)
    assert any("Program Version 6.0.1" in text for text in texts)
    assert any("Program Version 6.1.0" in text for text in texts)


def test_orca_output_auto_detection_prefers_orca_for_shared_out_suffix() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "local" / "H2_sp_orca.out"
    batch = AutoParser(str(path), parser_detection="auto", n_jobs=1, only_last_frame=True)

    assert len(batch) == 1
    assert batch[0].detected_format_id == "orcaout"
    assert batch[0].qm_software == "ORCA"


@pytest.mark.parametrize(
    "feature, minimum",
    [
        ("gradient", 5),
        ("optimization", 3),
        ("vibrations", 3),
        ("solvation", 3),
        ("excited_state", 7),
        ("nmr", 3),
        ("polarizability", 3),
    ],
)
def test_orca_output_fixture_feature_counts(feature: str, minimum: int) -> None:
    assert sum(1 for _entry in _entries_with_feature(feature)) >= minimum


@pytest.mark.parametrize("entry", _fixture_entries(), ids=lambda entry: entry["path"])
def test_orca_output_structured_parse_contract(entry: dict[str, Any]) -> None:
    path = _fixture_path(entry)
    batch = AutoParser(str(path), parser_detection="orcaout", n_jobs=1)

    assert len(batch) == 1
    file_model = batch[0]
    expected = entry["structured_expectation"]
    assert len(file_model) >= expected["min_frames"]
    assert file_model.qm_software == "ORCA"
    assert file_model.qm_software_version.startswith(expected["version_prefix"])

    last_frame = file_model[-1]
    assert last_frame.qm_software == "ORCA"
    assert last_frame.qm_software_version.startswith(expected["version_prefix"])
    if expected.get("normal_terminated"):
        assert file_model.status is not None
        assert file_model.status.normal_terminated is True
        assert last_frame.status is not None
        assert last_frame.status.normal_terminated is True
    if energy := expected.get("final_single_point_energy_hartree"):
        assert last_frame.energies is not None
        assert last_frame.energies.total_energy is not None
        assert last_frame.energies.total_energy.to("hartree").magnitude == pytest.approx(energy)
    if expected.get("has_forces"):
        assert last_frame.forces is not None
    if expected.get("has_vibrations"):
        assert last_frame.vibrations is not None
        assert len(last_frame.vibrations) > 0
    if expected.get("has_charge_spin_populations"):
        assert last_frame.charge_spin_populations is not None
    if expected.get("has_polarizability"):
        assert last_frame.polarizability is not None
    if expected.get("has_optimization_status"):
        assert last_frame.geometry_optimization_status is not None
    if expected.get("has_solvation"):
        assert file_model.solvent is not None or last_frame.solvent is not None
    if expected.get("has_electronic_states"):
        assert last_frame.electronic_states is not None
        assert len(last_frame.electronic_states) > 0
