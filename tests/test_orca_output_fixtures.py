import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from molop.io import AutoParser
from molop.io.logic.orca.log.frame_parsers._orca_extractors import (
    extract_orca_atomic_masses,
    extract_orca_energies,
    extract_orca_forces,
    extract_orca_vibrations,
)
from molop.io.logic.orca.log.frame_parsers.ORCALogFileFrameParser import (
    ORCALogFileFrameParserMemory,
)
from molop.io.logic.orca.log.locators import locate_orca_job_frames, locate_orca_jobs
from molop.io.logic.orca.log.parsers._orca_log_shared import extract_orca_status
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory
from molop.unit import atom_ureg


ORCA_OUTPUT_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "orca" / "output_files"
MANIFEST_PATH = ORCA_OUTPUT_FIXTURE_DIR / "manifest.json"

pytestmark = pytest.mark.orca_output_corpus

if not MANIFEST_PATH.is_file():
    pytest.skip(
        "optional ORCA output fixture corpus manifest is not available",
        allow_module_level=True,
    )


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


def test_orca_closed_shell_population_schemes_are_all_structured() -> None:
    source = (ORCA_OUTPUT_FIXTURE_DIR / "local" / "h2o_orca_v5_charges.out").read_text()
    populations = ORCALogFileParserMemory().parse(source)[-1].charge_spin_populations

    assert populations is not None
    assert populations.population_names == [
        "mulliken_charges",
        "lowdin_charges",
        "hirshfeld_charges",
        "hirshfeld_spins",
    ]
    assert all(len(series.values) == 3 for _name, series in populations.population_items())


def test_orca_open_shell_mulliken_and_lowdin_charge_spin_tables_are_structured() -> None:
    source = (ORCA_OUTPUT_FIXTURE_DIR / "cclib" / "basicORCA6.0" / "dvb_rocis.out").read_text()
    populations = ORCALogFileParserMemory().parse(source)[-1].charge_spin_populations

    assert populations is not None
    assert len(populations["mulliken_charges"].values) == 20
    assert len(populations["mulliken_spins"].values) == 20
    assert len(populations["lowdin_charges"].values) == 20
    assert len(populations["lowdin_spins"].values) == 20


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


def test_orca_output_metadata_result_is_model_ready() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "local" / "H2_sp_orca.out"
    parser = ORCALogFileParserMemory()

    metadata = parser._parse_segment_metadata_result(
        path.read_text(encoding="utf-8", errors="replace")
    ).model_data()

    assert metadata["qm_software"] == "ORCA"
    assert metadata["qm_software_version"] == "4.1.1"
    assert metadata["input_file_name"] == "H2_sp.inp"
    assert metadata["model_chemistry"].method_family == "DFT"
    assert metadata["model_chemistry"].functional == "PBE"
    assert metadata["model_chemistry"].basis_set == "def2-SVP"
    assert metadata["task_requests"][0].task_type == "sp"
    assert metadata["charge"] == 0
    assert metadata["multiplicity"] == 1
    assert metadata["status"].normal_terminated is True
    assert metadata["status"].scf_converged is True
    assert metadata["running_time"].to("second").magnitude > 0


def test_orca_parses_masses_from_cartesian_au_coordinates() -> None:
    source = (ORCA_OUTPUT_FIXTURE_DIR / "local" / "h2o_orca_v5_charges.out").read_text()
    frame = ORCALogFileParserMemory().parse(source)[-1]

    assert frame.atomic_masses_source == "orca_cartesian_au_mass"
    assert frame.atomic_masses is not None
    np.testing.assert_allclose(frame.atomic_masses.m_as("amu"), [15.999, 1.008, 1.008])


def test_orca_atomic_mass_extractor_uses_latest_table_and_supports_d_notation() -> None:
    source = """
----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     9.000    0.000000    0.000000    0.000000
----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
  NO LB      ZA    FRAG     MASS         X           Y           Z
   0 H     1.0000    0     1.007825D+00    0.000000    0.000000    0.000000
   1 O     8.0000    0     15.994915    1.000000    0.000000    0.000000
 next section
 """

    masses = extract_orca_atomic_masses(source, expected_atom_count=2)

    assert masses is not None
    assert masses.units == atom_ureg.amu
    np.testing.assert_allclose(masses.magnitude, [1.007825, 15.994915])


def test_orca_atomic_mass_extractor_rejects_missing_or_invalid_tables() -> None:
    header = """
----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
"""
    non_contiguous = (
        header
        + """
   0 H     1.0000    0     1.008    0.000000    0.000000    0.000000
   2 H     1.0000    0     1.008    1.000000    0.000000    0.000000
"""
    )
    one_row = (
        header
        + """
   0 H     1.0000    0     1.008    0.000000    0.000000    0.000000
"""
    )

    assert extract_orca_atomic_masses("no coordinate table") is None
    assert extract_orca_atomic_masses(header) is None
    assert extract_orca_atomic_masses(non_contiguous) is None
    assert extract_orca_atomic_masses(one_row, expected_atom_count=2) is None


def test_orca_log_state_machine_rejects_unexpected_phase(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = ORCALogFileFrameParserMemory()

    monkeypatch.setattr(parser, "_run_structure_phase", lambda _text, _result: object())

    with pytest.raises(AssertionError, match="Unexpected ORCA log frame parse phase"):
        parser._parse_block_to_result("")


def test_orca_log_metadata_state_machine_rejects_unexpected_phase(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = ORCALogFileParserMemory()

    monkeypatch.setattr(parser, "_run_software_metadata_phase", lambda _context, _result: object())

    with pytest.raises(AssertionError, match="Unexpected ORCA log metadata parse phase"):
        parser._parse_segment_metadata_result("")


def test_orca_cartesian_gradient_is_exported_as_elementwise_negative_force() -> None:
    text = """
------------------
CARTESIAN GRADIENT
------------------

   1   H   :   -0.125000000    0.250000000   -0.500000000
   2   O   :    1.250000000   -2.500000000    5.000000000
"""

    forces = extract_orca_forces(text, num_atoms=2)

    assert forces is not None
    assert str(forces.units) == "hartree / bohr"
    np.testing.assert_allclose(
        forces.magnitude,
        [[0.125, -0.25, 0.5], [-1.25, 2.5, -5.0]],
        rtol=0.0,
        atol=0.0,
    )


def test_orca_mp2_gradient_norm_without_rows_does_not_fabricate_forces() -> None:
    text = "NORM OF THE MP2 GRADIENT:  0.063637\n"

    assert extract_orca_forces(text, num_atoms=6) is None


def test_orca_mp2_gradient_rows_are_exported_as_forces() -> None:
    text = """
The final MP2 gradient
  0:  -0.00971201  -0.00773534  -0.02473580
  1:  -0.00329997   0.00499343   0.00012475

NORM OF THE MP2 GRADIENT:  0.030000
"""

    forces = extract_orca_forces(text, num_atoms=2)

    assert forces is not None
    np.testing.assert_allclose(
        forces.magnitude,
        [[0.00971201, 0.00773534, 0.02473580], [0.00329997, -0.00499343, -0.00012475]],
        rtol=0.0,
        atol=0.0,
    )


def test_orca_vibration_trimming_preserves_printed_mode_indices() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "local" / "H2O_hess_orca.out"

    vibrations = extract_orca_vibrations(
        path.read_text(encoding="utf-8", errors="replace"),
        num_atoms=3,
    )

    assert vibrations is not None
    assert vibrations.mode_indices == [6, 7, 8]
    assert vibrations.frequencies.magnitude.tolist() == pytest.approx([1567.61, 3467.70, 3651.46])


def test_orca_final_geometry_does_not_inherit_previous_frame_forces() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "local" / "opt_orca.out"
    file_model = AutoParser(str(path), parser_detection="orcaout", n_jobs=1)[0]

    assert any(frame.forces is not None for frame in file_model[:-1])
    assert file_model[-1].forces is None


def test_orca_mp2_correlation_component_is_not_treated_as_total_energy() -> None:
    energies = extract_orca_energies(
        """
Total Energy       :        -74.90000000000000 Eh
E(MP2)                                     ...     -0.035800372
FINAL SINGLE POINT ENERGY                  -74.935800372000
"""
    )

    assert energies is not None
    assert energies.mp2_energy is None
    assert energies.reference_energy.to("hartree").magnitude == pytest.approx(-74.9)
    assert energies.electronic_energy.to("hartree").magnitude == pytest.approx(-74.935800372)


def test_orca_energy_extraction_keeps_ccsd_and_ccsd_t_totals_distinct() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "cclib" / "basicORCA6.0" / "water_ccsd_t.out"
    energies = extract_orca_energies(path.read_text(encoding="utf-8", errors="replace"))

    assert energies is not None
    assert energies.ccsd_energy.to("hartree").magnitude == pytest.approx(-75.013487814)
    assert energies.ccsd_t_energy.to("hartree").magnitude == pytest.approx(-75.013556306)
    assert energies.total_energy == energies.electronic_energy


def test_orca_coupled_cluster_energy_block_is_typed_as_ccsd_total() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "cclib" / "basicORCA6.0" / "water_ccsd.out"
    energies = extract_orca_energies(path.read_text(encoding="utf-8", errors="replace"))

    assert energies is not None
    assert energies.ccsd_energy.to("hartree").magnitude == pytest.approx(-75.013487814)
    assert energies.ccsd_t_energy is None
    assert energies.total_energy == energies.electronic_energy


def test_orca_energy_source_semantics_are_opt_in() -> None:
    text = """
Total Energy       :        -74.90000000000000 Eh
MP2 TOTAL ENERGY:           -74.93580037200000 Eh
MP2 CORRELATION ENERGY :     -0.03580037200000 Eh
E(MP3) =                    -74.94100000000000
EC(MP3) =                    -0.04100000000000
E3 =                         -0.00519962800000
FINAL SINGLE POINT ENERGY   -74.94100000000000
"""

    plain = extract_orca_energies(text)
    captured = extract_orca_energies(text, capture_source_evidence=True)

    assert plain is not None
    assert plain.observations == []
    assert captured is not None
    assert {
        (observation.method, observation.quantity_semantics, observation.source_label)
        for observation in captured.observations
    } >= {
        ("reference", "total_energy", "Total Energy"),
        ("MP2", "total_energy", "MP2 TOTAL ENERGY"),
        ("MP2", "correlation_correction", "MP2 CORRELATION ENERGY"),
        ("MP3", "total_energy", "E(MP3)"),
        ("MP3", "correlation_correction", "EC(MP3)"),
        ("MP3", "component", "E3"),
        ("electronic", "total_energy", "FINAL SINGLE POINT ENERGY"),
    }


def test_orca_frame_source_metadata_is_opt_in() -> None:
    path = ORCA_OUTPUT_FIXTURE_DIR / "local" / "opt_orca.out"
    source = path.read_text(encoding="utf-8", errors="replace")
    jobs = locate_orca_jobs(source)
    block = locate_orca_job_frames(source, jobs[0])[0].text(source)

    plain = ORCALogFileFrameParserMemory().parse(block)
    captured = ORCALogFileFrameParserMemory(capture_source_evidence=True).parse(block)

    assert plain.energies is not None
    assert plain.energies.observations == []
    assert plain.coordinate_source is None
    assert plain.coordinate_provenance is None
    assert plain.coordinate_decimal_places is None
    assert plain.force_source_field is None
    assert plain.force_transformation is None
    assert plain.forces is not None
    assert plain.forces_axis_order == ("atom", "cartesian")
    assert plain.forces_atom_order == "source"
    assert plain.forces_orientation == "source"
    assert plain.geometry_optimization_status is not None
    assert plain.geometry_optimization_status.source_converged is None
    assert plain.geometry_optimization_status.source_labels is None

    assert captured.energies is not None
    assert captured.energies.observations
    assert captured.coordinate_source == "observed"
    assert captured.coordinate_provenance is not None
    assert captured.coordinate_decimal_places is not None
    assert captured.coordinate_decimal_places > 0
    assert captured.force_source_field == "gradient"
    assert captured.force_transformation is not None
    assert captured.forces_axis_order == ("atom", "cartesian")
    assert captured.forces_atom_order == "source"
    assert captured.forces_orientation == "source"
    assert captured.geometry_optimization_status is not None
    assert captured.geometry_optimization_status.source_converged
    assert captured.geometry_optimization_status.source_labels
    assert captured.geometry_optimization_status.rms_force is not None
    assert str(captured.geometry_optimization_status.rms_force.units) == "hartree / bohr"
    assert captured.geometry_optimization_status.rms_displacement is not None
    assert str(captured.geometry_optimization_status.rms_displacement.units) == "bohr"


@pytest.mark.parametrize(
    ("text", "expected_normal", "expected_scf"),
    [
        ("****ORCA TERMINATED NORMALLY****", True, None),
        ("SCF CONVERGED AFTER 8 CYCLES\n****ORCA TERMINATED NORMALLY****", True, True),
        ("SCF NOT CONVERGED AFTER 125 CYCLES\nORCA TERMINATED ABNORMALLY", False, False),
        (
            "SCF CONVERGED AFTER 8 CYCLES\nSCF DID NOT CONVERGE\nORCA TERMINATED ABNORMALLY",
            False,
            False,
        ),
        (
            "SCF DID NOT CONVERGE\nSCF CONVERGED AFTER 9 CYCLES\nORCA TERMINATED ABNORMALLY",
            False,
            True,
        ),
    ],
)
def test_orca_status_uses_explicit_scf_evidence(
    text: str, expected_normal: bool, expected_scf: bool | None
) -> None:
    status = extract_orca_status(text)

    assert status is not None
    assert status.normal_terminated is expected_normal
    assert status.scf_converged is expected_scf


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
    batch = AutoParser(
        str(path),
        parser_detection="orcaout",
        n_jobs=1,
        capture_source_evidence=True,
    )

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
        assert last_frame.segment_index is not None
        last_segment = next(
            segment
            for segment in file_model.source_segments
            if segment.segment_index == last_frame.segment_index
        )
        assert last_segment.termination_status is True
        assert last_segment.parse_presence["termination_status"] == "parsed"
        if last_frame.status is not None:
            assert last_frame.status.normal_terminated is None
    if energy := expected.get("final_single_point_energy_hartree"):
        assert last_frame.energies is not None
        assert last_frame.energies.total_energy is not None
        assert last_frame.energies.total_energy.to("hartree").magnitude == pytest.approx(energy)
    if expected.get("has_forces"):
        assert any(frame.forces is not None for frame in file_model)
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
