import json
from pathlib import Path

import numpy as np
import pytest

from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
from molop.io.logic.gaussian.log.frame_parsers._g16_extractors import (
    ParseState,
    extract_populations_from_state,
)
from molop.io.logic.gaussian.log.frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.locators import (
    locate_g16_section_frames,
    locate_g16_sections,
)
from molop.io.logic.gaussian.log.models.G16LogFile import G16LogFileMemory
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import (
    extract_g16_atomic_masses,
)
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.unit import atom_ureg


FIXTURES = [
    Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log",
    Path(__file__).resolve().parent / "test_files" / "g16log" / "1-INT1-Sp.log",
]

ARCHIVE_ONLY_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
CCSD_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "CH3-ccsd-sp.log"
OPEN_SHELL_NPA_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "dsgdb9nsd_130472-1.log"
)
ATOMIC_MASS_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "dsgdb9nsd_000180-9-.log"
)


def _last_frame_block(fixture: Path) -> str:
    file_content = fixture.read_text()
    frame_blocks = [
        frame
        for section in locate_g16_sections(file_content)
        for frame in locate_g16_section_frames(file_content, section)
    ]
    return frame_blocks[-1].text(file_content)


def test_g16log_locator_lifecycle_retains_a_zero_frame_segment() -> None:
    source = "Gaussian 16:\n Error termination via Lnk1e.\n"

    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(source)

    assert len(parsed) == 0
    assert parsed.parse_completeness == "partial"
    assert len(parsed.source_segments) == 1
    segment = parsed.source_segments[0]
    assert segment.frame_count == 0
    assert segment.captured_frame_indices == []
    assert segment.parse_presence["geometry"] == "absent_in_source"
    assert segment.parse_completeness == "partial"


def test_g16log_zero_frame_segment_retains_portable_protocol_evidence() -> None:
    source = (
        FIXTURES[0]
        .read_text()
        .replace("Input orientation:", "Input geometry omitted:")
        .replace("Standard orientation:", "Standard geometry omitted:")
    )

    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(source)

    assert len(parsed) == 0
    assert all(segment.frame_count == 0 for segment in parsed.source_segments)
    segment = next(segment for segment in parsed.source_segments if segment.protocol is not None)
    assert segment.protocol is not None
    assert segment.protocol["method_family"]
    assert segment.task_requests
    dump = parsed.to_unitless_dump_with_unit_keys(exclude_none=True)
    dumped_segment = next(
        segment for segment in dump["source_segments"] if segment.get("protocol") is not None
    )
    assert isinstance(dumped_segment["protocol"], dict)
    assert isinstance(dumped_segment["task_requests"], list)
    json.dumps(dumped_segment)


def test_g16log_state_machine_frame_parser_matches_file_parser_frames():
    for fixture in FIXTURES:
        block = _last_frame_block(fixture)
        direct_frame = G16LogFileFrameParserMemory().parse(block)
        file_frame = G16LogFileParserMemory().parse(fixture.read_text())[-1]

        assert direct_frame.atoms == file_frame.atoms
        assert direct_frame.coords is not None
        assert file_frame.coords is not None
        assert direct_frame.standard_coords is not None or file_frame.standard_coords is None
        assert direct_frame.energies is not None
        assert file_frame.energies is not None
        assert direct_frame.energies.reference_energy is not None
        assert file_frame.energies.reference_energy is not None
        assert direct_frame.energies.reference_energy.to("hartree").m == pytest.approx(
            file_frame.energies.reference_energy.to("hartree").m
        )
        assert bool(direct_frame.vibrations) == bool(file_frame.vibrations)
        assert bool(direct_frame.molecular_orbitals) == bool(file_frame.molecular_orbitals)
        assert bool(direct_frame.charge_spin_populations) == bool(
            file_frame.charge_spin_populations
        )
        assert bool(direct_frame.polarizability) == bool(file_frame.polarizability)
        assert bool(direct_frame.hessian is not None) == bool(file_frame.hessian is not None)
        assert bool(direct_frame.forces is not None) == bool(file_frame.forces is not None)
        assert bool(direct_frame.geometry_optimization_status) == bool(
            file_frame.geometry_optimization_status
        )


def test_g16log_parses_atmwgt_masses_for_every_geometry_frame() -> None:
    parsed = G16LogFileParserMemory().parse(ATOMIC_MASS_FIXTURE.read_text())

    assert len(parsed.frames) == 18
    expected_masses = [
        12.0,
        12.0,
        12.0,
        12.0,
        14.003074,
        15.9949146,
        *([1.007825] * 8),
    ]
    for frame in parsed.frames:
        assert frame.atomic_masses_source == "gaussian_atmwgt"
        assert frame.atomic_masses is not None
        np.testing.assert_allclose(frame.atomic_masses.m_as("amu"), expected_masses)


def test_g16log_only_last_frame_keeps_masses_on_frame_only() -> None:
    parsed = G16LogFileParserMemory(only_last_frame=True).parse(ATOMIC_MASS_FIXTURE.read_text())

    assert len(parsed.frames) == 1
    frame = parsed.frames[0]
    assert frame.atomic_masses_source == "gaussian_atmwgt"
    assert frame.atomic_masses is not None
    assert frame.atomic_masses.shape == (14,)
    assert "atomic_masses" not in parsed.model_dump()


def test_g16log_backfills_thermochemistry_masses_to_matching_frames() -> None:
    parsed = G16LogFileParserMemory().parse(FIXTURES[0].read_text())

    expected_masses = [
        12.0,
        12.0,
        12.0,
        12.0,
        1.00783,
        12.0,
        1.00783,
        12.0,
        1.00783,
        1.00783,
        14.00307,
        12.0,
        1.00783,
        1.00783,
        1.00783,
        12.0,
        1.00783,
        1.00783,
        1.00783,
    ]
    assert parsed.frames
    for frame in parsed.frames:
        assert frame.atomic_masses_source == "gaussian_thermochemistry"
        assert frame.atomic_masses is not None
        np.testing.assert_allclose(frame.atomic_masses.m_as("amu"), expected_masses)
    assert "atomic_masses" not in parsed.model_dump()
    assert "atomic_masses_source" not in parsed.model_dump()


def test_extract_g16_atomic_masses_uses_latest_atmwgt_group_and_supports_d_notation() -> None:
    source = """
 Isotopes and Nuclear Properties:
 AtmWgt=  99.0000000
 Isotopes and Nuclear Properties:
 AtmWgt=  12.0000000  1.007825D+00
 - Thermochemistry -
 Atom     1 has atomic number  6 and mass  12.00000
 Atom     2 has atomic number  1 and mass   1.00783
 """

    masses, source_name = extract_g16_atomic_masses(source, expected_atom_count=2)

    assert source_name == "gaussian_atmwgt"
    assert masses is not None
    np.testing.assert_allclose(masses.m_as("amu"), [12.0, 1.007825])


def test_extract_g16_atomic_masses_skips_invalid_thermochemistry_group() -> None:
    source = """
 - Thermochemistry -
 Atom     1 has atomic number  6 and mass  12.00000
 Atom     2 has atomic number  1 and mass   1.00783
 - Thermochemistry -
 Atom     1 has atomic number  6 and mass  13.00000
 Atom     3 has atomic number  1 and mass   1.00783
 """

    masses, source_name = extract_g16_atomic_masses(source, expected_atom_count=2)

    assert source_name == "gaussian_thermochemistry"
    assert masses is not None
    np.testing.assert_allclose(masses.m_as("amu"), [12.0, 1.00783])

    assert extract_g16_atomic_masses("AtmWgt= 12.0000000", expected_atom_count=2) == (None, None)


def _g16_mass_frame(
    masses: list[float] | None,
    source: str | None,
) -> G16LogFileFrameMemory:
    return G16LogFileFrameMemory(
        atoms=[6, 1],
        coords=np.zeros((2, 3)) * atom_ureg.angstrom,
        atomic_masses=None if masses is None else np.asarray(masses) * atom_ureg.amu,
        atomic_masses_source=source,
    )


def test_g16log_mass_backfill_prefers_atmwgt_over_thermochemistry() -> None:
    chem_file = G16LogFileMemory()
    thermochemistry_frame = _g16_mass_frame([12.0, 1.00783], "gaussian_thermochemistry")
    missing_frame = _g16_mass_frame(None, None)
    atmwgt_frame = _g16_mass_frame([12.0, 1.007825], "gaussian_atmwgt")
    for frame in (thermochemistry_frame, missing_frame, atmwgt_frame):
        chem_file.append(frame)

    G16LogFileParserMemory()._update_file_metadata_from_frames(chem_file, {})

    assert thermochemistry_frame.atomic_masses_source == "gaussian_thermochemistry"
    assert missing_frame.atomic_masses_source == "gaussian_atmwgt"
    assert missing_frame.atomic_masses is not None
    np.testing.assert_allclose(missing_frame.atomic_masses.m_as("amu"), [12.0, 1.007825])


def test_g16log_mass_backfill_does_not_resolve_conflicting_equal_priority_masses() -> None:
    chem_file = G16LogFileMemory()
    first_frame = _g16_mass_frame([12.0, 1.007825], "gaussian_atmwgt")
    second_frame = _g16_mass_frame([13.0, 1.007825], "gaussian_atmwgt")
    missing_frame = _g16_mass_frame(None, None)
    for frame in (first_frame, second_frame, missing_frame):
        chem_file.append(frame)

    G16LogFileParserMemory()._update_file_metadata_from_frames(chem_file, {})

    assert missing_frame.atomic_masses is None
    assert missing_frame.atomic_masses_source is None


def test_g16log_open_shell_mulliken_and_npa_populations_are_both_retained() -> None:
    parsed = G16LogFileParserMemory().parse(OPEN_SHELL_NPA_FIXTURE.read_text())
    populations = parsed[-1].charge_spin_populations

    assert populations is not None
    assert len(populations["mulliken_charges"].values) == 11
    assert len(populations["mulliken_spins"].values) == 11
    assert len(populations["npa_charges"].values) == 11
    assert sum(populations["mulliken_charges"].values) == pytest.approx(0.0, abs=1.0e-5)
    assert sum(populations["mulliken_spins"].values) == pytest.approx(1.0, abs=1.0e-5)


def test_g16log_esp_population_uses_extensible_series() -> None:
    source = """
 Population analysis using the SCF Density.
 Mulliken charges:
     1  C   -0.100000
     2  H    0.100000
 Sum of Mulliken charges =   0.00000
 N-N= 1.0 E-N=-2.0 KE= 1.0
 ESP charges:
     1  C   -0.250000
     2  H    0.250000
 Sum of ESP charges =   0.00000
"""

    payload = extract_populations_from_state(ParseState(source))
    populations = payload["charge_spin_populations"]

    assert populations["populations"]["mulliken_charges"]["values"] == [-0.1, 0.1]
    assert populations["populations"]["esp_charges"]["values"] == [-0.25, 0.25]


def test_g16log_state_machine_result_is_canonical_model_data():
    block = _last_frame_block(FIXTURES[0])
    result = G16LogFileFrameParserMemory()._parse_block_to_result(block)
    model_data = result.model_data()

    assert isinstance(result, ModelParseResult)
    assert model_data["qm_software"] == "Gaussian"
    assert "component_tree" not in model_data
    assert "_component_tree" not in model_data
    assert "component_tree_input" not in model_data
    assert model_data["energies"]["reference_energy"] is not None
    assert model_data["vibrations"]["frequencies"] is not None
    assert model_data["thermal_informations"]["ZPVE"] is not None

    mutated = result.model_data()
    mutated["qm_software"] = "Mutated"
    assert result.model_data()["qm_software"] == "Gaussian"

    frame = G16LogFileFrameMemory.model_validate({"frame_content": block, **model_data})
    assert frame.energies is not None
    assert frame.vibrations is not None
    assert frame.component_tree is not None


def test_g16log_state_machine_uses_archive_energies_without_inventing_live_status():
    frame = G16LogFileFrameParserMemory().parse(_last_frame_block(ARCHIVE_ONLY_FIXTURE))

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert frame.energies.reference_energy.to("hartree").m == pytest.approx(-1.1330046)
    assert frame.energies.ccsd_energy is not None
    assert frame.energies.ccsd_energy.to("hartree").m == pytest.approx(-1.1726356)
    assert frame.status is None
    assert frame.temperature is None
    assert frame.pressure is None


def test_g16log_state_machine_preserves_live_reference_energy_over_archive_value():
    frame = G16LogFileFrameParserMemory().parse(_last_frame_block(CCSD_FIXTURE))
    parsed_frame = G16LogFileParserMemory().parse(CCSD_FIXTURE.read_text())[-1]

    assert frame.energies is not None
    assert parsed_frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert parsed_frame.energies.reference_energy is not None
    assert frame.energies.reference_energy.to("hartree").m == pytest.approx(
        parsed_frame.energies.reference_energy.to("hartree").m
    )


def test_g16log_state_machine_parses_external_energy_and_status() -> None:
    block = _last_frame_block(FIXTURES[0])
    matched = g16_log_patterns.SCF_ENERGY_AND_FUNCTIONAL.search(block)
    assert matched is not None
    block = (
        block[: matched.start()]
        + " External calculation of energy and first derivatives.\n"
        + ' Running external command "calculator input.json R"\n'
        + " Energy=    -123.456789     NIter=   0.\n"
        + block[matched.end() :]
    )

    frame = G16LogFileFrameParserMemory(capture_source_evidence=True).parse(block)

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert frame.energies.reference_energy.to("hartree").m == pytest.approx(-123.456789)
    assert frame.status is not None
    assert frame.status.scf_converged is True
    assert frame.is_error is False
    assert any(
        observation.source_label == "Gaussian External Energy"
        for observation in frame.energies.observations
    )


def test_g16log_external_archive_frame_uses_archive_energy_as_success_evidence() -> None:
    block = _last_frame_block(ARCHIVE_ONLY_FIXTURE)
    archive_start = block.find("1\\1\\GINC")
    assert archive_start >= 0
    functional_start = block.find("\\RCCSD-FC\\", archive_start)
    assert functional_start >= 0
    block = (
        block[:functional_start]
        + "\\RExternal='calculator input.json'\\"
        + block[functional_start + len("\\RCCSD-FC\\") :]
    )

    frame = G16LogFileFrameParserMemory().parse(block)

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert frame.status is not None
    assert frame.status.scf_converged is True
    assert frame.is_error is False


def test_g16log_component_tree_is_model_derived_view():
    frame = G16LogFileFrameParserMemory().parse(_last_frame_block(FIXTURES[0]))
    tree = frame.component_tree

    assert tree is not None
    assert tree.source_text == ""
    assert tree.validate_contracts() == []

    component_names = set(tree.component_names(include_synthetic=True))
    assert "l502.cycle" in component_names
    assert "l601.popanal" in component_names
    assert "l716.freq" in component_names
    assert "l716.thermochemistry" in component_names
    assert "l716.thermochemistry.mass" in component_names
    assert "l9999.archive" in component_names
    assert all(
        component.raw_text == "" for component in tree.iter_components(include_synthetic=True)
    )


def test_g16log_file_parser_frames_expose_model_derived_component_tree():
    parsed = G16LogFileParserMemory().parse(FIXTURES[0].read_text())

    assert parsed[0].component_tree is not None
    assert "l1.header" in parsed[0].component_tree.component_names()
    assert parsed[0].component_tree.source_text == ""


def test_g16log_frame_model_no_longer_accepts_component_tree_input_field():
    assert "component_tree_input" not in G16LogFileFrameMemory.model_fields


def test_g16log_state_machine_rejects_unexpected_phase(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = G16LogFileFrameParserMemory()

    monkeypatch.setattr(parser, "_run_header_phase", lambda _block, _result: object())

    with pytest.raises(AssertionError, match="Unexpected G16 frame parse phase"):
        parser._parse_block_to_result("")


def test_g16log_segment_metadata_state_machine_rejects_unexpected_phase(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parser = G16LogFileParserMemory()

    monkeypatch.setattr(parser, "_run_route_metadata_phase", lambda _context, _result: object())

    with pytest.raises(AssertionError, match="Unexpected G16 metadata parse phase"):
        parser._parse_segment_metadata_result("")
