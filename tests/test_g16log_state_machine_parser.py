import json
from pathlib import Path

import pytest

from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
from molop.io.logic.gaussian.log.frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.locators import (
    locate_g16_section_frames,
    locate_g16_sections,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURES = [
    Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log",
    Path(__file__).resolve().parent / "test_files" / "g16log" / "1-INT1-Sp.log",
]

ARCHIVE_ONLY_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
CCSD_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "CH3-ccsd-sp.log"


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
