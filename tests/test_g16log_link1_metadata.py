from pathlib import Path

import pytest

from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.logic.gaussian.log.locators import locate_g16_sections
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import (
    extract_g16_termination_status,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import (
    G16LogFileParserMemory,
    G16LogFileParserMixin,
)


CCSD_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
DFT_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"
GEN_INHERITANCE_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "dsgdb9nsd_000696-4.log"
)
MINIMUM_FREQ_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"
)
TS_FREQ_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "5-TS2-Opt.log"


def test_g16log_file_parser_uses_shared_parse_lifecycle() -> None:
    assert "_parse" not in G16LogFileParserMixin.__dict__


def test_g16log_artifact_metadata_does_not_run_segment_parser(monkeypatch) -> None:
    def fail_segment_parse(*_args, **_kwargs):
        raise AssertionError("artifact metadata must not parse segment metadata")

    monkeypatch.setattr(
        G16LogFileParserMixin,
        "_parse_segment_metadata_result",
        fail_segment_parse,
    )

    metadata = G16LogFileParserMemory()._parse_artifact_metadata(DFT_FIXTURE.read_text())

    assert metadata == {
        "qm_software": "Gaussian",
        "qm_software_version": "ES64L-G16RevA.03",
    }


def test_g16log_shared_lifecycle_preserves_gen_basis_inheritance() -> None:
    parsed = G16LogFileParserMemory().parse(GEN_INHERITANCE_FIXTURE.read_text())

    assert parsed[-2].basis_set == "6-311+G(d,p)"
    assert parsed[-1].basis_set == parsed[-2].basis_set


def _split_raw_gaussian_blocks(file_content: str) -> tuple[str, list[str]]:
    split_token = "Input orientation:"
    if split_token not in file_content:
        split_token = "Standard orientation:"
    parts = file_content.split(split_token)
    header = parts[0]
    frames = [f"{split_token}{fragment}" for fragment in parts[1:]]
    return header, frames


def test_g16log_link1_sections_propagate_section_metadata_to_later_frames() -> None:
    ccsd_header, ccsd_frames = _split_raw_gaussian_blocks(CCSD_FIXTURE.read_text())
    dft_header, dft_frames = _split_raw_gaussian_blocks(DFT_FIXTURE.read_text())

    combined = "\n".join(
        [
            ccsd_header + ccsd_frames[0],
            "Link1:  Proceeding to internal job step number  2.\n"
            + dft_header
            + dft_frames[0]
            + dft_frames[1],
        ]
    )

    parsed = G16LogFileParserMemory().parse(combined)

    assert len(parsed) == 3

    assert "ccsd/aug-cc-pvtz" in parsed[0].keywords.lower()
    assert parsed[0].title_card.strip() == "Title Card Required"

    for frame in parsed[1:]:
        assert "b3lyp/6-31g(d)" in frame.keywords.lower()
        assert "ccsd/aug-cc-pvtz" not in frame.keywords.lower()
        assert frame.title_card.strip() == "opt and freq in PhCl/MeOH=20/1"


def test_g16log_link1_segments_export_distinct_protocol_evidence() -> None:
    ccsd_header, ccsd_frames = _split_raw_gaussian_blocks(CCSD_FIXTURE.read_text())
    dft_header, dft_frames = _split_raw_gaussian_blocks(DFT_FIXTURE.read_text())
    combined = "\n".join(
        [
            ccsd_header + ccsd_frames[0],
            "Link1:  Proceeding to internal job step number  2.\n"
            + dft_header
            + dft_frames[0]
            + dft_frames[1],
        ]
    )

    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(combined)

    assert [frame.file_frame_index for frame in parsed.frames] == [0, 1, 2]
    first, second = [segment for segment in parsed.source_segments if segment.protocol is not None]
    assert first.protocol is not None
    assert first.protocol["method_family"] == "CCSD"
    assert second.protocol is not None
    assert second.protocol["method_family"] == "DFT"
    assert [request["task_type"] for request in second.task_requests] == ["opt", "freq"]


def test_g16log_segment_metadata_result_is_canonical_model_data() -> None:
    parser = G16LogFileParserMemory()
    result = parser._parse_segment_metadata_result(DFT_FIXTURE.read_text())
    metadata = result.model_data()

    assert isinstance(result, ModelParseResult)
    assert "_file_content" not in parser.__dict__
    assert metadata["qm_software"] == "Gaussian"
    assert "component_tree" not in metadata
    assert "_component_tree" not in metadata
    assert "frame_content" not in metadata
    assert "file_content" not in metadata
    assert "b3lyp/6-31g(d)" in metadata["keywords"].lower()
    assert metadata["title_card"].strip() == "opt and freq in PhCl/MeOH=20/1"
    assert metadata["charge"] == 0
    assert metadata["multiplicity"] == 1

    mutated = result.model_data()
    mutated["qm_software"] = "Mutated"
    assert result.model_data()["qm_software"] == "Gaussian"

    file_model = parser._chem_file.model_validate(metadata)
    assert file_model.qm_software == "Gaussian"
    assert file_model.method == "DFT"


@pytest.mark.parametrize("fixture", [MINIMUM_FREQ_FIXTURE, TS_FREQ_FIXTURE])
def test_g16log_generated_freq_link1_scans_full_segment_for_scf_status(
    fixture: Path,
) -> None:
    source = fixture.read_text()
    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(source)
    segment = next(segment for segment in parsed.source_segments if segment.task_types == ["freq"])
    segment_content = locate_g16_sections(source)[segment.segment_index].text(source)

    assert "Geom=AllCheck" in segment_content
    assert "SCF Done" in segment_content
    status = (
        G16LogFileParserMemory()
        ._parse_segment_metadata_result(segment_content)
        .model_data()["status"]
    )
    assert status.scf_converged is True
    assert segment.parse_presence["scf_status"] == "parsed"
    assert segment.scf_status is True


def test_g16log_generated_freq_section_inherits_previous_calculation_configuration() -> None:
    parsed = G16LogFileParserMemory(capture_source_evidence=True).parse(DFT_FIXTURE.read_text())

    previous_frame = parsed[0]
    frequency_frame = parsed[-1]
    assert [task.task_type for task in frequency_frame.task_requests] == ["freq"]
    assert frequency_frame.keywords != previous_frame.keywords
    assert frequency_frame.functional == "B3LYP"
    assert frequency_frame.model_chemistry.dispersion_correction == "GD3BJ"
    assert frequency_frame.model_chemistry.solvation_model == "smd"
    assert frequency_frame.model_chemistry.solvent == "generic"
    assert frequency_frame.solvent == previous_frame.solvent

    frequency_protocol = parsed.source_segments[-1].protocol
    assert frequency_protocol is not None
    assert frequency_protocol["dispersion_correction"] == "GD3BJ"
    assert frequency_protocol["solvation_model"] == "smd"


def test_g16log_generated_freq_section_inherits_configuration_when_only_last_frame_is_selected() -> (
    None
):
    parsed = G16LogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(DFT_FIXTURE.read_text())

    assert len(parsed) == 1
    assert parsed[0].task_requests[0].task_type == "freq"
    assert parsed[0].model_chemistry.dispersion_correction == "GD3BJ"
    assert parsed[0].functional == "B3LYP"


def test_g16log_entering_link1_sections_propagate_section_metadata_to_later_frames() -> None:
    ccsd_header, ccsd_frames = _split_raw_gaussian_blocks(CCSD_FIXTURE.read_text())
    dft_header, dft_frames = _split_raw_gaussian_blocks(DFT_FIXTURE.read_text())

    combined = "\n".join(
        [
            ccsd_header + ccsd_frames[0],
            dft_header + dft_frames[0] + dft_frames[1],
        ]
    )

    parsed = G16LogFileParserMemory().parse(combined)

    assert len(parsed) == 3
    assert "ccsd/aug-cc-pvtz" in parsed[0].keywords.lower()

    for frame in parsed[1:]:
        assert "b3lyp/6-31g(d)" in frame.keywords.lower()
        assert "ccsd/aug-cc-pvtz" not in frame.keywords.lower()
        assert frame.title_card.strip() == "opt and freq in PhCl/MeOH=20/1"


def test_gaussian_normal_termination_does_not_fabricate_scf_convergence() -> None:
    status = extract_g16_termination_status(
        TextParseContext(" Normal termination of Gaussian 16 at Mon Jan  1 00:00:00 2024.\n")
    )

    assert status is not None
    assert status.normal_terminated is True
    assert status.scf_converged is None


def test_gaussian_external_energy_is_successful_electronic_calculation_evidence() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " External calculation of energy and first derivatives.\n"
            ' Running external command "calculator input.json R"\n'
            '         output file      "/tmp/Gau-1.EOu"\n'
            " Energy=    -858.897871     NIter=   0.\n"
            " Normal termination of Gaussian 16 at Mon Jan  1 00:00:00 2024.\n"
        )
    )

    assert status is not None
    assert status.normal_terminated is True
    assert status.scf_converged is True


def test_gaussian_isolated_energy_line_does_not_fabricate_scf_convergence() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " Energy=    -858.897871     NIter=   0.\n"
            " Normal termination of Gaussian 16 at Mon Jan  1 00:00:00 2024.\n"
        )
    )

    assert status is not None
    assert status.normal_terminated is True
    assert status.scf_converged is None


def test_gaussian_external_archive_energy_is_success_evidence() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " 1\\1\\GINC-HOST\\FOpt\\RExternal='calculator input.json'\\ZDO\\Molecule\\User"
            "\\01-Jan-2024\\0\\\\# opt external('calculator input.json')\\\\Title\\\\0,1"
            "\\H,0.,0.,0.\\\\Version=ES64L-G16RevA.03\\HF=-1.23456789\\RMSD=0.000e+00\\@\n"
        )
    )

    assert status is not None
    assert status.normal_terminated is None
    assert status.scf_converged is True


def test_gaussian_non_external_archive_energy_does_not_fabricate_scf_convergence() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " 1\\1\\GINC-HOST\\FOpt\\RB3LYP\\6-31G(d)\\Molecule\\User\\01-Jan-2024\\0"
            "\\\\# opt b3lyp/6-31g(d)\\\\Title\\\\0,1\\H,0.,0.,0.\\\\Version="
            "ES64L-G16RevA.03\\HF=-1.23456789\\RMSD=0.000e+00\\@\n"
        )
    )

    assert status is None


def test_gaussian_explicit_scf_failure_is_preserved_independently_of_termination() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " Convergence failure -- run terminated.\n"
            " Error termination via Lnk1e in /tmp/l502.exe at Mon Jan  1 00:00:00 2024.\n"
        )
    )

    assert status is not None
    assert status.normal_terminated is False
    assert status.scf_converged is False


def test_gaussian_later_scf_failure_overrides_external_energy_success() -> None:
    status = extract_g16_termination_status(
        TextParseContext(
            " External calculation of energy, first and second derivatives.\n"
            ' Running external command "calculator input.json R"\n'
            " Energy=    -858.897871     NIter=   0.\n"
            " SCF failed to converge\n"
            " Error termination via Lnk1e in /tmp/l502.exe at Mon Jan  1 00:00:00 2024.\n"
        )
    )

    assert status is not None
    assert status.normal_terminated is False
    assert status.scf_converged is False
