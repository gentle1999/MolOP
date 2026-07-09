from pathlib import Path

from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


CCSD_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
DFT_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"


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


def test_g16log_file_metadata_result_is_canonical_model_data() -> None:
    parser = G16LogFileParserMemory()
    result = parser._parse_metadata_result(DFT_FIXTURE.read_text())
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
