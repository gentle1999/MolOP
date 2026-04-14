from pathlib import Path

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
