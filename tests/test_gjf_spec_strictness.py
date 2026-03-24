import pytest

from molop.io.logic.qminput_frame_parsers.GJFFileFrameParser import GJFFileFrameParserMemory
from molop.io.logic.qminput_parsers.GJFFileParser import GJFFileParserMemory


def test_gjf_title_card_is_required() -> None:
    parser = GJFFileFrameParserMemory()
    block = """%chk=t.chk
#p hf/3-21g


0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
"""
    with pytest.raises(ValueError, match="title card"):
        parser.parse(block)


def test_gjf_title_card_cannot_exceed_five_lines() -> None:
    parser = GJFFileFrameParserMemory()
    block = """%chk=t.chk
#p hf/3-21g

line1
line2
line3
line4
line5
line6

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
"""
    with pytest.raises(ValueError, match="cannot exceed 5 lines"):
        parser.parse(block)


def test_link1_requires_blank_line_before_separator() -> None:
    parser = GJFFileParserMemory()
    content = """%chk=a.chk
#p hf/3-21g

title

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
--Link1--
%chk=b.chk
#p hf/3-21g

title2

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
"""
    with pytest.raises(ValueError, match="preceded by a blank line"):
        parser.parse(content)


def test_link1_with_blank_line_is_accepted() -> None:
    parser = GJFFileParserMemory()
    content = """%chk=a.chk
#p hf/3-21g

title

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7

--Link1--
%chk=b.chk
#p hf/3-21g

title2

0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
"""
    parsed = parser.parse(content)
    assert len(parsed) == 2


def test_geom_allcheck_allows_missing_title_and_molecule_sections() -> None:
    parser = GJFFileFrameParserMemory()
    block = """%chk=t.chk
#p hf/3-21g geom=allcheck

FreezeAll
"""
    frame = parser.parse(block)
    assert frame.title_card.title_card == ""
    assert frame.molecule_specifications.molecule_fragments == []
    assert "FreezeAll" in frame.additional_sections


def test_geom_checkpoint_allows_missing_molecule_with_title() -> None:
    parser = GJFFileFrameParserMemory()
    block = """%chk=t.chk
#p hf/3-21g geom=checkpoint

checkpoint title

FreezeAll
"""
    frame = parser.parse(block)
    assert frame.title_card.title_card.strip() == "checkpoint title"
    assert frame.molecule_specifications.molecule_fragments == []
    assert "FreezeAll" in frame.additional_sections
