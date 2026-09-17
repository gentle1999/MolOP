from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from molop.io.base_models.DataClasses import Comment, CommentContainer
from molop.io.logic.coords.frame_parsers.SDFFileFrameParser import SDFFileFrameParserMemory
from molop.io.logic.coords.frame_parsers.XYZFileFrameParser import XYZFileFrameParserMemory
from molop.io.logic.coords.models.XYZFile import XYZFileMemory
from molop.io.logic.gaussian.input.frame_parsers.GJFFileFrameParser import (
    GJFFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.input.frame_parsers.ORCAInpFileFrameParser import (
    ORCAInpFileFrameParserMemory,
)


def test_comment_container_is_ordered_and_list_like() -> None:
    comments = CommentContainer()
    comments.add("first", marker="#")
    comments.ensure(Comment(text="second", kind="title"))

    assert len(comments) == 2
    assert comments.texts == ["first", "second"]
    assert comments[0].marker == "#"
    title_comment = comments.first(kind="title")
    assert title_comment is not None
    assert title_comment.text == "second"


def test_file_and_frame_comments_have_one_access_pattern() -> None:
    parsed_frame = XYZFileFrameParserMemory().parse("2\nframe note\nH 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    parsed_file = XYZFileMemory.model_validate(
        {"comments": {"items": [{"text": "file note", "source_format": "xyz"}]}}
    )
    parsed_file.append(parsed_frame)

    assert parsed_file.comments.texts == ["file note"]
    assert parsed_file[0].comments.texts == ["frame note"]
    assert [comment.text for comment in parsed_file.iter_comments(scope="file")] == ["file note"]
    assert [comment.text for comment in parsed_file.iter_comments(scope="frame")] == ["frame note"]
    assert [comment.text for comment in parsed_file.iter_comments()] == [
        "file note",
        "frame note",
    ]


def test_format_specific_comment_fields_project_to_common_container() -> None:
    orca = ORCAInpFileFrameParserMemory().parse("! HF\n# ORCA note\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n")
    gaussian = GJFFileFrameParserMemory().parse(
        "#p hf/sto-3g\n\ngjf title\n\n0 1\nH 0.0 0.0 0.0\n\n"
    )

    assert orca.comments.texts == ["ORCA note"]
    assert orca.comments[0].marker == "#"
    assert gaussian.comments.texts == ["gjf title"]
    assert gaussian.comments[0].kind == "title"


def test_g16log_title_card_projects_to_file_and_frame_comments() -> None:
    source = Path("tests/test_files/g16log/H2O.log").read_text(encoding="utf-8")
    parsed = G16LogFileParserMemory().parse(source)

    expected = "opt and freq in PhCl/MeOH=20/1"
    assert parsed.title_card == expected
    assert parsed.comments.texts == [expected]
    assert parsed.comments[0].kind == "title"
    assert all(frame.comments.texts == [expected] for frame in parsed.frames)


def test_sdf_title_is_available_as_a_frame_comment() -> None:
    molecule = Chem.MolFromSmiles("CO")
    assert molecule is not None
    molecule.SetProp("_Name", "methanol")
    frame = SDFFileFrameParserMemory().parse(Chem.MolToMolBlock(molecule))

    assert frame.comments.texts == ["methanol"]
    assert frame.comments[0].kind == "title"
