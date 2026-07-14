# pyright: reportPrivateUsage=false

from rdkit import Chem

from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.coords.frame_parsers.SDFFileFrameParser import SDFFileFrameParserMemory
from molop.io.logic.coords.frame_parsers.SMIFileFrameParser import SMIFileFrameParserMemory
from molop.io.logic.coords.frame_parsers.XYZFileFrameParser import XYZFileFrameParserMemory
from molop.io.logic.coords.parsers.SDFFileParser import SDFFileParserMemory
from molop.io.logic.coords.parsers.SMIFileParser import SMIFileParserMemory
from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserMemory


def test_xyz_frame_parser_builds_model_parse_result() -> None:
    block = "\n".join(
        [
            "3",
            "charge -1 multiplicity 2",
            "O 0.000000 0.000000 0.000000",
            "H 0.957200 0.000000 0.000000",
            "H -0.239987 0.927297 0.000000",
        ]
    )

    result = XYZFileFrameParserMemory()._parse_block_to_result(block)
    data = result.model_data()

    assert isinstance(result, ModelParseResult)
    assert data["comment"] == "charge -1 multiplicity 2"
    assert data["atoms"] == [8, 1, 1]
    assert data["charge"] == -1
    assert data["multiplicity"] == 2
    assert data["coords"].m.shape == (3, 3)


def test_sdf_frame_parser_builds_model_parse_result() -> None:
    rdmol = Chem.MolFromSmiles("CCO")
    block = Chem.MolToMolBlock(rdmol)

    result = SDFFileFrameParserMemory()._parse_block_to_result(block)
    data = result.model_data()

    assert isinstance(result, ModelParseResult)
    assert data["atoms"] == [6, 6, 8]
    assert data["charge"] == 0
    assert data["multiplicity"] == 1
    assert data["coords"].m.shape == (3, 3)
    assert len(data["bonds"]) == 2


def test_smi_frame_parser_builds_model_parse_result_from_first_token() -> None:
    result = SMIFileFrameParserMemory()._parse_block_to_result("CCO ethanol")
    data = result.model_data()

    assert isinstance(result, ModelParseResult)
    assert data["atoms"] == [6, 6, 8]
    assert data["charge"] == 0
    assert data["multiplicity"] == 1
    assert data["coords"].m.shape == (3, 3)
    assert len(data["bonds"]) == 2


def test_coords_file_parsers_need_no_empty_metadata_hook() -> None:
    parsers = [
        XYZFileParserMemory(),
        SDFFileParserMemory(),
        SMIFileParserMemory(),
    ]

    for parser in parsers:
        assert parser._parse_artifact_metadata("") is None
        assert not hasattr(parser, "_parse_metadata")
        assert not hasattr(parser, "_split_file")
