from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import pytest
from rdkit import Chem

from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.base_models.ChemFileFrame import BaseChemFileFrame
from molop.io.codec_registry import Registry
from molop.io.logic.coords.parsers.SDFFileParser import (
    SDFFileParserDisk,
    SDFFileParserMemory,
)
from molop.io.logic.coords.parsers.SDFFileParser import (
    register as register_sdf,
)
from molop.io.logic.coords.parsers.SMIFileParser import (
    SMIFileParserDisk,
    SMIFileParserMemory,
)
from molop.io.logic.coords.parsers.SMIFileParser import (
    register as register_smi,
)
from molop.io.logic.coords.parsers.XYZFileParser import (
    XYZFileParserDisk,
    XYZFileParserMemory,
)
from molop.io.logic.coords.parsers.XYZFileParser import (
    register as register_xyz,
)
from molop.io.logic.gaussian.input.parsers.GJFFileParser import (
    GJFFileParserDisk,
    GJFFileParserMemory,
)
from molop.io.logic.gaussian.input.parsers.GJFFileParser import (
    register as register_gjf,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import (
    G16LogFileParserDisk,
    G16LogFileParserMemory,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import (
    register as register_g16log,
)
from molop.io.logic.orca.input.parsers.ORCAInpFileParser import (
    ORCAInpFileParserDisk,
    ORCAInpFileParserMemory,
)
from molop.io.logic.orca.input.parsers.ORCAInpFileParser import (
    register as register_orcainp,
)
from molop.io.logic.orca.log.parsers.ORCALogFileParser import (
    ORCALogFileParserDisk,
    ORCALogFileParserMemory,
)
from molop.io.logic.orca.log.parsers.ORCALogFileParser import (
    register as register_orcaout,
)


FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files"


def _sdf_source() -> str:
    molecule = Chem.MolFromSmiles("C")
    assert molecule is not None
    return f"{Chem.MolToMolBlock(molecule)}\n$$$$\n"


FORMAT_CASES: tuple[
    tuple[
        type[Any],
        type[Any],
        Callable[[Registry], None],
        str,
        str,
        str,
    ],
    ...,
] = (
    (
        XYZFileParserMemory,
        XYZFileParserDisk,
        register_xyz,
        "xyz",
        ".xyz",
        "1\nprobe\nH 0.0 0.0 0.0\n",
    ),
    (
        SDFFileParserMemory,
        SDFFileParserDisk,
        register_sdf,
        "sdf",
        ".mol",
        _sdf_source(),
    ),
    (
        SMIFileParserMemory,
        SMIFileParserDisk,
        register_smi,
        "smi",
        ".txt",
        "C methane\n",
    ),
    (
        GJFFileParserMemory,
        GJFFileParserDisk,
        register_gjf,
        "gjf",
        ".com",
        "#p hf/3-21g\n\nprobe\n\n0 1\nH 0.0 0.0 0.0\n\n",
    ),
    (
        ORCAInpFileParserMemory,
        ORCAInpFileParserDisk,
        register_orcainp,
        "orcainp",
        ".inp",
        "! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n",
    ),
    (
        G16LogFileParserMemory,
        G16LogFileParserDisk,
        register_g16log,
        "g16log",
        ".out",
        (FIXTURE_ROOT / "g16log" / "H2O.log").read_text(),
    ),
    (
        ORCALogFileParserMemory,
        ORCALogFileParserDisk,
        register_orcaout,
        "orcaout",
        ".log",
        (FIXTURE_ROOT / "orca" / "opt_orca.out").read_text(),
    ),
)


@pytest.mark.parametrize(
    ("memory_parser", "_disk_parser", "_register", "expected", "_suffix", "source"),
    FORMAT_CASES,
    ids=[case[3] for case in FORMAT_CASES],
)
def test_memory_parsers_set_canonical_file_source_format(
    memory_parser: type[Any],
    _disk_parser: type[Any],
    _register: Callable[[Registry], None],
    expected: str,
    _suffix: str,
    source: str,
) -> None:
    parsed = memory_parser().parse(source)

    assert memory_parser.format_id == expected
    assert parsed.source_format == expected
    assert parsed.model_dump()["source_format"] == expected
    assert all("source_format" not in frame.model_dump() for frame in parsed.frames)


@pytest.mark.parametrize(
    ("_memory_parser", "disk_parser", "_register", "expected", "suffix", "source"),
    FORMAT_CASES,
    ids=[case[3] for case in FORMAT_CASES],
)
def test_disk_parsers_set_source_format_independently_of_file_extension(
    tmp_path: Path,
    _memory_parser: type[Any],
    disk_parser: type[Any],
    _register: Callable[[Registry], None],
    expected: str,
    suffix: str,
    source: str,
) -> None:
    path = tmp_path / f"source{suffix}"
    path.write_text(source)

    parsed = disk_parser().parse(str(path))

    assert disk_parser.format_id == expected
    assert parsed.source_format == expected
    assert Path(parsed.file_path).suffix == suffix
    assert all("source_format" not in frame.model_dump() for frame in parsed.frames)


@pytest.mark.parametrize(
    ("_memory_parser", "disk_parser", "register", "expected", "suffix", "_source"),
    FORMAT_CASES,
    ids=[case[3] for case in FORMAT_CASES],
)
def test_builtin_register_uses_the_parser_format_id_classvar(
    monkeypatch: pytest.MonkeyPatch,
    _memory_parser: type[Any],
    disk_parser: type[Any],
    register: Callable[[Registry], None],
    expected: str,
    suffix: str,
    _source: str,
) -> None:
    sentinel = f"test-{expected}"
    monkeypatch.setattr(disk_parser, "format_id", sentinel)
    registry = Registry(autoload_defaults=False)

    register(registry)
    reader = registry.select_reader(f"source{suffix}", hint_format=sentinel)[0]
    concrete_reader = cast(Any, reader)._factory()

    assert reader.format_id == sentinel
    assert concrete_reader.format_id == sentinel
    assert concrete_reader.parser_cls is disk_parser
    assert concrete_reader.parser_cls.format_id == sentinel


def test_source_format_is_optional_for_direct_file_models_and_absent_from_frames() -> None:
    assert BaseChemFile().source_format is None
    assert "source_format" not in BaseChemFileFrame.model_fields


def test_parser_lifecycle_rejects_a_missing_source_format_classvar() -> None:
    class MissingFormatXYZParser(XYZFileParserMemory):
        format_id = None

    with pytest.raises(ValueError, match="format_id must be a non-empty normalized identifier"):
        MissingFormatXYZParser().parse("1\nprobe\nH 0.0 0.0 0.0\n")
