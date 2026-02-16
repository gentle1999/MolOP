import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read(rel_path: str) -> str:
    return (ROOT / rel_path).read_text(encoding="utf-8")


def _extract_overload_blocks(text: str, function_name: str) -> list[str]:
    pattern = re.compile(
        rf"@overload\n"
        rf"(\s*def\s+{function_name}\([\s\S]*?\)\s*->\s*[\s\S]*?:\n"
        rf"\s*\"\"\"[\s\S]*?\"\"\"\n"
        rf"\s*\.\.\.)",
        re.MULTILINE,
    )
    return pattern.findall(text)


def _extract_docstring(block: str) -> str:
    match = re.search(r'"""([\s\S]*?)"""', block)
    assert match is not None
    return match.group(1)


def _extract_args_names(docstring: str) -> list[str]:
    section = re.search(
        r"\n\s*Args:\n(?P<body>[\s\S]*?)(?:\n\s*(?:Format-specific Args|Returns|Source):|\Z)",
        docstring,
    )
    if section is None:
        return []
    return re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*):", section.group("body"), re.MULTILINE)


def _assert_no_duplicate_doc_args(block: str) -> None:
    names = _extract_args_names(_extract_docstring(block))
    assert len(names) == len(set(names)), f"duplicate Args entries: {names}"


def _assert_sorted_literal_overloads(blocks: list[str], arg_name: str) -> None:
    format_ids: list[str] = []
    for block in blocks:
        match = re.search(rf"{arg_name}:\s*Literal\[\"([^\"]+)\"\]", block)
        if match:
            format_ids.append(match.group(1))
    assert format_ids == sorted(format_ids)


def test_io_stub_overloads_are_sorted_and_have_google_args() -> None:
    for rel_path in (
        "src/molop/io/base_models/_format_transform.pyi",
        "src/molop/io/_batch_format_transform.pyi",
    ):
        blocks = _extract_overload_blocks(_read(rel_path), "format_transform")
        assert blocks
        _assert_sorted_literal_overloads(blocks, "format")

        literal_blocks = [b for b in blocks if 'format: Literal["' in b]
        assert literal_blocks
        for block in literal_blocks:
            doc = _extract_docstring(block)
            assert "Args:" in doc
            _assert_no_duplicate_doc_args(block)


def test_cli_stub_transform_overloads_include_literal_to_args_and_source() -> None:
    blocks = _extract_overload_blocks(_read("src/molop/cli/app.pyi"), "transform")
    assert blocks
    _assert_sorted_literal_overloads(blocks, "to")

    literal_blocks = [b for b in blocks if 'to: Literal["' in b]
    assert literal_blocks
    for block in literal_blocks:
        assert 'to: Literal["' in block
        doc = _extract_docstring(block)
        assert "Args:" in doc
        assert "Source:" in doc
        _assert_no_duplicate_doc_args(block)


def test_stub_generators_check_mode_is_deterministic() -> None:
    commands = [
        [
            "uv",
            "run",
            "python",
            "scripts/generate_chemfile_format_transform_stubs.py",
            "--check",
        ],
        [
            "uv",
            "run",
            "python",
            "scripts/generate_cli_transform_stubs.py",
            "--check",
        ],
    ]
    for command in commands:
        subprocess.run(command, cwd=ROOT, check=True)
