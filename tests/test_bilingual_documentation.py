"""Check language-independent facts across all paired documentation.

Prose still needs editorial review. These checks retain translatable messages
while detecting divergent APIs, defaults, calculation inputs and rendered scope.
"""

from __future__ import annotations

import ast
import json
import re
import shlex
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from textwrap import dedent


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_LEVELS = {
    "yes": "supported",
    "是": "supported",
    "supported": "supported",
    "支持": "supported",
    "已支持": "supported",
    "no": "unsupported",
    "否": "unsupported",
    "unsupported": "unsupported",
    "不支持": "unsupported",
    "partial": "partial",
    "部分": "partial",
    "部分支持": "partial",
    "unverified": "unverified",
    "未验证": "unverified",
    "not applicable": "not applicable",
    "不适用": "not applicable",
    "raw": "raw",
    "raw/unsupported": "raw/unsupported",
    "raw/不支持": "raw/unsupported",
    "not guaranteed": "not guaranteed",
    "不保证": "not guaranteed",
    "indirect": "indirect",
    "间接": "indirect",
}


def _notebook_pairs() -> Iterator[tuple[str, dict, dict]]:
    english_root = ROOT / "docs/en"
    chinese_root = ROOT / "docs/zh"
    english = {path.relative_to(english_root) for path in english_root.rglob("*.ipynb")}
    chinese = {path.relative_to(chinese_root) for path in chinese_root.rglob("*.ipynb")}
    assert english == chinese
    for relative in sorted(english):
        yield (
            str(relative),
            json.loads((english_root / relative).read_text(encoding="utf-8")),
            json.loads((chinese_root / relative).read_text(encoding="utf-8")),
        )


def _markdown_pairs() -> Iterator[tuple[str, str, str]]:
    for english in sorted((ROOT / "docs/en").rglob("*.md")):
        relative = english.relative_to(ROOT / "docs/en")
        chinese = ROOT / "docs/zh" / relative
        yield (
            str(relative),
            english.read_text(encoding="utf-8"),
            chinese.read_text(encoding="utf-8"),
        )
    yield (
        "README",
        (ROOT / "README.md").read_text(encoding="utf-8"),
        (ROOT / "README.zh.md").read_text(encoding="utf-8"),
    )
    for name, english, chinese in _notebook_pairs():
        assert [cell["cell_type"] for cell in english["cells"]] == [
            cell["cell_type"] for cell in chinese["cells"]
        ], name
        for index, (first, second) in enumerate(
            zip(english["cells"], chinese["cells"], strict=True)
        ):
            if first["cell_type"] == "markdown":
                yield f"{name}#cell-{index}", "".join(first["source"]), "".join(second["source"])


def _markdown_parts(text: str) -> Iterator[tuple[str | None, str]]:
    """Recognize fence lengths so Markdown examples can contain shorter fences."""
    fence: str | None = None
    language: str | None = None
    content: list[str] = []
    for line in text.splitlines():
        if fence is None:
            match = re.match(r"^\s*(`{3,}|~{3,})([^\s]*)\s*$", line)
            if match:
                yield None, "\n".join(content)
                fence, language = match.groups()
                content = []
            else:
                content.append(line)
        elif re.fullmatch(rf"\s*{re.escape(fence[0])}{{{len(fence)},}}\s*", line):
            yield language, dedent("\n".join(content))
            fence = None
            language = None
            content = []
        else:
            content.append(line)
    assert fence is None, "Unclosed Markdown fence"
    yield None, "\n".join(content)


def _fenced_blocks(text: str) -> Iterator[tuple[str, str]]:
    for language, content in _markdown_parts(text):
        if language is not None:
            yield language, content


def _without_fences(text: str) -> str:
    # Keep only prose for links and generated-reference directives, not examples of Markdown.
    return "\n".join(content for language, content in _markdown_parts(text) if language is None)


def _link_targets(text: str) -> Counter[str]:
    text = _without_fences(text)
    targets = re.findall(r"\]\(([^\s)]+)(?:\s+\"[^\"\n]*\")?\)", text)
    targets.extend(re.findall(r"^\[[^]]+\]:\s*(\S+)", text, re.M))
    # README links select different locales of the same deployed page.
    targets = [
        target.replace(
            "https://gentle1999.github.io/MolOP/en/", "https://gentle1999.github.io/MolOP/"
        )
        for target in targets
    ]
    return Counter(targets)


def _api_directives(text: str) -> list[str]:
    lines = _without_fences(text).splitlines()
    blocks: list[str] = []
    for index, line in enumerate(lines):
        if not line.startswith("::: "):
            continue
        block = [line]
        for following in lines[index + 1 :]:
            if following and not following.startswith((" ", "\t")):
                break
            block.append(following)
        blocks.append("\n".join(block).rstrip())
    return blocks


def _tables(text: str) -> list[list[list[str]]]:
    tables: list[list[list[str]]] = []
    current: list[list[str]] = []
    for line in _without_fences(text).splitlines():
        if line.strip().startswith("|"):
            current.append([cell.strip() for cell in re.split(r"(?<!\\)\|", line)[1:-1]])
        elif current:
            tables.append(current)
            current = []
    if current:
        tables.append(current)
    return tables


class _LocalizedMessages(ast.NodeTransformer):
    """Localize display/error prose; preserve computation and data-valued strings."""

    def visit_Expr(self, node: ast.Expr) -> ast.AST | None:
        if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            return None  # Docstrings have no effect on the demonstrated calculation.
        return self.generic_visit(node)

    @staticmethod
    def _message(node: ast.expr) -> ast.expr:
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return ast.Constant(value="<localized message>")
        if isinstance(node, ast.JoinedStr):
            node.values = [
                ast.Constant(value="<localized message>")
                if isinstance(part, ast.Constant) and isinstance(part.value, str)
                else part
                for part in node.values
            ]
        return node

    def visit_Call(self, node: ast.Call) -> ast.Call:
        self.generic_visit(node)
        if isinstance(node.func, ast.Name) and (
            node.func.id == "print" or node.func.id.endswith(("Error", "Exception"))
        ):
            node.args = [self._message(argument) for argument in node.args]
        for keyword in node.keywords:
            if keyword.arg == "missing_reason":
                keyword.value = self._message(keyword.value)
        return node


def _python_signature(code: str) -> str:
    return ast.dump(_LocalizedMessages().visit(ast.parse(dedent(code))))


def test_bilingual_links_and_generated_api_scope_match() -> None:
    for name, english, chinese in _markdown_pairs():
        assert _link_targets(english) == _link_targets(chinese), name
        assert _api_directives(english) == _api_directives(chinese), name
        reference_links = r"\]\[([^]\n]+)\]"
        assert Counter(re.findall(reference_links, english)) == Counter(
            re.findall(reference_links, chinese)
        ), name


def test_bilingual_tables_keep_rows_identifiers_and_support_levels() -> None:
    for name, english, chinese in _markdown_pairs():
        first_tables, second_tables = _tables(english), _tables(chinese)
        assert len(first_tables) == len(second_tables), name
        for table_index, (first_table, second_table) in enumerate(
            zip(first_tables, second_tables, strict=True)
        ):
            assert [len(row) for row in first_table] == [len(row) for row in second_table], name
            default_columns = {
                index
                for index, heading in enumerate(first_table[0])
                if heading.lower() == "default"
            }
            assert default_columns == {
                index for index, heading in enumerate(second_table[0]) if heading == "默认值"
            }, (name, table_index)
            for row_index, (first, second) in enumerate(
                zip(first_table, second_table, strict=True)
            ):
                location = (name, table_index, row_index)
                assert re.findall(r"`([^`]+)`", first[0]) == re.findall(r"`([^`]+)`", second[0]), (
                    location
                )
                for column, (left, right) in enumerate(zip(first, second, strict=True)):
                    left_level = SUPPORT_LEVELS.get(left.lower())
                    right_level = SUPPORT_LEVELS.get(right.lower())
                    if left_level is not None or right_level is not None:
                        assert left_level == right_level, (*location, left, right)
                    if column in default_columns and row_index > 1:
                        assert re.findall(r"`([^`]+)`", left) == re.findall(r"`([^`]+)`", right), (
                            location
                        )


def test_bilingual_examples_keep_operations_and_calculation_parameters() -> None:
    for name, english, chinese in _markdown_pairs():
        first_blocks = [
            block for block in _fenced_blocks(english) if block[0] in {"python", "bash"}
        ]
        second_blocks = [
            block for block in _fenced_blocks(chinese) if block[0] in {"python", "bash"}
        ]
        assert len(first_blocks) == len(second_blocks), name
        for index, (first, second) in enumerate(zip(first_blocks, second_blocks, strict=True)):
            assert first[0] == second[0], (name, index)
            if first[0] == "python":
                assert _python_signature(first[1]) == _python_signature(second[1]), (name, index)
            else:
                assert shlex.split(first[1], comments=True) == shlex.split(
                    second[1], comments=True
                ), (
                    name,
                    index,
                )


def test_bilingual_notebooks_keep_cell_structure_and_calculation_parameters() -> None:
    for name, english, chinese in _notebook_pairs():
        assert [cell["cell_type"] for cell in english["cells"]] == [
            cell["cell_type"] for cell in chinese["cells"]
        ], name
        for index, (first, second) in enumerate(
            zip(english["cells"], chinese["cells"], strict=True)
        ):
            if first["cell_type"] == "code":
                assert _python_signature("".join(first["source"])) == _python_signature(
                    "".join(second["source"])
                ), (name, index)
