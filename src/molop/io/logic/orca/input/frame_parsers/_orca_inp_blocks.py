from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from molop.io.logic.orca.common import ORCABlock, ORCABlockLine
from molop.io.logic.orca.input.frame_parsers._orca_inp_tokens import (
    strip_orca_inline_comment,
)


_NESTED_BLOCK_KEYWORDS_BY_BLOCK = {
    "compound": {"new_step"},
    "geom": {
        "connectfragments",
        "constraints",
        "constrainfragments",
        "hybrid_hess",
        "inhess",
        "hess_internal",
        "modify_internal",
        "potentials",
        "scan",
        "ts_active_atoms",
        "ts_mode",
    },
    "mrci": {"newblock", "refs"},
}


@dataclass(slots=True)
class LineSpan:
    start: int
    end: int
    line: str


@dataclass(slots=True)
class ORCABlockSpan:
    block: ORCABlock
    line_start: int
    line_end: int


def build_line_spans(block: str) -> list[LineSpan]:
    spans: list[LineSpan] = []
    cursor = 0
    for line in block.splitlines(keepends=True):
        next_cursor = cursor + len(line)
        spans.append(LineSpan(start=cursor, end=next_cursor, line=line))
        cursor = next_cursor
    if not spans and block:
        spans.append(LineSpan(start=0, end=len(block), line=block))
    return spans


def is_coords_open_line(line: str) -> bool:
    lower = line.lower()
    return lower.startswith("coords") and (len(lower) == 6 or lower[6].isspace())


def _parse_percent_header(line: str) -> tuple[str, str] | None:
    stripped = line.strip()
    if not stripped.startswith("%"):
        return None
    rest = stripped[1:].strip()
    if not rest:
        return None
    tokens = rest.split(maxsplit=1)
    return tokens[0].lower(), tokens[1].strip() if len(tokens) > 1 else ""


def _is_top_level_directive(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return False
    lowered = stripped.lower()
    return stripped.startswith(("!", "*", "%")) or lowered.startswith("$new_job")


def _first_code_token(line: str) -> str:
    content = strip_orca_inline_comment(line).strip()
    if not content:
        return ""
    return content.split(maxsplit=1)[0].lower()


def _code_tokens(line: str) -> list[str]:
    return strip_orca_inline_comment(line).strip().lower().split()


def _opens_nested_block(parent_name: str, line: str) -> bool:
    tokens = _code_tokens(line)
    if not tokens:
        return False
    return (
        tokens[0] in _NESTED_BLOCK_KEYWORDS_BY_BLOCK.get(parent_name.lower(), set())
        and "end" not in tokens[1:]
    )


def _closes_nested_block(parent_name: str, line: str) -> bool:
    return parent_name.lower() == "compound" and _first_code_token(line) == "step_end"


def _block_body_from_inline(inline: str) -> tuple[list[ORCABlockLine], bool]:
    if not inline:
        return [], False
    tokens = inline.split()
    if tokens and tokens[-1].lower() == "end":
        body = " ".join(tokens[:-1]).strip()
        return ([ORCABlockLine(text=body)] if body else []), True
    return [ORCABlockLine(text=inline)], False


def extract_block_spans(block: str) -> list[ORCABlockSpan]:
    spans = build_line_spans(block)
    result: list[ORCABlockSpan] = []
    idx = 0
    while idx < len(spans):
        header = _parse_percent_header(spans[idx].line)
        if header is None:
            idx += 1
            continue

        name, inline = header
        body_lines, inline_closed = _block_body_from_inline(inline)
        end_idx = idx + 1

        nested_depth = (
            1 if inline and not inline_closed and _opens_nested_block(name, inline) else 0
        )
        if not inline_closed:
            scan_idx = idx + 1
            while scan_idx < len(spans):
                candidate = spans[scan_idx].line.strip()
                lowered = candidate.lower()
                if is_coords_open_line(lowered):
                    nested_depth += 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if _opens_nested_block(name, candidate):
                    nested_depth += 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if _closes_nested_block(name, candidate) and nested_depth > 0:
                    nested_depth -= 1
                    body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                    scan_idx += 1
                    continue
                if lowered == "end":
                    if nested_depth > 0:
                        nested_depth -= 1
                        body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                        scan_idx += 1
                        continue
                    end_idx = scan_idx + 1
                    break
                if nested_depth == 0 and _is_top_level_directive(candidate):
                    end_idx = scan_idx
                    break
                body_lines.append(ORCABlockLine(text=spans[scan_idx].line.rstrip("\r\n")))
                scan_idx += 1
            else:
                end_idx = len(spans)

        result.append(
            ORCABlockSpan(
                block=ORCABlock(
                    name=name,
                    lines=[line for line in body_lines if line.text.strip()],
                    raw_header=spans[idx].line.rstrip("\r\n"),
                    raw_text="".join(span.line for span in spans[idx:end_idx]).rstrip("\r\n"),
                ),
                line_start=idx,
                line_end=end_idx,
            )
        )
        idx = max(end_idx, idx + 1)
    return result


def line_in_spans(line_idx: int, spans: Sequence[tuple[int, int]]) -> bool:
    return any(start <= line_idx < end for start, end in spans)
