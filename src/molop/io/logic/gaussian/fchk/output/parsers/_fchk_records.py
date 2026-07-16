from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from molop.io.codec_exceptions import FormatMismatchError


FCHKScalar: TypeAlias = int | float | str | bool
FCHKValue: TypeAlias = FCHKScalar | list[int] | list[float] | list[str] | list[bool]


@dataclass(frozen=True, slots=True)
class FCHKHeader:
    title: str
    job_type: str
    method: str
    basis_set: str


@dataclass(frozen=True, slots=True)
class FCHKRecord:
    label: str
    data_type: str
    value: FCHKValue
    count: int | None = None


def _to_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_fchk_header(file_content: str) -> FCHKHeader:
    lines = file_content.splitlines()
    if len(lines) < 3:
        raise FormatMismatchError("Not a Gaussian formatted checkpoint: missing header.")
    tokens = lines[1].split()
    if len(tokens) < 2:
        raise FormatMismatchError("Not a Gaussian formatted checkpoint: invalid method header.")
    return FCHKHeader(
        title=lines[0].rstrip(),
        job_type=tokens[0],
        method=tokens[1],
        basis_set=tokens[2] if len(tokens) >= 3 else "",
    )


def ensure_fchk_content(file_content: str) -> None:
    parse_fchk_header(file_content)
    prefix_lines = file_content.splitlines()[:80]
    labels = {line[:40].rstrip() for line in prefix_lines if len(line) >= 44}
    required = {"Number of atoms", "Atomic numbers", "Current cartesian coordinates"}
    if not required.issubset(labels):
        raise FormatMismatchError(
            "Not a Gaussian formatted checkpoint: missing required structure records."
        )


def _parse_scalar(data_type: str, value: str) -> FCHKScalar:
    if data_type == "I":
        return int(value)
    if data_type == "R":
        return _to_float(value)
    if data_type == "L":
        return value.strip().upper().startswith(("T", "1"))
    return value.rstrip()


def _parse_array_tokens(data_type: str, tokens: list[str]) -> FCHKValue:
    if data_type == "I":
        return [int(token) for token in tokens]
    if data_type == "R":
        return [_to_float(token) for token in tokens]
    if data_type == "L":
        return [token.upper().startswith(("T", "1")) for token in tokens]
    return tokens


def parse_fchk_records(
    file_content: str,
    *,
    wanted_labels: set[str] | frozenset[str] | None = None,
) -> dict[str, FCHKRecord]:
    lines = file_content.splitlines()
    records: dict[str, FCHKRecord] = {}
    index = 2
    while index < len(lines):
        line = lines[index]
        index += 1
        if len(line) < 44:
            continue
        data_type = line[43]
        if data_type not in {"I", "R", "C", "L", "H"}:
            continue
        label = line[:40].rstrip()
        if not label:
            continue
        remainder = line[44:].strip()
        if not remainder.startswith("N="):
            if wanted_labels is None or label in wanted_labels:
                records[label] = FCHKRecord(
                    label=label,
                    data_type=data_type,
                    value=_parse_scalar(data_type, remainder),
                )
            continue

        try:
            count = int(remainder.split("=", 1)[1])
        except ValueError as exc:
            raise FormatMismatchError(f"Invalid fchk array count for {label!r}.") from exc

        keep = wanted_labels is None or label in wanted_labels
        if data_type in {"C", "H"}:
            chunks: list[str] = []
            while len(chunks) < count and index < len(lines):
                data_line = lines[index]
                index += 1
                padded_length = max(12, ((len(data_line) + 11) // 12) * 12)
                padded = data_line.ljust(padded_length)
                chunks.extend(padded[offset : offset + 12] for offset in range(0, len(padded), 12))
            if len(chunks) < count:
                raise FormatMismatchError(f"Truncated fchk character array {label!r}.")
            if keep:
                records[label] = FCHKRecord(
                    label=label,
                    data_type=data_type,
                    value="".join(chunks[:count]).rstrip(),
                    count=count,
                )
            continue

        tokens: list[str] = []
        while len(tokens) < count and index < len(lines):
            tokens.extend(lines[index].split())
            index += 1
        if len(tokens) < count:
            raise FormatMismatchError(f"Truncated fchk numeric array {label!r}.")
        if keep:
            records[label] = FCHKRecord(
                label=label,
                data_type=data_type,
                value=_parse_array_tokens(data_type, tokens[:count]),
                count=count,
            )
    return records


__all__ = [
    "FCHKHeader",
    "FCHKRecord",
    "FCHKScalar",
    "FCHKValue",
    "ensure_fchk_content",
    "parse_fchk_header",
    "parse_fchk_records",
]
