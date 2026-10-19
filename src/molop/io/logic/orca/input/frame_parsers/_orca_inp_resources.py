"""ORCA requested resources and print settings.

Basis: ORCA-RESOURCE, ORCA-OUTPUT. See
``docs/en/reference/formats/orcainp.md`` for verified
official manuals and tutorials. For field mappings and MolOP conventions, see
``docs/en/reference/data_model.md``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from molop.io.logic.orca.common import ORCABlock, ORCAOutputPrintSetting
from molop.io.logic.orca.input.frame_parsers._orca_inp_patterns import orca_inp_patterns
from molop.io.logic.orca.input.frame_parsers._orca_inp_tokens import (
    parse_orca_float,
    parse_orca_int,
    split_orca_key_value,
    strip_orca_inline_comment,
)
from molop.unit import atom_ureg


def render_resource_blocks(blocks: Sequence[ORCABlock]) -> str:
    resources: list[str] = []
    for block in blocks:
        if block.raw_text:
            resources.append(block.raw_text)
            continue
        header = block.raw_header or f"%{block.name}"
        body = block.body_text()
        resources.append(header if not body else f"{header}\n{body}")
    return "\n".join(resources)


def parse_request_num_cpu(blocks: Sequence[ORCABlock]) -> int | None:
    for block in blocks:
        if block.name.lower() != "pal":
            continue
        for line in block.lines:
            key, value = split_orca_key_value(line.text)
            if key.lower() != "nprocs":
                continue
            parsed = parse_orca_int(value.rstrip(";"))
            if parsed is not None:
                return parsed
    return None


def parse_request_memory(blocks: Sequence[ORCABlock]) -> Any | None:
    for block in blocks:
        if block.name.lower() != "maxcore":
            continue
        values: list[str] = []
        header_tokens = strip_orca_inline_comment(block.raw_header).split()
        if len(header_tokens) > 1:
            values.extend(header_tokens[1:])
        for line in block.lines:
            values.extend(strip_orca_inline_comment(line.text).split())
        for value in values:
            parsed = parse_orca_float(value)
            if parsed is not None:
                return parsed * atom_ureg.megabyte
    return None


def parse_output_print_settings(blocks: Sequence[ORCABlock]) -> list[ORCAOutputPrintSetting]:
    settings: list[ORCAOutputPrintSetting] = []
    for block in blocks:
        if block.name.lower() != "output":
            continue
        for line in block.lines:
            matched = orca_inp_patterns.OUTPUT_PRINT_SETTING.match(
                strip_orca_inline_comment(line.text).strip()
            )
            if matched is None:
                continue
            settings.append(
                ORCAOutputPrintSetting(
                    target=matched.group("target").strip(),
                    value=matched.group("value").strip(),
                )
            )
    return settings
