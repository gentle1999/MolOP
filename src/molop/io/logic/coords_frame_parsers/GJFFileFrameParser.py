"""
Author: TMJ
Date: 2025-10-09 15:23:26
LastEditors: TMJ
LastEditTime: 2026-02-04 15:11:28
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, cast

import numpy as np
from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.coords_frame_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.patterns.G16Patterns import g16_input_patterns
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()


class GJFFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        block = typed_self._block
        metadata: dict[str, Any] = {}
        if matches := g16_input_patterns.OPTIONS.match_content(block):
            options = "\n".join([f"{match[0]}={match[1]}" for match in matches])
            metadata.update({"options": options})
        if indexes := g16_input_patterns.ROUTE.locate_content(block):
            start_start, start_end, end_start, end_end = indexes
            if g16_input_patterns.ROUTE.content_pattern_compiled:  # noqa: SIM102
                if matches := g16_input_patterns.ROUTE.content_pattern_compiled.search(
                    block[start_start:end_end]
                ):
                    route = matches.groups()[0].strip()
                    metadata.update({"keywords": route})
            block = block[end_end:]
        if matches := g16_input_patterns.TITLE.match_content(block):
            title = matches[0][0]
            metadata.update({"title_card": title})
        if matches := g16_input_patterns.CHARGE_MULTIPLICITY.match_content(block):
            charge, multiplicity = matches[0]
            metadata.update({"charge": charge, "multiplicity": multiplicity})
        if g16_input_patterns.ATOMS.content_pattern_compiled:
            while match := g16_input_patterns.ATOMS.content_pattern_compiled.search(block):
                block = block[match.end() :]
            suffix = block.strip()
            metadata.update({"suffix": suffix})
        charge, multiplicity = None, None
        atoms: list[int] = []
        coords: np.ndarray = np.empty((0, 3), dtype=np.float32)
        if matches := g16_input_patterns.CHARGE_MULTIPLICITY.match_content(block):
            charge, multiplicity = matches[0]
        if matches := g16_input_patterns.ATOMS.match_content(typed_self._block):
            atoms = [pt.GetAtomicNumber(row[0]) for row in matches]
            coords = np.array(
                [(float(row[1]), float(row[2]), float(row[3])) for row in matches],
                dtype=np.float32,
            )
        return (
            metadata
            | (
                {
                    "atoms": atoms,
                    "coords": coords * atom_ureg.angstrom,
                }
                if atoms
                else {}
            )
            | (
                {
                    "charge": int(charge),
                    "multiplicity": int(multiplicity),
                }
                if charge and multiplicity
                else {}
            )
        )


class GJFFileFrameParserMemory(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameMemory]):
    _file_frame_class_ = GJFFileFrameMemory


class GJFFileFrameParserDisk(GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameDisk]):
    _file_frame_class_ = GJFFileFrameDisk
