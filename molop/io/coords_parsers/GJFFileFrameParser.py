"""
Author: TMJ
Date: 2025-10-09 15:23:26
LastEditors: TMJ
LastEditTime: 2025-11-21 16:55:39
Description: 请填写简介
"""

from typing import TYPE_CHECKING, Any, Mapping, Protocol

import numpy as np
from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.coords_models.GJFFileFrame import GJFFileFrameDisk, GJFFileFrameMemory
from molop.io.patterns.G16Patterns import g16_input_patterns
from molop.unit import atom_ureg

pt = Chem.GetPeriodicTable()


class GJFFileFrameParserProtocol(Protocol):
    _block: str
    only_extract_structure: bool
    _file_frame_class_: type[GJFFileFrameMemory] | type[GJFFileFrameDisk]


if TYPE_CHECKING:

    class _GJFFileFrameParserProtocol(GJFFileFrameParserProtocol): ...
else:

    class _GJFFileFrameParserProtocol(object): ...


class GJFFileFrameParserMixin(_GJFFileFrameParserProtocol):
    def _parse_frame(self) -> Mapping[str, Any]:
        block = self._block
        if matches := g16_input_patterns.CHARGE_MULTIPLICITY.match_content(block):
            charge, multiplicity = matches[0]
        if matches := g16_input_patterns.ATOMS.match_content(self._block):
            atoms = [pt.GetAtomicNumber(row[0]) for row in matches]
            coords = np.array(
                [(float(row[1]), float(row[2]), float(row[3])) for row in matches],
                dtype=np.float32,
            )
        return {
            "atoms": atoms,
            "coords": coords * atom_ureg.angstrom,
            "charge": int(charge),
            "multiplicity": int(multiplicity),
        }


class GJFFileFrameParserMemory(
    GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameMemory]
):
    _file_frame_class_ = GJFFileFrameMemory


class GJFFileFrameParserDisk(
    GJFFileFrameParserMixin, BaseFrameParser[GJFFileFrameDisk]
):
    _file_frame_class_ = GJFFileFrameDisk
