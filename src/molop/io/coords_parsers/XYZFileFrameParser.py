"""
Author: TMJ
Date: 2025-10-09 15:33:21
LastEditors: TMJ
LastEditTime: 2026-02-04 15:16:09
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, cast

import numpy as np
from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.patterns.XYZPatterns import xyz_patterns
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()


class XYZFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        charge = 0
        multiplicity = 1
        lines = typed_self._block.splitlines()
        atom_num = int(lines[0].strip())
        comment = lines[1].strip()
        if matches := xyz_patterns.CHARGE_MULTIPLICITY.match_content(comment):
            charge = int(matches[0][0])
            multiplicity = int(matches[0][1])
        if matches := xyz_patterns.ATOMS.match_content(typed_self._block):
            assert len(matches) == atom_num, "Atom number does not match."
            atoms = [pt.GetAtomicNumber(row[0]) for row in matches]
            coords = np.array(
                [(float(row[1]), float(row[2]), float(row[3])) for row in matches],
                dtype=np.float32,
            )
            return {
                "comment": comment,
                "atoms": atoms,
                "coords": coords * atom_ureg.angstrom,
                "charge": charge,
                "multiplicity": multiplicity,
            }
        raise ValueError("No valid atom coordinates found.")


class XYZFileFrameParserMemory(XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameMemory]):
    _file_frame_class_ = XYZFileFrameMemory


class XYZFileFrameParserDisk(XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameDisk]):
    _file_frame_class_ = XYZFileFrameDisk
