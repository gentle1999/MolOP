"""
Author: TMJ
Date: 2025-07-29 22:54:56
LastEditors: TMJ
LastEditTime: 2025-09-17 20:02:05
Description: 请填写简介
"""

from typing import Any, Mapping, Protocol

import numpy as np
from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.coords_models.XYZFileFrame import XYZFileFrameDisk, XYZFileFrameMemory
from molop.io.patterns.XYZPatterns import XYZPatterns
from molop.unit import atom_ureg

pt = Chem.GetPeriodicTable()

xyz_patterns = XYZPatterns()


class XYZFileFrameParser(Protocol):
    _block: str
    only_extract_structure: bool
    _file_frame_class_: type[XYZFileFrameMemory] | type[XYZFileFrameDisk]

    def _parse_frame(self) -> Mapping[str, Any]: ...


class XYZFileFrameParserMixin:

    def _parse_frame(self: XYZFileFrameParser) -> Mapping[str, Any]:
        lines = self._block.splitlines()
        atom_num = int(lines[0].strip())
        comment = lines[1].strip()
        if matches := xyz_patterns.ATOMS.match_content(self._block):
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
                "charge": 0,
                "multiplicity": 1,
            }
        raise ValueError("No valid atom coordinates found.")


class XYZFileFrameParserMemory(
    XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameMemory]
):
    _file_frame_class_ = XYZFileFrameMemory


class XYZFileFrameParserDisk(
    XYZFileFrameParserMixin, BaseFrameParser[XYZFileFrameDisk]
):
    _file_frame_class_ = XYZFileFrameDisk
