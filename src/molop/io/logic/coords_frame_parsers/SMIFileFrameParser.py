"""
Author: TMJ
Date: 2025-12-14 23:30:32
LastEditors: TMJ
LastEditTime: 2026-02-04 15:15:05
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, cast

from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.coords_frame_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory
from molop.structure.StructureTransformation import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg


class SMIFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        rdmol = Chem.MolFromSmiles(typed_self._block)
        Compute2DCoords(rdmol)
        coords = rdmol.GetConformer().GetPositions()
        formal_charges = get_formal_charges(rdmol)
        formal_num_radicals = get_formal_num_radicals(rdmol)
        return {
            "atoms": [cast(Chem.Atom, a).GetAtomicNum() for a in rdmol.GetAtoms()],
            "coords": coords * atom_ureg.angstrom,
            "charge": sum(formal_charges),
            "multiplicity": sum(formal_num_radicals) + 1,
            "bonds": get_bond_pairs(rdmol),
            "formal_charges": formal_charges,
            "formal_num_radicals": formal_num_radicals,
        }


class SMIFileFrameParserMemory(SMIFileFrameParserMixin, BaseFrameParser[SMIFileFrameMemory]):
    _file_frame_class_ = SMIFileFrameMemory


class SMIFileFrameParserDisk(SMIFileFrameParserMixin, BaseFrameParser[SMIFileFrameDisk]):
    _file_frame_class_ = SMIFileFrameDisk
