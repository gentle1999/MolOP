"""
Author: TMJ
Date: 2025-07-30 10:30:03
LastEditors: TMJ
LastEditTime: 2026-02-04 15:13:04
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, cast

from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.coords_frame_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.structure.StructureTransformation import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg


class SDFFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        typed_self = cast(_HasParseMethod, self)
        suppl = Chem.SDMolSupplier()
        suppl.SetData(typed_self._block, removeHs=False, sanitize=False)
        fake_mol: Chem.Mol = next(suppl)
        formal_charges = get_formal_charges(fake_mol)
        formal_num_radicals = get_formal_num_radicals(fake_mol)
        return {
            "atoms": [cast(Chem.Atom, a).GetAtomicNum() for a in fake_mol.GetAtoms()],
            "coords": fake_mol.GetConformer().GetPositions() * atom_ureg.angstrom,
            "charge": sum(formal_charges),
            "multiplicity": sum(formal_num_radicals) + 1,
            "bonds": get_bond_pairs(fake_mol),
            "formal_charges": formal_charges,
            "formal_num_radicals": formal_num_radicals,
        }


class SDFFileFrameParserMemory(SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameMemory]):
    _file_frame_class_ = SDFFileFrameMemory


class SDFFileFrameParserDisk(SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameDisk]):
    _file_frame_class_ = SDFFileFrameDisk
