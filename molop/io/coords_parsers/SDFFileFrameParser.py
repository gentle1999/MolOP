"""
Author: TMJ
Date: 2025-07-30 10:30:03
LastEditors: TMJ
LastEditTime: 2025-08-21 00:06:43
Description: 请填写简介
"""

from typing import Any, Mapping, Protocol

from rdkit import Chem

from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.coords_models.SDFFileFrame import SDFFileFrameDisk, SDFFileFrameMemory
from molop.structure.StructureTransformation import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg


class SDFFileFrameParser(Protocol):
    _block: str
    only_extract_structure: bool
    _file_frame_class_: type[SDFFileFrameMemory] | type[SDFFileFrameDisk]

    def _parse_frame(self) -> Mapping[str, Any]: ...


class SDFFileFrameParserMixin:

    def _parse_frame(self: SDFFileFrameParser) -> Mapping[str, Any]:
        fake_mol = Chem.MolFromMolBlock(self._block, removeHs=False, sanitize=False)
        formal_charges = get_formal_charges(fake_mol)
        formal_num_radicals = get_formal_num_radicals(fake_mol)
        return {
            "atoms": [a.GetAtomicNum() for a in fake_mol.GetAtoms()],
            "coords": fake_mol.GetConformer().GetPositions() * atom_ureg.angstrom,
            "charge": sum(formal_charges),
            "multiplicity": sum(formal_num_radicals) + 1,
            "bonds": get_bond_pairs(fake_mol),
            "formal_charges": formal_charges,
            "formal_num_radicals": formal_num_radicals,
        }


class SDFFileFrameParserMemory(
    SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameMemory]
):
    _file_frame_class_ = SDFFileFrameMemory


class SDFFileFrameParserDisk(
    SDFFileFrameParserMixin, BaseFrameParser[SDFFileFrameDisk]
):
    _file_frame_class_ = SDFFileFrameDisk
