"""
Author: TMJ
Date: 2025-12-14 23:30:32
LastEditors: TMJ
LastEditTime: 2025-12-14 23:34:21
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Protocol

from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.coords_models.SMIFileFrame import SMIFileFrameDisk, SMIFileFrameMemory
from molop.structure.StructureTransformation import (
    get_bond_pairs,
    get_formal_charges,
    get_formal_num_radicals,
)
from molop.unit import atom_ureg


class SMIFileFrameParserProtocol(Protocol):
    _block: str
    only_extract_structure: bool
    _file_frame_class_: type[SMIFileFrameMemory] | type[SMIFileFrameDisk]


if TYPE_CHECKING:

    class _SMIFileFrameParserProtocol(SMIFileFrameParserProtocol): ...
else:

    class _SMIFileFrameParserProtocol: ...


class SMIFileFrameParserMixin(_SMIFileFrameParserProtocol):
    def _parse_frame(self) -> Mapping[str, Any]:
        rdmol = Chem.MolFromSmiles(self._block)
        Compute2DCoords(rdmol)
        coords = rdmol.GetConformer().GetPositions()
        formal_charges = get_formal_charges(rdmol)
        formal_num_radicals = get_formal_num_radicals(rdmol)
        return {
            "atoms": [a.GetAtomicNum() for a in rdmol.GetAtoms()],
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
