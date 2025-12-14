"""
Author: TMJ
Date: 2025-07-29 22:14:37
LastEditors: TMJ
LastEditTime: 2025-12-14 21:34:43
Description: 请填写简介
"""

import os
from abc import abstractmethod
from typing import Any, Generic, Optional, Sequence

from pydantic import Field, PrivateAttr

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFile import ChemFile
from molop.io.base_models.ChemFileFrame import ChemFileFrame
from molop.io.base_models.FrameParser import FrameParser


class BaseFileParser(
    BaseDataClassWithUnit, Generic[ChemFile, ChemFileFrame, FrameParser]
):
    # add by subclass
    _frame_parser: type[FrameParser] = PrivateAttr()
    _chem_file: type[ChemFile] = PrivateAttr()
    _file_path: Optional[str] = PrivateAttr(default=None)

    forced_charge: Optional[int] = Field(
        default=None,
        description="The forced charge of the molecule.",
        exclude=True,
        repr=False,
    )
    forced_multiplicity: Optional[int] = Field(
        default=None,
        description="The forced multiplicity of the molecule.",
        ge=1,
        exclude=True,
        repr=False,
    )
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    only_last_frame: bool = Field(default=False, exclude=True, repr=False)

    @abstractmethod
    def _split_file(self, file_content: str) -> Sequence[str]:
        """Split the file content into frames."""
        raise NotImplementedError

    @abstractmethod
    def _parse_metadata(self, file_content: str) -> Optional[dict[str, Any]]:
        """Parse the metadata of the file. Return a dictionary as additional information in chem file frame."""
        raise NotImplementedError

    def _parse_frame(
        self, frame_content: str, *, additional_data: Optional[dict[str, Any]] = None
    ) -> ChemFileFrame:
        """Parse a single frame."""
        frame_parser = self._frame_parser()
        return frame_parser.parse(frame_content, additional_data=additional_data)

    def _parse(
        self,
        file_content: str,
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
    ) -> ChemFile:
        """Parse the file content."""
        metadata: dict[str, Any] = {
            "file_path": self._file_path,
            "file_content": file_content,
        }
        if mata := self._parse_metadata(file_content):
            metadata.update(mata)
        frame_contents = self._split_file(file_content)
        if self.only_last_frame:
            frame_contents = [frame_contents[-1]]
        final_charge = total_charge or self.forced_charge
        final_multiplicity = total_multiplicity or self.forced_multiplicity
        if final_charge:
            metadata["charge"] = final_charge
        if final_multiplicity:
            metadata["multiplicity"] = final_multiplicity
        _chem_file = self._chem_file.model_validate(metadata)
        for frame_content in frame_contents:
            frame = self._parse_frame(frame_content, additional_data=metadata)
            if final_charge:
                frame.charge = final_charge
            if final_multiplicity:
                frame.multiplicity = final_multiplicity
            _chem_file.append(frame)
        return _chem_file

    def _add_default_units(self) -> None: ...

    def to_summary_dict(self, **kwargs) -> dict[tuple[str, str], Any]:
        return {
            ("FileParser", "forced_charge"): self.forced_charge,
            ("FileParser", "forced_multiplicity"): self.forced_multiplicity,
            ("FileParser", "only_extract_structure"): self.only_extract_structure,
            ("FileParser", "only_last_frame"): self.only_last_frame,
        }


class BaseFileParserMemory(BaseFileParser[ChemFile, ChemFileFrame, FrameParser]):
    def parse(
        self,
        file_content: str,
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
        release_file_content: bool = False,
    ) -> ChemFile:
        _chem_file = self._parse(file_content, total_charge, total_multiplicity)
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file


class BaseFileParserDisk(BaseFileParser[ChemFile, ChemFileFrame, FrameParser]):
    _allowed_formats_: Sequence[str] = PrivateAttr(default_factory=tuple)

    @classmethod
    def _check_file_path(cls, file_path: str) -> str:
        """
        Validate the file path.

        Parameters:
            file_path (str): The file path to validate.

        Returns:
            str: The validated file path.

        Raises:
            ValueError: If the file path does not exist or is not a file.
        """
        if not os.path.exists(file_path):
            raise ValueError(f"File {file_path} does not exist.")
        if not os.path.isfile(file_path):
            raise ValueError(f"Path {file_path} is not a file.")
        for fmt in cls._allowed_formats_.default:  # type: ignore
            if file_path.endswith(fmt):
                return file_path
        raise ValueError(
            f"File {file_path} has an invalid format. "
            f"Allowed formats: {cls._allowed_formats_.default}."  # type: ignore
        )

    def parse(
        self,
        file_path: str,
        total_charge: Optional[int] = None,
        total_multiplicity: Optional[int] = None,
        release_file_content: bool = False,
    ) -> ChemFile:
        with open(file_path, "r") as f:
            file_content = f.read()
        self._file_path = self._check_file_path(os.path.abspath(file_path))
        _chem_file = self._parse(file_content, total_charge, total_multiplicity)
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file
