"""
Author: TMJ
Date: 2025-08-01 16:13:58
LastEditors: TMJ
LastEditTime: 2026-05-11 14:17:37
Description: 请填写简介
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal, Protocol, cast

import numpy as np
from pint.facets.plain import PlainQuantity

from molop.io.base_models.DataClasses import ImplicitSolvation, Status
from molop.io.base_models.FileParser import (
    BaseFileParserDisk,
    BaseFileParserMemory,
    _HasFileParseMethod,
)
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.QM_frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserDisk,
    G16LogFileFrameParserMemory,
)
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParserV2 import (
    G16LogFileFrameParserV2Disk,
    G16LogFileFrameParserV2Memory,
)
from molop.io.logic.QM_models.G16LogFile import BaseCalcFile, G16LogFileDisk, G16LogFileMemory
from molop.io.logic.QM_parsers._g16log_archive_tail import parse_archive_tail
from molop.io.patterns.G16Patterns import MolOPPattern, g16_log_patterns
from molop.unit import atom_ureg
from molop.utils.functions import find_rigid_transform


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


SPLIT_PATTERN = "Input orientation:"
SPLIT_PATTERN_2 = "Standard orientation:"
LINK1_SECTION_PATTERN = re.compile(
    r"^\s*(?:Link1:\s+Proceeding to internal job step number\s+\d+\.|Entering Link 1 = .*?)\s*$",
    re.MULTILINE,
)
_GAUSSIAN_PROBE_BYTES = 20000
_GAUSSIAN_FINGERPRINTS = (
    "Entering Gaussian System",
    "Gaussian 16:",
    "This is part of the Gaussian(R)",
    "Gaussian, Inc.",
)


def _ensure_gaussian_output(file_content: str) -> None:
    if not any(fingerprint in file_content for fingerprint in _GAUSSIAN_FINGERPRINTS):
        raise FormatMismatchError("Not a Gaussian output file.")


class _HasParseMethod(Protocol):
    forced_charge: int | None = None
    forced_multiplicity: int | None = None
    only_extract_structure: bool = False
    only_last_frame: bool = False

    _chem_file: type[BaseCalcFile]

    def _parse_frame(
        self, frame_content: str, additional_data: dict[str, Any]
    ) -> G16LogFileFrameDisk | G16LogFileFrameMemory: ...


class _HasMetadataFinalizeMethod(Protocol):
    def _update_file_metadata_from_frames(
        self,
        chem_file: Any,
        metadata: dict[str, Any],
    ) -> None: ...


class G16LogFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        _ensure_gaussian_output(file_content[:_GAUSSIAN_PROBE_BYTES])

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
        raw_file_content = file_content
        self._file_content = file_content
        if version := self._parse_version():
            metadata["qm_software_version"] = version
        if options := self._parse_options():
            metadata["options"] = options
        if route := self._parse_keywords():
            metadata["keywords"] = route
        if title := self._parse_title():
            metadata["title_card"] = title
        if charge_multiplicity := self._parse_charge_multiplicity():
            metadata["charge"], metadata["multiplicity"] = charge_multiplicity
        transform_matrix = self._parse_standard_orientation_transformation_matrix()
        if transform_matrix is not None:
            metadata["standard_orientation_transformation_matrix"] = transform_matrix
        if cast(_HasFileParseMethod, self).only_extract_structure:
            return metadata
        if solvent := self._parse_solvent():
            metadata["solvent"] = solvent
        temperature, pressure = self._parse_temperature_and_pressure(raw_file_content)
        if temperature:
            metadata["temperature"] = temperature
        if pressure:
            metadata["pressure"] = pressure
        if tail := parse_archive_tail(file_content):
            metadata.update(tail[0])
        if running_time := self._parse_running_time():
            metadata["running_time"] = running_time
        if status := self._parse_termination_status():
            metadata["status"] = status
        return metadata

    def _split_sections(self, file_content: str) -> list[str]:
        matches = list(LINK1_SECTION_PATTERN.finditer(file_content))
        if not matches:
            return [file_content]
        sections: list[str] = []
        # first_start = matches[0].start()
        # if first_start > 0:
        #     sections.append(file_content[:first_start])
        for idx, match in enumerate(matches):
            next_start = matches[idx + 1].start() if idx + 1 < len(matches) else len(file_content)
            sections.append(file_content[match.start() : next_start])
        return [section for section in sections if section.strip()]

    def _split_section_frames(self, section_content: str) -> list[str]:
        split_ = SPLIT_PATTERN
        if SPLIT_PATTERN in section_content:
            split_ = SPLIT_PATTERN
        elif SPLIT_PATTERN_2 in section_content:
            split_ = SPLIT_PATTERN_2
        else:
            return []
        fragments = section_content.split(split_)
        header = fragments[0]
        frame_contents: list[str] = []
        for idx, fragment in enumerate(fragments[1:]):
            frame_block = f"{split_}\n{fragment}"
            if idx == 0 and header.strip():
                frame_block = f"{header}{frame_block}"
            frame_contents.append(frame_block)
        return frame_contents

    def _split_file(self, file_content: str) -> Sequence[str]:
        sections = self._split_sections(file_content)
        return [frame for section in sections for frame in self._split_section_frames(section)]

    @staticmethod
    def _first_frame_value(frames: Sequence[Any], field: str) -> Any:
        for frame in frames:
            value = getattr(frame, field, None)
            if value is not None:
                return value
        return None

    @staticmethod
    def _last_frame_value(frames: Sequence[Any], field: str) -> Any:
        for frame in reversed(frames):
            value = getattr(frame, field, None)
            if value is not None:
                return value
        return None

    def _update_file_metadata_from_frames(self, chem_file: Any, metadata: dict[str, Any]) -> None:
        base_parser = cast(
            _HasMetadataFinalizeMethod,
            super(),
        )
        base_parser._update_file_metadata_from_frames(chem_file, metadata)
        frames = list(getattr(chem_file, "frames", []))
        if not frames:
            return

        for field in ("solvent", "temperature", "pressure"):
            if metadata.get(field) is not None:
                continue
            value = self._first_frame_value(frames, field)
            if value is not None:
                metadata[field] = value
                setattr(chem_file, field, value)

        for field in ("status", "geometry_optimization_status"):
            value = self._last_frame_value(frames, field)
            if value is not None:
                metadata[field] = value
                setattr(chem_file, field, value)

    # override the _parse method
    def _parse(
        self,
        source: str,
        source_type: Literal["file_path", "string"] = "file_path",
        *,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
    ) -> Any:
        self._file_path: str | None
        metadata_base: dict[str, Any]
        if source_type == "file_path":
            self._file_path = source
            with open(source) as f:
                file_content = f.read()
            metadata_base = {"file_path": source, "file_content": file_content}
        elif source_type == "string":
            self._file_path = None
            file_content = source
            metadata_base = {"file_content": file_content}
        else:
            raise ValueError(f"Invalid source_type: {source_type}")
        self._quick_check_file_format(file_content)

        typed_self = cast(_HasParseMethod, self)
        final_charge = total_charge if total_charge is not None else typed_self.forced_charge
        final_multiplicity = (
            total_multiplicity if total_multiplicity is not None else typed_self.forced_multiplicity
        )
        if final_charge is not None:
            metadata_base["charge"] = final_charge
        if final_multiplicity is not None:
            metadata_base["multiplicity"] = final_multiplicity
        if typed_self.only_last_frame:
            sections = self._split_sections(file_content)
            if sections:
                section = sections[-1]
            else:
                raise ValueError("No section found in the file content.")
            metadata = self._parse_metadata(section)
            file_metadata = metadata | metadata_base
            _chem_file = typed_self._chem_file.model_validate(file_metadata)
            frame_contents = self._split_section_frames(section)
            if frame_contents:
                last_frame_content = frame_contents[-1]
            else:
                raise ValueError("No frame found in the section.")
            frame = typed_self._parse_frame(last_frame_content, additional_data=file_metadata)
            if final_charge is not None:
                frame.charge = final_charge
            if final_multiplicity is not None:
                frame.multiplicity = final_multiplicity
            _chem_file.append(frame)
        else:
            sections = self._split_sections(file_content)
            metadata_list = [self._parse_metadata(section) for section in sections]
            if not metadata_list:
                raise ValueError("No metadata found in the file content.")
            file_metadata = metadata_list[0] | metadata_base
            running_times = [
                running_time
                for metadata in metadata_list
                if (running_time := metadata.get("running_time")) is not None
            ]
            if running_times:
                file_metadata["running_time"] = sum(running_times[1:], running_times[0])
            _chem_file = typed_self._chem_file.model_validate(file_metadata)
            for metadata, section in zip(metadata_list, sections, strict=True):
                section_metadata = metadata | metadata_base
                for frame_content in self._split_section_frames(section):
                    frame = typed_self._parse_frame(frame_content, additional_data=section_metadata)
                    if final_charge is not None:
                        frame.charge = final_charge
                    if final_multiplicity is not None:
                        frame.multiplicity = final_multiplicity
                    if (
                        frame.basis_set.lower() == "gen"
                        and _chem_file
                        and _chem_file[-1].basis_set.lower() != "gen"
                    ):
                        frame.basis_set = _chem_file[-1].basis_set
                    _chem_file.append(frame)

        self._update_file_metadata_from_frames(_chem_file, file_metadata)
        return _chem_file

    # -------------------------------------------
    # Details of the parsing process for metadata
    # -------------------------------------------

    def _parse_version(self) -> str | None:
        focus_content, self._file_content = g16_log_patterns.VERSION.split_content(
            self._file_content
        )
        if matches := g16_log_patterns.VERSION.get_matches(focus_content):
            return matches[0][0]
        return None

    def _parse_options(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.OPTIONS.split_content(
            self._file_content
        )
        if matches := g16_log_patterns.OPTIONS.get_matches(focus_content):
            options = "\n".join([f"{match[0]}={match[1]}" for match in matches])
            self._file_content = continued_content
            return options
        return None

    def _parse_keywords(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.KEYWORDS.split_content(
            self._file_content
        )
        if len(keyword_lines := focus_content.splitlines()) >= 3:
            self._file_content = continued_content
            return "\n".join(keyword_lines[1:-1]).replace("\n ", "")
        return None

    def _parse_title(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.TITLE.split_content(self._file_content)
        if len(title_lines := focus_content.splitlines()) >= 3:
            self._file_content = continued_content
            return "\n".join(title_lines[1:-1]).replace("\n ", "")
        return None

    def _parse_charge_multiplicity(self) -> tuple[int, int] | tuple[None, None]:
        if match := g16_log_patterns.CHARGE_MULTIPLICITY.match_content(self._file_content):
            charge, multiplicity = int(match[0][0]), int(match[0][1])
            return charge, multiplicity
        return None, None

    def _parse_standard_orientation_transformation_matrix(self) -> np.ndarray | None:
        try:
            coords = self._parse_coordinates(g16_log_patterns.INITIAL_INPUT_COORDS)
            if coords is None:
                return None
            standard_coords = self._parse_coordinates(g16_log_patterns.STANDARD_COORDS)
            if standard_coords is None:
                return None
            return find_rigid_transform(coords, standard_coords)
        except Exception as e:
            raise e

    def _parse_coordinates(self, pattern: MolOPPattern) -> np.ndarray | None:
        focus_content, self._file_content = pattern.split_content(self._file_content)
        if matches := pattern.match_content(focus_content):
            return np.array(
                [list(map(float, match[1:])) for match in matches],
            )
        return None

    def _parse_solvent(self) -> ImplicitSolvation | None:
        solvent_dict: dict[str, Any] = {}
        focus_content, self._file_content = g16_log_patterns.SOLVENT_PARAMETERS.split_content(
            self._file_content
        )
        if matches := g16_log_patterns.SOLVENT_MODEL.get_matches(focus_content):
            solvent_dict["solvent_model"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_ATOM_RADII.get_matches(focus_content):
            solvent_dict["atomic_radii"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_TYPE.get_matches(focus_content):
            solvent_dict["solvent"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_EPS.get_matches(focus_content):
            solvent_dict["solvent_epsilon"] = float(matches[0][0])
        if matches := g16_log_patterns.SOLVENT_EPS_INF.get_matches(focus_content):
            solvent_dict["solvent_epsilon_infinite"] = float(matches[0][-1])
        if solvent_dict:
            return ImplicitSolvation.model_validate(solvent_dict)
        return None

    def _parse_temperature_and_pressure(
        self,
        content: str | None = None,
    ) -> tuple[PlainQuantity, PlainQuantity] | tuple[None, None]:
        source_content = self._file_content if content is None else content
        index = source_content.find(" - Thermochemistry -")
        if index == -1:
            return None, None
        focus_content = source_content[index:]
        # self._file_content = self._file_content[index:]
        if match := g16_log_patterns.TEMPEREATURE_PRESSURE.match_content(focus_content):
            temperature, pressure = (
                float(match[0][0]) * atom_ureg.K,
                float(match[0][1]) * atom_ureg.atm,
            )
            return temperature, pressure
        return None, None

    def _parse_running_time(self) -> PlainQuantity | None:
        if matches := g16_log_patterns.JOB_TIME.match_content(self._file_content):
            total_seconds = 0.0
            for match in matches:
                days, hours, minutes, seconds = [float(i) for i in match[1:]]
                total_seconds += (days * 24 + hours) * 3600 + minutes * 60 + seconds
            return total_seconds * atom_ureg.second
        return None

    def _parse_termination_status(self) -> Status | None:
        status_dict: dict[str, Any] = {}
        if matches := g16_log_patterns.TERMINATION_STATUS.match_content(self._file_content):
            status = matches[-1][0]
            if status == "Normal":
                status_dict["normal_terminated"] = True
                status_dict["scf_converged"] = True
            else:
                status_dict["normal_terminated"] = False
        if status_dict:
            return Status.model_validate(status_dict)
        return None


class G16LogFileParserMemory(
    G16LogFileParserMixin,
    BaseFileParserMemory[G16LogFileMemory, G16LogFileFrameMemory, G16LogFileFrameParserMemory],
):
    _chem_file = G16LogFileMemory
    _frame_parser = G16LogFileFrameParserMemory


class G16LogFileParserDisk(
    G16LogFileParserMixin,
    BaseFileParserDisk[G16LogFileDisk, G16LogFileFrameDisk, G16LogFileFrameParserDisk],
):
    allowed_formats = (".log", ".g16", ".gal", ".out", ".irc", "gau")
    _chem_file = G16LogFileDisk
    _frame_parser = G16LogFileFrameParserDisk


class G16LogFileParserV2Memory(
    G16LogFileParserMixin,
    BaseFileParserMemory[G16LogFileMemory, G16LogFileFrameMemory, G16LogFileFrameParserV2Memory],
):
    _chem_file = G16LogFileMemory
    _frame_parser = G16LogFileFrameParserV2Memory


class G16LogFileParserV2Disk(
    G16LogFileParserMixin,
    BaseFileParserDisk[G16LogFileDisk, G16LogFileFrameDisk, G16LogFileFrameParserV2Disk],
):
    allowed_formats = (".log", ".g16", ".gal", ".out", ".irc", "gau")
    _chem_file = G16LogFileDisk
    _frame_parser = G16LogFileFrameParserV2Disk


def register(registry: Registry) -> None:
    """Register this file parser as a reader codec.

    Called by lazy activation via `molop.io.codecs.catalog`.
    """

    from typing import cast

    from molop.io.codecs._shared.reader_helpers import (
        ParserDiskReader,
        ReaderCodec,
        StructureLevel,
        extensions_for_parser,
    )

    extensions = frozenset(extensions_for_parser(G16LogFileParserV2Disk))
    priority = 100

    @registry.reader_factory(format_id="g16log", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="g16log",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=G16LogFileParserV2Disk,
                priority=priority,
            ),
        )
