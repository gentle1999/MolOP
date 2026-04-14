"""
Author: TMJ
Date: 2025-08-01 16:13:58
LastEditors: TMJ
LastEditTime: 2026-02-05 19:54:25
Description: 请填写简介
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Protocol, cast

import numpy as np
from pint.facets.plain import PlainQuantity

from molop.io.base_models.DataClasses import ImplicitSolvation, Status
from molop.io.base_models.FileParser import (
    BaseFileParserDisk,
    BaseFileParserMemory,
    _HasFileParseMethod,
)
from molop.io.logic.QM_frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParserV2 import (
    G16LogFileFrameParserV2Disk,
    G16LogFileFrameParserV2Memory,
)
from molop.io.logic.QM_models.G16LogFile import G16LogFileDisk, G16LogFileMemory
from molop.io.logic.QM_parsers._g16log_archive_tail import (
    parse_archive_tail,
)
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
SECTION_METADATA_KEYS = (
    "qm_software",
    "qm_software_version",
    "options",
    "keywords",
    "title_card",
    "charge",
    "multiplicity",
    "standard_orientation_transformation_matrix",
    "solvent",
    "temperature",
    "pressure",
    "job_type",
    "functional",
    "basis_set",
)


class _HasG16ParserState(Protocol):
    _frame_parser: Any
    _chem_file: Any
    _file_path: str | None
    only_last_frame: bool
    forced_charge: int | None
    forced_multiplicity: int | None


class G16LogFileParserMixin:
    @staticmethod
    def _has_metadata_value(value: Any) -> bool:
        if value is None:
            return False
        if isinstance(value, str):
            return value.strip() != ""
        return True

    def _filter_section_metadata(self, metadata: dict[str, Any] | None) -> dict[str, Any]:
        if metadata is None:
            return {}
        return {
            key: metadata[key]
            for key in SECTION_METADATA_KEYS
            if key in metadata and self._has_metadata_value(metadata[key])
        }

    def _split_sections(self, file_content: str) -> list[str]:
        matches = list(LINK1_SECTION_PATTERN.finditer(file_content))
        if not matches:
            return [file_content]
        sections: list[str] = []
        first_start = matches[0].start()
        if first_start > 0:
            sections.append(file_content[:first_start])
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

    def _resolve_frame_additional_data(
        self,
        *,
        section_metadata: dict[str, Any],
        previous_metadata: dict[str, Any] | None,
    ) -> dict[str, Any]:
        resolved: dict[str, Any] = {}
        for key in SECTION_METADATA_KEYS:
            if key in section_metadata and self._has_metadata_value(section_metadata[key]):
                resolved[key] = section_metadata[key]
            elif (
                previous_metadata is not None
                and key in previous_metadata
                and self._has_metadata_value(previous_metadata[key])
            ):
                resolved[key] = previous_metadata[key]
        return resolved

    def _derive_file_metadata_from_frames(self, frames: Sequence[Any]) -> dict[str, Any]:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
        if not frames:
            return metadata

        first_frame = frames[0]
        last_frame = frames[-1]

        def first_nonempty(key: str) -> Any:
            for frame in frames:
                value = getattr(frame, key, None)
                if self._has_metadata_value(value):
                    return value
            return None

        for key in SECTION_METADATA_KEYS:
            value = getattr(first_frame, key, None)
            if self._has_metadata_value(value):
                metadata[key] = value

        for key in ("running_time", "status"):
            value = getattr(last_frame, key, None)
            if self._has_metadata_value(value):
                metadata[key] = value

        for key in ("solvent", "temperature", "pressure"):
            value = first_nonempty(key)
            if self._has_metadata_value(value):
                metadata[key] = value

        return metadata

    def _parse_metadata(self, file_content: str) -> dict[str, Any] | None:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
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
        temperature, pressure = self._parse_temperature_and_pressure()
        if temperature:
            metadata["temperature"] = temperature
        if pressure:
            metadata["pressure"] = pressure
        if tail := self._parse_tail():
            metadata.update(tail)
        if running_time := self._parse_running_time():
            metadata["running_time"] = running_time
        if status := self._parse_termination_status():
            metadata["status"] = status
        return metadata

    def _parse_section_metadata(self, section_content: str) -> dict[str, Any]:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
        self._file_content = section_content
        if version := self._parse_version():
            metadata["qm_software_version"] = version
        if options := self._parse_options():
            metadata["options"] = options
        if route := self._parse_keywords():
            metadata["keywords"] = route
        if title := self._parse_title():
            metadata["title_card"] = title
        charge_multiplicity = self._parse_charge_multiplicity()
        if charge_multiplicity != (None, None):
            metadata["charge"], metadata["multiplicity"] = charge_multiplicity
        return metadata

    def _split_file(self, file_content: str) -> Sequence[str]:
        frame_contents: list[str] = []
        for section in self._split_sections(file_content):
            frame_contents.extend(self._split_section_frames(section))
        return frame_contents

    def _parse_frame(self, frame_content: str, *, additional_data: dict[str, Any] | None = None):
        typed_self = cast(_HasG16ParserState, self)
        frame_parser = typed_self._frame_parser()
        return cast(Any, frame_parser.parse(frame_content, additional_data=additional_data))

    def _parse(
        self,
        file_content: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
    ):
        typed_self = cast(_HasG16ParserState, self)

        frame_entries: list[tuple[str, dict[str, Any]]] = []
        for section in self._split_sections(file_content):
            section_metadata = self._filter_section_metadata(self._parse_section_metadata(section))
            for frame_content in self._split_section_frames(section):
                frame_entries.append((frame_content, section_metadata))

        if typed_self.only_last_frame and frame_entries:
            frame_entries = [frame_entries[-1]]

        final_charge = total_charge or typed_self.forced_charge
        final_multiplicity = total_multiplicity or typed_self.forced_multiplicity

        parsed_frames: list[Any] = []
        previous_metadata: dict[str, Any] | None = None

        for frame_content, section_metadata in frame_entries:
            resolved_additional_data = self._resolve_frame_additional_data(
                section_metadata=section_metadata,
                previous_metadata=previous_metadata,
            )
            frame = self._parse_frame(frame_content, additional_data=resolved_additional_data)
            previous_metadata = {
                key: getattr(frame, key)
                for key in SECTION_METADATA_KEYS
                if hasattr(frame, key) and self._has_metadata_value(getattr(frame, key))
            }
            if final_charge:
                frame.charge = final_charge
            if final_multiplicity:
                frame.multiplicity = final_multiplicity
            parsed_frames.append(frame)

        metadata: dict[str, Any] = {
            "file_path": typed_self._file_path,
            "file_content": file_content,
        }
        metadata.update(self._derive_file_metadata_from_frames(parsed_frames))
        if final_charge:
            metadata["charge"] = final_charge
        if final_multiplicity:
            metadata["multiplicity"] = final_multiplicity

        chem_file = typed_self._chem_file.model_validate(metadata)
        for frame in parsed_frames:
            chem_file.append(frame)
        return cast(Any, chem_file)

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
    ) -> tuple[PlainQuantity, PlainQuantity] | tuple[None, None]:
        index = self._file_content.find(" - Thermochemistry -")
        if index == -1:
            return None, None
        focus_content = self._file_content[index:]
        # self._file_content = self._file_content[index:]
        if match := g16_log_patterns.TEMPEREATURE_PRESSURE.match_content(focus_content):
            temperature, pressure = (
                float(match[0][0]) * atom_ureg.K,
                float(match[0][1]) * atom_ureg.atm,
            )
            return temperature, pressure
        return None, None

    def _parse_tail(self) -> dict[str, Any]:
        tail_dict, remaining_content = parse_archive_tail(self._file_content, include_coords=False)
        self._file_content = remaining_content
        return tail_dict

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
    BaseFileParserMemory[G16LogFileMemory, G16LogFileFrameMemory, G16LogFileFrameParserV2Memory],
):
    _chem_file = G16LogFileMemory
    _frame_parser = G16LogFileFrameParserV2Memory


class G16LogFileParserDisk(
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

    extensions = frozenset(extensions_for_parser(G16LogFileParserDisk))
    priority = 100

    @registry.reader_factory(format_id="g16log", extensions=extensions, priority=priority)
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id="g16log",
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=G16LogFileParserDisk,
                priority=priority,
            ),
        )
