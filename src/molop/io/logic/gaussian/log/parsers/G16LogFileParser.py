"""
Author: TMJ
Date: 2025-08-01 16:13:58
LastEditors: TMJ
LastEditTime: 2026-05-11 14:17:37
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, Literal, Protocol, cast

from molop.io.base_models.FileParser import (
    BaseFileParserDisk,
    BaseFileParserMemory,
    _HasFileParseMethod,
)
from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.logic.gaussian.input.GaussianRouteParsing import parse_gaussian_route_semantic
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.gaussian.log.frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserDisk,
    G16LogFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.models.G16LogFile import (
    BaseCalcFile,
    G16LogFileDisk,
    G16LogFileMemory,
)
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import (
    ensure_g16_output_content,
    extract_g16_charge_multiplicity,
    extract_g16_keywords,
    extract_g16_options,
    extract_g16_running_time,
    extract_g16_solvent,
    extract_g16_standard_orientation_transformation_matrix,
    extract_g16_temperature_and_pressure,
    extract_g16_termination_status,
    extract_g16_title,
    extract_g16_version,
    first_frame_value,
    last_frame_value,
    split_g16_section_frames,
    split_g16_sections,
)
from molop.io.logic.gaussian.log.parsers._g16log_archive_tail import parse_archive_tail


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


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


class G16MetadataParsePhase(Enum):
    """Explicit stages for Gaussian output file metadata parsing."""

    ROUTE = auto()
    STRUCTURE = auto()
    CONDITIONS = auto()
    ARCHIVE = auto()
    TIMING_STATUS = auto()
    DONE = auto()


class G16LogFileParserMixin:
    @classmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        ensure_g16_output_content(file_content)

    def _run_route_metadata_phase(
        self, context: TextParseContext, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        if version := extract_g16_version(context):
            result.set("qm_software_version", version)
        if options := extract_g16_options(context):
            result.set("options", options)
        if route := extract_g16_keywords(context):
            result.set("keywords", route)
            result.set("semantic_route", parse_gaussian_route_semantic(route))
        if title := extract_g16_title(context):
            result.set("title_card", title)
        if charge_multiplicity := extract_g16_charge_multiplicity(context):
            charge, multiplicity = charge_multiplicity
            result.set("charge", charge)
            result.set("multiplicity", multiplicity)
        return G16MetadataParsePhase.STRUCTURE

    def _run_structure_metadata_phase(
        self, context: TextParseContext, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        transform_matrix = extract_g16_standard_orientation_transformation_matrix(context)
        if transform_matrix is not None:
            result.set("standard_orientation_transformation_matrix", transform_matrix)
        if cast(_HasFileParseMethod, self).only_extract_structure:
            return G16MetadataParsePhase.DONE
        return G16MetadataParsePhase.CONDITIONS

    def _run_conditions_metadata_phase(
        self, file_content: str, context: TextParseContext, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        if solvent := extract_g16_solvent(context):
            result.set("solvent", solvent)
        temperature, pressure = extract_g16_temperature_and_pressure(file_content)
        if temperature:
            result.set("temperature", temperature)
        if pressure:
            result.set("pressure", pressure)
        return G16MetadataParsePhase.ARCHIVE

    def _run_archive_metadata_phase(
        self, file_content: str, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        if tail := parse_archive_tail(file_content):
            result.update(tail[0])
        return G16MetadataParsePhase.TIMING_STATUS

    def _run_timing_status_metadata_phase(
        self, context: TextParseContext, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        if running_time := extract_g16_running_time(context):
            result.set("running_time", running_time)
        if status := extract_g16_termination_status(context):
            result.set("status", status)
        return G16MetadataParsePhase.DONE

    def _parse_metadata_result(self, file_content: str) -> ModelParseResult:
        context = TextParseContext(file_content)
        result = ModelParseResult({"qm_software": "Gaussian"})
        phase = G16MetadataParsePhase.ROUTE
        while phase is not G16MetadataParsePhase.DONE:
            if phase is G16MetadataParsePhase.ROUTE:
                phase = self._run_route_metadata_phase(context, result)
            elif phase is G16MetadataParsePhase.STRUCTURE:
                phase = self._run_structure_metadata_phase(context, result)
            elif phase is G16MetadataParsePhase.CONDITIONS:
                phase = self._run_conditions_metadata_phase(file_content, context, result)
            elif phase is G16MetadataParsePhase.ARCHIVE:
                phase = self._run_archive_metadata_phase(file_content, result)
            elif phase is G16MetadataParsePhase.TIMING_STATUS:
                phase = self._run_timing_status_metadata_phase(context, result)
            else:
                raise AssertionError(f"Unexpected G16 metadata parse phase: {phase!r}")
        return result

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        return self._parse_metadata_result(file_content).model_data()

    def _split_sections(self, file_content: str) -> list[str]:
        return split_g16_sections(file_content)

    def _split_section_frames(self, section_content: str) -> list[str]:
        return split_g16_section_frames(section_content)

    def _split_file(self, file_content: str) -> Sequence[str]:
        sections = self._split_sections(file_content)
        return [frame for section in sections for frame in self._split_section_frames(section)]

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
            value = first_frame_value(frames, field)
            if value is not None:
                metadata[field] = value
                setattr(chem_file, field, value)

        for field in ("status", "geometry_optimization_status"):
            value = last_frame_value(frames, field)
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
