"""
Author: TMJ
Date: 2025-08-01 16:13:58
LastEditors: TMJ
LastEditTime: 2026-07-14 14:34:59
Description: 请填写简介
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from enum import Enum, auto
from typing import TYPE_CHECKING, Any, ClassVar, Protocol, cast

from molop.io.base_models.FileParser import (
    BaseFileParserDisk,
    BaseFileParserMemory,
    _HasFileParseMethod,
)
from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.base_models.source import LocatedSourceSegment
from molop.io.logic.gaussian.input.GaussianRoute import (
    _with_dispersion_suffix,
    build_gaussian_model_chemistry,
    build_gaussian_task_requests,
)
from molop.io.logic.gaussian.input.GaussianRouteParsing import parse_gaussian_route_semantic
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.gaussian.log.frame_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserDisk,
    G16LogFileFrameParserMemory,
)
from molop.io.logic.gaussian.log.locators import (
    locate_g16_section_frames,
    locate_g16_sections,
)
from molop.io.logic.gaussian.log.models.G16LogFile import G16LogFileDisk, G16LogFileMemory
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import (
    ensure_g16_output_content,
    extract_g16_artifact_version,
    extract_g16_atomic_masses,
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
)
from molop.io.logic.gaussian.log.parsers._g16log_archive_tail import parse_archive_tail


if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


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
    format_id: ClassVar[str] = "g16log"
    _parse_unselected_segment_metadata = True
    _assess_segment_calculation_status = True
    _calculation_configuration_fields: ClassVar[tuple[str, ...]] = (
        "options",
        "method",
        "basis_set",
        "functional",
        "solvent",
        "resources_raw",
        "request_num_cpu",
        "request_memory",
        "resource_request",
    )

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
            archive_metadata = tail[0]
            result.update(archive_metadata)
            if archive_keywords := archive_metadata.get("keywords"):
                result.set(
                    "semantic_route",
                    parse_gaussian_route_semantic(archive_keywords),
                )
        return G16MetadataParsePhase.TIMING_STATUS

    def _run_timing_status_metadata_phase(
        self, segment_content: str, result: ModelParseResult
    ) -> G16MetadataParsePhase:
        context = TextParseContext(segment_content)
        if running_time := extract_g16_running_time(context):
            result.set("running_time", running_time)
        if status := extract_g16_termination_status(context):
            result.set("status", status)
        return G16MetadataParsePhase.DONE

    def _parse_segment_metadata_result(self, segment_content: str) -> ModelParseResult:
        context = TextParseContext(segment_content)
        result = ModelParseResult({"qm_software": "Gaussian"})
        phase = G16MetadataParsePhase.ROUTE
        while phase is not G16MetadataParsePhase.DONE:
            if phase is G16MetadataParsePhase.ROUTE:
                phase = self._run_route_metadata_phase(context, result)
            elif phase is G16MetadataParsePhase.STRUCTURE:
                phase = self._run_structure_metadata_phase(context, result)
            elif phase is G16MetadataParsePhase.CONDITIONS:
                phase = self._run_conditions_metadata_phase(segment_content, context, result)
            elif phase is G16MetadataParsePhase.ARCHIVE:
                phase = self._run_archive_metadata_phase(segment_content, result)
            elif phase is G16MetadataParsePhase.TIMING_STATUS:
                phase = self._run_timing_status_metadata_phase(segment_content, result)
            else:
                raise AssertionError(f"Unexpected G16 metadata parse phase: {phase!r}")
        return result

    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        return tuple(
            LocatedSourceSegment(
                segment=section,
                frames=locate_g16_section_frames(file_content, section),
            )
            for section in locate_g16_sections(file_content)
        )

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any]:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
        if version := extract_g16_artifact_version(file_content):
            metadata["qm_software_version"] = version
        if cast(_HasFileParseMethod, self).only_last_frame:
            atomic_masses, atomic_masses_source = extract_g16_atomic_masses(file_content)
            if atomic_masses is not None:
                metadata["atomic_masses"] = atomic_masses
                metadata["atomic_masses_source"] = atomic_masses_source
        return metadata

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        _ = artifact_metadata
        metadata = self._parse_segment_metadata_result(segment_content).model_data()
        metadata.pop("qm_software", None)
        metadata.pop("qm_software_version", None)
        semantic_route = metadata.get("semantic_route")
        if semantic_route is not None:
            metadata["model_chemistry"] = build_gaussian_model_chemistry(
                semantic_route,
                keywords=str(metadata.get("keywords") or ""),
                legacy_method=str(metadata.get("method") or ""),
                legacy_basis_set=str(metadata.get("basis_set") or ""),
                legacy_functional=str(metadata.get("functional") or ""),
            )
            metadata["task_requests"] = build_gaussian_task_requests(semantic_route)
        atomic_masses, atomic_masses_source = extract_g16_atomic_masses(segment_content)
        if atomic_masses is not None:
            metadata["atomic_masses"] = atomic_masses
            metadata["atomic_masses_source"] = atomic_masses_source
        return {key: value for key, value in metadata.items() if value is not None}

    @staticmethod
    def _is_missing_configuration_value(value: Any) -> bool:
        if value is None:
            return True
        if isinstance(value, str):
            return not value.strip()
        if isinstance(value, (list, tuple, set, dict)):
            return not value
        return False

    @classmethod
    def _inherit_missing_model_fields(cls, target: Any, source: Any) -> None:
        if target is None or source is None:
            return
        for field_name in getattr(type(target), "model_fields", {}):
            target_value = getattr(target, field_name, None)
            source_value = getattr(source, field_name, None)
            if cls._is_missing_configuration_value(target_value):
                setattr(target, field_name, deepcopy(source_value))
            elif isinstance(target_value, dict) and isinstance(source_value, Mapping):
                merged = dict(deepcopy(source_value))
                for key, value in target_value.items():
                    if key not in merged or not cls._is_missing_configuration_value(merged[key]):
                        merged[key] = value
                setattr(target, field_name, merged)

    @staticmethod
    def _is_frequency_only_segment(metadata: Mapping[str, Any]) -> bool:
        requests = metadata.get("task_requests")
        if not isinstance(requests, Sequence) or isinstance(requests, str):
            return False
        task_types: set[str] = set()
        for request in requests:
            if isinstance(request, Mapping):
                enabled = request.get("enabled", True)
                task_type = request.get("task_type")
            else:
                enabled = getattr(request, "enabled", True)
                task_type = getattr(request, "task_type", None)
            if enabled and isinstance(task_type, str):
                task_types.add(task_type)
        return task_types == {"freq"}

    @staticmethod
    def _is_scrf_checkpoint_route(route: Any) -> bool:
        option_maps = getattr(route, "option_maps", {})
        scrf_option = option_maps.get("scrf") if isinstance(option_maps, Mapping) else None
        scalar_value = getattr(scrf_option, "scalar_value", None)
        if isinstance(scalar_value, str) and scalar_value.lower() in {"check", "restart"}:
            return True
        solvation_model = getattr(route, "solvation_model", None)
        return isinstance(solvation_model, str) and solvation_model.replace(" ", "").lower() in {
            "scrf=check",
            "scrf=restart",
        }

    @classmethod
    def _inherit_frequency_segment_configuration(
        cls,
        current: dict[str, Any],
        previous: Mapping[str, Any],
    ) -> None:
        for field_name in cls._calculation_configuration_fields:
            current_value = current.get(field_name)
            previous_value = previous.get(field_name)
            if cls._is_missing_configuration_value(current_value) and previous_value is not None:
                current[field_name] = deepcopy(previous_value)

        current_route = current.get("semantic_route")
        previous_route = previous.get("semantic_route")
        if current_route is None or previous_route is None:
            if cls._is_missing_configuration_value(current.get("model_chemistry")):
                current["model_chemistry"] = deepcopy(previous.get("model_chemistry"))
            return

        merged_route = current_route.model_copy(deep=True)
        cls._inherit_missing_model_fields(
            merged_route.model_chemistry,
            getattr(previous_route, "model_chemistry", None),
        )

        if cls._is_missing_configuration_value(merged_route.empirical_dispersion):
            merged_route.empirical_dispersion = getattr(
                previous_route,
                "empirical_dispersion",
                None,
            )

        previous_scrf = getattr(previous_route, "scrf_options", None)
        if previous_scrf is not None and getattr(previous_scrf, "enabled", False):
            if cls._is_scrf_checkpoint_route(merged_route) or not getattr(
                merged_route.scrf_options, "enabled", False
            ):
                merged_route.scrf_options = previous_scrf.model_copy(deep=True)
                merged_route.solvation_model = getattr(
                    previous_route,
                    "solvation_model",
                    merged_route.solvation_model,
                )
            else:
                cls._inherit_missing_model_fields(merged_route.scrf_options, previous_scrf)

        previous_route_modifiers = getattr(previous_route, "route_modifiers", [])
        for modifier in previous_route_modifiers:
            if (
                modifier in {"em", "empiricaldispersion"}
                and modifier not in merged_route.route_modifiers
            ):
                merged_route.route_modifiers.append(modifier)

        previous_capabilities = getattr(previous_route, "capabilities", [])
        for capability in ("Dispersion", "Solvation"):
            if capability in previous_capabilities and capability not in merged_route.capabilities:
                merged_route.capabilities.append(capability)

        previous_option_maps = getattr(previous_route, "option_maps", {})
        for option_name in ("em", "empiricaldispersion"):
            if option_name not in merged_route.option_maps and option_name in previous_option_maps:
                merged_route.option_maps[option_name] = deepcopy(previous_option_maps[option_name])

        current["semantic_route"] = merged_route
        current["model_chemistry"] = build_gaussian_model_chemistry(
            merged_route,
            keywords=str(current.get("keywords") or ""),
            legacy_method=str(current.get("method") or ""),
            legacy_basis_set=str(current.get("basis_set") or ""),
            legacy_functional=str(current.get("functional") or ""),
        )
        current["task_requests"] = build_gaussian_task_requests(merged_route)

    @classmethod
    def _supplement_frame_calculation_configuration(
        cls,
        frame: Any,
        section_metadata: Mapping[str, Any],
    ) -> None:
        """Complete frame configuration from the resolved metadata of its section."""

        for field_name in cls._calculation_configuration_fields:
            section_value = section_metadata.get(field_name)
            if cls._is_missing_configuration_value(section_value):
                continue
            if cls._is_missing_configuration_value(getattr(frame, field_name, None)):
                setattr(frame, field_name, deepcopy(section_value))

        section_route = section_metadata.get("semantic_route")
        frame_route = getattr(frame, "semantic_route", None)
        if section_route is not None and frame_route is None:
            frame.semantic_route = deepcopy(section_route)
            frame_route = frame.semantic_route
        if section_route is not None and frame_route is not None:
            cls._inherit_missing_model_fields(frame_route, section_route)
            cls._inherit_missing_model_fields(
                getattr(frame_route, "model_chemistry", None),
                getattr(section_route, "model_chemistry", None),
            )
            if cls._is_missing_configuration_value(
                getattr(frame_route, "empirical_dispersion", None)
            ):
                frame_route.empirical_dispersion = deepcopy(
                    getattr(section_route, "empirical_dispersion", None)
                )
            cls._inherit_missing_model_fields(
                getattr(frame_route, "scrf_options", None),
                getattr(section_route, "scrf_options", None),
            )

        section_model = section_metadata.get("model_chemistry")
        frame_model = getattr(frame, "model_chemistry", None)
        if section_model is not None and frame_model is None:
            frame.model_chemistry = deepcopy(section_model)
            frame_model = frame.model_chemistry
        elif section_model is not None and frame_model is not None:
            cls._inherit_missing_model_fields(frame_model, section_model)

        if frame_model is None:
            return

        section_route_model = getattr(section_route, "model_chemistry", None)
        section_functional = getattr(section_model, "functional", None)
        if section_functional is None:
            section_functional = getattr(section_route_model, "functional", None)
        section_dispersion = getattr(section_model, "dispersion_correction", None)
        if section_dispersion is None:
            section_dispersion = getattr(section_route, "empirical_dispersion", None)
        canonical_section_functional = _with_dispersion_suffix(
            section_functional,
            section_dispersion,
        )
        frame_functional = getattr(frame_model, "functional", None) or getattr(
            frame,
            "functional",
            None,
        )
        canonical_frame_functional = _with_dispersion_suffix(
            frame_functional,
            section_dispersion,
        )
        if canonical_section_functional and (
            cls._is_missing_configuration_value(frame_functional)
            or (
                canonical_frame_functional is not None
                and canonical_frame_functional.upper() == canonical_section_functional.upper()
            )
        ):
            frame_model.functional = canonical_section_functional.upper()

        backfill = getattr(frame, "backfill_common_qm_containers_from_legacy", None)
        if callable(backfill):
            backfill()
        project = getattr(frame, "project_common_qm_fields", None)
        if callable(project):
            project()

    def _postprocess_segment_metadata(
        self,
        segment_metadata: Sequence[Mapping[str, Any]],
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> Sequence[Mapping[str, Any]]:
        _ = artifact_metadata
        resolved_metadata = [dict(metadata) for metadata in segment_metadata]
        for segment_index in range(1, len(resolved_metadata)):
            current = resolved_metadata[segment_index]
            if not self._is_frequency_only_segment(current):
                continue
            self._inherit_frequency_segment_configuration(
                current,
                resolved_metadata[segment_index - 1],
            )
        return resolved_metadata

    def _prepare_file_metadata(
        self,
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Sequence[Mapping[str, Any]],
    ) -> dict[str, Any]:
        metadata = dict(artifact_metadata)
        if segment_metadata:
            metadata.update(segment_metadata[0])
        metadata.pop("atomic_masses", None)
        metadata.pop("atomic_masses_source", None)
        running_times = [
            running_time
            for segment in segment_metadata
            if (running_time := segment.get("running_time")) is not None
        ]
        if running_times:
            metadata["running_time"] = sum(running_times[1:], running_times[0])
        for segment in reversed(segment_metadata):
            if (status := segment.get("status")) is not None:
                metadata["status"] = status
                break
        return metadata

    def _postprocess_parsed_frame(
        self,
        frame: G16LogFileFrameDisk | G16LogFileFrameMemory,
        chem_file: Any,
        *,
        segment_index: int | None,
        segment_frame_index: int,
        segment_frame_count: int,
    ) -> None:
        _ = segment_index, segment_frame_index, segment_frame_count
        if (
            frame.basis_set.lower() == "gen"
            and chem_file
            and chem_file[-1].basis_set.lower() != "gen"
        ):
            frame.basis_set = chem_file[-1].basis_set
            frame.model_chemistry.basis_set = frame.basis_set

    def _postprocess_frame_configuration(
        self,
        frame: G16LogFileFrameDisk | G16LogFileFrameMemory,
        segment_metadata: Mapping[str, Any],
        *,
        segment_index: int | None,
        segment_frame_index: int,
        segment_frame_count: int,
    ) -> None:
        _ = segment_index, segment_frame_index, segment_frame_count
        self._supplement_frame_calculation_configuration(frame, segment_metadata)

    def _source_frame_role(
        self,
        frame: G16LogFileFrameDisk | G16LogFileFrameMemory,
        *,
        frame_index: int,
        frame_count: int,
    ) -> str:
        _ = self
        task_types = {task.task_type for task in frame.task_requests if task.enabled}
        if "sp" in task_types and "opt" not in task_types and frame_count == 1:
            return "single_point"
        if frame_count == 1 or frame_index == frame_count - 1:
            return "terminal"
        if frame_index == 0:
            return "initial"
        return "intermediate"

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

        if metadata.get("status") is None:
            value = last_frame_value(frames, "status")
            if value is not None:
                metadata["status"] = value
                chem_file.status = value

        for field in ("geometry_optimization_status",):
            value = last_frame_value(frames, field)
            if value is not None:
                metadata[field] = value
                setattr(chem_file, field, value)

        mass_source_priority = {
            "gaussian_thermochemistry": 1,
            "gaussian_atmwgt": 2,
        }
        mass_references: dict[tuple[int, ...], tuple[Any, str | None]] = {}
        ambiguous_mass_keys: set[tuple[int, ...]] = set()
        for frame in frames:
            if frame.atomic_masses is None:
                continue
            atom_key = tuple(frame.atoms)
            reference = mass_references.get(atom_key)
            if reference is None:
                mass_references[atom_key] = (
                    frame.atomic_masses,
                    frame.atomic_masses_source,
                )
                continue
            reference_source = reference[1]
            current_source = frame.atomic_masses_source
            reference_priority = mass_source_priority.get(reference_source or "", 0)
            current_priority = mass_source_priority.get(current_source or "", 0)
            if current_priority > reference_priority:
                mass_references[atom_key] = (
                    frame.atomic_masses,
                    current_source,
                )
                ambiguous_mass_keys.discard(atom_key)
            elif current_priority == reference_priority and (
                reference[0].m_as("amu").tolist() != frame.atomic_masses.m_as("amu").tolist()
            ):
                ambiguous_mass_keys.add(atom_key)

        for frame in frames:
            if frame.atomic_masses is not None:
                continue
            atom_key = tuple(frame.atoms)
            if atom_key in ambiguous_mass_keys or atom_key not in mass_references:
                continue
            atomic_masses, atomic_masses_source = mass_references[atom_key]
            frame.atomic_masses = atomic_masses
            frame.atomic_masses_source = atomic_masses_source


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

    @registry.reader_factory(
        format_id=G16LogFileParserDisk.format_id,
        extensions=extensions,
        priority=priority,
    )
    def _factory() -> ReaderCodec:
        return cast(
            ReaderCodec,
            ParserDiskReader(
                format_id=G16LogFileParserDisk.format_id,
                extensions=extensions,
                level=StructureLevel.COORDS,
                parser_cls=G16LogFileParserDisk,
                priority=priority,
                memory_parser_cls=G16LogFileParserMemory,
            ),
        )
