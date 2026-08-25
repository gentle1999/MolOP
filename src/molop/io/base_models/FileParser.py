"""
Author: TMJ
Date: 2025-07-29 22:14:37
LastEditors: TMJ
LastEditTime: 2026-05-11 14:36:16
Description: 请填写简介
"""

import os
from abc import abstractmethod
from collections.abc import Mapping, Sequence
from dataclasses import asdict, replace
from enum import Enum
from hashlib import sha256
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any, ClassVar, Generic, Literal, Protocol, TypeVar, cast

from molgr.config import CONFIG as MOLGR_CONFIG
from pydantic import Field, PrivateAttr

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.base_models.ChemFileFrame import BaseChemFileFrame
from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.base_models.source import (
    DecodedSource,
    LocatedSourceSegment,
    LocatedTextBlock,
    ParseCompleteness,
    ParseDiagnostic,
    ParsePresence,
    ParserProvenance,
    SourceSegmentEvidence,
    SourceSpan,
    canonical_json_sha256,
    normalize_parser_line_endings,
)
from molop.io.base_models.summary import SummaryDict, summary_column
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.codec_types import ParseOptions


FileT = TypeVar("FileT", bound=BaseChemFile)
FrameT = TypeVar("FrameT", bound=BaseChemFileFrame)
FrameParserT = TypeVar("FrameParserT", bound=BaseFrameParser[Any])


class _HasFileParseMethod(Protocol):
    _file_content: str
    only_extract_structure: bool


class BaseFileParser(BaseDataClassWithUnit, Generic[FileT, FrameT, FrameParserT]):
    # add by subclass
    _frame_parser: type[FrameParserT] = PrivateAttr()
    _chem_file: type[FileT] = PrivateAttr()
    _file_path: str | None = PrivateAttr(default=None)
    _segment_scoped_frame_metadata: ClassVar[frozenset[str]] = frozenset({"status", "running_time"})
    _assess_segment_calculation_status: ClassVar[bool] = False
    format_id: ClassVar[str] = ""
    parser_id: ClassVar[str | None] = None
    parser_version: ClassVar[str | None] = None

    forced_charge: int | None = Field(
        default=None,
        description="The forced charge of the molecule.",
        exclude=True,
        repr=False,
    )
    forced_multiplicity: int | None = Field(
        default=None,
        description="The forced multiplicity of the molecule.",
        ge=1,
        exclude=True,
        repr=False,
    )
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    only_last_frame: bool = Field(default=False, exclude=True, repr=False)
    capture_source_evidence: bool = Field(default=False, exclude=True, repr=False)
    source_encoding: str = Field(default="utf-8", min_length=1, exclude=True, repr=False)
    parse_options: ParseOptions | None = Field(default=None, exclude=True, repr=False)

    @staticmethod
    def _distribution_version(distribution: str) -> str:
        try:
            return distribution_version(distribution)
        except PackageNotFoundError:
            return "unknown"

    def _snapshot_parser_provenance(
        self,
        *,
        charge_override: int | None,
        multiplicity_override: int | None,
        source_encoding: str,
        parse_options: ParseOptions | None = None,
    ) -> ParserProvenance:
        """Capture versions and all parser-relevant global configuration by value."""

        options = (parse_options or ParseOptions()).resolved()
        molgr_config = asdict(MOLGR_CONFIG)
        molgr_config["interface"]["reconstruction_failure_policy"] = (
            options.reconstruction_failure_policy
        )

        effective_config: dict[str, Any] = {
            "parser": {
                "total_charge_override": charge_override,
                "total_multiplicity_override": multiplicity_override,
                "only_extract_structure": options.only_extract_structure,
                "only_last_frame": options.only_last_frame,
                "capture_source_evidence": options.capture_source_evidence,
                "source_encoding": source_encoding,
            },
            "molop": {
                "force_unit_transform": options.force_unit_transform,
                "graph_reconstruction_backend": options.graph_reconstruction_backend,
                "reconstruction_failure_policy": options.reconstruction_failure_policy,
                "make_dative_bonds": options.make_dative_bonds,
                "make_stereochemistry": options.make_stereochemistry,
            },
            "molgr": molgr_config,
        }
        molop_version = self._distribution_version("molop")
        parser_class = type(self)
        parser_id = self.parser_id or (f"{parser_class.__module__}.{parser_class.__qualname__}")
        parser_version = self.parser_version or molop_version
        return ParserProvenance(
            parser_id=parser_id,
            parser_version=parser_version,
            molop_version=molop_version,
            molgr_version=self._distribution_version("molgr"),
            rdkit_version=self._distribution_version("rdkit"),
            effective_config=effective_config,
            effective_config_sha256=canonical_json_sha256(effective_config),
        )

    def _load_source(
        self,
        source: str,
        source_type: Literal["file_path", "string"],
    ) -> DecodedSource:
        """Load and strictly decode one parser source without losing byte offsets."""

        if source_type == "file_path":
            self._file_path = source
            return DecodedSource.from_bytes(Path(source).read_bytes(), self.source_encoding)
        if source_type == "string":
            self._file_path = None
            return DecodedSource.from_text(source, self.source_encoding)
        raise ValueError(f"Invalid source_type: {source_type}")

    @staticmethod
    def _source_artifact_identity(source: DecodedSource) -> dict[str, str | int]:
        """Return format-independent identity fields for one decoded artifact."""

        return {
            "artifact_sha256": sha256(source.raw_bytes).hexdigest(),
            "artifact_size_bytes": len(source.raw_bytes),
            "source_encoding": source.encoding,
        }

    @staticmethod
    def _source_block_identity(
        source: DecodedSource,
        block: LocatedTextBlock,
    ) -> tuple[SourceSpan, str]:
        """Return the exact span and raw-byte digest for one located text block."""

        source.prepare_blocks((block,))
        span = source.span(block)
        return span, source.block_sha256(span)

    @abstractmethod
    def _locate_segments(self, file_content: str) -> Sequence[LocatedSourceSegment]:
        """Return every exact source segment and its exact parser frame blocks."""

        raise NotImplementedError

    def _parse_artifact_metadata(self, file_content: str) -> dict[str, Any] | None:
        """Parse metadata that applies to every segment in the artifact."""

        _ = file_content
        return None

    def _parse_segment_metadata(
        self,
        segment_content: str,
        *,
        artifact_metadata: Mapping[str, Any],
    ) -> dict[str, Any] | None:
        """Parse metadata local to one located source segment."""

        _ = segment_content, artifact_metadata
        return None

    def _validated_located_segments(
        self,
        file_content: str,
    ) -> tuple[LocatedSourceSegment, ...]:
        """Validate the format locator contract before any source slice is parsed."""

        located = self._locate_segments(file_content)
        if located is None:
            raise TypeError(
                f"{type(self).__name__}._locate_segments() must return a non-empty sequence, "
                "not None"
            )
        segments = tuple(located)
        if not segments:
            raise ValueError("The file content contains no locator-provided source segments.")

        previous_segment_end = 0
        for segment_index, located_segment in enumerate(segments):
            if not isinstance(located_segment, LocatedSourceSegment):
                raise TypeError(
                    f"{type(self).__name__}._locate_segments() item {segment_index} must be "
                    "LocatedSourceSegment"
                )
            segment = located_segment.segment
            if segment.end_char > len(file_content):
                raise ValueError(f"Located source segment {segment_index} exceeds the source text")
            if segment_index and segment.start_char < previous_segment_end:
                raise ValueError("Located source segments must be ordered and non-overlapping")
            previous_segment_end = segment.end_char

            previous_frame_end = segment.start_char
            for frame_index, frame in enumerate(located_segment.frames):
                if frame.end_char > len(file_content):
                    raise ValueError(
                        f"Located frame {segment_index}:{frame_index} exceeds the source text"
                    )
                if frame.start_char < previous_frame_end:
                    raise ValueError(
                        "Located source frames must be ordered and non-overlapping within a segment"
                    )
                if file_content[previous_frame_end : frame.start_char].strip():
                    raise ValueError(
                        "Located source frames must cover every non-whitespace character "
                        "in a segment that contains frames"
                    )
                previous_frame_end = frame.end_char
            if (
                located_segment.frames
                and file_content[previous_frame_end : segment.end_char].strip()
            ):
                raise ValueError(
                    "Located source frames must cover every non-whitespace character "
                    "in a segment that contains frames"
                )
        return segments

    def _prepare_file_metadata(
        self,
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Sequence[Mapping[str, Any]],
    ) -> dict[str, Any]:
        """Compose model metadata before parsing the selected frames."""

        metadata = dict(artifact_metadata)
        if segment_metadata:
            metadata.update(segment_metadata[0])
        return metadata

    def _postprocess_parsed_frame(
        self,
        frame: FrameT,
        chem_file: FileT,
        *,
        segment_index: int | None,
        segment_frame_index: int,
        segment_frame_count: int,
    ) -> None:
        """Apply format-owned adjustments before a parsed frame is appended."""

        _ = frame, chem_file, segment_index, segment_frame_index, segment_frame_count

    def _source_frame_role(
        self,
        frame: FrameT,
        *,
        frame_index: int,
        frame_count: int,
    ) -> str | None:
        """Return an optional format-owned role for one located frame."""

        _ = frame, frame_index, frame_count
        return None

    def _source_frame_fields(
        self,
        frame: FrameT,
        *,
        segment_index: int,
        segment_frame_index: int,
        segment_frame_count: int,
    ) -> Mapping[str, Any]:
        """Return optional model fields derived from a located frame's source context.

        This is part of the source-evidence lifecycle, so it must remain a
        source-only operation.  In particular, do not read ``frame.rdmol``
        here: coordinate-only calculation frames expose that property lazily,
        and reading it would reconstruct a molecular graph while evidence is
        being attached.  Topology is intentionally left unassessed until a
        graph-dependent operation explicitly requests it.
        """

        role = self._source_frame_role(
            frame,
            frame_index=segment_frame_index,
            frame_count=segment_frame_count,
        )
        fields: dict[str, Any] = {} if role is None else {"frame_role": role}
        if not hasattr(frame, "geometry_optimization_status"):
            return fields

        presence = dict(getattr(frame, "parse_presence", {}))
        diagnostics = list(getattr(frame, "parse_diagnostics", []))

        if self.only_extract_structure:
            fields["parse_presence"] = presence
            if diagnostics:
                fields["parse_diagnostics"] = diagnostics
            fields["parse_completeness"] = self._parse_completeness(presence, diagnostics)
            return fields

        task_types = {
            task_type
            for request in getattr(frame, "task_requests", ())
            if getattr(request, "enabled", True)
            and (task_type := getattr(request, "task_type", None))
        }
        optimization_status = getattr(frame, "geometry_optimization_status", None)
        tasks_known = bool(task_types)

        atoms = getattr(frame, "atoms", ())
        coords = getattr(frame, "coords", ())
        try:
            geometry_parsed = bool(atoms) and len(coords) == len(atoms)
        except TypeError:
            geometry_parsed = False
        presence["geometry"] = (
            ParsePresence.PARSED if geometry_parsed else ParsePresence.PARSE_FAILED
        )
        if not geometry_parsed:
            diagnostics.append(
                ParseDiagnostic(
                    code="MOL.PARSE.GEOMETRY_INCOMPLETE",
                    severity="error",
                    scope="frame",
                    message="The parsed atom and coordinate payloads do not form a complete geometry.",
                    segment_index=segment_index,
                    segment_frame_index=segment_frame_index,
                    field="geometry",
                )
            )

        def optional_presence(value: Any, requested: bool | None) -> ParsePresence:
            if value is not None:
                return ParsePresence.PARSED
            if requested is True:
                return ParsePresence.ABSENT_IN_SOURCE
            if requested is False:
                return ParsePresence.NOT_REQUESTED
            return ParsePresence.UNSUPPORTED

        presence["energy"] = optional_presence(
            getattr(frame, "energies", None),
            True if tasks_known else None,
        )
        force_requested = bool(task_types & {"opt", "force", "gradient"}) if tasks_known else None
        frequency_requested = "freq" in task_types if tasks_known else None
        presence["forces"] = optional_presence(
            getattr(frame, "forces", None),
            force_requested,
        )
        presence["hessian"] = optional_presence(
            getattr(frame, "hessian", None),
            frequency_requested,
        )
        presence["vibrations"] = optional_presence(
            getattr(frame, "vibrations", None),
            frequency_requested,
        )
        presence["thermochemistry"] = optional_presence(
            getattr(frame, "thermal_informations", None),
            frequency_requested,
        )

        optimization_requested = "opt" in task_types
        single_point_requested = "sp" in task_types and not optimization_requested
        if optimization_status is not None:
            presence["optimization_status"] = ParsePresence.PARSED
            if single_point_requested:
                diagnostics.append(
                    ParseDiagnostic(
                        code="MOL.CALC.UNEXPECTED_OPTIMIZATION",
                        severity="error",
                        scope="frame",
                        message=(
                            "Optimization evidence was parsed although the segment did not "
                            "request an optimization task."
                        ),
                        segment_index=segment_index,
                        segment_frame_index=segment_frame_index,
                        field="optimization_status",
                    )
                )
        elif optimization_requested:
            presence["optimization_status"] = ParsePresence.ABSENT_IN_SOURCE
        elif tasks_known:
            presence["optimization_status"] = ParsePresence.NOT_REQUESTED
        else:
            presence["optimization_status"] = ParsePresence.UNSUPPORTED

        status = getattr(frame, "status", None)
        scf_converged = getattr(status, "scf_converged", None)
        presence["scf_status"] = optional_presence(
            scf_converged,
            True if tasks_known else None,
        )
        fields["parse_presence"] = presence
        if diagnostics:
            fields["parse_diagnostics"] = diagnostics
        fields["parse_completeness"] = self._parse_completeness(presence, diagnostics)
        return fields

    @staticmethod
    def _parse_completeness(
        presence: Mapping[str, ParsePresence],
        diagnostics: Sequence[ParseDiagnostic],
    ) -> ParseCompleteness:
        if not presence:
            return ParseCompleteness.NOT_ASSESSED
        if any(
            state in {ParsePresence.PARSE_FAILED, ParsePresence.UNSUPPORTED}
            for state in presence.values()
        ):
            return ParseCompleteness.PARTIAL
        if any(
            diagnostic.severity == "error" and diagnostic.code.startswith("MOL.PARSE.")
            for diagnostic in diagnostics
        ):
            return ParseCompleteness.PARTIAL
        return ParseCompleteness.COMPLETE

    @staticmethod
    def _explicit_model_field(value: Any, field: str) -> tuple[bool, Any]:
        if isinstance(value, Mapping):
            return field in value, value.get(field)
        fields_set: set[str] = getattr(value, "model_fields_set", set())
        return field in fields_set, getattr(value, field, None)

    @staticmethod
    def _task_types(metadata: Mapping[str, Any]) -> list[str]:
        task_types: list[str] = []
        for request in metadata.get("task_requests", ()):
            if isinstance(request, Mapping):
                enabled = request.get("enabled", True)
                task_type = request.get("task_type")
            else:
                enabled = getattr(request, "enabled", True)
                task_type = getattr(request, "task_type", None)
            if enabled and isinstance(task_type, str) and task_type not in task_types:
                task_types.append(task_type)
        semantic_route = metadata.get("semantic_route")
        for task_type in getattr(semantic_route, "job_types", ()):
            if isinstance(task_type, str) and task_type not in task_types:
                task_types.append(task_type)
        return task_types

    @classmethod
    def _portable_evidence_value(cls, value: Any) -> Any:
        """Project common metadata models into JSON-portable evidence values."""

        model_dump = getattr(value, "model_dump", None)
        if callable(model_dump):
            value = model_dump(mode="json", exclude_none=True)
        if isinstance(value, Mapping):
            return {
                str(key): cls._portable_evidence_value(item)
                for key, item in value.items()
                if item is not None
            }
        if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
            return [cls._portable_evidence_value(item) for item in value]
        if isinstance(value, Enum):
            return value.value
        return value

    @classmethod
    def _segment_protocol(cls, metadata: Mapping[str, Any]) -> dict[str, Any] | None:
        model_chemistry = metadata.get("model_chemistry")
        if model_chemistry is None:
            return None
        protocol = cls._portable_evidence_value(model_chemistry)
        if not isinstance(protocol, dict):
            raise TypeError("segment model_chemistry must project to a mapping")
        return protocol

    @classmethod
    def _segment_task_requests(cls, metadata: Mapping[str, Any]) -> list[dict[str, Any]]:
        projected = cls._portable_evidence_value(metadata.get("task_requests", ()))
        if not isinstance(projected, list) or not all(
            isinstance(request, dict) for request in projected
        ):
            raise TypeError("segment task_requests must project to a list of mappings")
        return projected

    def _build_segment_evidence(
        self,
        *,
        segment_index: int,
        source_span: SourceSpan,
        source_block_sha256: str,
        frame_count: int,
        captured_frame_indices: list[int],
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Mapping[str, Any],
    ) -> SourceSegmentEvidence:
        status = segment_metadata.get("status")
        termination_set, normal_terminated = self._explicit_model_field(status, "normal_terminated")
        scf_set, scf_converged = self._explicit_model_field(status, "scf_converged")
        if not self._assess_segment_calculation_status:
            termination_presence = None
            scf_presence = None
            normal_terminated = None
            scf_converged = None
        elif self.only_extract_structure:
            termination_presence = ParsePresence.NOT_REQUESTED
            scf_presence = ParsePresence.NOT_REQUESTED
            normal_terminated = None
            scf_converged = None
        else:
            termination_presence = (
                ParsePresence.PARSED if termination_set else ParsePresence.ABSENT_IN_SOURCE
            )
            scf_presence = ParsePresence.PARSED if scf_set else ParsePresence.ABSENT_IN_SOURCE

        diagnostics: list[ParseDiagnostic] = []
        if termination_presence is ParsePresence.ABSENT_IN_SOURCE:
            diagnostics.append(
                ParseDiagnostic(
                    code="MOL.PARSE.TERMINATION_ABSENT",
                    severity="warning",
                    scope="segment",
                    message="No calculation termination marker was found in the source segment.",
                    segment_index=segment_index,
                    field="termination_status",
                )
            )
        elif normal_terminated is False:
            diagnostics.append(
                ParseDiagnostic(
                    code="MOL.CALC.ABNORMAL_TERMINATION",
                    severity="error",
                    scope="segment",
                    message="The source segment contains abnormal termination evidence.",
                    segment_index=segment_index,
                    field="termination_status",
                )
            )
        if scf_converged is False:
            diagnostics.append(
                ParseDiagnostic(
                    code="MOL.CALC.SCF_NOT_CONVERGED",
                    severity="error",
                    scope="segment",
                    message="The last SCF evidence in the source segment reports non-convergence.",
                    segment_index=segment_index,
                    field="scf_status",
                )
            )

        qm_software = segment_metadata.get("qm_software") or artifact_metadata.get("qm_software")
        qm_software_version = segment_metadata.get("qm_software_version") or artifact_metadata.get(
            "qm_software_version"
        )
        parse_presence = {}
        if termination_presence is not None and scf_presence is not None:
            parse_presence = {
                "termination_status": termination_presence,
                "scf_status": scf_presence,
            }
        if self._assess_segment_calculation_status and frame_count == 0:
            parse_presence["geometry"] = ParsePresence.ABSENT_IN_SOURCE
            diagnostics.append(
                ParseDiagnostic(
                    code="MOL.PARSE.SEGMENT_FRAMES_ABSENT",
                    severity="error",
                    scope="segment",
                    message="The calculation segment contains no locator-provided geometry frame.",
                    segment_index=segment_index,
                    field="geometry",
                )
            )
        return SourceSegmentEvidence(
            segment_index=segment_index,
            source_span=source_span,
            source_block_sha256=source_block_sha256,
            frame_count=frame_count,
            captured_frame_indices=captured_frame_indices,
            qm_software=qm_software or None,
            qm_software_version=qm_software_version or None,
            protocol=self._segment_protocol(segment_metadata),
            task_requests=self._segment_task_requests(segment_metadata),
            task_types=self._task_types(segment_metadata),
            termination_status=normal_terminated,
            scf_status=scf_converged,
            parse_presence=parse_presence,
            parse_completeness=self._parse_completeness(parse_presence, diagnostics),
            diagnostics=diagnostics,
        )

    def _frame_metadata(
        self,
        artifact_metadata: Mapping[str, Any],
        segment_metadata: Mapping[str, Any],
        context_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        metadata = dict(artifact_metadata)
        metadata.update(segment_metadata)
        metadata.update(context_metadata)
        for field in self._segment_scoped_frame_metadata:
            metadata.pop(field, None)
        if self.parse_options is not None:
            metadata["topology_reconstruction_backend"] = (
                self.parse_options.graph_reconstruction_backend
            )
            metadata["topology_reconstruction_failure_policy"] = (
                self.parse_options.reconstruction_failure_policy
            )
            metadata["topology_make_dative_bonds"] = self.parse_options.make_dative_bonds
            metadata["topology_make_stereochemistry"] = self.parse_options.make_stereochemistry
        return metadata

    @classmethod
    @abstractmethod
    def _quick_check_file_format(cls, file_content: str) -> None:
        """Run a cheap format fingerprint check after reading content and before full parsing."""
        raise NotImplementedError

    def _parse_frame(
        self, frame_content: str, *, additional_data: dict[str, Any] | None = None
    ) -> FrameT:
        """Parse a single frame."""
        frame_parser = self._frame_parser(
            only_extract_structure=self.only_extract_structure,
            capture_source_evidence=self.capture_source_evidence,
            parse_options=self.parse_options,
        )
        return cast(FrameT, frame_parser.parse(frame_content, additional_data=additional_data))

    def _apply_forced_charge_and_multiplicity(
        self,
        frame: FrameT,
        *,
        charge: int | None,
        multiplicity: int | None,
    ) -> None:
        if charge is not None:
            frame.charge = charge
        if multiplicity is not None:
            frame.multiplicity = multiplicity

    def _attach_artifact_identity(self, source: DecodedSource, chem_file: FileT) -> None:
        for field, value in self._source_artifact_identity(source).items():
            setattr(chem_file, field, value)

    @staticmethod
    def _update_file_parse_completeness(chem_file: FileT) -> None:
        states = [
            *(segment.parse_completeness for segment in chem_file.source_segments),
            *(frame.parse_completeness for frame in chem_file.frames),
        ]
        assessed = [state for state in states if state is not ParseCompleteness.NOT_ASSESSED]
        if not assessed:
            chem_file.parse_completeness = ParseCompleteness.NOT_ASSESSED
        elif any(state is ParseCompleteness.PARTIAL for state in assessed):
            chem_file.parse_completeness = ParseCompleteness.PARTIAL
        else:
            chem_file.parse_completeness = ParseCompleteness.COMPLETE

    def _attach_located_frame_source(
        self,
        source: DecodedSource,
        frame: FrameT,
        block: LocatedTextBlock,
        *,
        segment_index: int,
        segment_frame_index: int,
        file_frame_index: int,
        segment_frame_count: int,
    ) -> None:
        span, block_hash = self._source_block_identity(source, block)
        frame.source_span = span
        frame.source_block_sha256 = block_hash
        frame.segment_index = segment_index
        frame.segment_frame_index = segment_frame_index
        frame.file_frame_index = file_frame_index
        for field, value in self._source_frame_fields(
            frame,
            segment_index=segment_index,
            segment_frame_index=segment_frame_index,
            segment_frame_count=segment_frame_count,
        ).items():
            setattr(frame, field, value)

    def _update_file_metadata_from_frames(
        self,
        chem_file: FileT,
        metadata: dict[str, Any],
    ) -> None:
        """Finalize file-level metadata after frames have been parsed."""
        if not chem_file:
            return
        if not chem_file.frames:
            return
        first_frame = chem_file.frames[0]
        if metadata.get("charge") is None:
            metadata["charge"] = first_frame.charge
            chem_file.charge = first_frame.charge
        if metadata.get("multiplicity") is None:
            metadata["multiplicity"] = first_frame.multiplicity
            chem_file.multiplicity = first_frame.multiplicity

    def _parse(
        self,
        source: str,
        source_type: Literal["file_path", "string"] = "file_path",
        *,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
    ) -> FileT:
        """Run the shared source, metadata, frame, and evidence lifecycle."""
        decoded_source = self._load_source(source, source_type)
        return self._parse_decoded_source(
            decoded_source,
            source_path=source if source_type == "file_path" else None,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )

    def _parse_decoded_source(
        self,
        decoded_source: DecodedSource,
        *,
        source_path: str | None = None,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
    ) -> FileT:
        """Parse one already-loaded source through the shared in-memory lifecycle."""

        self._file_path = source_path
        file_content = decoded_source.text
        parser_file_content = normalize_parser_line_endings(file_content)
        source_format = self.format_id
        if not source_format or source_format != source_format.strip().lower():
            raise ValueError(
                f"{type(self).__name__}.format_id must be a non-empty normalized identifier"
            )
        context_metadata: dict[str, Any] = {
            "file_content": file_content,
            "source_format": source_format,
        }
        if source_path is not None:
            context_metadata["file_path"] = source_path
        self._quick_check_file_format(parser_file_content)

        final_charge = total_charge if total_charge is not None else self.forced_charge
        final_multiplicity = (
            total_multiplicity if total_multiplicity is not None else self.forced_multiplicity
        )
        self.parse_options = replace(
            self.parse_options
            or ParseOptions(
                only_extract_structure=self.only_extract_structure,
                only_last_frame=self.only_last_frame,
                capture_source_evidence=self.capture_source_evidence,
                source_encoding=decoded_source.encoding,
            ),
            total_charge=final_charge,
            total_multiplicity=final_multiplicity,
            only_extract_structure=self.only_extract_structure,
            only_last_frame=self.only_last_frame,
            capture_source_evidence=self.capture_source_evidence,
            source_encoding=decoded_source.encoding,
        ).resolved()
        parser_provenance = self._snapshot_parser_provenance(
            charge_override=final_charge,
            multiplicity_override=final_multiplicity,
            source_encoding=decoded_source.encoding,
            parse_options=self.parse_options,
        )
        if final_charge is not None:
            context_metadata["charge"] = final_charge
        if final_multiplicity is not None:
            context_metadata["multiplicity"] = final_multiplicity

        artifact_metadata = dict(self._parse_artifact_metadata(parser_file_content) or {})
        located_segments = self._validated_located_segments(file_content)
        segment_layouts = [
            (
                segment_index,
                located_segment,
                tuple(enumerate(located_segment.frames)),
            )
            for segment_index, located_segment in enumerate(located_segments)
        ]
        file_frame_indices = {
            (segment_index, segment_frame_index): file_frame_index
            for file_frame_index, (segment_index, segment_frame_index) in enumerate(
                (segment_index, segment_frame_index)
                for segment_index, _, located_frame_blocks in segment_layouts
                for segment_frame_index, _ in located_frame_blocks
            )
        }

        if self.only_last_frame:
            last_layout_with_frames = next(
                (layout for layout in reversed(segment_layouts) if layout[2]),
                None,
            )
            if last_layout_with_frames is None:
                segment_index, segment, located_frame_blocks = segment_layouts[-1]
                selected_source_layouts = [(segment_index, segment, located_frame_blocks)]
            else:
                segment_index, segment, located_frame_blocks = last_layout_with_frames
                selected_source_layouts = [(segment_index, segment, located_frame_blocks[-1:])]
        else:
            selected_source_layouts = segment_layouts

        selected_layouts = [
            (
                segment_index,
                located_segment,
                selected_source_frames,
                dict(
                    self._parse_segment_metadata(
                        normalize_parser_line_endings(located_segment.segment.text(file_content)),
                        artifact_metadata=artifact_metadata,
                    )
                    or {}
                ),
            )
            for segment_index, located_segment, selected_source_frames in selected_source_layouts
        ]

        metadata = self._prepare_file_metadata(
            artifact_metadata,
            [segment_metadata for _, _, _, segment_metadata in selected_layouts],
        )
        metadata.update(context_metadata)
        metadata["parser_provenance"] = parser_provenance
        _chem_file = self._chem_file.model_validate(
            metadata,
            context={"force_unit_transform": self.parse_options.force_unit_transform},
        )
        if self.capture_source_evidence:
            decoded_source.prepare_blocks(
                tuple(
                    block
                    for _, located_segment, selected_source_frames, _ in selected_layouts
                    for block in (
                        located_segment.segment,
                        *(frame_block for _, frame_block in selected_source_frames),
                    )
                )
            )
            self._attach_artifact_identity(decoded_source, _chem_file)
            selected_frame_keys = {
                (segment_index, frame_index)
                for segment_index, _, selected_source_frames, _ in selected_layouts
                for frame_index, _ in selected_source_frames
            }
            all_frame_keys = {
                (segment_index, frame_index)
                for segment_index, _, located_frame_blocks in segment_layouts
                for frame_index, _ in located_frame_blocks
            }
            selected_segment_indices = {
                segment_index for segment_index, _, _, _ in selected_layouts
            }
            all_segment_indices = {segment_index for segment_index, _, _ in segment_layouts}
            _chem_file.source_complete = (
                selected_segment_indices == all_segment_indices
                and selected_frame_keys == all_frame_keys
            )

        for (
            segment_index,
            located_segment,
            selected_source_frames,
            segment_metadata,
        ) in selected_layouts:
            all_frame_count = len(located_segment.frames)
            frame_metadata = self._frame_metadata(
                artifact_metadata,
                segment_metadata,
                context_metadata,
            )
            if self.capture_source_evidence:
                segment_span, segment_hash = self._source_block_identity(
                    decoded_source,
                    located_segment.segment,
                )
                segment_evidence = self._build_segment_evidence(
                    segment_index=segment_index,
                    source_span=segment_span,
                    source_block_sha256=segment_hash,
                    frame_count=all_frame_count,
                    captured_frame_indices=[
                        frame_index for frame_index, _ in selected_source_frames
                    ],
                    artifact_metadata=artifact_metadata,
                    segment_metadata=segment_metadata,
                )
                _chem_file.source_segments.append(segment_evidence)
                _chem_file.source_diagnostics.extend(segment_evidence.diagnostics)
            for segment_frame_index, frame_block in selected_source_frames:
                exact_frame_content = frame_block.text(file_content)
                frame = self._parse_frame(
                    normalize_parser_line_endings(exact_frame_content),
                    additional_data=frame_metadata,
                )
                frame.frame_content = exact_frame_content
                self._apply_forced_charge_and_multiplicity(
                    frame,
                    charge=final_charge,
                    multiplicity=final_multiplicity,
                )
                self._postprocess_parsed_frame(
                    frame,
                    _chem_file,
                    segment_index=segment_index,
                    segment_frame_index=segment_frame_index,
                    segment_frame_count=all_frame_count,
                )
                if self.capture_source_evidence:
                    self._attach_located_frame_source(
                        decoded_source,
                        frame,
                        frame_block,
                        segment_index=segment_index,
                        segment_frame_index=segment_frame_index,
                        file_frame_index=file_frame_indices[(segment_index, segment_frame_index)],
                        segment_frame_count=all_frame_count,
                    )
                _chem_file.append(frame)

        self._update_file_metadata_from_frames(_chem_file, metadata)
        if self.capture_source_evidence:
            self._update_file_parse_completeness(_chem_file)
        return _chem_file

    def to_summary_dict(self, **kwargs) -> SummaryDict:
        return {
            summary_column("FileParser", "forced_charge"): self.forced_charge,
            summary_column("FileParser", "forced_multiplicity"): self.forced_multiplicity,
            summary_column("FileParser", "only_extract_structure"): self.only_extract_structure,
            summary_column("FileParser", "only_last_frame"): self.only_last_frame,
            summary_column("FileParser", "capture_source_evidence"): (self.capture_source_evidence),
            summary_column("FileParser", "source_encoding"): self.source_encoding,
        }


class BaseFileParserMemory(BaseFileParser[FileT, FrameT, FrameParserT]):
    def parse(
        self,
        file_content: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        _chem_file = self._parse(
            source=file_content,
            source_type="string",
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file

    def parse_bytes(
        self,
        raw_bytes: bytes,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        """Parse already-loaded source bytes without a filesystem round trip."""

        return self.parse_decoded_source(
            DecodedSource.from_bytes(raw_bytes, self.source_encoding),
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            release_file_content=release_file_content,
        )

    def parse_decoded_source(
        self,
        source: DecodedSource,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        """Parse a decoded in-memory source while preserving exact byte offsets."""

        _chem_file = self._parse_decoded_source(
            source,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file


class BaseFileParserDisk(BaseFileParser[FileT, FrameT, FrameParserT]):
    allowed_formats: ClassVar[tuple[str, ...]] = ()
    probe_bytes: ClassVar[int] = 20_000

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
        for fmt in cls.allowed_formats:
            if file_path.endswith(fmt):
                return file_path
        raise ValueError(
            f"File {file_path} has an invalid format. Allowed formats: {cls.allowed_formats}."
        )

    def parse(
        self,
        file_path: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        file_path = self._check_file_path(os.path.abspath(file_path))
        _chem_file = self._parse(
            source=file_path,
            source_type="file_path",
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file

    def parse_bytes(
        self,
        raw_bytes: bytes,
        *,
        file_path: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        """Parse loaded bytes while retaining a disk source identity."""

        return self.parse_decoded_source(
            DecodedSource.from_bytes(raw_bytes, self.source_encoding),
            file_path=file_path,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
            release_file_content=release_file_content,
        )

    def parse_decoded_source(
        self,
        source: DecodedSource,
        *,
        file_path: str,
        total_charge: int | None = None,
        total_multiplicity: int | None = None,
        release_file_content: bool = False,
    ) -> FileT:
        """Parse a decoded source without reopening its filesystem path."""

        normalized_path = os.path.abspath(file_path)
        _chem_file = self._parse_decoded_source(
            source,
            source_path=normalized_path,
            total_charge=total_charge,
            total_multiplicity=total_multiplicity,
        )
        if release_file_content:
            _chem_file.release_file_content()
        return _chem_file

    @classmethod
    def probe_file_format(cls, file_path: str | Path) -> bool:
        """Return whether the file prefix matches this parser's own quick check."""

        try:
            with open(file_path, encoding="utf-8", errors="ignore") as f:
                file_content = f.read(cls.probe_bytes)
            cls._quick_check_file_format(file_content)
        except FormatMismatchError:
            return False
        return True
