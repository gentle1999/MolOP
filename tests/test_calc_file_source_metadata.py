from molop.io.base_models.ChemFile import BaseCalcFile, BaseChemFile
from molop.io.base_models.source import (
    ParseCompleteness,
    ParseDiagnostic,
    ParsePresence,
    SourceSegmentEvidence,
    SourceSpan,
)


def test_file_source_metadata_is_generic_sparse_and_dumpable() -> None:
    source_fields = {
        "artifact_sha256",
        "artifact_size_bytes",
        "source_encoding",
        "source_segments",
        "source_complete",
        "source_diagnostics",
        "parse_completeness",
    }

    assert source_fields <= set(BaseChemFile.model_fields)
    assert source_fields.isdisjoint(BaseChemFile().model_dump())
    assert source_fields.isdisjoint(BaseCalcFile().model_dump())

    artifact_sha256 = "a" * 64
    source_span = SourceSpan(
        start_byte=0,
        end_byte=128,
        start_char=0,
        end_char=128,
        start_line=1,
        end_line=2,
    )
    source_segments = [
        SourceSegmentEvidence(
            segment_index=0,
            source_span=source_span,
            source_block_sha256="b" * 64,
            frame_count=1,
            captured_frame_indices=[0],
            parse_presence={"energy": ParsePresence.PARSE_FAILED},
            parse_completeness=ParseCompleteness.PARTIAL,
        )
    ]
    source_diagnostics = [
        ParseDiagnostic(
            code="MOL.PARSE.PARTIAL",
            severity="warning",
            scope="artifact",
            message="Only part of the requested scientific payload was parsed.",
        )
    ]
    dumped = BaseChemFile(
        artifact_sha256=artifact_sha256,
        artifact_size_bytes=128,
        source_encoding="utf-8",
        source_segments=source_segments,
        source_complete=False,
        source_diagnostics=source_diagnostics,
        parse_completeness=ParseCompleteness.PARTIAL,
    ).model_dump()

    assert {field: dumped[field] for field in source_fields} == {
        "artifact_sha256": artifact_sha256,
        "artifact_size_bytes": 128,
        "source_encoding": "utf-8",
        "source_segments": [segment.model_dump() for segment in source_segments],
        "source_complete": False,
        "source_diagnostics": [diagnostic.model_dump() for diagnostic in source_diagnostics],
        "parse_completeness": ParseCompleteness.PARTIAL,
    }
