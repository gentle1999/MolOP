import pytest
from pydantic import ValidationError

from molop.io.base_models.source import SourceSpan


def _span(start: int, end: int) -> SourceSpan:
    return SourceSpan(
        start_byte=start,
        end_byte=end,
        start_char=start,
        end_char=end,
        start_line=start + 1,
        end_line=end + 1,
    )


def test_evidence_value_models_are_frozen_and_forbid_extra_fields() -> None:
    span = _span(0, 4)

    with pytest.raises(ValidationError, match="frozen"):
        span.start_byte = 1

    with pytest.raises(ValidationError, match="extra"):
        SourceSpan.model_validate({**span.model_dump(), "unexpected": True})


def test_source_span_rejects_empty_half_open_ranges() -> None:
    with pytest.raises(ValidationError, match="end_byte must be greater than start_byte"):
        SourceSpan(
            start_byte=2,
            end_byte=2,
            start_char=2,
            end_char=2,
            start_line=1,
            end_line=1,
        )
