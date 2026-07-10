from __future__ import annotations

from typing import Any

from molop.io.logic.orca.common import project_orca_printed_input_metadata
from molop.io.logic.orca.input.frame_parsers.ORCAInpFileFrameParser import (
    parse_orca_input_frame_result,
)


def parse_orca_input_metadata(input_text: str) -> dict[str, Any]:
    if not input_text.strip():
        return {}
    result = parse_orca_input_frame_result(input_text)
    return project_orca_printed_input_metadata(result.model_data())
