from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, TypeAlias


FrameSelector: TypeAlias = int | Sequence[int] | Literal["all"]


def normalize_frame_selector(
    frame_selector: FrameSelector,
    frame_count: int,
    *,
    parameter_name: str = "frame",
    validate_range: bool = False,
) -> list[int]:
    if isinstance(frame_selector, str):
        if frame_selector == "all":
            return list(range(frame_count))
        raise ValueError(f'{parameter_name} must be an integer, a sequence of integers, or "all"')

    raw_frame_ids = [frame_selector] if isinstance(frame_selector, int) else list(frame_selector)
    frame_ids: list[int] = []
    for raw_frame_id in raw_frame_ids:
        if not isinstance(raw_frame_id, int):
            raise TypeError(f"{parameter_name} sequence must contain only integers")
        frame_id = raw_frame_id if raw_frame_id >= 0 else frame_count + raw_frame_id
        if validate_range and not 0 <= frame_id < frame_count:
            raise IndexError(f"Frame index {raw_frame_id} is out of range")
        frame_ids.append(frame_id)
    return frame_ids
