from typing import Literal


def parse_frame_selection(spec: str) -> int | list[int] | Literal["all"]:
    spec = spec.strip()
    if not spec:
        raise ValueError("Frame selection string cannot be empty.")

    if spec.lower() == "all":
        return "all"

    if "," in spec:
        try:
            return [int(s.strip()) for s in spec.split(",")]
        except ValueError as e:
            raise ValueError(
                f"Invalid frame selection: '{spec}'. Expected 'all', an integer, or a comma-separated list of integers."
            ) from e

    try:
        return int(spec)
    except ValueError as e:
        raise ValueError(
            f"Invalid frame selection: '{spec}'. Expected 'all', an integer, or a comma-separated list of integers."
        ) from e
