"""Exact source locators for Gaussian Link1 sections and geometry frames."""

from molop.io.base_models.source import LocatedTextBlock
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns


INPUT_ORIENTATION_MARKER = "Input orientation:"
STANDARD_ORIENTATION_MARKER = "Standard orientation:"


def locate_g16_sections(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate Gaussian Link1 sections without rebuilding their source text."""

    matches = g16_log_patterns.LINK1_SECTION.find_matches(file_content)
    if not matches:
        return (LocatedTextBlock(0, len(file_content)),)

    starts = [0]
    for match_index, matched in enumerate(matches):
        label = matched.group().strip()
        if match_index == 0 and label.startswith("Entering Link 1"):
            continue
        if matched.start() > starts[-1]:
            starts.append(matched.start())

    return tuple(
        LocatedTextBlock(start, starts[index + 1] if index + 1 < len(starts) else len(file_content))
        for index, start in enumerate(starts)
    )


def _marker_positions(content: str, marker: str) -> list[int]:
    positions: list[int] = []
    start = 0
    while (position := content.find(marker, start)) >= 0:
        positions.append(content.rfind("\n", 0, position) + 1)
        start = position + len(marker)
    return positions


def locate_g16_section_frames(
    file_content: str,
    section: LocatedTextBlock,
) -> tuple[LocatedTextBlock, ...]:
    """Locate geometry frames within one Link1 section as exact source slices."""

    section_content = section.text(file_content)
    if INPUT_ORIENTATION_MARKER in section_content:
        marker = INPUT_ORIENTATION_MARKER
    elif STANDARD_ORIENTATION_MARKER in section_content:
        marker = STANDARD_ORIENTATION_MARKER
    else:
        return ()

    marker_positions = _marker_positions(section_content, marker)
    if not marker_positions:
        return ()

    starts = [section.start_char]
    starts.extend(section.start_char + position for position in marker_positions[1:])
    return tuple(
        LocatedTextBlock(start, starts[index + 1] if index + 1 < len(starts) else section.end_char)
        for index, start in enumerate(starts)
    )


__all__ = [
    "INPUT_ORIENTATION_MARKER",
    "LocatedTextBlock",
    "STANDARD_ORIENTATION_MARKER",
    "locate_g16_section_frames",
    "locate_g16_sections",
]
