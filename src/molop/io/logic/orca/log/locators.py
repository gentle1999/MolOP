"""Exact source locators for ORCA jobs and coordinate-bearing frames."""

from molop.io.base_models.source import LocatedTextBlock
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns


def _line_start(text: str, position: int) -> int:
    return text.rfind("\n", 0, position) + 1


def _job_starts(file_content: str) -> list[int]:
    numbered_starts = [
        _line_start(file_content, matched.start())
        for matched in orca_log_patterns.JOB_NUMBER_LINE.find_matches(file_content)
        if int(matched.group("number")) > 1
    ]
    if numbered_starts:
        return [0, *numbered_starts]

    banner_matches = orca_log_patterns.BANNER.find_matches(file_content)
    return [
        0,
        *(
            _line_start(file_content, matched.start())
            for matched in banner_matches[1:]
            if matched.start() > 0
        ),
    ]


def locate_orca_jobs(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate independent ORCA jobs without rebuilding or copying their text."""

    if not file_content:
        return ()
    starts = list(dict.fromkeys(_job_starts(file_content)))
    return tuple(
        LocatedTextBlock(start, starts[index + 1] if index + 1 < len(starts) else len(file_content))
        for index, start in enumerate(starts)
    )


def locate_orca_job_frames(
    file_content: str,
    job: LocatedTextBlock,
) -> tuple[LocatedTextBlock, ...]:
    """Locate the exact parser blocks associated with ORCA coordinate prints."""

    job_content = job.text(file_content)
    matches = orca_log_patterns.COORD_HEADER.find_matches(job_content)
    if not matches:
        return ()
    starts = [job.start_char]
    starts.extend(job.start_char + matched.start() for matched in matches[1:])
    candidates = tuple(
        LocatedTextBlock(start, starts[index + 1] if index + 1 < len(starts) else job.end_char)
        for index, start in enumerate(starts)
    )
    return tuple(
        frame for frame in candidates if _contains_parseable_coordinates(frame.text(file_content))
    )


def _contains_parseable_coordinates(frame_content: str) -> bool:
    headers = orca_log_patterns.COORD_HEADER.find_matches(frame_content)
    if not headers:
        return False
    return any(
        orca_log_patterns.COORD_ROW.match(line) is not None
        for line in frame_content[headers[-1].end() :].splitlines()
    )


__all__ = ["locate_orca_job_frames", "locate_orca_jobs"]
