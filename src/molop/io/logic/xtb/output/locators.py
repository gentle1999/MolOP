"""Exact source locators for concatenated xTB command-line runs."""

from molop.io.base_models.source import LocatedSourceSegment, LocatedTextBlock
from molop.io.logic.xtb.output.parsers._xtb_output_extractors import iter_xtb_version_matches


def locate_xtb_runs(file_content: str) -> tuple[LocatedSourceSegment, ...]:
    if not file_content:
        return ()

    version_matches = tuple(iter_xtb_version_matches(file_content))
    if not version_matches:
        return ()

    starts = [0]
    starts.extend(
        file_content.rfind("\n", 0, matched.start()) + 1 for matched in version_matches[1:]
    )
    starts = list(dict.fromkeys(starts))
    return tuple(
        LocatedSourceSegment(
            segment=(
                block := LocatedTextBlock(
                    start,
                    starts[index + 1] if index + 1 < len(starts) else len(file_content),
                )
            ),
            frames=(block,),
        )
        for index, start in enumerate(starts)
    )


__all__ = ["locate_xtb_runs"]
