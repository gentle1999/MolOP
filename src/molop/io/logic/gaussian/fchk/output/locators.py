from molop.io.base_models.source import LocatedSourceSegment, LocatedTextBlock


def locate_fchk_content(file_content: str) -> tuple[LocatedSourceSegment, ...]:
    if not file_content:
        return ()
    block = LocatedTextBlock(0, len(file_content))
    return (LocatedSourceSegment(segment=block, frames=(block,)),)


__all__ = ["locate_fchk_content"]
