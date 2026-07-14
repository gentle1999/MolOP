from __future__ import annotations

from molop.io.base_models.source import LocatedTextBlock
from molop.io.codec_exceptions import FormatMismatchError


def _line_offsets(file_content: str) -> tuple[list[str], list[int]]:
    lines = file_content.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    return lines, offsets


def locate_xyz_frames(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate XYZ frames without normalizing their original line endings."""

    lines, offsets = _line_offsets(file_content)
    anchor = 0
    frames: list[LocatedTextBlock] = []
    while anchor < len(lines):
        try:
            num_atoms = int(lines[anchor].strip())
        except ValueError:
            anchor += 1
            continue
        if num_atoms <= 0:
            raise FormatMismatchError("Not an XYZ file: atom count must be positive.")
        frame_end = anchor + num_atoms + 2
        if frame_end > len(lines):
            raise FormatMismatchError("Not an XYZ file: incomplete frame.")
        frames.append(LocatedTextBlock(offsets[anchor], offsets[frame_end]))
        anchor = frame_end
    if not frames:
        raise FormatMismatchError("Not an XYZ file: no atom-count header found.")
    return tuple(frames)


def locate_sdf_frames(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate original SDF records, or one un-delimited MOL block."""

    lines, offsets = _line_offsets(file_content)
    frames: list[LocatedTextBlock] = []
    start_line = 0
    for line_index, line in enumerate(lines):
        if line.strip() != "$$$$":
            continue
        end_line = line_index + 1
        if file_content[offsets[start_line] : offsets[end_line]].strip():
            frames.append(LocatedTextBlock(offsets[start_line], offsets[end_line]))
        start_line = end_line
    if start_line < len(lines) and file_content[offsets[start_line] :].strip():
        frames.append(LocatedTextBlock(offsets[start_line], len(file_content)))
    if not frames and file_content.strip():
        frames.append(LocatedTextBlock(0, len(file_content)))
    if not frames:
        raise FormatMismatchError("Not an SDF/MOL file: no non-empty record found.")
    return tuple(frames)


def locate_smi_frames(file_content: str) -> tuple[LocatedTextBlock, ...]:
    """Locate non-empty SMILES records in the original text."""

    lines, offsets = _line_offsets(file_content)
    frames = tuple(
        LocatedTextBlock(offsets[index], offsets[index + 1])
        for index, line in enumerate(lines)
        if line.strip()
    )
    if not frames:
        raise FormatMismatchError("Not a SMILES file: empty file.")
    return frames


__all__ = [
    "locate_sdf_frames",
    "locate_smi_frames",
    "locate_xyz_frames",
]
