"""Animated molecular rendering backed by :mod:`rdkit_dof`."""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any, Literal, Protocol, TypeVar, cast

from rdkit import Chem


AnimationFormat = Literal["gif", "svg"]
FrameValueT = TypeVar("FrameValueT")


class _HasRDKitMol(Protocol):
    @property
    def rdmol(self) -> Chem.Mol | None: ...


def _select_rendered_frame_values(
    values: Sequence[FrameValueT],
    rendered_indices: Sequence[int],
    original_count: int,
    parameter_name: str,
) -> list[FrameValueT]:
    if len(values) != original_count:
        raise ValueError(f"{parameter_name} must have the same length as molecules")
    return [values[index] for index in rendered_indices]


def render_molecule_animation(
    molecules: Sequence[_HasRDKitMol],
    *,
    image_format: AnimationFormat = "gif",
    file_path: os.PathLike[str] | str | None = None,
    duration: int | Sequence[int] = 200,
    loop: int = 0,
    legends: Sequence[str | None] | None = None,
    **kwargs: Any,
) -> Any:
    """Render a sequence of MolOP molecules as an animated GIF or SVG.

    Every renderable molecule becomes one animation frame. Molecules without a
    usable RDKit graph are skipped; rendering fails only when none are usable.
    """

    normalized_format = image_format.lower()
    if normalized_format not in {"gif", "svg"}:
        raise ValueError("image_format must be 'gif' or 'svg'")
    if not molecules:
        raise ValueError("molecules must contain at least one molecule")
    if legends is not None and len(legends) != len(molecules):
        raise ValueError("legends must have the same length as molecules")

    rdmols: list[Chem.Mol] = []
    rendered_indices: list[int] = []
    for index, molecule in enumerate(molecules):
        rdmol = molecule.rdmol
        if rdmol is None:
            continue
        rdmols.append(rdmol)
        rendered_indices.append(index)
    if not rdmols:
        raise ValueError("No animation frames have a usable RDKit molecule")

    rendered_legends = (
        _select_rendered_frame_values(legends, rendered_indices, len(molecules), "legends")
        if legends is not None
        else None
    )
    rendered_duration: int | Sequence[int] = duration
    if not isinstance(duration, int):
        rendered_duration = _select_rendered_frame_values(
            duration,
            rendered_indices,
            len(molecules),
            "duration",
        )

    rendered_kwargs = dict(kwargs)
    for parameter_name in ("highlightAtomLists", "highlightBondLists"):
        values = rendered_kwargs.get(parameter_name)
        if values is not None:
            rendered_kwargs[parameter_name] = _select_rendered_frame_values(
                cast(Sequence[Any], values),
                rendered_indices,
                len(molecules),
                parameter_name,
            )

    filename = os.fspath(file_path) if file_path is not None else None
    if normalized_format == "gif":
        from rdkit_dof import MolsToDofGif

        return MolsToDofGif(
            rdmols,
            legends=rendered_legends,
            duration=rendered_duration,
            loop=loop,
            filename=filename,
            **rendered_kwargs,
        )

    from rdkit_dof import MolsToDofSvgAnimation

    return MolsToDofSvgAnimation(
        rdmols,
        legends=rendered_legends,
        duration=rendered_duration,
        loop=loop,
        filename=filename,
        **rendered_kwargs,
    )
