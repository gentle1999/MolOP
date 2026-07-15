from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, cast

import numpy as np
from rdkit import Chem

from molop.io.logic.orca.common import (
    ORCABlock,
    ORCABlockLine,
    ORCAGeometry,
    ORCAGeometryAtom,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_blocks import (
    extract_block_spans,
)
from molop.unit import atom_ureg


def build_orca_geometry_from_frame_payload(
    payload: Mapping[str, Any],
) -> ORCAGeometry | None:
    geometry = payload.get("geometry")
    if isinstance(geometry, ORCAGeometry):
        return geometry
    if isinstance(geometry, Mapping):
        return ORCAGeometry.model_validate(geometry)

    atoms = payload.get("atoms")
    coords = payload.get("coords")
    if not isinstance(atoms, list) or not atoms or coords is None:
        return None

    if hasattr(coords, "to"):
        coord_values = np.asarray(cast(Any, coords).to(atom_ureg.angstrom).magnitude)
    elif hasattr(coords, "magnitude"):
        coord_values = np.asarray(cast(Any, coords).magnitude)
    else:
        coord_values = np.asarray(coords)
    if coord_values.shape != (len(atoms), 3):
        return None

    periodic_table = Chem.GetPeriodicTable()
    geometry_atoms: list[ORCAGeometryAtom] = []
    for atom, (x, y, z) in zip(atoms, coord_values, strict=True):
        if isinstance(atom, str):
            symbol = atom
            atomic_number = periodic_table.GetAtomicNumber(symbol)
        else:
            atomic_number = int(atom)
            symbol = periodic_table.GetElementSymbol(atomic_number)
        geometry_atoms.append(
            ORCAGeometryAtom(
                symbol=symbol,
                atomic_number=atomic_number,
                x=float(x),
                y=float(y),
                z=float(z),
            )
        )

    return ORCAGeometry(
        items=geometry_atoms,
        charge=int(payload.get("charge", 0) or 0),
        multiplicity=int(payload.get("multiplicity", 1) or 1),
        source="unknown",
    )


def adapt_orca_writer_payload(data: Mapping[str, Any]) -> dict[str, Any]:
    payload = dict(data)
    geometry = build_orca_geometry_from_frame_payload(payload)
    if geometry is not None:
        payload["geometry"] = geometry

    # Keywords from another program are not valid ORCA input keywords. Cross-format
    # callers must provide them explicitly to the renderer.
    qm_software = str(payload.get("qm_software", "")).strip().lower()
    if qm_software != "orca":
        payload["keywords"] = ""
        payload["comment_lines"] = []
        payload["keyword_lines"] = []
        payload["blocks"] = []
        payload["trailing_lines"] = []
    return payload


def render_orca_input_frame(
    *,
    fallback_keywords: str,
    fallback_keyword_lines: Sequence[Any],
    fallback_blocks: Sequence[ORCABlock],
    fallback_comments: Sequence[Any],
    fallback_trailing_lines: Sequence[str],
    geometry: ORCAGeometry | None,
    keywords: str | Sequence[str] | None = None,
    nprocs: int | None = None,
    maxcore: int | None = None,
    blocks: str | Mapping[str, Any] | Sequence[ORCABlock] | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    coordinate_decimal_places: int = 10,
) -> str:
    keyword_lines = _resolve_keyword_lines(keywords, fallback_keywords, fallback_keyword_lines)
    if not keyword_lines:
        raise ValueError("ORCA input rendering requires at least one keyword line.")
    if geometry is None:
        raise ValueError("ORCA input rendering requires molecular geometry.")
    if nprocs is not None and nprocs <= 0:
        raise ValueError("nprocs must be greater than zero.")
    if maxcore is not None and maxcore <= 0:
        raise ValueError("maxcore must be greater than zero.")
    if not 0 <= coordinate_decimal_places <= 18:
        raise ValueError("coordinate_decimal_places must be between 0 and 18.")

    resolved_blocks = [
        block
        for block in _merge_blocks(list(fallback_blocks), blocks)
        if block.name.lower() != "coords"
    ]
    if nprocs is not None:
        resolved_blocks = [block for block in resolved_blocks if block.name.lower() != "pal"]
    if maxcore is not None:
        resolved_blocks = [block for block in resolved_blocks if block.name.lower() != "maxcore"]

    parts: list[str] = []
    parts.extend(f"# {line.text}" for line in fallback_comments if line.text.strip())
    parts.extend(f"! {line}" for line in keyword_lines)
    if nprocs is not None:
        parts.append(f"%pal\n  nprocs {nprocs}\nend")
    if maxcore is not None:
        parts.append(f"%maxcore {maxcore}")
    parts.extend(_render_block(block) for block in resolved_blocks)
    parts.append(
        _render_geometry(
            geometry,
            charge=charge,
            multiplicity=multiplicity,
            decimal_places=coordinate_decimal_places,
        )
    )
    parts.extend(line.rstrip("\r\n") for line in fallback_trailing_lines if line.strip())
    return "\n\n".join(part for part in parts if part.strip())


def _resolve_keyword_lines(
    keywords: str | Sequence[str] | None,
    fallback_keywords: str,
    fallback_keyword_lines: Sequence[Any],
) -> list[str]:
    if keywords is None:
        lines = [str(line.text) for line in fallback_keyword_lines if line.text.strip()]
        if not lines:
            lines = fallback_keywords.splitlines()
    elif isinstance(keywords, str):
        lines = keywords.splitlines()
    else:
        lines = list(keywords)
    return [line.strip().removeprefix("!").strip() for line in lines if line.strip()]


def _merge_blocks(
    fallback_blocks: list[ORCABlock],
    overrides: str | Mapping[str, Any] | Sequence[ORCABlock] | None,
) -> list[ORCABlock]:
    if overrides is None:
        return fallback_blocks
    override_blocks = _coerce_blocks(overrides)
    override_names = {block.name.lower() for block in override_blocks}
    return [
        block for block in fallback_blocks if block.name.lower() not in override_names
    ] + override_blocks


def _coerce_blocks(
    value: str | Mapping[str, Any] | Sequence[ORCABlock],
) -> list[ORCABlock]:
    if isinstance(value, str):
        spans = extract_block_spans(value)
        if not spans and value.strip():
            raise ValueError("blocks must contain one or more valid ORCA % blocks.")
        return [span.block for span in spans]
    if isinstance(value, Mapping):
        return [_block_from_value(str(name), body) for name, body in value.items()]
    return [
        block if isinstance(block, ORCABlock) else ORCABlock.model_validate(block)
        for block in value
    ]


def _block_from_value(name: str, body: Any) -> ORCABlock:
    if isinstance(body, Mapping):
        lines = [
            ORCABlockLine(text=f"  {key} {_render_scalar(value)}") for key, value in body.items()
        ]
    elif isinstance(body, str):
        lines = [ORCABlockLine(text=line) for line in body.splitlines() if line.strip()]
    elif isinstance(body, Sequence):
        lines = [ORCABlockLine(text=str(line)) for line in body]
    else:
        lines = [ORCABlockLine(text=str(body))]
    return ORCABlock(name=name.removeprefix("%"), lines=lines)


def _render_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def _render_block(block: ORCABlock) -> str:
    if block.raw_text.strip():
        return block.raw_text.strip("\r\n")
    body = "\n".join(line.text.rstrip("\r\n") for line in block.lines)
    if not body:
        return f"%{block.name}"
    return f"%{block.name}\n{body}\nend"


def _render_geometry(
    geometry: ORCAGeometry,
    *,
    charge: int | None,
    multiplicity: int | None,
    decimal_places: int,
) -> str:
    final_charge = geometry.charge if charge is None else charge
    final_multiplicity = geometry.multiplicity if multiplicity is None else multiplicity
    if final_multiplicity <= 0:
        raise ValueError("multiplicity must be greater than zero.")

    if geometry.ctype in {"xyzfile", "gzmtfile", "pdbfile"}:
        if not geometry.external_path:
            raise ValueError(f"{geometry.ctype} geometry requires external_path.")
        return f"* {geometry.ctype} {final_charge} {final_multiplicity} {geometry.external_path}"
    if geometry.ctype not in {"xyz", "cart", "cartesian"}:
        raise NotImplementedError(f"ORCA writer does not support {geometry.ctype!r} geometry yet.")

    coordinate_lines: list[str] = []
    width = decimal_places + 6
    for atom in geometry:
        if atom.x is None or atom.y is None or atom.z is None:
            raise ValueError(f"ORCA Cartesian atom {atom.symbol!r} has incomplete coordinates.")
        symbol = atom.symbol + (":" if atom.is_ghost else "")
        if atom.fragment_id is not None:
            symbol = f"{symbol}({atom.fragment_id})"
        prefix = "$" if atom.frozen else ""
        line = (
            f"{symbol:<2} {prefix}{atom.x:{width}.{decimal_places}f} "
            f"{atom.y:{width}.{decimal_places}f} {atom.z:{width}.{decimal_places}f}"
        )
        if atom.basis_overrides:
            override_tokens = [
                " ".join([override.kind, *override.tokens]) for override in atom.basis_overrides
            ]
            line = f"{line} {' '.join(override_tokens)}"
        coordinate_lines.append(line)
    for point_charge in geometry.point_charges:
        coordinate_lines.append(
            "Q  "
            f"{point_charge['charge']:.6f} "
            f"{point_charge['x']:{width}.{decimal_places}f} "
            f"{point_charge['y']:{width}.{decimal_places}f} "
            f"{point_charge['z']:{width}.{decimal_places}f}"
        )

    header = f"* xyz {final_charge} {final_multiplicity}"
    body = "\n".join(coordinate_lines)
    return f"{header}\n{body}\n*" if body else f"{header}\n*"


__all__ = [
    "adapt_orca_writer_payload",
    "build_orca_geometry_from_frame_payload",
    "render_orca_input_frame",
]
