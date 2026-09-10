"""File export helpers for inferred transition-state endpoints."""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from rdkit import Chem

from molop.utils.types import RdMol


EndpointFormat = Literal["xyz", "sdf"]
XYZSerializer = Callable[[RdMol], str]


def export_ts_endpoints(
    output_dir: os.PathLike[str] | str,
    pre_rdmol: RdMol,
    post_rdmol: RdMol,
    *,
    prefix: str,
    format: EndpointFormat = "xyz",
    xyz_serializer: XYZSerializer,
) -> tuple[Path, Path]:
    """Write inferred pre- and post-TS endpoint molecules to disk."""

    normalized_format = format.lower()
    if normalized_format not in {"xyz", "sdf"}:
        raise ValueError(f"Unsupported endpoint format: {format!r}. Use 'xyz' or 'sdf'.")

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    pre_path = destination / f"{prefix}_pre.{normalized_format}"
    post_path = destination / f"{prefix}_post.{normalized_format}"

    if normalized_format == "xyz":
        pre_path.write_text(xyz_serializer(pre_rdmol), encoding="utf-8")
        post_path.write_text(xyz_serializer(post_rdmol), encoding="utf-8")
    else:
        for endpoint_name, rdmol, path in (
            ("pre-TS endpoint candidate", pre_rdmol, pre_path),
            ("post-TS endpoint candidate", post_rdmol, post_path),
        ):
            endpoint_rdmol = Chem.Mol(rdmol)
            endpoint_rdmol.SetProp("_Name", endpoint_name)
            writer = Chem.SDWriter(str(path))
            try:
                writer.write(endpoint_rdmol)
            finally:
                writer.close()
    return pre_path, post_path


__all__ = ["EndpointFormat", "XYZSerializer", "export_ts_endpoints"]
