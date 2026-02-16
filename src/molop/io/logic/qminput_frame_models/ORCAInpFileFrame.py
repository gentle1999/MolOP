"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-16 12:18:32
Description: ORCA input frame models
"""

from __future__ import annotations

import re
from typing import cast

from pydantic import Field, PrivateAttr, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin


class ORCAInpFileFrameMixin:
    orca_raw_preamble: str = Field(default="", description="Raw preamble (preservation-only)")
    orca_raw_postamble: str = Field(default="", description="Raw postamble (preservation-only)")

    def _render(
        self,
        keywords: str | None = None,
        resources_raw: str | None = None,
        nprocs: int | None = None,
        maxcore: int | None = None,
        blocks: str | None = None,
        use_raw_geometry: bool = True,
        **kwargs,
    ) -> str:
        _ = kwargs

        has_overrides = (
            keywords is not None
            or resources_raw is not None
            or nprocs is not None
            or maxcore is not None
            or blocks is not None
            or not use_raw_geometry
        )

        typed_self = cast(BaseQMInputFrame, self)

        preamble_parts: list[str] = []

        if has_overrides:
            keywords_to_use = keywords if keywords is not None else typed_self.keywords
            if keywords_to_use:
                keyword_lines = [
                    line.strip() for line in keywords_to_use.splitlines() if line.strip()
                ]
                preamble_parts.extend(
                    line if line.startswith("!") else f"! {line}" for line in keyword_lines
                )

            resources_raw_to_use = (
                resources_raw if resources_raw is not None else typed_self.resources_raw
            )
            if resources_raw_to_use.strip():
                preamble_parts.append(resources_raw_to_use.strip("\n"))
            if nprocs is not None:
                preamble_parts.append(f"%pal\n  nprocs {nprocs}\nend")
            if maxcore is not None:
                preamble_parts.append(f"%maxcore {maxcore}")
            if blocks is not None and blocks.strip():
                preamble_parts.append(blocks.strip("\n"))
            preamble = "\n".join(preamble_parts)
        else:
            preamble = self.orca_raw_preamble.strip("\n")

        coords_ctype = getattr(self, "_orca_coords_ctype", None)
        external_path = getattr(self, "_orca_external_path", None)

        if isinstance(coords_ctype, str) and coords_ctype.lower() in {"xyzfile", "gzmtfile"}:
            coords_type = coords_ctype.lower()
            charge = int(getattr(typed_self, "charge", 0))
            multiplicity = int(getattr(typed_self, "multiplicity", 1))
            suffix = f" {external_path}" if isinstance(external_path, str) and external_path else ""
            geometry = f"* {coords_type} {charge} {multiplicity}{suffix}".rstrip()
        else:
            header = f"* xyz {typed_self.charge} {typed_self.multiplicity}"
            lines: list[str] = []
            if typed_self.atoms:
                coords_m = typed_self.coords.m
                lines.extend(
                    f"{atom} {x:.10f} {y:.10f} {z:.10f}"
                    for atom, (x, y, z) in zip(typed_self.atom_symbols, coords_m, strict=False)
                )
            point_charges = getattr(self, "_orca_point_charges", None)
            if isinstance(point_charges, list) and point_charges:
                for pc in point_charges:
                    try:
                        q = float(pc["charge"])
                        x = float(pc["x"])
                        y = float(pc["y"])
                        z = float(pc["z"])
                        lines.append(f"Q {q:.2f} {x:.6f} {y:.6f} {z:.6f}")
                    except (KeyError, ValueError, TypeError):
                        continue
            body = "\n".join(lines)
            geometry = f"{header}\n{body}\n*" if body else f"{header}\n*"

        postamble = self.orca_raw_postamble.strip("\n")
        rendered = f"{preamble}\n{geometry}" if preamble else geometry
        if postamble:
            return f"{rendered}\n{postamble}"
        return rendered

    @model_validator(mode="after")
    def _fill_qm_input_metadata(self) -> Self:
        # Populate BaseQMInputFrame fields best-effort from preserved raw input.
        typed_self = cast(BaseQMInputFrame, self)

        if not getattr(typed_self, "qm_software", ""):
            typed_self.qm_software = "ORCA"

        # Preserve resource directives as raw text (no normalization).
        if not getattr(typed_self, "resources_raw", ""):
            resources_blocks: list[str] = []
            for part_name in ["orca_raw_preamble", "orca_raw_postamble"]:
                part_text = cast(str, getattr(typed_self, part_name, "") or "")
                if part_text:
                    lines = part_text.splitlines()
                    i = 0
                    while i < len(lines):
                        line = lines[i]
                        stripped = line.strip()
                        lower = stripped.lower()
                        if lower.startswith("%pal"):
                            block = [line]
                            if not lower.endswith("end"):
                                j = i + 1
                                while j < len(lines):
                                    block.append(lines[j])
                                    if lines[j].strip().lower() == "end":
                                        i = j
                                        break
                                    j += 1
                            resources_blocks.append("\n".join(block))
                        elif lower.startswith("%maxcore"):
                            resources_blocks.append(line)
                        i += 1
            if resources_blocks:
                typed_self.resources_raw = "\n".join(resources_blocks)

        # Very conservative heuristics from !-keywords.
        # Only fill fields when we have a single high-confidence match.
        kw = getattr(typed_self, "keywords", "") or ""
        if kw:
            tokens = [t for t in re.split(r"\s+", kw.strip()) if t]
            upper_tokens = [t.upper() for t in tokens]

            def _pick_one(matches: list[str]) -> str | None:
                uniq = list(dict.fromkeys(matches))
                if len(uniq) == 1:
                    return uniq[0]
                return None

            if not getattr(typed_self, "basis_set", ""):
                basis_hits: list[str] = []
                for t in tokens:
                    if re.match(r"^(?:ma-)?def2-[A-Za-z0-9+\-]+$", t, flags=re.I) or re.match(
                        r"^(?:aug-)?cc-pV[A-Za-z0-9]+Z(?:-PP)?$", t, flags=re.I
                    ):
                        basis_hits.append(t)
                picked = _pick_one(basis_hits)
                if picked:
                    typed_self.basis_set = picked

            if not getattr(typed_self, "functional", "") and not getattr(typed_self, "method", ""):
                # High-precision functional allowlist
                functional_allow = {
                    "B3LYP",
                    "PBE",
                    "PBE0",
                    "BP86",
                    "TPSSh".upper(),
                    "M06",
                    "M06-2X",
                    "WB97X-D",
                    "R2SCAN",
                }
                func_hits = [
                    tokens[i] for i, ut in enumerate(upper_tokens) if ut in functional_allow
                ]
                picked_func = _pick_one(func_hits)
                if picked_func:
                    typed_self.functional = picked_func
                    typed_self.method = "DFT"

            if not getattr(typed_self, "method", ""):
                method_allow = {
                    "HF",
                    "MP2",
                    "CCSD",
                    "CCSD(T)",
                    "DLPNO-CCSD(T)",
                }
                meth_hits = [tokens[i] for i, ut in enumerate(upper_tokens) if ut in method_allow]
                picked_meth = _pick_one(meth_hits)
                if picked_meth:
                    typed_self.method = picked_meth

        return self


class ORCAInpFileFrameMemory(
    MemoryStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameMemory"]
):
    _orca_coords_ctype: str | None = PrivateAttr(default=None)
    _orca_external_path: str | None = PrivateAttr(default=None)
    _orca_point_charges: list[dict[str, float]] | None = PrivateAttr(default=None)


class ORCAInpFileFrameDisk(
    DiskStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameDisk"]
):
    _orca_coords_ctype: str | None = PrivateAttr(default=None)
    _orca_external_path: str | None = PrivateAttr(default=None)
    _orca_point_charges: list[dict[str, float]] | None = PrivateAttr(default=None)
