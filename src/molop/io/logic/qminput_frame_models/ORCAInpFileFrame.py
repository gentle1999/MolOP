"""
Author: TMJ
Date: 2026-02-10 00:00:00
LastEditors: TMJ
LastEditTime: 2026-02-10 21:20:17
Description: ORCA input frame models
"""

from __future__ import annotations

import re
from typing import Any, cast

from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin


class ORCAInpFileFrameMixin:
    orca_raw_preamble: str | None = Field(
        default=None, description="Raw text before geometry block"
    )
    orca_raw_geometry: str | None = Field(default=None, description="Raw geometry block text")
    orca_raw_postamble: str | None = Field(
        default=None, description="Raw text after geometry block"
    )

    orca_fragments: list[list[int]] | None = Field(
        default=None,
        description="Best-effort ORCA fragment markers for round-trip",
    )
    orca_ghost_markers: list[int] | None = Field(
        default=None,
        description="Best-effort ORCA ghost atom markers for round-trip",
    )
    orca_dummy_markers: list[int] | None = Field(
        default=None,
        description="Best-effort ORCA dummy atom markers for round-trip",
    )
    orca_point_charges: list[dict[str, Any]] | None = Field(
        default=None,
        description="Best-effort ORCA point-charge placeholders for round-trip",
    )
    orca_freeze_markers: list[int] | None = Field(
        default=None,
        description="Best-effort ORCA coordinate freeze markers for round-trip",
    )
    orca_isotope_tokens: dict[int, str] | None = Field(
        default=None,
        description="Best-effort ORCA isotope tokens keyed by atom index",
    )
    orca_nuclear_charge_tokens: dict[int, str] | None = Field(
        default=None,
        description="Best-effort ORCA nuclear charge tokens keyed by atom index",
    )

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

        if not has_overrides and self.orca_raw_geometry is not None:
            parts = [self.orca_raw_preamble, self.orca_raw_geometry, self.orca_raw_postamble]
            return "".join(part for part in parts if part)

        typed_self = cast(BaseQMInputFrame, self)

        preamble_parts: list[str] = []

        keywords_to_use = keywords if keywords is not None else typed_self.keywords
        if keywords_to_use:
            keyword_lines = [line.strip() for line in keywords_to_use.splitlines() if line.strip()]
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

        if self.orca_raw_geometry is not None and use_raw_geometry:
            geometry = self.orca_raw_geometry
        else:
            header = f"* xyz {typed_self.charge} {typed_self.multiplicity}"
            body = "\n".join(
                f"{atom} {x:.10f} {y:.10f} {z:.10f}"
                for atom, (x, y, z) in zip(
                    typed_self.atom_symbols, typed_self.coords.m, strict=True
                )
            )
            geometry = f"{header}\n{body}\n*" if body else f"{header}\n*"

        preamble = "\n".join(preamble_parts)
        if preamble:
            return f"{preamble}\n{geometry}"
        return geometry

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
                if part_text := getattr(typed_self, part_name, None):
                    lines = cast(str, part_text).splitlines()
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
): ...


class ORCAInpFileFrameDisk(
    DiskStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameDisk"]
): ...
