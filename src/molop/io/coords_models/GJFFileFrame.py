"""
Author: TMJ
Date: 2025-07-28 23:09:31
LastEditors: TMJ
LastEditTime: 2026-02-04 11:47:50
Description: 请填写简介
"""

import os
from typing import cast

from pydantic import Field

from molop.io.base_models.ChemFileFrame import BaseCoordsFrame, _HasCoords, _HasKeywords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.patterns.G16Patterns import options_parser, route_section_parser
from molop.structure.FormatConverter import rdmol_to_gjf_connectivity


class GJFFileFrameMixin:
    options: str = Field("")
    title_card: str = Field("")
    suffix: str = Field("")

    def _render(
        self,
        options: str | None = None,
        route: str | None = None,
        title_card: str | None = None,
        suffix: str | None = None,
        template: str | None = None,
        use_link1: bool = False,
        chk: str | None = None,
        old_chk: str | None = None,
        **kwargs,
    ) -> str:
        if template:
            if os.path.exists(template) and os.path.isfile(template):
                from molop.io import GJFFileParserDisk

                template_file = GJFFileParserDisk().parse(template)
                if use_link1 and len(template_file) > 1:
                    return "\n--Link1--\n".join(
                        [
                            self._render_one_frame(
                                options=frame.options,
                                route=cast(_HasKeywords, frame).keywords,
                                title_card=frame.title_card,
                                suffix=frame.suffix,
                                chk=f"chk_{frame.frame_id + 1}_{chk}" if chk else None,
                                old_chk=f"chk_{frame.frame_id}_{old_chk}"
                                if old_chk and frame.frame_id > 0
                                else None,
                                link1_mode=True,
                                **kwargs,
                            )
                            for frame in template_file
                        ]
                    )
                options = options or template_file[0].options
                route = route or cast(_HasKeywords, template_file[0]).keywords
                title_card = title_card or template_file[0].title_card
                suffix = suffix or template_file[0].suffix
                return self._render_one_frame(
                    options=options,
                    route=route,
                    title_card=title_card,
                    suffix=suffix,
                    chk=chk,
                    old_chk=old_chk,
                    **kwargs,
                )

            else:
                raise FileNotFoundError(f"template file {template} not found or not a file")
        else:
            return self._render_one_frame(
                options=options,
                route=route,
                title_card=title_card,
                suffix=suffix,
                chk=chk,
                old_chk=old_chk,
                **kwargs,
            )

    def _render_one_frame(
        self,
        options: str | None = None,
        route: str | None = None,
        title_card: str | None = None,
        suffix: str | None = None,
        chk: str | None = None,
        old_chk: str | None = None,
        link1_mode: bool = False,
    ) -> str:
        typed_self = cast(_HasCoords, self)
        options_to_use = options or self.options
        route_to_use = route or cast(_HasKeywords, self).keywords
        title_card_to_use = title_card or self.title_card
        suffix_to_use = suffix or self.suffix
        _options = options_parser(options_to_use)
        if chk:
            _options[r"%chk"] = chk
        if old_chk:
            _options[r"%oldchk"] = old_chk
        options_lines = (
            "\n".join([f"{key}={val}" for key, val in _options.items()]) + "\n" if _options else ""
        )
        print_title_card = True
        print_electronic_state = True
        print_coords = True
        if link1_mode:
            route_params, dieze_tag = route_section_parser(route_to_use)
            if isinstance(route_params.get("geom"), str):
                if route_params.get("geom") in ("check", "checkpoint"):
                    print_coords = False
                if route_params.get("geom") in ("allcheck", "allcheckpoint"):
                    print_electronic_state = False
                    print_title_card = False
                    print_coords = False
            elif isinstance(route_params.get("geom"), dict):
                if any(key in ("check", "checkpoint") for key in route_params.get("geom", {})):
                    print_coords = False
                if any(
                    key in ("allcheck", "allcheckpoint") for key in route_params.get("geom", {})
                ):
                    print_electronic_state = False
                    print_title_card = False
                    print_coords = False
        route_to_use = route_to_use or "#"
        return (
            options_lines
            + f"{route_to_use}\n\n"
            + (f"{title_card_to_use or 'title'}\n\n" if print_title_card else "")
            + (f"{typed_self.charge} {typed_self.multiplicity}\n" if print_electronic_state else "")
            + (
                "\n".join(
                    [
                        f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                        for atom, (x, y, z) in zip(
                            typed_self.atom_symbols, typed_self.coords.m, strict=True
                        )
                    ]
                )
                + "\n\n"
                + self._generate_connectivity()
                if print_coords
                else ""
            )
            + suffix_to_use
            + "\n\n"
        )

    def _generate_connectivity(self) -> str:
        typed_self = cast(_HasCoords, self)
        if typed_self.rdmol:
            return rdmol_to_gjf_connectivity(typed_self.rdmol) + "\n\n"
        else:
            return ""


class GJFFileFrameMemory(
    MemoryStorageMixin, GJFFileFrameMixin, BaseCoordsFrame["GJFFileFrameMemory"]
): ...


class GJFFileFrameDisk(DiskStorageMixin, GJFFileFrameMixin, BaseCoordsFrame["GJFFileFrameDisk"]):
    def _render(
        self,
        options: str | None = None,
        route: str | None = None,
        title_card: str | None = None,
        suffix: str | None = None,
        template: str | None = None,
        use_link1: bool = False,
        chk: bool | str | None = False,
        old_chk: bool | str | None = False,
        **kwargs,
    ) -> str:
        valid_chk, valid_old_chk = None, None
        if chk:
            if isinstance(chk, str):
                valid_chk = chk
            elif isinstance(chk, bool):
                valid_chk = f"{self.pure_filename}.chk"
            else:
                raise ValueError(f"chk must be a bool or str, but got {chk}")
        if old_chk:
            if isinstance(old_chk, str):
                valid_old_chk = old_chk
            elif isinstance(old_chk, bool):
                valid_old_chk = f"{self.pure_filename}.chk"
            else:
                raise ValueError(f"old_chk must be a bool or str, but got {old_chk}")
        return super()._render(
            options=options,
            route=route,
            title_card=title_card,
            suffix=suffix,
            template=template,
            use_link1=use_link1,
            chk=valid_chk,
            old_chk=valid_old_chk,
            **kwargs,
        )
