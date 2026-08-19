"""
Author: TMJ
Date: 2025-07-28 23:09:31
LastEditors: TMJ
LastEditTime: 2026-04-05 17:29:39
Description: 请填写简介
"""

from collections.abc import Mapping
from typing import Any, Literal, cast

from pydantic import Field, model_validator

from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFAtomSpecification as GJFAtomSpecification,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFGICLine as GJFGICLine,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFGICOption as GJFGICOption,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFGICSection as GJFGICSection,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFLink0 as GJFLink0,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFLink0Commands as GJFLink0Commands,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFModRedundantLine as GJFModRedundantLine,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFModRedundantSection as GJFModRedundantSection,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFMoleculeSpecifications as GJFMoleculeSpecifications,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFMoleculeSpecificationsFragment as GJFMoleculeSpecificationsFragment,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFNBOSection as GJFNBOSection,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFRouteSection as GJFRouteSection,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFSectionParsingDiagnostic as GJFSectionParsingDiagnostic,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFTitleCard as GJFTitleCard,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    GJFUnknownSection as GJFUnknownSection,
)
from molop.io.logic.gaussian.input.GaussianInput import (
    adapt_gjf_writer_payload,
    normalize_gjf_frame_geometry_payload,
    resolve_gjf_render_parts,
)
from molop.io.logic.gaussian.input.GaussianRoute import (
    populate_gaussian_legacy_qm_fields_from_semantic,
)
from molop.structure.FormatConverter import rdmol_to_gjf_connectivity
from molop.utils.progressbar import NativeReconstructionConcurrencyError


class GJFFileFrameMixin:
    @classmethod
    def adapt_writer_payload(cls, data: dict[str, Any]) -> dict[str, Any]:
        return adapt_gjf_writer_payload(data)

    @model_validator(mode="before")
    @classmethod
    def normalize_frame_payload(cls, data: Any):
        if not isinstance(data, Mapping):
            return data

        return normalize_gjf_frame_geometry_payload(data)

    link0_commands: GJFLink0Commands = Field(
        default=GJFLink0Commands(), description="Link0 commands, define the calculation options"
    )
    route_section: GJFRouteSection = Field(
        default=GJFRouteSection(), description="Route section, define the calculation keywords"
    )
    title_card: GJFTitleCard = Field(default=GJFTitleCard(), description="Title card")
    molecule_specifications: GJFMoleculeSpecifications = Field(
        default=GJFMoleculeSpecifications(), description="Molecule specifications"
    )
    additional_sections: str = Field(default="", description="Additional sections")
    parsed_additional_sections: list[
        GJFGICSection | GJFModRedundantSection | GJFNBOSection | GJFUnknownSection
    ] = Field(default_factory=list, description="Parsed additional sections")
    additional_section_diagnostics: list[GJFSectionParsingDiagnostic] = Field(
        default_factory=list, description="Diagnostics for additional section parsing"
    )

    def _render(
        self,
        link0_commands: str | GJFLink0Commands | dict[str, str | None] | None = None,
        route_section: str | GJFRouteSection | None = None,
        title_card: str | GJFTitleCard | None = None,
        molecule_specifications: str | GJFMoleculeSpecifications | None = None,
        additional_sections: str | None = None,
        parsed_additional_sections: list[
            GJFGICSection | GJFModRedundantSection | GJFNBOSection | GJFUnknownSection
        ]
        | None = None,
        chk: str | bool | None = None,
        old_chk: bool | str | None = None,
        coords_type: Literal["cartesian", "internal", "auto"] = "auto",
        add_gjf_connectivity: bool = False,
        **kwargs,
    ) -> str:
        """
        Render the current GJF frame as Gaussian input text.

        The method assembles the output in standard Gaussian input order:
        1. Link0 commands
        2. Route section
        3. Title card
        4. Molecule specifications
        5. Additional sections

        Callers may temporarily override any of these parts for a single render
        without mutating the instance first. Any argument left as ``None`` falls
        back to the corresponding attribute on ``self``.

        Parameters
        ----------
        link0_commands : str | GJFLink0Commands | dict[str, str | None] | None, optional
            Link0 content, i.e. the resource and job-control directives at the
            beginning of the file such as ``%mem``, ``%nprocshared``, and
            ``%chk``.

            Actual behavior:
            - ``str``: parsed as Link0 syntax into ``GJFLink0Commands``;
            - ``dict``: converted into a Link0 keyword list;
            - ``GJFLink0Commands``: used directly;
            - ``None``: falls back to ``self.link0_commands``.

        route_section : str | GJFRouteSection | None, optional
            Gaussian route section, i.e. the keyword line beginning with ``#``,
            for example ``#p b3lyp/6-31g(d) opt freq``.

            Actual behavior:
            - ``str``: converted to ``GJFRouteSection`` and normalized there;
            - ``GJFRouteSection``: used directly;
            - ``None``: falls back to ``self.route_section``.

        title_card : str | GJFTitleCard | None, optional
            Gaussian title section.

            Actual behavior:
            - ``str``: converted to ``GJFTitleCard``;
            - ``GJFTitleCard``: used directly;
            - ``None``: falls back to ``self.title_card``;
            - if the final title is empty, it falls back to ``self.pure_filename``
              or, if unavailable, ``"title"``.

        molecule_specifications : str | GJFMoleculeSpecifications | None, optional
            Molecule specification block, i.e. the charge/multiplicity line plus
            the coordinate block that follows it.

            Actual behavior:
            - ``str``: parsed as Gaussian molecule-block syntax into
              ``GJFMoleculeSpecifications``;
            - ``GJFMoleculeSpecifications``: used directly;
            - ``None``: falls back to ``self.molecule_specifications``.

        additional_sections : str | None, optional
            Raw additional-section text appended after the main molecule block,
            such as GIC, ModRedundant, NBO, or other Gaussian extra-input blocks.

            Actual behavior:
            - whenever this argument is explicitly provided, it has top priority;
            - even an empty string overrides the structured additional-section
              sources.

        parsed_additional_sections : list[...] | None, optional
            Structured additional sections. Each section is rendered via its own
            ``_render()`` and then joined with blank lines.

            Actual behavior:
            - only used when ``additional_sections is None``;
            - when provided, it takes priority over
              ``self.parsed_additional_sections``;
            - sections whose rendered text is blank are skipped.

            Additional-section priority is:
            1. ``additional_sections``
            2. ``parsed_additional_sections``
            3. ``self.parsed_additional_sections``
            4. ``self.additional_sections``
            5. empty string

        chk : str | bool | None, optional
            Whether to append an extra ``%chk`` Link0 keyword for this render.

            Actual behavior:
            - ``None`` / ``False``: do not append anything;
            - ``True``: append ``%chk=<basename>.chk``;
            - ``str``: append ``%chk=<given string>``.

            ``<basename>`` is taken from ``self.pure_filename`` when available,
            otherwise from the final title card. This is an append operation, not
            a replacement, and existing ``%chk`` entries are not deduplicated.

        old_chk : str | bool | None, optional
            Whether to append an extra ``%oldchk`` Link0 keyword for this render,
            typically to reference a previous checkpoint file.

            Actual behavior matches ``chk``:
            - ``None`` / ``False``: do not append anything;
            - ``True``: append ``%oldchk=<basename>.chk``;
            - ``str``: append ``%oldchk=<given string>``.

            This is also an append operation and does not replace or deduplicate
            existing ``%oldchk`` entries.

        coords_type : {"cartesian", "internal", "auto"}, default "auto"
            Controls the preferred output form of the coordinate block.

            Actual behavior in the current implementation:
            - ``"auto"``: preserve the stored representation. Each atom
              specification renders itself according to the structure of its own
              ``coords_part``: Cartesian coordinates, internal coordinates, or
              raw text.
            - ``"cartesian"``: require the molecule block to be rendered in
              Cartesian coordinates.
              - If all atoms in a fragment already store Cartesian coordinates,
                they are rendered directly in their current order.
              - If a fragment stores internal coordinates, the fragment is first
                converted through ``frag.coords()`` and then reconstructed as
                Cartesian atom lines.
              - If the fragment contains dummy atoms or ghost atoms, the
                conversion is rejected with ``ValueError`` because the current
                ``coords()`` / ``to_XYZ_block()`` path only returns coordinates
                for real atoms and cannot be aligned back to the original atom
                list safely.
            - ``"internal"``: require the molecule block to be rendered in
              internal-coordinate form.
              - If all atoms in a fragment already store internal coordinates,
                they are rendered directly in their current order.
              - If a fragment stores Cartesian coordinates, the fragment is
                converted with ``InternalCoords.from_cartesian_coords()`` and
                then reconstructed as Gaussian Z-matrix-style atom lines.
              - If the fragment contains dummy atoms or ghost atoms, the
                conversion is rejected with ``ValueError`` because the current
                Cartesian -> internal conversion path operates on the real-atom
                symbol/coordinate sequence returned by ``self.symbols()`` and
                ``self.coords()``, so it cannot be mapped back to the original
                mixed atom list safely.

            In other words, ``coords_type`` is now an active rendering control,
            but its supported scope is intentionally limited:
            - internal -> Cartesian is supported;
            - Cartesian -> internal is supported for ordinary real-atom
              fragments;
            - fragments containing dummy or ghost atoms are not converted
              implicitly in either direction.

        add_gjf_connectivity : bool, default False
            Whether to append a Gaussian connectivity section after the molecule
            coordinate block.

            Actual behavior:
            - ``False``: render only charge/multiplicity plus coordinates;
            - ``True``: append the connectivity information generated by
              ``self.connectivity()``.

        **kwargs
            Extra keyword arguments forwarded to downstream section, fragment,
            and atom ``_render()`` methods. This method consumes only a small
            subset directly; the rest are interpreted by child renderers if they
            recognize them.

        Returns
        -------
        str
            Complete Gaussian input text. The returned string always ends with
            two trailing newlines.

        Notes
        -----
        - This is a render-time override interface: supplied arguments affect
          only the current output and are not meant to be permanent field
          updates.
        - ``chk`` / ``old_chk`` append keywords only to a temporary Link0 copy
          used for this render.
        """
        render_parts = resolve_gjf_render_parts(
            link0_commands=link0_commands,
            route_section=route_section,
            title_card=title_card,
            molecule_specifications=molecule_specifications,
            fallback_link0_commands=self.link0_commands,
            fallback_route_section=self.route_section,
            fallback_title_card=self.title_card,
            fallback_molecule_specifications=self.molecule_specifications,
        )
        link0_commands_to_use = render_parts.link0_commands
        route_section_to_use = render_parts.route_section
        title_card_to_use = render_parts.title_card
        molecule_specifications_to_use = render_parts.molecule_specifications

        if title_card_to_use.title_card == "":
            title_card_to_use = title_card_to_use.model_copy(deep=True)
            title_card_to_use.title_card = getattr(self, "pure_filename", None) or "title"
        if additional_sections is not None:
            additional_sections_to_use = additional_sections
        elif parsed_additional_sections is not None:
            additional_sections_to_use = "\n\n".join(
                section._render(**kwargs)
                for section in parsed_additional_sections
                if section._render(**kwargs).strip()
            )
        elif self.parsed_additional_sections:
            additional_sections_to_use = "\n\n".join(
                section._render(**kwargs)
                for section in self.parsed_additional_sections
                if section._render(**kwargs).strip()
            )
        elif self.additional_sections:
            additional_sections_to_use = self.additional_sections
        else:
            additional_sections_to_use = ""

        chk_basename = getattr(self, "pure_filename", None) or title_card_to_use.title_card

        if chk or old_chk:
            link0_commands_to_use = link0_commands_to_use.model_copy(deep=True)
        if chk:
            link0_commands_to_use.add_link0_keyword(
                "chk", chk if isinstance(chk, str) else f"{chk_basename}.chk"
            )
        if old_chk:
            link0_commands_to_use.add_link0_keyword(
                "oldchk",
                old_chk if isinstance(old_chk, str) else f"{chk_basename}.chk",
            )
        render_kwargs = dict(kwargs)
        connectivity_text: str | None = None
        if add_gjf_connectivity:
            typed_self = cast(BaseQMInputFrame, self)
            topology_status = getattr(self, "topology_reconstruction_status", None)
            if topology_status not in {"failed", "suspicious_fallback"}:
                try:
                    rdmol = typed_self.rdmol
                    if rdmol is not None:
                        connectivity_text = rdmol_to_gjf_connectivity(rdmol)
                except NativeReconstructionConcurrencyError:
                    raise
                except Exception:
                    pass
            if connectivity_text is None:
                connectivity_text = ""
        return (
            link0_commands_to_use._render(**render_kwargs)
            + route_section_to_use._render(**render_kwargs)
            + title_card_to_use._render(**render_kwargs)
            + molecule_specifications_to_use._render(
                coords_type=coords_type,
                add_gjf_connectivity=add_gjf_connectivity,
                connectivity_text=connectivity_text,
                **render_kwargs,
            )
            + additional_sections_to_use
            + "\n\n"
        )

    @model_validator(mode="after")
    def set_properties(self):
        typed_self = cast(BaseQMInputFrame, self)
        semantic_route = self.route_section.semantic_route
        typed_self.qm_software = "Gaussian"
        typed_self.qm_software_version = "Any"
        typed_self.keywords = self.route_section.route
        typed_self.resources_raw = self.link0_commands._render()
        typed_self.request_num_cpu = self.link0_commands.cpu_request()
        typed_self.request_memory = self.link0_commands.memory_request()
        typed_self.charge = self.molecule_specifications.total_charge
        typed_self.multiplicity = self.molecule_specifications.spin_multiplicity
        typed_self.atoms = self.molecule_specifications.atomic_numbers()
        typed_self.coords = self.molecule_specifications.coords()
        populate_gaussian_legacy_qm_fields_from_semantic(typed_self, semantic_route)
        return self


class GJFFileFrameMemory(
    MemoryStorageMixin, GJFFileFrameMixin, BaseQMInputFrame["GJFFileFrameMemory"]
): ...


class GJFFileFrameDisk(
    DiskStorageMixin, GJFFileFrameMixin, BaseQMInputFrame["GJFFileFrameDisk"]
): ...
