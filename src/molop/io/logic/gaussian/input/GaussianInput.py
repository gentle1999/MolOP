from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, ClassVar, Literal, cast

import numpy as np
from molgr.interface import xyz_to_rdmol
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field, model_validator
from rdkit import Chem

from molop.config import molopconfig, moloplogger
from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.DataClasses import AtomInInternalCoords, InternalCoords
from molop.io.logic.gaussian.input.GaussianLink0 import (
    GJFLink0 as GJFLink0,
)
from molop.io.logic.gaussian.input.GaussianLink0 import (
    GJFLink0Commands as GJFLink0Commands,
)
from molop.io.logic.gaussian.input.GaussianRoute import GaussianRouteSemantic
from molop.structure.FormatConverter import rdmol_to_gjf_connectivity
from molop.structure.GeometryTransformation import merge_mols_directly
from molop.unit import atom_ureg
from molop.utils.progressbar import (
    NativeReconstructionConcurrencyError,
    native_reconstruction_guard,
)
from molop.utils.types import PintArrayNx3


pt = Chem.GetPeriodicTable()


def _gaussian_element_symbol_from_label(element_label: str) -> str:
    if len(element_label) >= 2 and element_label[0].isupper() and element_label[1].isdigit():
        return element_label[0]
    if (
        len(element_label) >= 3
        and element_label[0].isupper()
        and element_label[1].islower()
        and element_label[2].isdigit()
    ):
        return element_label[:2]
    return element_label


class GJFRouteSection(BaseDataClassWithUnit):
    route: str = Field(default="", description="Route section")
    semantic_route: GaussianRouteSemantic = Field(
        default_factory=GaussianRouteSemantic, description="Structured Gaussian route semantics"
    )

    @classmethod
    def from_str(cls, data: str) -> GJFRouteSection:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_route_section,
        )

        return parse_gjf_route_section(data)

    @model_validator(mode="before")
    @classmethod
    def normalize_route_payload(cls, data: Any):
        if not isinstance(data, Mapping):
            return data
        payload = dict(data)
        semantic = payload.get("semantic_route")
        if not payload.get("route") and isinstance(semantic, GaussianRouteSemantic):
            payload["route"] = semantic.raw_route or semantic.normalized_route
        return payload

    @model_validator(mode="after")
    def validate_route(self):
        if self.route and not self.route.startswith("#"):
            self.route = f"# {self.route}".replace(" \n", "").replace("\n", "")
        return self

    def _render(self, **kwargs) -> str:
        return self.route + "\n\n"

    def to_dict(self) -> dict[str, Any]:
        return self.semantic_route.to_route_dict()

    def dieze_tag(self) -> str | None:
        return self.semantic_route.dieze_tag


class GJFTitleCard(BaseDataClassWithUnit):
    # follow G16 manual (https://gaussian.com/input/)
    INVALID_CHARS: ClassVar[list[str]] = ["@", "#", "!", "-", "_", r"\\"]

    title_card: str = Field(default="", description="Title card")

    @classmethod
    def from_str(cls, data: str) -> GJFTitleCard:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_title_card,
        )

        return parse_gjf_title_card(data)

    @model_validator(mode="after")
    def validate_title_card(self):
        for char in self.INVALID_CHARS:
            self.title_card = self.title_card.replace(char, " ").strip("\n")
        if self.title_card:
            line_count = len(self.title_card.splitlines())
            if line_count > 5:
                raise ValueError("Title card cannot exceed 5 lines")
        return self

    def _render(self, **kwargs) -> str:
        return self.title_card + "\n\n"


class GJFAtomSpecification(BaseDataClassWithUnit):
    element_label: str = Field(description="Element label")
    atom_type: str | None = Field(default=None, description="Atom type")
    charge: float | None = Field(default=None, description="Atom charge")
    params: dict[str, str] = Field(default={}, description="Atom parameters")
    frozen_tag: str | None = Field(default=None, description="Frozen tag")
    coords_part: str = Field(
        description="Atom coordinates part, may be cartesian or internal coordinates"
    )

    @classmethod
    def from_str(cls, data: str) -> GJFAtomSpecification:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_atom_specification,
        )

        return parse_gjf_atom_specification(data)

    def _render(self, **kwargs) -> str:
        atom_part = (
            f"{self.element_label}"
            + (f"-{self.atom_type}" if self.atom_type is not None else "")
            + (f"-{self.charge}" if self.charge is not None else "")
            + (
                "(" + ",".join([f"{k}={v}" for k, v in self.params.items()]) + ")"
                if self.params
                else ""
            )
            + (f" {self.frozen_tag}" if self.frozen_tag is not None else "")
        )
        if self.is_cartesian_coords():
            coord_part = "".join([f"{float(coord):15.6f}" for coord in self.coords_part.split()])
        elif self.is_internal_coords():
            values = self.coords_part.split()
            format_flag: int | None = None
            if len(values) % 2 == 1 and len(values) >= 3 and values[-1] in {"0", "1"}:
                format_flag = int(values[-1])
                values = values[:-1]

            pair_count = len(values) // 2
            pairs: list[tuple[int, float]] = []
            for pair_idx in range(pair_count):
                idx_token = values[pair_idx * 2]
                value_token = values[pair_idx * 2 + 1]
                try:
                    idx_val = int(idx_token)
                except Exception:
                    idx_val = 0
                try:
                    number_val = float(value_token)
                except Exception:
                    number_val = 0.0
                pairs.append((idx_val, number_val))

            coord_part = "".join([f"{idx:5d}{value:15.6f}" for (idx, value) in pairs])
            if format_flag is not None:
                coord_part += f"{format_flag:5d}"
        else:
            coord_part = self.coords_part
        return f"{atom_part:6s}{coord_part}\n"

    def get_fragment_id(self) -> int:
        if self.params is None:
            return 0
        return int(self.params.get("fragment", "0"))

    def is_cartesian_coords(self) -> bool:
        return len(self.coords_part.split()) == 3

    def is_internal_coords(self) -> bool:
        values = self.coords_part.split()
        if len(values) < 2:
            return False
        if len(values) == 3:
            return False
        if len(values) % 2 == 0:
            return True
        if len(values) >= 3 and values[-1] in {"0", "1"}:
            return len(values[:-1]) % 2 == 0
        return False

    @property
    def symbol(self) -> str:
        return _gaussian_element_symbol_from_label(self.element_label)

    @property
    def cartesian_coords(self) -> NumpyQuantity:
        if not self.is_cartesian_coords():
            raise ValueError("Atom coordinates must be cartesian coordinates")
        return np.array([float(coord) for coord in self.coords_part.split()]) * atom_ureg.angstrom

    @property
    def internal_coords(self) -> tuple[int, int, int, PlainQuantity, PlainQuantity, PlainQuantity]:
        if not self.is_internal_coords():
            raise ValueError("Atom coordinates must be internal coordinates")
        raw_values = self.coords_part.split()
        format_flag = (
            int(raw_values[-1]) if len(raw_values) % 2 == 1 and raw_values[-1] in {"0", "1"} else 0
        )
        values = (
            raw_values[:-1] if format_flag in {0, 1} and len(raw_values) % 2 == 1 else raw_values
        )
        indices: list[int] = []
        for i in (0, 2, 4):
            try:
                indices.append(int(values[i]))
            except Exception:
                indices.append(0)
        distance: PlainQuantity
        angle: PlainQuantity
        dihedral: PlainQuantity
        try:
            distance = cast(PlainQuantity, float(values[1]) * atom_ureg.angstrom)
        except Exception:
            distance = cast(PlainQuantity, 0 * atom_ureg.angstrom)
        try:
            angle = cast(PlainQuantity, float(values[3]) * atom_ureg.degree)
        except Exception:
            angle = cast(PlainQuantity, 0 * atom_ureg.degree)
        try:
            dihedral = cast(PlainQuantity, float(values[5]) * atom_ureg.degree)
        except Exception:
            dihedral = cast(PlainQuantity, 0 * atom_ureg.degree)

        pair_count = len(values) // 2

        if distance.to(atom_ureg.angstrom).m <= 0:
            raise ValueError("Z-matrix bond length must be positive")
        if pair_count >= 2 and not 0 < angle.to(atom_ureg.degree).m < 180:
            raise ValueError("Z-matrix bond angle must be between 0 and 180 degrees")
        if pair_count >= 3 and format_flag == 1 and not 0 < dihedral.to(atom_ureg.degree).m < 180:
            raise ValueError("Alternate Z-matrix second angle must be between 0 and 180 degrees")
        return indices[0], indices[1], indices[2], distance, angle, dihedral

    @property
    def is_dummy(self) -> bool:
        return self.element_label.upper() == "X"

    @property
    def is_ghost(self) -> bool:
        if self.atom_type is None:
            return False
        return self.atom_type.upper() == "BQ"

    def to_atom_in_internal_coordinates_coords(self) -> AtomInInternalCoords:
        if not self.coords_part.strip():
            return AtomInInternalCoords(
                symbol=self.element_label,
                is_dummy=self.is_dummy,
                is_ghost=self.is_ghost,
            )

        internal_coords = self.internal_coords
        values = self.coords_part.split()
        zmat_format = int(values[-1]) if len(values) % 2 == 1 and values[-1] in {"0", "1"} else 0

        def _to_zero_based(index: int) -> int:
            return max(index - 1, 0) if index > 0 else 0

        return AtomInInternalCoords(
            symbol=self.element_label,
            distance_to_index=_to_zero_based(internal_coords[0]),
            distance=internal_coords[3],
            angle_to_index=_to_zero_based(internal_coords[1]),
            angle=internal_coords[4],
            dihedral_to_index=_to_zero_based(internal_coords[2]),
            dihedral=internal_coords[5],
            zmat_format=zmat_format,
            is_dummy=self.is_dummy,
            is_ghost=self.is_ghost,
        )

    @property
    def atomic_number(self) -> int:
        return pt.GetAtomicNumber(self.symbol)


class GJFMoleculeSpecificationsFragment(BaseDataClassWithUnit):
    fragment_id: int = Field(default=0, description="Fragment ID")
    total_charge: int = Field(default=0, description="Total charge")
    spin_multiplicity: int = Field(default=1, description="Spin multiplicity")
    atom_specifications: list[GJFAtomSpecification] = Field(
        default_factory=list, description="Atom specifications"
    )

    def _all_atoms_cartesian(self) -> bool:
        return all(
            (not atom_spec.coords_part.strip()) or atom_spec.is_cartesian_coords()
            for atom_spec in self.atom_specifications
        )

    def _all_atoms_internal(self) -> bool:
        return all(
            (not atom_spec.coords_part.strip()) or atom_spec.is_internal_coords()
            for atom_spec in self.atom_specifications
        )

    def _render_cartesian(self, **kwargs) -> str:
        _ = kwargs
        if self._all_atoms_cartesian():
            return "".join([atom_spec._render(**kwargs) for atom_spec in self.atom_specifications])

        if any(atom_spec.is_dummy or atom_spec.is_ghost for atom_spec in self.atom_specifications):
            raise ValueError(
                "coords_type='cartesian' cannot convert internal coordinates when dummy or ghost atoms are present"
            )

        coords = self.coords().to(atom_ureg.angstrom).magnitude
        rendered_lines: list[str] = []
        coord_idx = 0
        for atom_spec in self.atom_specifications:
            if not atom_spec.coords_part.strip():
                rendered_lines.append(atom_spec._render(**kwargs))
                continue
            x, y, z = coords[coord_idx]
            coord_idx += 1
            rendered_atom = atom_spec.model_copy(update={"coords_part": f"{x:.6f} {y:.6f} {z:.6f}"})
            rendered_lines.append(rendered_atom._render(**kwargs))
        return "".join(rendered_lines)

    def _render_internal(self, **kwargs) -> str:
        _ = kwargs
        if self._all_atoms_internal():
            return "".join([atom_spec._render(**kwargs) for atom_spec in self.atom_specifications])

        if any(atom_spec.is_dummy or atom_spec.is_ghost for atom_spec in self.atom_specifications):
            raise ValueError(
                "coords_type='internal' cannot convert cartesian coordinates when dummy or ghost atoms are present"
            )

        internal = self.to_internal_coords()
        rendered_lines: list[str] = []
        internal_idx = 0
        for atom_spec in self.atom_specifications:
            if not atom_spec.coords_part.strip():
                rendered_lines.append(atom_spec._render(**kwargs))
                continue

            atom_internal = internal.atoms[internal_idx]
            internal_idx += 1

            tokens: list[str] = []
            if internal_idx >= 2:
                tokens.extend(
                    [
                        str(atom_internal.distance_to_index + 1),
                        f"{atom_internal.distance.to(atom_ureg.angstrom).magnitude:.6f}",
                    ]
                )
            if internal_idx >= 3:
                tokens.extend(
                    [
                        str(atom_internal.angle_to_index + 1),
                        f"{atom_internal.angle.to(atom_ureg.degree).magnitude:.6f}",
                    ]
                )
            if internal_idx >= 4:
                tokens.extend(
                    [
                        str(atom_internal.dihedral_to_index + 1),
                        f"{atom_internal.dihedral.to(atom_ureg.degree).magnitude:.6f}",
                    ]
                )
                if atom_internal.zmat_format in {0, 1}:
                    tokens.append(str(atom_internal.zmat_format))

            rendered_atom = atom_spec.model_copy(update={"coords_part": " ".join(tokens)})
            rendered_lines.append(rendered_atom._render(**kwargs))
        return "".join(rendered_lines)

    def _render(self, **kwargs) -> str:
        coords_type = kwargs.get("coords_type", "auto")
        if coords_type == "cartesian":
            return self._render_cartesian(**kwargs)
        if coords_type == "internal":
            return self._render_internal(**kwargs)
        return "".join([atom_spec._render(**kwargs) for atom_spec in self.atom_specifications])

    def to_internal_coords(self) -> InternalCoords:
        return InternalCoords(
            items=[
                atom_spec.to_atom_in_internal_coordinates_coords()
                for atom_spec in self.atom_specifications
            ]
        )

    def fragment_molecule(self) -> Chem.rdchem.Mol | None:
        try:
            with native_reconstruction_guard():
                return xyz_to_rdmol(
                    self.to_XYZ_block(),
                    total_charge=self.total_charge,
                    spin_multiplicity=self.spin_multiplicity,
                    backend=molopconfig.graph_reconstruction_backend,
                    make_dative_bonds=molopconfig.make_dative_bonds,
                    make_stereochemistry=molopconfig.make_stereochemistry,
                )
        except NativeReconstructionConcurrencyError:
            # A concurrency boundary violation must reach the caller.  Turning
            # it into an empty fragment would silently produce an incomplete
            # Gaussian graph and hide an unsafe scheduling decision.
            raise
        except Exception as e:
            moloplogger.error(f"{e}")
            return None

    def symbols(self) -> list[str]:
        return [
            atom_spec.symbol
            for atom_spec in self.atom_specifications
            if not atom_spec.is_dummy and not atom_spec.is_ghost
        ]

    def atomic_numbers(self) -> list[int]:
        return [
            atom_spec.atomic_number
            for atom_spec in self.atom_specifications
            if not atom_spec.is_dummy and not atom_spec.is_ghost
        ]

    def coords(self) -> PintArrayNx3:
        if self.atom_specifications[0].is_cartesian_coords():
            return (
                np.array(
                    [
                        atom_spec.cartesian_coords.magnitude
                        for atom_spec in self.atom_specifications
                        if not atom_spec.is_dummy and not atom_spec.is_ghost
                    ]
                )
                * atom_ureg.angstrom
            )
        else:
            return self.to_internal_coords().to_cartesian_coords()

    def to_XYZ_block(self) -> str:
        symbols = self.symbols()
        coords = self.coords().magnitude
        return f"{len(symbols)}\n\n" + "\n".join(
            [
                f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
                for symbol, (x, y, z) in zip(symbols, coords, strict=True)
            ]
        )


class GJFMoleculeSpecifications(BaseDataClassWithUnit):
    total_charge: int = Field(default=0, description="Total charge")
    spin_multiplicity: int = Field(default=1, description="Spin multiplicity")

    molecule_fragments: list[GJFMoleculeSpecificationsFragment] = Field(
        default_factory=list, description="Molecule fragments"
    )

    @classmethod
    def from_str(cls, data: str) -> GJFMoleculeSpecifications:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_molecule_specifications,
        )

        return parse_gjf_molecule_specifications(data)

    def _render(
        self,
        add_gjf_connectivity: bool = False,
        connectivity_text: str | None = None,
        **kwargs: Any,
    ) -> str:
        fragment_kwargs = dict(kwargs)

        if len(self.molecule_fragments) == 1:
            spin_multiplicity = f"{self.total_charge} {self.spin_multiplicity}\n"
            molecule_fragment_part = self.molecule_fragments[0]._render(**fragment_kwargs) + "\n"
            connectivity_part = (
                f"{connectivity_text if connectivity_text is not None else self.connectivity()}\n\n"
                if add_gjf_connectivity
                else ""
            )
            return f"{spin_multiplicity}{molecule_fragment_part}{connectivity_part}"
        if len(self.molecule_fragments) > 1:
            spin_multiplicity = (
                f"{self.total_charge} {self.spin_multiplicity}"
                + "".join(
                    f" {frag.total_charge} {frag.spin_multiplicity}"
                    for frag in self.molecule_fragments
                )
                + "\n"
            )
            molecule_fragment_part = (
                "".join([frag._render(**fragment_kwargs) for frag in self.molecule_fragments])
                + "\n"
            )
            connectivity_part = (
                f"{connectivity_text if connectivity_text is not None else self.connectivity()}\n\n"
                if add_gjf_connectivity
                else ""
            )
            return f"{spin_multiplicity}{molecule_fragment_part}{connectivity_part}"
        return ""

    def symbols(self) -> list[str]:
        return [
            atom_spec.symbol
            for frag in self.molecule_fragments
            for atom_spec in frag.atom_specifications
        ]

    def atomic_numbers(self) -> list[int]:
        return [
            atom_spec.atomic_number
            for frag in self.molecule_fragments
            for atom_spec in frag.atom_specifications
        ]

    def coords(self) -> PintArrayNx3:
        if len(self.molecule_fragments) == 0:
            return np.zeros((0, 3)) * atom_ureg.angstrom
        return cast(
            NumpyQuantity,
            np.concatenate([frag.coords() for frag in self.molecule_fragments], axis=0),
        )

    def to_XYZ_block(self) -> str:
        symbols = self.symbols()
        coords = self.coords().magnitude
        return f"{len(symbols)}\n\n" + "\n".join(
            [
                f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
                for symbol, (x, y, z) in zip(symbols, coords, strict=True)
            ]
        )

    @property
    def rdmol_fragments(self) -> list[Chem.rdchem.Mol]:
        return [
            rdmol
            for frag in self.molecule_fragments
            if (rdmol := frag.fragment_molecule()) is not None
        ]

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
        return merge_mols_directly(self.rdmol_fragments)

    def connectivity(self) -> str:
        return rdmol_to_gjf_connectivity(self.rdmol)


def build_gjf_molecule_specifications_from_frame_payload(
    payload: Mapping[str, Any],
) -> GJFMoleculeSpecifications | None:
    atoms = payload.get("atoms")
    coords = payload.get("coords")
    if not isinstance(atoms, list) or not atoms:
        return None

    coord_values: list[list[float]] | None = None
    if isinstance(coords, np.ndarray):
        coord_values = coords.tolist()
    elif coords is not None and hasattr(coords, "magnitude"):
        coord_values = np.asarray(cast(Any, coords).magnitude).tolist()
    elif isinstance(coords, list):
        coord_values = coords

    if not isinstance(coord_values, list) or len(coord_values) != len(atoms):
        return None

    atom_specifications: list[GJFAtomSpecification] = []
    for atom, row in zip(atoms, coord_values, strict=True):
        if not isinstance(row, list) or len(row) < 3:
            return None
        try:
            x, y, z = float(row[0]), float(row[1]), float(row[2])
        except Exception:
            return None

        if isinstance(atom, str):
            symbol = atom
        else:
            try:
                symbol = pt.GetElementSymbol(int(atom))
            except Exception:
                return None
        atom_specifications.append(
            GJFAtomSpecification(
                element_label=symbol,
                coords_part=f"{x:.6f} {y:.6f} {z:.6f}",
            )
        )

    total_charge = int(payload.get("charge", 0) or 0)
    spin_multiplicity = int(payload.get("multiplicity", 1) or 1)
    return GJFMoleculeSpecifications(
        total_charge=total_charge,
        spin_multiplicity=spin_multiplicity,
        molecule_fragments=[
            GJFMoleculeSpecificationsFragment(
                fragment_id=0,
                total_charge=total_charge,
                spin_multiplicity=spin_multiplicity,
                atom_specifications=atom_specifications,
            )
        ],
    )


def normalize_gjf_frame_geometry_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    if normalized.get("molecule_specifications") is not None:
        return normalized

    built = build_gjf_molecule_specifications_from_frame_payload(normalized)
    if built is None:
        return normalized

    normalized["molecule_specifications"] = built
    coords = normalized.get("coords")
    if isinstance(coords, np.ndarray):
        normalized["coords"] = coords * atom_ureg.angstrom
    elif isinstance(coords, list):
        normalized["coords"] = np.asarray(coords, dtype=float) * atom_ureg.angstrom
    return normalized


def adapt_gjf_writer_payload(data: Mapping[str, Any]) -> dict[str, Any]:
    payload = dict(data)

    link0_commands = payload.get("link0_commands")
    if isinstance(link0_commands, str):
        payload["link0_commands"] = GJFLink0Commands.from_str(link0_commands)
    elif link0_commands is None:
        resources_raw = payload.get("resources_raw")
        if isinstance(resources_raw, str) and resources_raw.strip():
            payload["link0_commands"] = GJFLink0Commands.from_str(resources_raw)

    route_section = payload.get("route_section")
    if isinstance(route_section, str):
        if route_section.strip():
            payload["route_section"] = GJFRouteSection.from_str(route_section)
        else:
            payload.pop("route_section", None)
    elif route_section is None:
        keywords = payload.get("keywords")
        if isinstance(keywords, str) and keywords.strip():
            payload["route_section"] = GJFRouteSection.from_str(keywords)

    title_card = payload.get("title_card")
    if isinstance(title_card, str):
        payload["title_card"] = GJFTitleCard.from_str(title_card)

    molecule_specifications = payload.get("molecule_specifications")
    if isinstance(molecule_specifications, str):
        if molecule_specifications.strip():
            payload["molecule_specifications"] = GJFMoleculeSpecifications.from_str(
                molecule_specifications
            )
        else:
            payload.pop("molecule_specifications", None)

    return payload


@dataclass(frozen=True, slots=True)
class GJFRenderParts:
    link0_commands: GJFLink0Commands
    route_section: GJFRouteSection
    title_card: GJFTitleCard
    molecule_specifications: GJFMoleculeSpecifications


def resolve_gjf_render_parts(
    *,
    link0_commands: str | GJFLink0Commands | dict[str, str | None] | None,
    route_section: str | GJFRouteSection | None,
    title_card: str | GJFTitleCard | None,
    molecule_specifications: str | GJFMoleculeSpecifications | None,
    fallback_link0_commands: GJFLink0Commands,
    fallback_route_section: GJFRouteSection,
    fallback_title_card: GJFTitleCard,
    fallback_molecule_specifications: GJFMoleculeSpecifications,
) -> GJFRenderParts:
    if isinstance(link0_commands, str):
        link0_commands_to_use = GJFLink0Commands.from_str(link0_commands)
    elif isinstance(link0_commands, GJFLink0Commands):
        link0_commands_to_use = link0_commands
    elif isinstance(link0_commands, dict):
        link0_commands_to_use = GJFLink0Commands.from_dict(link0_commands)
    elif link0_commands is None:
        link0_commands_to_use = fallback_link0_commands
    else:
        raise ValueError(f"Invalid type for link0_commands: {type(link0_commands)}")

    if isinstance(route_section, str):
        route_section_to_use = GJFRouteSection.from_str(route_section)
    elif isinstance(route_section, GJFRouteSection):
        route_section_to_use = route_section
    elif route_section is None:
        route_section_to_use = fallback_route_section
    else:
        raise ValueError(f"Invalid type for route_section: {type(route_section)}")

    if isinstance(title_card, str):
        title_card_to_use = GJFTitleCard.from_str(title_card)
    elif isinstance(title_card, GJFTitleCard):
        title_card_to_use = title_card
    elif title_card is None:
        title_card_to_use = fallback_title_card
    else:
        raise ValueError(f"Invalid type for title_card: {type(title_card)}")

    if isinstance(molecule_specifications, str):
        molecule_specifications_to_use = GJFMoleculeSpecifications.from_str(molecule_specifications)
    elif isinstance(molecule_specifications, GJFMoleculeSpecifications):
        molecule_specifications_to_use = molecule_specifications
    elif molecule_specifications is None:
        molecule_specifications_to_use = fallback_molecule_specifications
    else:
        raise ValueError(
            f"Invalid type for molecule_specifications: {type(molecule_specifications)}"
        )

    return GJFRenderParts(
        link0_commands=link0_commands_to_use,
        route_section=route_section_to_use,
        title_card=title_card_to_use,
        molecule_specifications=molecule_specifications_to_use,
    )


class GJFSectionParsingDiagnostic(BaseDataClassWithUnit):
    section_type: str = Field(description="Detected additional section type")
    message: str = Field(description="Diagnostic message")
    section_index: int = Field(description="0-based additional section index")


class GJFGICOption(BaseDataClassWithUnit):
    raw: str = Field(description="Raw option token")
    key: str = Field(description="Normalized option key")
    value: str | None = Field(default=None, description="Optional option value")
    is_flag: bool = Field(default=True, description="Whether the option is a bare flag")


class GJFGICLine(BaseDataClassWithUnit):
    raw_line: str = Field(description="Original GIC line")
    label: str | None = Field(default=None, description="Optional coordinate label")
    label_options: list[str] = Field(default_factory=list, description="Label options")
    expression: str = Field(default="", description="Raw GIC expression")
    is_standalone_option: bool = Field(default=False, description="Standalone global option")
    expression_kind: str = Field(default="raw", description="Parsed expression kind")
    function_name: str | None = Field(default=None, description="Function-like expression name")
    function_args: list[str] = Field(
        default_factory=list, description="Function-like expression args"
    )
    standalone_keyword: str | None = Field(default=None, description="Standalone option keyword")
    standalone_args: list[str] = Field(default_factory=list, description="Standalone option args")
    parsed_label_options: list[GJFGICOption] = Field(
        default_factory=list, description="Structured label options"
    )
    standalone_action: str | None = Field(default=None, description="Normalized standalone action")
    standalone_target: str | None = Field(default=None, description="Standalone target identifier")
    normalized_state: str | None = Field(default=None, description="Normalized semantic state")
    option_values: dict[str, str] = Field(
        default_factory=dict, description="Normalized key/value options"
    )

    @classmethod
    def from_str(cls, data: str) -> GJFGICLine:
        from molop.io.logic.gaussian.input.GaussianInputParsing import parse_gjf_gic_line

        return parse_gjf_gic_line(data)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        if self.is_standalone_option:
            tokens = []
            if self.standalone_keyword:
                tokens.append(self.standalone_keyword)
            tokens.extend(self.standalone_args)
            return " ".join(tokens) if tokens else self.expression

        if self.function_name is not None:
            expression = f"{self.function_name}({', '.join(self.function_args)})"
        else:
            expression = self.expression

        if self.label is None:
            return expression

        option_tokens = [opt.raw for opt in self.parsed_label_options] or self.label_options
        if option_tokens:
            return f"{self.label}({', '.join(option_tokens)})={expression}"
        return f"{self.label}={expression}"


class GJFGICSection(BaseDataClassWithUnit):
    section_type: Literal["gic"] = Field(default="gic", description="Additional section type")
    raw: str = Field(description="Raw GIC section")
    lines: list[GJFGICLine] = Field(default_factory=list, description="Parsed GIC lines")

    @classmethod
    def from_str(cls, data: str) -> GJFGICSection:
        from molop.io.logic.gaussian.input.GaussianInputParsing import parse_gjf_gic_section

        return parse_gjf_gic_section(data)

    def _render(self, **kwargs) -> str:
        if self.lines:
            return "\n".join(line._render(**kwargs) for line in self.lines)
        _ = kwargs
        return self.raw


class GJFModRedundantLine(BaseDataClassWithUnit):
    raw_line: str = Field(description="Original ModRedundant line")
    coordinate_type: str | None = Field(default=None, description="Coordinate descriptor")
    atom_refs: list[str] = Field(default_factory=list, description="Atom references")
    action: str | None = Field(default=None, description="ModRedundant action")
    parameters: list[str] = Field(default_factory=list, description="Remaining parameters")

    @classmethod
    def from_str(cls, data: str) -> GJFModRedundantLine:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_modredundant_line,
        )

        return parse_gjf_modredundant_line(data)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        tokens = []
        if self.coordinate_type is not None:
            tokens.append(self.coordinate_type)
        tokens.extend(self.atom_refs)
        if self.action is not None:
            tokens.append(self.action)
        tokens.extend(self.parameters)
        return " ".join(tokens) if tokens else self.raw_line


class GJFModRedundantSection(BaseDataClassWithUnit):
    section_type: Literal["modredundant"] = Field(
        default="modredundant", description="Additional section type"
    )
    raw: str = Field(description="Raw ModRedundant section")
    lines: list[GJFModRedundantLine] = Field(
        default_factory=list, description="Parsed ModRedundant lines"
    )

    @classmethod
    def from_str(cls, data: str) -> GJFModRedundantSection:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_modredundant_section,
        )

        return parse_gjf_modredundant_section(data)

    def _render(self, **kwargs) -> str:
        if self.lines:
            return "\n".join(line._render(**kwargs) for line in self.lines)
        _ = kwargs
        return self.raw


class GJFUnknownSection(BaseDataClassWithUnit):
    section_type: Literal["unknown"] = Field(
        default="unknown", description="Additional section type"
    )
    raw: str = Field(description="Raw section text")

    @classmethod
    def from_str(cls, data: str) -> GJFUnknownSection:
        return cls(raw=data)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        return self.raw


class GJFNBOSection(BaseDataClassWithUnit):
    section_type: Literal["nbo"] = Field(default="nbo", description="Additional section type")
    raw: str = Field(description="Raw NBO section")
    header: str = Field(default="$nbo", description="NBO section header")
    commands: list[str] = Field(default_factory=list, description="NBO command lines")
    footer: str = Field(default="$end", description="NBO section footer")

    @classmethod
    def from_str(cls, data: str) -> GJFNBOSection:
        from molop.io.logic.gaussian.input.GaussianInputParsing import parse_gjf_nbo_section

        return parse_gjf_nbo_section(data)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        parts = [self.header]
        parts.extend(self.commands)
        parts.append(self.footer)
        return "\n".join(parts)
