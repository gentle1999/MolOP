"""
Author: TMJ
Date: 2025-07-28 23:09:31
LastEditors: TMJ
LastEditTime: 2026-04-05 17:29:39
Description: 请填写简介
"""

import re
from collections.abc import Mapping
from typing import Any, ClassVar, Literal, cast

import numpy as np
from molgr.interface import xyz_to_rdmol
from pint.facets.numpy.quantity import NumpyQuantity
from pint.facets.plain import PlainQuantity
from pydantic import Field, model_validator
from rdkit import Chem

from molop.config import molopconfig, moloplogger
from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.DataClasses import AtomInInternalCoords, InternalCoords
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.io.logic.gaussian_route_models import (
    GaussianRouteSemantic,
    parse_gaussian_route_semantic,
)
from molop.io.patterns.G16Patterns import options_parser
from molop.structure.FormatConverter import rdmol_to_gjf_connectivity
from molop.structure.GeometryTransformation import merge_mols_directly
from molop.unit import atom_ureg


pt = Chem.GetPeriodicTable()


class GJFLink0(BaseDataClassWithUnit):
    key: str = Field(description="Link0 keyword key")
    value: str | None = Field(description="Link0 keyword value")

    def _render(self, **kwargs) -> str:
        return f"%{self.key}={self.value}\n" if self.value else f"%{self.key}\n"

    def memory_request(self) -> PlainQuantity | None:
        if self.key.lower() == "mem":
            if self.value and self.value.endswith("B"):
                return atom_ureg.Quantity(self.value)
            if self.value and self.value.endswith("W"):
                return atom_ureg.Quantity(self.value.replace("W", "B")) * 8

    def cpu_request(self) -> int | None:
        if self.key.lower() == "cpu" and self.value:
            return len(self.value.split(","))
        if (self.key.lower() == "nproc" or self.key.lower() == "nprocshared") and self.value:
            return int(self.value)


class GJFLink0Commands(BaseDataClassWithUnit):
    link0_keywords: list[GJFLink0] = Field(default_factory=list, description="Link0 keywords")

    @classmethod
    def from_dict(cls, data: dict[str, str | None]) -> "GJFLink0Commands":
        link0_keywords: list[GJFLink0] = []
        for key, value in data.items():
            link0_keywords.append(GJFLink0(key=key, value=value))
        return cls(link0_keywords=link0_keywords)

    @classmethod
    def from_str(cls, data: str) -> "GJFLink0Commands":
        return cls.from_dict(options_parser(data))

    def _render(self, **kwargs) -> str:
        return "".join([link0._render() for link0 in self.link0_keywords])

    def memory_request(self) -> PlainQuantity:
        for link0 in self.link0_keywords:
            if (memory_request := link0.memory_request()) is not None:
                return memory_request
        return atom_ureg.Quantity("800MB")

    def cpu_request(self) -> int:
        for link0 in self.link0_keywords:
            if (cpu_request := link0.cpu_request()) is not None:
                return cpu_request
        return 1

    def add_link0_keyword(self, key: str, value: str | None = None) -> None:
        self.link0_keywords.append(GJFLink0(key=key, value=value))


class GJFRouteSection(BaseDataClassWithUnit):
    route: str = Field(default="", description="Route section")
    semantic_route: GaussianRouteSemantic = Field(
        default_factory=GaussianRouteSemantic, description="Structured Gaussian route semantics"
    )

    @classmethod
    def from_str(cls, data: str) -> "GJFRouteSection":
        return cls(route=data)

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
        if not self.route.startswith("#"):
            self.route = f"# {self.route}".replace(" \n", "").replace("\n", "")
        self.semantic_route = parse_gaussian_route_semantic(self.route)
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
    def from_str(cls, data: str) -> "GJFTitleCard":
        return cls(title_card=data)

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
    def _normalize_free_separators(cls, data: str) -> str:
        normalized_chars: list[str] = []
        depth = 0
        for ch in data:
            if ch in "([{":
                depth += 1
                normalized_chars.append(ch)
                continue
            if ch in ")]}":
                depth = max(depth - 1, 0)
                normalized_chars.append(ch)
                continue
            if depth == 0 and ch in {",", "/", "\t"}:
                normalized_chars.append(" ")
                continue
            normalized_chars.append(ch)
        return "".join(normalized_chars)

    @classmethod
    def from_str(cls, data: str) -> "GJFAtomSpecification":
        normalized = cls._normalize_free_separators(data)
        parts = normalized.split()
        if len(parts) < 2:
            return cls(element_label=parts[0], coords_part="")
        if re.match(r"^[-+]?\d+(\.\d+)?$", parts[1]):
            coords_part = " ".join(parts[1:])
            frozen_tag = None
        else:
            frozen_tag = parts[1]
            coords_part = " ".join(parts[2:])
        atom_info = parts[0]
        matches = re.match(
            r"^([A-Z][A-Za-z0-9]*)(?:-([A-Za-z][A-Za-z0-9]*)(?:-(-?\d+(?:\.\d+)?))?)?(?:\(([A-Za-z][A-Za-z0-9]*=\w+(?:,[A-Za-z][A-Za-z0-9]*=\w+)*)\))?$",
            atom_info,
        )
        if not matches:
            raise ValueError("Atom specification must match the pattern")
        element_label, atom_type, charge, params_str = matches.groups()
        if params_str:
            params = {
                sec.split("=")[0].lower(): sec.split("=")[1].lower()
                for sec in params_str.split(",")
            }
        else:
            params = {}
        return cls(
            element_label=element_label,
            atom_type=atom_type,
            charge=float(charge) if charge else None,
            params=params,
            frozen_tag=frozen_tag,
            coords_part=coords_part,
        )

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
        if len(values) % 2 == 0:
            return True
        if len(values) >= 3 and values[-1] in {"0", "1"}:
            return len(values[:-1]) % 2 == 0
        return False

    @property
    def symbol(self) -> str:
        if matches := re.match(r"([A-Z][a-z]?)\d+", self.element_label):
            return matches.group(1)
        return self.element_label

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
        res = []
        for i in (0, 2, 4):
            try:
                res.append(int(values[i]))
            except Exception:
                res.append(0)
        try:
            res.append(float(values[1]) * atom_ureg.angstrom)
        except Exception:
            res.append(0 * atom_ureg.angstrom)
        for i in (3, 5):
            try:
                res.append(float(values[i]) * atom_ureg.degree)
            except Exception:
                res.append(0 * atom_ureg.degree)

        pair_count = len(values) // 2

        if res[3].to(atom_ureg.angstrom).m <= 0:
            raise ValueError("Z-matrix bond length must be positive")
        if pair_count >= 2 and not 0 < res[4].to(atom_ureg.degree).m < 180:
            raise ValueError("Z-matrix bond angle must be between 0 and 180 degrees")
        if pair_count >= 3 and format_flag == 1 and not 0 < res[5].to(atom_ureg.degree).m < 180:
            raise ValueError("Alternate Z-matrix second angle must be between 0 and 180 degrees")
        return tuple(res)

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
            atoms=[
                atom_spec.to_atom_in_internal_coordinates_coords()
                for atom_spec in self.atom_specifications
            ]
        )

    def fragment_molecule(self) -> Chem.rdchem.Mol | None:
        try:
            return xyz_to_rdmol(
                self.to_XYZ_block(),
                total_charge=self.total_charge,
                spin_multiplicity=self.spin_multiplicity,
                backend=molopconfig.graph_reconstruction_backend,
                make_dative_bonds=molopconfig.make_dative_bonds,
            )
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

    def coords(self) -> NumpyQuantity:
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

    _zmat_variable_assignment_pattern: ClassVar[re.Pattern[str]] = re.compile(
        r"^(?P<name>[A-Za-z][A-Za-z0-9_]*)\s+(?P<value>[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)$"
    )

    @classmethod
    def _normalize_free_separators(cls, line: str) -> str:
        return line.replace("\t", " ").replace(",", " ").replace("/", " ")

    @classmethod
    def _strip_comment_content(cls, line: str) -> str:
        return line.split("!", 1)[0].rstrip()

    @classmethod
    def _is_zmat_variable_label(cls, line: str) -> bool:
        return line.strip().lower() in {"variables:", "constants:"}

    @classmethod
    def _parse_zmat_variable_assignment(cls, line: str) -> tuple[str, str] | None:
        normalized = " ".join(cls._normalize_free_separators(line).split())
        match = cls._zmat_variable_assignment_pattern.match(normalized)
        if match is None:
            return None
        return match.group("name"), match.group("value")

    @classmethod
    def _split_atom_and_variable_lines(cls, lines: list[str]) -> tuple[list[str], dict[str, str]]:
        atom_lines: list[str] = []
        variable_map: dict[str, str] = {}
        in_variable_block = False
        blank_seen_after_atoms = False

        for line in lines:
            stripped = line.strip()
            if not stripped:
                if atom_lines:
                    blank_seen_after_atoms = True
                continue

            if cls._is_zmat_variable_label(stripped):
                in_variable_block = True
                continue

            assignment = cls._parse_zmat_variable_assignment(stripped)
            if assignment is not None and (in_variable_block or blank_seen_after_atoms):
                var_name, var_value = assignment
                variable_map[var_name] = var_value
                in_variable_block = True
                continue

            atom_lines.append(stripped)

        return atom_lines, variable_map

    @classmethod
    def _replace_zmat_variables_in_line(cls, line: str, variable_map: Mapping[str, str]) -> str:
        normalized_line = cls._normalize_free_separators(line)
        if not variable_map:
            return " ".join(normalized_line.split())
        tokens = normalized_line.split()
        if len(tokens) < 2:
            return " ".join(tokens)

        replaced: list[str] = [tokens[0]]
        for token in tokens[1:]:
            if token in variable_map:
                replaced.append(variable_map[token])
                continue
            if len(token) >= 2 and token[0] in "+-" and token[1:] in variable_map:
                replaced.append(f"{token[0]}{variable_map[token[1:]]}")
                continue
            replaced.append(token)
        return " ".join(replaced)

    @classmethod
    def _validate_internal_coordinate_references(
        cls, atom_specifications: list[GJFAtomSpecification]
    ) -> None:
        for atom_index, atom_spec in enumerate(atom_specifications, start=1):
            if not atom_spec.coords_part.strip() or not atom_spec.is_internal_coords():
                continue

            values = atom_spec.coords_part.split()
            if len(values) % 2 == 1 and len(values) >= 3 and values[-1] in {"0", "1"}:
                values = values[:-1]
            pair_count = len(values) // 2
            refs = [int(values[i]) for i in range(0, len(values), 2)]

            for ref_pos, ref in enumerate(refs, start=1):
                if ref <= 0:
                    raise ValueError(
                        f"Z-matrix reference {ref_pos} on atom {atom_index} must be a positive 1-based index"
                    )
                if ref >= atom_index:
                    raise ValueError(
                        f"Z-matrix reference {ref_pos} on atom {atom_index} must refer to a previously defined atom"
                    )

            if pair_count >= 2 and refs[0] == refs[1]:
                raise ValueError(
                    f"Z-matrix distance and angle references on atom {atom_index} must be different"
                )
            if pair_count >= 3 and len(set(refs[:3])) < 3:
                raise ValueError(
                    f"Z-matrix distance/angle/third references on atom {atom_index} must be distinct"
                )

    @classmethod
    def from_str(cls, data: str) -> "GJFMoleculeSpecifications":
        lines = [cls._strip_comment_content(line) for line in data.splitlines()]
        if "".join(lines).strip() == "":
            return cls()
        first_non_blank_idx = next((idx for idx, line in enumerate(lines) if line.strip()), None)
        if first_non_blank_idx is None:
            return cls()

        body_lines = lines[first_non_blank_idx + 1 :]
        if len(body_lines) < 1:
            raise ValueError(
                "Molecule specifications must have at least 2 lines(charge and spin multiplicity line"
                " & coordinates line)"
            )
        electron_config = cls._normalize_free_separators(lines[first_non_blank_idx]).split()
        if len(electron_config) < 2 or len(electron_config) % 2 != 0:
            raise ValueError("charge and spin multiplicity line must have even number of values")
        atom_lines, variable_map = cls._split_atom_and_variable_lines(body_lines)
        normalized_atom_lines = [
            cls._replace_zmat_variables_in_line(line, variable_map) for line in atom_lines
        ]

        if len(electron_config) == 2:
            total_charge, spin_multiplicity = map(int, electron_config)
            atom_specifications = [
                GJFAtomSpecification.from_str(line) for line in normalized_atom_lines
            ]
            cls._validate_internal_coordinate_references(atom_specifications)
            molecule_fragments = [
                GJFMoleculeSpecificationsFragment(
                    total_charge=total_charge,
                    spin_multiplicity=spin_multiplicity,
                    atom_specifications=atom_specifications,
                )
            ]
        else:
            total_charge, spin_multiplicity = map(int, electron_config[0:2])
            fragment_charge_spin_multiplicity: dict[int, tuple[int, int]] = {}
            for fragment_id, i in enumerate(range(2, len(electron_config), 2)):
                fragment_charge_spin_multiplicity[fragment_id] = (
                    int(electron_config[i]),
                    int(electron_config[i + 1]),
                )
            atom_specifications = [
                GJFAtomSpecification.from_str(line) for line in normalized_atom_lines
            ]
            cls._validate_internal_coordinate_references(atom_specifications)
            fragment_ids = [atom_spec.get_fragment_id() for atom_spec in atom_specifications]
            declared_fragment_count = len(fragment_charge_spin_multiplicity)

            if any(fragment_id <= 0 for fragment_id in fragment_ids):
                raise ValueError(
                    "Multi-fragment molecule specifications require every atom to declare Fragment=n"
                )

            expected_fragment_ids = set(range(1, declared_fragment_count + 1))
            actual_fragment_ids = set(fragment_ids)
            if actual_fragment_ids != expected_fragment_ids:
                raise ValueError(
                    "Fragment assignments must be contiguous and match declared fragment charge/spin pairs"
                )

            molecule_fragments = [
                GJFMoleculeSpecificationsFragment(
                    fragment_id=fragment_id,
                    total_charge=sub_charge,
                    spin_multiplicity=sub_spin_multiplicity,
                    atom_specifications=[
                        atom_specification
                        for atom_specification in atom_specifications
                        if atom_specification.params.get("fragment", "0") == str(fragment_id + 1)
                    ],
                )
                for fragment_id, (
                    sub_charge,
                    sub_spin_multiplicity,
                ) in fragment_charge_spin_multiplicity.items()
            ]
            if any(len(fragment.atom_specifications) == 0 for fragment in molecule_fragments):
                raise ValueError(
                    "Each declared fragment charge/spin pair must correspond to at least one atom"
                )
        return cls(
            total_charge=total_charge,
            spin_multiplicity=spin_multiplicity,
            molecule_fragments=molecule_fragments,
        )

    def _render(self, **kwargs) -> str:
        add_gjf_connectivity = kwargs.get("add_gjf_connectivity", False)

        if len(self.molecule_fragments) == 1:
            spin_multiplicity = f"{self.total_charge} {self.spin_multiplicity}\n"
            molecule_fragment_part = self.molecule_fragments[0]._render(**kwargs) + "\n"
            connectivity_part = f"{self.connectivity()}\n\n" if add_gjf_connectivity else ""
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
                "".join([frag._render(**kwargs) for frag in self.molecule_fragments]) + "\n"
            )
            connectivity_part = f"{self.connectivity()}\n\n" if add_gjf_connectivity else ""
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

    def coords(self) -> NumpyQuantity:
        if len(self.molecule_fragments) == 0:
            return np.array([]) * atom_ureg.angstrom
        return cast(
            NumpyQuantity,
            np.concatenate([frag.coords() for frag in self.molecule_fragments], axis=1),
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
        return [frag.fragment_molecule() for frag in self.molecule_fragments]

    @property
    def rdmol(self) -> Chem.rdchem.Mol:
        return merge_mols_directly(self.rdmol_fragments)

    def connectivity(self) -> str:
        return rdmol_to_gjf_connectivity(self.rdmol)


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

    _gic_standalone_options: ClassVar[set[str]] = {
        "freezeall",
        "unfreezeall",
        "removeall",
    }

    @classmethod
    def _split_top_level_args(cls, raw: str) -> list[str]:
        if not raw.strip():
            return []
        args: list[str] = []
        current: list[str] = []
        depth = 0
        for ch in raw:
            if ch == "," and depth == 0:
                token = "".join(current).strip()
                if token:
                    args.append(token)
                current = []
                continue
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth = max(depth - 1, 0)
            current.append(ch)
        token = "".join(current).strip()
        if token:
            args.append(token)
        return args

    @classmethod
    def _parse_function_expression(cls, expression: str) -> tuple[str | None, list[str]]:
        match = re.fullmatch(r"(?P<name>[A-Za-z][A-Za-z0-9_]*)\((?P<args>.*)\)", expression.strip())
        if match is None:
            return None, []
        return match.group("name"), cls._split_top_level_args(match.group("args"))

    @classmethod
    def _parse_option_token(cls, token: str) -> GJFGICOption:
        raw = token.strip()
        if "=" in raw:
            key, value = raw.split("=", 1)
            return GJFGICOption(
                raw=raw,
                key=key.strip().lower(),
                value=value.strip(),
                is_flag=False,
            )
        return GJFGICOption(raw=raw, key=raw.lower(), value=None, is_flag=True)

    @classmethod
    def _derive_state_from_options(
        cls, options: list[GJFGICOption], standalone_action: str | None = None
    ) -> tuple[str | None, dict[str, str]]:
        option_values = {opt.key: opt.value for opt in options if opt.value is not None}
        normalized_state = None

        if standalone_action is not None:
            action = standalone_action.lower()
            if action in {"freeze", "frozen"}:
                normalized_state = "frozen"
            elif action in {"remove", "inactive", "kill", "removeall"}:
                normalized_state = "inactive"
            elif action in {"active", "activate", "modify"}:
                normalized_state = "active"
            elif action in {"printonly", "print-only"}:
                normalized_state = "print-only"

        if normalized_state is None:
            option_keys = {opt.key for opt in options}
            if option_keys & {"freeze", "frozen"}:
                normalized_state = "frozen"
            elif option_keys & {"inactive", "remove", "kill", "removeall"}:
                normalized_state = "inactive"
            elif option_keys & {"active", "modify"}:
                normalized_state = "active"
            elif option_keys & {"printonly", "print-only"}:
                normalized_state = "print-only"

        return normalized_state, option_values

    @classmethod
    def _split_top_level_assignment(cls, line: str) -> tuple[str, str] | None:
        depth = 0
        for idx, ch in enumerate(line):
            if ch in "([{":
                depth += 1
                continue
            if ch in ")]}":
                depth = max(depth - 1, 0)
                continue
            if ch == "=" and depth == 0:
                return line[:idx].strip(), line[idx + 1 :].strip()
        return None

    @classmethod
    def from_str(cls, data: str) -> "GJFGICLine":
        line = data.strip()
        if not line:
            raise ValueError("GIC line cannot be empty")
        line = line.split("!", 1)[0].strip()
        lower_line = line.lower()
        if lower_line in cls._gic_standalone_options or re.match(r"^atom\s+\S+", lower_line):
            parts = line.split()
            standalone_action = parts[2].lower() if len(parts) > 2 else parts[0].lower()
            normalized_state, option_values = cls._derive_state_from_options([], standalone_action)
            return cls(
                raw_line=data,
                expression=line,
                is_standalone_option=True,
                expression_kind="standalone",
                standalone_keyword=parts[0],
                standalone_args=parts[1:],
                standalone_target=parts[1] if len(parts) > 1 else None,
                standalone_action=standalone_action,
                normalized_state=normalized_state,
                option_values=option_values,
            )

        assignment = cls._split_top_level_assignment(line)
        if assignment is None:
            function_name, function_args = cls._parse_function_expression(line)
            return cls(
                raw_line=data,
                expression=line,
                expression_kind="function" if function_name else "raw",
                function_name=function_name,
                function_args=function_args,
            )

        left, right = assignment
        left_match = re.match(
            r"^(?P<label>[A-Za-z][A-Za-z0-9]*)(?:\s*[\(\[\{](?P<opts>[^\)\]\}]*)[\)\]\}])?$",
            left,
        )
        if left_match is None:
            function_name, function_args = cls._parse_function_expression(right)
            return cls(
                raw_line=data,
                expression=line,
                expression_kind="function" if function_name else "raw",
                function_name=function_name,
                function_args=function_args,
            )

        opts = left_match.group("opts")
        function_name, function_args = cls._parse_function_expression(right)
        parsed_options = (
            [cls._parse_option_token(opt) for opt in cls._split_top_level_args(opts)]
            if opts
            else []
        )
        normalized_state, option_values = cls._derive_state_from_options(parsed_options)
        return cls(
            raw_line=data,
            label=left_match.group("label"),
            label_options=cls._split_top_level_args(opts) if opts else [],
            expression=right,
            expression_kind="function" if function_name else "assignment",
            function_name=function_name,
            function_args=function_args,
            parsed_label_options=parsed_options,
            normalized_state=normalized_state,
            option_values=option_values,
        )

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
    def looks_like(cls, data: str) -> bool:
        lines = [line.strip() for line in data.splitlines() if line.strip()]
        if not lines:
            return False
        gic_head_tokens = (
            "r(",
            "bond(",
            "stretch(",
            "a(",
            "angle(",
            "bend(",
            "d(",
            "dihedral(",
            "torsion(",
            "l(",
            "linear(",
            "linearbend(",
            "x(",
            "y(",
            "z(",
            "cartesian(",
            "cart(",
            "dotdiff(",
            "xcntr(",
            "ycntr(",
            "zcntr(",
            "freezeall",
            "unfreezeall",
            "removeall",
            "atom ",
        )
        for line in lines:
            lowered = line.lower()
            if lowered.startswith(gic_head_tokens) or "=" in lowered:
                return True
        return False

    @classmethod
    def from_str(cls, data: str) -> "GJFGICSection":
        parsed_lines: list[GJFGICLine] = []
        for line in data.splitlines():
            if not line.strip():
                continue
            parsed_lines.append(GJFGICLine.from_str(line))
        return cls(raw=data, lines=parsed_lines)

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
    def from_str(cls, data: str) -> "GJFModRedundantLine":
        line = data.strip()
        if not line:
            raise ValueError("ModRedundant line cannot be empty")
        parts = line.split()
        if len(parts) <= 1:
            return cls(raw_line=data, coordinate_type=parts[0])

        coordinate_type = parts[0]
        atom_refs: list[str] = []
        idx = 1
        while idx < len(parts) and re.match(r"^(\*|-?\d+)$", parts[idx]):
            atom_refs.append(parts[idx])
            idx += 1
        action = parts[idx] if idx < len(parts) else None
        parameters = parts[idx + 1 :] if idx + 1 < len(parts) else []
        return cls(
            raw_line=data,
            coordinate_type=coordinate_type,
            atom_refs=atom_refs,
            action=action,
            parameters=parameters,
        )

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
    def looks_like(cls, data: str) -> bool:
        lines = [line.strip() for line in data.splitlines() if line.strip()]
        if not lines:
            return False
        modredundant_actions = {"a", "f", "b", "k", "r", "d", "h", "s"}
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            if not re.match(r"^[A-Za-z]$", parts[0]):
                continue
            if any(part.lower() in modredundant_actions for part in parts[1:]):
                return True
        return False

    @classmethod
    def from_str(cls, data: str) -> "GJFModRedundantSection":
        parsed_lines: list[GJFModRedundantLine] = []
        for line in data.splitlines():
            if not line.strip():
                continue
            parsed_lines.append(GJFModRedundantLine.from_str(line))
        return cls(raw=data, lines=parsed_lines)

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
    def from_str(cls, data: str) -> "GJFUnknownSection":
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
    def looks_like(cls, data: str) -> bool:
        lines = [line.strip() for line in data.splitlines() if line.strip()]
        if len(lines) < 2:
            return False
        return lines[0].lower().startswith("$nbo") and lines[-1].lower() == "$end"

    @classmethod
    def from_str(cls, data: str) -> "GJFNBOSection":
        lines = [line.rstrip() for line in data.splitlines() if line.strip()]
        if not lines:
            raise ValueError("NBO section cannot be empty")
        header = lines[0]
        footer = lines[-1] if len(lines) > 1 else "$end"
        commands = lines[1:-1] if len(lines) > 2 else []
        return cls(raw=data, header=header, commands=commands, footer=footer)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        parts = [self.header]
        parts.extend(self.commands)
        parts.append(self.footer)
        return "\n".join(parts)


class GJFFileFrameMixin:
    @classmethod
    def _build_molecule_specifications_from_payload(
        cls, payload: dict[str, Any]
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

        atom_lines: list[str] = []
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
            atom_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")

        total_charge = int(payload.get("charge", 0) or 0)
        spin_multiplicity = int(payload.get("multiplicity", 1) or 1)
        return GJFMoleculeSpecifications.from_str(
            f"{total_charge} {spin_multiplicity}\n" + "\n".join(atom_lines)
        )

    @model_validator(mode="before")
    @classmethod
    def normalize_frame_payload(cls, data: Any):
        if not isinstance(data, Mapping):
            return data

        payload = dict(data)

        if isinstance(payload.get("link0_commands"), str):
            payload["link0_commands"] = GJFLink0Commands.from_str(payload["link0_commands"])
        elif payload.get("link0_commands") is None and isinstance(
            payload.get("resources_raw"), str
        ):
            payload["link0_commands"] = GJFLink0Commands.from_str(payload["resources_raw"])

        if isinstance(payload.get("route_section"), str):
            payload["route_section"] = GJFRouteSection.from_str(payload["route_section"])
        elif payload.get("route_section") is None and isinstance(payload.get("keywords"), str):
            payload["route_section"] = GJFRouteSection.from_str(payload["keywords"])

        if isinstance(payload.get("title_card"), str):
            payload["title_card"] = GJFTitleCard.from_str(payload["title_card"])

        if isinstance(payload.get("molecule_specifications"), str):
            payload["molecule_specifications"] = GJFMoleculeSpecifications.from_str(
                payload["molecule_specifications"]
            )
        elif payload.get("molecule_specifications") is None:
            built = cls._build_molecule_specifications_from_payload(payload)
            if built is not None:
                payload["molecule_specifications"] = built

        return payload

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

    @classmethod
    def split_additional_sections(cls, data: str) -> list[str]:
        raw = data.strip("\n")
        if not raw.strip():
            return []
        return [section for section in re.split(r"\n\s*\n", raw) if section.strip()]

    @classmethod
    def parse_additional_sections(
        cls, data: str
    ) -> tuple[
        list[GJFGICSection | GJFModRedundantSection | GJFNBOSection | GJFUnknownSection],
        list[GJFSectionParsingDiagnostic],
    ]:
        parsed_sections: list[
            GJFGICSection | GJFModRedundantSection | GJFNBOSection | GJFUnknownSection
        ] = []
        diagnostics: list[GJFSectionParsingDiagnostic] = []
        for section_index, section_raw in enumerate(cls.split_additional_sections(data)):
            looks_like_gic = GJFGICSection.looks_like(section_raw)
            looks_like_modredundant = GJFModRedundantSection.looks_like(section_raw)
            looks_like_nbo = GJFNBOSection.looks_like(section_raw)

            if looks_like_nbo:
                try:
                    parsed_sections.append(GJFNBOSection.from_str(section_raw))
                except Exception as exc:
                    parsed_sections.append(GJFUnknownSection.from_str(section_raw))
                    diagnostics.append(
                        GJFSectionParsingDiagnostic(
                            section_type="nbo",
                            message=f"Failed to parse as NBO section: {exc}",
                            section_index=section_index,
                        )
                    )
                continue

            if looks_like_gic and looks_like_modredundant:
                parsed_sections.append(GJFUnknownSection.from_str(section_raw))
                diagnostics.append(
                    GJFSectionParsingDiagnostic(
                        section_type="mixed-additional-section",
                        message=(
                            "A single additional section mixes GIC and ModRedundant syntax, "
                            "which Gaussian does not allow"
                        ),
                        section_index=section_index,
                    )
                )
                continue

            if looks_like_gic:
                try:
                    parsed_sections.append(GJFGICSection.from_str(section_raw))
                except Exception as exc:
                    parsed_sections.append(GJFUnknownSection.from_str(section_raw))
                    diagnostics.append(
                        GJFSectionParsingDiagnostic(
                            section_type="gic",
                            message=f"Failed to parse as GIC section: {exc}",
                            section_index=section_index,
                        )
                    )
                continue

            if looks_like_modredundant:
                try:
                    parsed_sections.append(GJFModRedundantSection.from_str(section_raw))
                except Exception as exc:
                    parsed_sections.append(GJFUnknownSection.from_str(section_raw))
                    diagnostics.append(
                        GJFSectionParsingDiagnostic(
                            section_type="modredundant",
                            message=f"Failed to parse as ModRedundant section: {exc}",
                            section_index=section_index,
                        )
                    )
                continue

            parsed_sections.append(GJFUnknownSection.from_str(section_raw))

        return parsed_sections, diagnostics

    def _render(
        self,
        link0_commands: str | GJFLink0Commands | dict[str, str | None] | None = None,
        route_section: str | GJFRouteSection | None = None,
        title_card: str | GJFTitleCard | None = None,
        molecule_specifications: str | GJFMoleculeSpecifications | None = None,
        additional_sections: str | None = None,
        parsed_additional_sections: list[GJFGICSection | GJFModRedundantSection | GJFUnknownSection]
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
        - ``chk`` / ``old_chk`` append keywords to the ``link0_commands_to_use``
          object used for this render. If that object is ``self.link0_commands``,
          the appended entries may therefore be visible on the instance-backed
          object as a side effect.
        """
        if isinstance(link0_commands, str):
            link0_commands_to_use = GJFLink0Commands.from_str(link0_commands)
        elif isinstance(link0_commands, GJFLink0Commands):
            link0_commands_to_use = link0_commands
        elif isinstance(link0_commands, dict):
            link0_commands_to_use = GJFLink0Commands.from_dict(link0_commands)
        elif link0_commands is None:
            link0_commands_to_use = self.link0_commands
        else:
            raise ValueError(f"Invalid type for link0_commands: {type(link0_commands)}")

        if isinstance(route_section, str):
            route_section_to_use = GJFRouteSection.from_str(route_section)
        elif isinstance(route_section, GJFRouteSection):
            route_section_to_use = route_section
        elif route_section is None:
            route_section_to_use = self.route_section
        else:
            raise ValueError(f"Invalid type for route_section: {type(route_section)}")

        if isinstance(title_card, str):
            title_card_to_use = GJFTitleCard.from_str(title_card)
        elif isinstance(title_card, GJFTitleCard):
            title_card_to_use = title_card
        elif title_card is None:
            title_card_to_use = self.title_card
        else:
            raise ValueError(f"Invalid type for title_card: {type(title_card)}")

        if isinstance(molecule_specifications, str):
            molecule_specifications_to_use = GJFMoleculeSpecifications.from_str(
                molecule_specifications
            )
        elif isinstance(molecule_specifications, GJFMoleculeSpecifications):
            molecule_specifications_to_use = molecule_specifications
        elif molecule_specifications is None:
            molecule_specifications_to_use = self.molecule_specifications
        else:
            raise ValueError(
                f"Invalid type for molecule_specifications: {type(molecule_specifications)}"
            )

        if title_card_to_use.title_card == "":
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

        if chk:
            link0_commands_to_use.add_link0_keyword(
                "chk", chk if isinstance(chk, str) else f"{chk_basename}.chk"
            )
        if old_chk:
            link0_commands_to_use.add_link0_keyword(
                "oldchk",
                old_chk if isinstance(old_chk, str) else f"{chk_basename}.chk",
            )
        return (
            link0_commands_to_use._render(**kwargs)
            + route_section_to_use._render(**kwargs)
            + title_card_to_use._render(**kwargs)
            + molecule_specifications_to_use._render(
                coords_type=coords_type, add_gjf_connectivity=add_gjf_connectivity, **kwargs
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
        if semantic_route.model_chemistry.method_family is not None:
            typed_self.method = semantic_route.model_chemistry.method_family
        if semantic_route.model_chemistry.basis_set is not None:
            typed_self.basis_set = semantic_route.model_chemistry.basis_set
        if semantic_route.model_chemistry.functional is not None:
            typed_self.functional = semantic_route.model_chemistry.functional
        typed_self.resources_raw = self.link0_commands._render()
        typed_self.request_num_cpu = self.link0_commands.cpu_request()
        typed_self.request_memory = self.link0_commands.memory_request()
        typed_self.charge = self.molecule_specifications.total_charge
        typed_self.multiplicity = self.molecule_specifications.spin_multiplicity
        typed_self.atoms = self.molecule_specifications.atomic_numbers()
        typed_self.coords = self.molecule_specifications.coords()
        return self

    @model_validator(mode="after")
    def set_additional_section_models(self):
        if self.parsed_additional_sections:
            return self
        parsed_sections, diagnostics = self.parse_additional_sections(self.additional_sections)
        self.parsed_additional_sections = parsed_sections
        self.additional_section_diagnostics = diagnostics
        return self


class GJFFileFrameMemory(
    MemoryStorageMixin, GJFFileFrameMixin, BaseQMInputFrame["GJFFileFrameMemory"]
): ...


class GJFFileFrameDisk(
    DiskStorageMixin, GJFFileFrameMixin, BaseQMInputFrame["GJFFileFrameDisk"]
): ...
