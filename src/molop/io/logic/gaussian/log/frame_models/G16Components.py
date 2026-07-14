from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, ClassVar, Protocol

import numpy as np
from pydantic import Field
from rdkit import Chem

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    Energies,
    GeometryOptimizationStatus,
    MolecularOrbitals,
    Polarizability,
    Status,
    ThermalInformations,
    TotalSpin,
    Vibrations,
)
from molop.unit import atom_ureg


RAW_MODEL_KEY = "__molop_model__"
RAW_DATA_KEY = "data"
RAW_QUANTITY_KEY = "__molop_quantity__"

pt = Chem.GetPeriodicTable()

_RAW_MODEL_REGISTRY: dict[str, type[BaseDataClassWithUnit]] = {}


def _register_raw_model(model_cls: type[BaseDataClassWithUnit]) -> type[BaseDataClassWithUnit]:
    _RAW_MODEL_REGISTRY[model_cls.__name__] = model_cls
    return model_cls


def _format_float(value: Any, digits: int = 6) -> str:
    value = _restore_raw_payload(value)
    if hasattr(value, "m"):
        value = value.m
    if isinstance(value, np.ndarray):
        value = float(value.reshape(-1)[0]) if value.size == 1 else float(value.flat[0])
    return f"{float(value):.{digits}f}"


def _format_job_time_seconds(seconds: float) -> str:
    days = int(seconds // 86400)
    seconds -= days * 86400
    hours = int(seconds // 3600)
    seconds -= hours * 3600
    minutes = int(seconds // 60)
    seconds -= minutes * 60
    return f"Job cpu time: {days} days {hours} hours {minutes} minutes {_format_float(seconds, 2)} seconds."


def _render_orientation_block(title: str, payload: Mapping[str, Any]) -> str:
    restored_payload = _restore_raw_payload(dict(payload))
    coords = restored_payload.get("coords")
    if coords is None:
        coords = restored_payload.get("standard_coords")
    atoms = restored_payload.get("atoms", [])
    if coords is None or not atoms:
        return ""
    coords_array = coords.m if hasattr(coords, "m") else np.asarray(coords)
    lines = [
        f"{title}:",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for idx, (atom, xyz) in enumerate(zip(atoms, coords_array, strict=False), start=1):
        lines.append(
            f" {idx:5d} {int(atom):11d} {0:11d}"
            f" {float(xyz[0]):13.6f} {float(xyz[1]):11.6f} {float(xyz[2]):11.6f}"
        )
    lines.append(" ---------------------------------------------------------------------")
    return "\n".join(lines)


def _render_frequency_payload(payload: Mapping[str, Any]) -> str:
    vibrations = payload.get("vibrations")
    if vibrations is None:
        return ""
    vib_map = _payload_mapping(_restore_raw_payload(vibrations))
    lines: list[str] = [
        " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering",
        " activities (A**4/AMU), depolarization ratios for plane and unpolarized",
        " incident light, reduced masses (AMU), force constants (mDyne/A),",
        " and normal coordinates:",
    ]
    frequencies = vib_map.get("frequencies")
    if frequencies is None or not len(frequencies):
        return ""
    atoms = payload.get("atoms") or []
    frequency_values = np.asarray(
        frequencies.m if hasattr(frequencies, "m") else frequencies
    ).reshape(-1)
    reduced_masses = vib_map.get("reduced_masses")
    reduced_mass_values = (
        np.asarray(reduced_masses.m if hasattr(reduced_masses, "m") else reduced_masses).reshape(-1)
        if reduced_masses is not None and len(reduced_masses)
        else None
    )
    force_constants = vib_map.get("force_constants")
    force_constant_values = (
        np.asarray(force_constants.m if hasattr(force_constants, "m") else force_constants).reshape(
            -1
        )
        if force_constants is not None and len(force_constants)
        else None
    )
    ir_intensities = vib_map.get("IR_intensities")
    ir_intensity_values = (
        np.asarray(ir_intensities.m if hasattr(ir_intensities, "m") else ir_intensities).reshape(-1)
        if ir_intensities is not None and len(ir_intensities)
        else None
    )
    vibration_modes = vib_map.get("vibration_modes")
    mode_arrays: list[np.ndarray[Any, Any]] = []
    if vibration_modes is not None and len(vibration_modes):
        restored_modes = vibration_modes.m if hasattr(vibration_modes, "m") else vibration_modes
        mode_arrays = [
            np.asarray(mode.m if hasattr(mode, "m") else mode) for mode in restored_modes
        ]

    for start in range(0, len(frequency_values), 3):
        stop = min(start + 3, len(frequency_values))
        chunk = frequency_values[start:stop]
        lines.append("")
        lines.append("".join(f"{mode_number:23d}" for mode_number in range(start + 1, stop + 1)))
        lines.append("".join(f"{'A':>23}" for _ in range(len(chunk))))
        lines.append(" Frequencies --" + "".join(f"{float(value):12.4f}" for value in chunk))
        if reduced_mass_values is not None:
            lines.append(
                " Red. masses --"
                + "".join(f"{float(value):12.4f}" for value in reduced_mass_values[start:stop])
            )
        if force_constant_values is not None:
            lines.append(
                " Frc consts  --"
                + "".join(f"{float(value):12.4f}" for value in force_constant_values[start:stop])
            )
        if ir_intensity_values is not None:
            lines.append(
                " IR Inten    --"
                + "".join(f"{float(value):12.4f}" for value in ir_intensity_values[start:stop])
            )
        lines.append("  Atom  AN" + "".join("      X      Y      Z" for _ in range(len(chunk))))
        if atoms and mode_arrays:
            for atom_index, atomic_number in enumerate(atoms, start=1):
                row = f"{atom_index:6d}{int(atomic_number):4d}"
                for mode_idx in range(start, stop):
                    if mode_idx < len(mode_arrays):
                        mode = np.asarray(mode_arrays[mode_idx])
                        vector = (
                            mode[atom_index - 1]
                            if mode.ndim > 1 and atom_index - 1 < len(mode)
                            else [0.0, 0.0, 0.0]
                        )
                        row += "".join(
                            f"{float(value):7.2f}" for value in np.asarray(vector).reshape(-1)[:3]
                        )
                    else:
                        row += f"{0.0:7.2f}{0.0:7.2f}{0.0:7.2f}"
                lines.append(row)
    lines.append(" -------------------")
    return "\n".join(lines)


for _model_cls in (
    Energies,
    TotalSpin,
    MolecularOrbitals,
    ChargeSpinPopulations,
    Polarizability,
    Vibrations,
    ThermalInformations,
    GeometryOptimizationStatus,
    Status,
):
    _register_raw_model(_model_cls)


@dataclass(frozen=True, slots=True)
class G16SyntheticChildSpec:
    component_cls: type[G16BaseComponent]
    payload: Mapping[str, Any]
    node_name: str | None = None
    include_in_aggregation: bool = False
    include_in_render: bool = False


def _payload_mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if hasattr(value, "model_dump"):
        return value.model_dump()
    return {}


if TYPE_CHECKING:
    pass


def _rawify_payload(value: Any) -> Any:
    if isinstance(value, BaseDataClassWithUnit):
        model_fields = type(value).model_fields
        return {
            RAW_MODEL_KEY: value.__class__.__name__,
            RAW_DATA_KEY: {
                field_name: _rawify_payload(getattr(value, field_name))
                for field_name in model_fields
                if getattr(value, field_name) is not None
            },
        }
    if hasattr(value, "m") and hasattr(value, "units"):
        magnitude = value.m.tolist() if isinstance(value.m, np.ndarray) else value.m
        return {
            RAW_QUANTITY_KEY: True,
            "magnitude": _rawify_payload(magnitude),
            "unit": str(value.units),
        }
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return [_rawify_payload(item) for item in value]
    if isinstance(value, tuple):
        return [_rawify_payload(item) for item in value]
    if isinstance(value, Mapping):
        return {key: _rawify_payload(item) for key, item in value.items()}
    return value


def _restore_raw_payload(value: Any) -> Any:
    if isinstance(value, list):
        return [_restore_raw_payload(item) for item in value]
    if isinstance(value, Mapping):
        if value.get(RAW_QUANTITY_KEY):
            magnitude = _restore_raw_payload(value["magnitude"])
            unit = atom_ureg.Unit(value["unit"])
            if isinstance(magnitude, list):
                return np.array(magnitude) * unit
            return magnitude * unit
        if RAW_MODEL_KEY in value:
            model_cls = _RAW_MODEL_REGISTRY[value[RAW_MODEL_KEY]]
            restored_data = {
                key: _restore_raw_payload(item) for key, item in value[RAW_DATA_KEY].items()
            }
            return model_cls.model_validate(restored_data)
        return {key: _restore_raw_payload(item) for key, item in value.items()}
    return value


def _has_meaningful_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    try:
        return len(value) > 0  # type: ignore[arg-type]
    except Exception:
        return True


class G16ComponentTreeProtocol(Protocol):
    only_extract_structure: bool


@dataclass(slots=True)
class G16BaseComponent:
    component_name: ClassVar[str] = "g16.component"
    gaussian_block: ClassVar[str] = ""
    parent_block: ClassVar[str] = "g16.frame"
    repeatable: ClassVar[bool] = True
    allowed_child_component_names: ClassVar[tuple[str, ...]] = ()
    required_child_component_names: ClassVar[tuple[str, ...]] = ()
    repeatable_child_component_names: ClassVar[tuple[str, ...]] = ()
    required_frame_fields: ClassVar[tuple[str, ...]] = ()
    optional_frame_fields: ClassVar[tuple[str, ...]] = ()

    span_start: int = 0
    span_end: int = 0
    raw_text: str = ""
    payload: dict[str, Any] = field(default_factory=dict)

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        return ()

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        return {}

    @classmethod
    def can_build_from_frame(cls, frame: Any) -> bool:
        return all(
            getattr(frame, field_name, None) is not None for field_name in cls.required_frame_fields
        )

    def _render_fakeg(self, **kwargs) -> str:
        return self.raw_text

    def render_fakeg(self, **kwargs) -> str:
        return self._render_fakeg(**kwargs)


class G16L1HeaderComponent(G16BaseComponent):
    component_name = "l1.header"
    gaussian_block = "l1.header"
    allowed_child_component_names = (
        "l1.options",
        "l1.keywords",
        "l101.title",
        "l101.charge_multiplicity",
    )
    repeatable = False
    optional_frame_fields = (
        "qm_software_version",
        "options",
        "keywords",
        "title_card",
        "charge",
        "multiplicity",
    )

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        if self.payload.get("options") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L1OptionsComponent,
                    payload={"options": self.payload["options"]},
                )
            )
        if self.payload.get("keywords") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L1KeywordsComponent,
                    payload={"keywords": self.payload["keywords"]},
                )
            )
        if self.payload.get("title_card") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L101TitleComponent,
                    payload={"title_card": self.payload["title_card"]},
                )
            )
        charge_mult_payload = {
            key: self.payload[key]
            for key in ("charge", "multiplicity")
            if self.payload.get(key) is not None
        }
        if charge_mult_payload:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L101ChargeMultiplicityComponent,
                    payload=charge_mult_payload,
                )
            )
        return tuple(specs)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        for key in (
            "qm_software_version",
            "options",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
            "running_time",
        ):
            if (value := getattr(frame, key, None)) is not None:
                payload[key] = value
        return payload

    def _render_fakeg(self, **kwargs) -> str:
        running_time = _restore_raw_payload(self.payload.get("running_time"))
        if running_time is None:
            return self.raw_text
        seconds = float(
            running_time.to("second").m if hasattr(running_time, "to") else running_time
        )
        return _format_job_time_seconds(seconds)


class G16L101TitleComponent(G16BaseComponent):
    component_name = "l101.title"
    gaussian_block = "l101.title"
    repeatable = False


class G16L1OptionsComponent(G16BaseComponent):
    component_name = "l1.options"
    gaussian_block = "l1.options"
    repeatable = False


class G16L1KeywordsComponent(G16BaseComponent):
    component_name = "l1.keywords"
    gaussian_block = "l1.keywords"
    repeatable = False


class G16L101ChargeMultiplicityComponent(G16BaseComponent):
    component_name = "l101.charge_multiplicity"
    gaussian_block = "l101.charge_multiplicity"
    repeatable = False

    def _render_fakeg(self, **kwargs) -> str:
        charge = self.payload.get("charge")
        multiplicity = self.payload.get("multiplicity")
        if charge is None or multiplicity is None:
            return self.raw_text
        return f" Charge = {int(charge):4d} Multiplicity = {int(multiplicity)}"


class G16L202OrientComponent(G16BaseComponent):
    component_name = "l202.orient"
    gaussian_block = "l202.orient"
    allowed_child_component_names = (
        "l202.orient.input",
        "l202.orient.standard",
        "l202.distmat",
        "l202.stoich",
    )
    required_frame_fields = ("atoms",)
    optional_frame_fields = ("coords", "standard_coords")

    def _render_fakeg(self, **kwargs) -> str:
        if "standard_coords" in self.payload:
            return _render_orientation_block("Standard orientation", self.payload)
        if "coords" in self.payload:
            return _render_orientation_block("Input orientation", self.payload)
        return self.raw_text

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        atoms = getattr(frame, "atoms", None)
        coords = getattr(frame, "coords", None)
        standard_coords = getattr(frame, "standard_coords", None)
        if atoms and coords is not None:
            payload["atoms"] = atoms
            payload["coords"] = coords
        if atoms and standard_coords is not None:
            payload["atoms"] = atoms
            payload["standard_coords"] = standard_coords
        return payload


class G16L202RotConstComponent(G16BaseComponent):
    component_name = "l202.rotconst"
    gaussian_block = "l202.rotconst"
    required_frame_fields = ("rotation_constants",)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (rotation_constants := getattr(frame, "rotation_constants", None)) is not None:
            return {"rotation_constants": rotation_constants}
        return {}


class G16L502CycleComponent(G16BaseComponent):
    component_name = "l502.cycle"
    gaussian_block = "l502.cycle"
    required_frame_fields = ("energies",)
    optional_frame_fields = ("total_spin",)

    def _render_fakeg(self, **kwargs) -> str:
        lines: list[str] = []
        energies = _restore_raw_payload(self.payload.get("energies"))
        if energies is not None and getattr(energies, "reference_energy", None) is not None:
            lines.append(
                f" SCF Done:  E(SCF) =  {_format_float(energies.reference_energy, 9)} A.U. after   1 cycles"
            )
        total_spin = _restore_raw_payload(self.payload.get("total_spin"))
        if total_spin is not None and getattr(total_spin, "spin_square", None) is not None:
            line = f" S**2 = {_format_float(total_spin.spin_square, 6)}"
            if getattr(total_spin, "spin_quantum_number", None) is not None:
                line += f"                 S = {_format_float(total_spin.spin_quantum_number, 6)}"
            lines.append(line)
        return "\n".join(lines) if lines else self.raw_text

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (energies := getattr(frame, "energies", None)) is not None:
            payload["energies"] = energies
        if (total_spin := getattr(frame, "total_spin", None)) is not None:
            payload["total_spin"] = total_spin
        return payload


class G16L601PopAnalComponent(G16BaseComponent):
    component_name = "l601.popanal"
    gaussian_block = "l601.popanal"
    allowed_child_component_names = (
        "l601.state",
        "l601.molecular_orbitals",
        "l601.charge_spin_populations",
        "l601.polarizability",
    )
    optional_frame_fields = (
        "molecular_orbitals",
        "charge_spin_populations",
        "polarizability",
    )

    _MO_SYMMETRY_MARKERS = ("Orbital symmetries:",)
    _ELECTRONIC_STATE_MARKERS = ("The electronic state is",)
    _MO_ENERGY_MARKERS = (
        "Alpha  occ. eigenvalues --",
        "Beta  occ. eigenvalues --",
        "Alpha virt. eigenvalues --",
        "Beta virt. eigenvalues --",
    )
    _POPULATION_SECTION_MARKERS = {
        "mulliken_spins": ("Mulliken charges and spin densities:",),
        "mulliken_charges": (
            "Mulliken charges:",
            "Mulliken atomic charges",
            "Mulliken charges and spin densities:",
        ),
        "apt_charges": ("APT charges:",),
        "lowdin_charges": ("Lowdin charges",),
    }
    _POPULATION_SECTION_END_MARKERS = {
        "apt_charges": ("Sum of APT charges",),
        "lowdin_charges": ("Sum of Lowdin charges",),
    }
    _POLAR_SECTION_MARKERS = {
        "electronic_spatial_extent": ("Electronic spatial extent (au):",),
        "dipole": ("Dipole moment (field-independent basis",),
        "quadrupole": ("Quadrupole moment (field-independent basis",),
        "traceless_quadrupole": ("Traceless Quadrupole moment (field-independent basis",),
        "octapole": ("Octapole moment (field-independent basis",),
        "hexadecapole": ("Hexadecapole moment (field-independent basis",),
    }
    _REMAINING_BLOCK_MARKERS = {
        "exact_polarizability": ("Exact polarizability:",),
        "hirshfeld": ("Hirshfeld charges, spin densities, dipoles, and CM5 charges",),
        "dipole_before_force": ("Dipole        =",),
        "polarizability_before_force": ("Polarizability=",),
    }

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (molecular_orbitals := getattr(frame, "molecular_orbitals", None)) is not None:
            payload["molecular_orbitals"] = molecular_orbitals
        if (charge_spin_populations := getattr(frame, "charge_spin_populations", None)) is not None:
            payload["charge_spin_populations"] = charge_spin_populations
        if (polarizability := getattr(frame, "polarizability", None)) is not None:
            payload["polarizability"] = polarizability
        return payload


class G16L716FreqComponent(G16BaseComponent):
    component_name = "l716.freq"
    gaussian_block = "l716.freq"
    allowed_child_component_names = (
        "l716.forceconstants",
        "l716.diagvib",
        "l716.irspectrum",
        "l716.vibration.mode",
    )

    repeatable_child_component_names = ("l716.vibration.mode",)
    required_frame_fields = ("vibrations",)

    def _render_fakeg(self, **kwargs) -> str:
        rendered = _render_frequency_payload(self.payload)
        return rendered or self.raw_text

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        vibrations = _restore_raw_payload(self.payload.get("vibrations"))
        if vibrations is None:
            return ()
        specs: list[G16SyntheticChildSpec] = []
        force_constants = getattr(vibrations, "force_constants", None)
        if force_constants is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ForceConstantsComponent,
                    payload={"force_constants": force_constants},
                )
            )
        vibration_modes = getattr(vibrations, "vibration_modes", None)
        if vibration_modes is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716DiagVibComponent,
                    payload={"vibration_modes": vibration_modes},
                )
            )
        ir_intensities = getattr(vibrations, "IR_intensities", None)
        if ir_intensities is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716IRSpectrumComponent,
                    payload={"IR_intensities": ir_intensities},
                )
            )
        for idx, vibration in enumerate(vibrations):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716VibrationModeComponent,
                    payload={
                        "mode_index": idx,
                        "is_imaginary": getattr(vibration, "is_imaginary", False),
                        **{
                            attr: value
                            for attr in (
                                "frequency",
                                "reduced_mass",
                                "force_constant",
                                "IR_intensity",
                                "vibration_mode",
                            )
                            if (value := getattr(vibration, attr, None)) is not None
                        },
                    },
                    node_name=f"l716.vibration.mode[{idx}]",
                )
            )
        return tuple(specs)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (vibrations := getattr(frame, "vibrations", None)) is not None:
            payload["vibrations"] = vibrations
        if (atoms := getattr(frame, "atoms", None)) is not None:
            payload["atoms"] = atoms
        return payload


class G16L716FreqChunkComponent(G16BaseComponent):
    component_name = "l716.freq.chunkx"
    gaussian_block = "l716.freq.chunkx"


class G16L716ForceConstantsComponent(G16BaseComponent):
    component_name = "l716.forceconstants"
    gaussian_block = "l716.forceconstants"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("force_constants"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Frc consts  -- {joined}"


class G16L716DiagVibComponent(G16BaseComponent):
    component_name = "l716.diagvib"
    gaussian_block = "l716.diagvib"

    def _render_fakeg(self, **kwargs) -> str:
        modes = _restore_raw_payload(self.payload.get("vibration_modes"))
        if modes is None or len(modes) == 0:
            return self.raw_text
        lines = [" Atom  AN      X      Y      Z"]
        for idx, mode in enumerate(modes[:3], start=1):
            arr = mode.m if hasattr(mode, "m") else mode
            vec = np.asarray(arr)
            first = vec[0] if vec.ndim > 1 else vec[:3]
            first = np.asarray(first).reshape(-1)
            padded = list(first[:3]) + [0.0] * max(0, 3 - len(first[:3]))
            lines.append(
                f" {idx:4d}   0 {float(padded[0]):7.2f} {float(padded[1]):7.2f} {float(padded[2]):7.2f}"
            )
        return "\n".join(lines)


class G16L716IRSpectrumComponent(G16BaseComponent):
    component_name = "l716.irspectrum"
    gaussian_block = "l716.irspectrum"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("IR_intensities"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" IR Inten    -- {joined}"


class G16L716VibrationModeComponent(G16BaseComponent):
    component_name = "l716.vibration.mode"
    gaussian_block = "l716.vibration.mode"

    def _render_fakeg(self, **kwargs) -> str:
        lines: list[str] = []
        if (mode_index := self.payload.get("mode_index")) is not None:
            lines.append(f"{int(mode_index) + 1:23d}")
            lines.append(f"{'A':>23}")
        if (value := _restore_raw_payload(self.payload.get("frequency"))) is not None:
            lines.append(f" Frequencies -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("reduced_mass"))) is not None:
            lines.append(f" Red. masses -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("force_constant"))) is not None:
            lines.append(f" Frc consts  -- {_format_float(value, 4)}")
        if (value := _restore_raw_payload(self.payload.get("IR_intensity"))) is not None:
            lines.append(f" IR Inten    -- {_format_float(value, 4)}")
        mode = _restore_raw_payload(self.payload.get("vibration_mode"))
        if mode is not None:
            arr = mode.m if hasattr(mode, "m") else mode
            vec = np.asarray(arr)
            first = vec[0] if vec.ndim > 1 else vec[:3]
            first = np.asarray(first).reshape(-1)
            padded = list(first[:3]) + [0.0] * max(0, 3 - len(first[:3]))
            lines.append(" Atom  AN      X      Y      Z")
            lines.append(
                f" {1:4d}   0 {float(padded[0]):7.2f} {float(padded[1]):7.2f} {float(padded[2]):7.2f}"
            )
        return "\n".join(lines) if lines else self.raw_text


class G16L716ThermochemistryComponent(G16BaseComponent):
    component_name = "l716.thermochemistry"
    gaussian_block = "l716.thermochemistry"
    allowed_child_component_names = (
        "l716.thermochemistry.temperature",
        "l716.thermochemistry.mass",
        "l716.thermochemistry.moi",
        "l716.thermochemistry.rotsymnum",
        "l716.thermochemistry.rottemp",
        "l716.thermochemistry.rotconsts",
        "l716.thermochemistry.vibtemp",
        "l716.thermochemistry.zpe",
        "l716.thermoprops",
        "l716.thermochemistry.energy",
        "l716.thermochemistry.enthalpy",
        "l716.thermochemistry.gibbs",
        "l716.thermochemistry.entropy",
        "l716.thermochemistry.heatcapacity",
    )

    required_frame_fields = ("thermal_informations",)
    optional_frame_fields = ("temperature", "pressure")

    def _render_fakeg(self, **kwargs) -> str:
        thermal = _restore_raw_payload(self.payload.get("thermal_informations"))
        temperature = _restore_raw_payload(self.payload.get("temperature"))
        pressure = _restore_raw_payload(self.payload.get("pressure"))
        if thermal is None and temperature is None and pressure is None:
            return self.raw_text
        lines: list[str] = [" - Thermochemistry -"]
        if temperature is not None or pressure is not None:
            temperature_text = (
                f" Temperature   {_format_float(temperature, 3)} Kelvin."
                if temperature is not None
                else ""
            )
            pressure_text = (
                f"  Pressure   {_format_float(pressure, 5)} Atm." if pressure is not None else ""
            )
            lines.append(f"{temperature_text}{pressure_text}".rstrip())
        if thermal is None:
            return "\n".join(lines) if lines else self.raw_text
        for key, label in (
            ("ZPVE", "Zero-point correction="),
            ("TCE", "Thermal correction to Energy="),
            ("TCH", "Thermal correction to Enthalpy="),
            ("TCG", "Thermal correction to Gibbs Free Energy="),
        ):
            value = getattr(thermal, key, None)
            if value is not None:
                lines.append(f"{label} {_format_float(value, 6)}")
        for key, label in (
            ("U_0", "Sum of electronic and zero-point Energies="),
            ("U_T", "Sum of electronic and thermal Energies="),
            ("H_T", "Sum of electronic and thermal Enthalpies="),
            ("G_T", "Sum of electronic and thermal Free Energies="),
        ):
            value = getattr(thermal, key, None)
            if value is not None:
                lines.append(f"{label} {_format_float(value, 6)}")
        entropy = getattr(thermal, "S", None)
        heat_capacity = getattr(thermal, "C_V", None)
        if entropy is not None or heat_capacity is not None:
            thermal_energy = getattr(thermal, "TCE", None)
            lines.append(" E (Thermal)             CV                S")
            lines.append(
                " Total"
                f"{_format_float(thermal_energy, 3) if thermal_energy is not None else '0.000':>13}"
                f"{_format_float(heat_capacity, 3) if heat_capacity is not None else '0.000':>18}"
                f"{_format_float(entropy, 3) if entropy is not None else '0.000':>17}"
            )
        lines.append(" Rotational           0.000000D+00    0.0000    0.0000")
        return "\n".join(lines) if lines else self.raw_text

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        temperature = _restore_raw_payload(self.payload.get("temperature"))
        pressure = _restore_raw_payload(self.payload.get("pressure"))
        if temperature is not None or pressure is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryTemperatureComponent,
                    payload={
                        key: value
                        for key, value in (
                            ("temperature", temperature),
                            ("pressure", pressure),
                        )
                        if value is not None
                    },
                )
            )
        thermal_info = _restore_raw_payload(self.payload.get("thermal_informations"))
        if thermal_info is not None:
            if getattr(thermal_info, "ZPVE", None) is not None:
                specs.append(
                    G16SyntheticChildSpec(
                        component_cls=G16L716ThermochemistryZPEComponent,
                        payload={"ZPVE": thermal_info.ZPVE},
                    )
                )
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermoPropsComponent,
                    payload={
                        attr: value
                        for attr in (
                            "U_0",
                            "U_T",
                            "TCE",
                            "H_T",
                            "TCH",
                            "G_T",
                            "TCG",
                            "S",
                            "C_V",
                        )
                        if (value := getattr(thermal_info, attr, None)) is not None
                    },
                )
            )
            for component_cls, attrs in (
                (G16L716ThermochemistryEnergyComponent, ("U_0", "U_T", "TCE")),
                (G16L716ThermochemistryEnthalpyComponent, ("H_T", "TCH")),
                (G16L716ThermochemistryGibbsComponent, ("G_T", "TCG")),
                (G16L716ThermochemistryEntropyComponent, ("S",)),
                (G16L716ThermochemistryHeatCapacityComponent, ("C_V",)),
            ):
                payload = {
                    attr: value
                    for attr in attrs
                    if (value := getattr(thermal_info, attr, None)) is not None
                }
                if payload:
                    specs.append(
                        G16SyntheticChildSpec(component_cls=component_cls, payload=payload)
                    )
        if thermal_info is not None and getattr(thermal_info, "molecular_mass", None) is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryMassComponent,
                    payload={"molecular_mass": thermal_info.molecular_mass},
                )
            )
        if (
            thermal_info is not None
            and getattr(thermal_info, "moments_of_inertia", None) is not None
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryMOIComponent,
                    payload={"moments_of_inertia": thermal_info.moments_of_inertia},
                )
            )
        if (
            thermal_info is not None
            and getattr(thermal_info, "rotational_symmetry_number", None) is not None
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryRotSymNumComponent,
                    payload={"rotational_symmetry_number": thermal_info.rotational_symmetry_number},
                )
            )
        if (
            thermal_info is not None
            and getattr(thermal_info, "rotational_temperatures", None) is not None
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryRotTempComponent,
                    payload={"rotational_temperatures": thermal_info.rotational_temperatures},
                )
            )
        if (
            thermal_info is not None
            and getattr(thermal_info, "rotational_constants", None) is not None
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryRotConstsComponent,
                    payload={"rotational_constants": thermal_info.rotational_constants},
                )
            )
        if (
            thermal_info is not None
            and getattr(thermal_info, "vibrational_temperatures", None) is not None
        ):
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716ThermochemistryVibTempComponent,
                    payload={
                        "vibrational_temperatures": thermal_info.vibrational_temperatures,
                        "vibrational_temperature_mode_indices": getattr(
                            thermal_info,
                            "vibrational_temperature_mode_indices",
                            None,
                        ),
                    },
                )
            )
        return tuple(specs)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if (thermal := getattr(frame, "thermal_informations", None)) is not None:
            payload["thermal_informations"] = thermal
        if (temperature := getattr(frame, "temperature", None)) is not None:
            payload["temperature"] = temperature
        if (pressure := getattr(frame, "pressure", None)) is not None:
            payload["pressure"] = pressure
        return payload


class G16L716ThermochemistryTemperatureComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.temperature"
    gaussian_block = "l716.thermochemistry.temperature"

    def _render_fakeg(self, **kwargs) -> str:
        temperature = _restore_raw_payload(self.payload.get("temperature"))
        pressure = _restore_raw_payload(self.payload.get("pressure"))
        if temperature is None and pressure is None:
            return self.raw_text
        if temperature is not None and pressure is not None:
            return (
                f" Temperature   {_format_float(temperature, 3)} Kelvin.  "
                f"Pressure   {_format_float(pressure, 5)} Atm."
            )
        if temperature is not None:
            return f" Temperature   {_format_float(temperature, 3)} Kelvin."
        return f" Pressure   {_format_float(pressure, 5)} Atm."


class G16L716ThermochemistryMassComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.mass"
    gaussian_block = "l716.thermochemistry.mass"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("molecular_mass"))
        if value is not None:
            return f" Molecular mass: {_format_float(value, 4)} amu."
        return self.raw_text


class G16L716ThermochemistryMOIComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.moi"
    gaussian_block = "l716.thermochemistry.moi"

    def _render_fakeg(self, **kwargs) -> str:
        values = _restore_raw_payload(self.payload.get("moments_of_inertia"))
        if values is None or not len(values):
            return self.raw_text
        arr = values.m if hasattr(values, "m") else values
        joined = " ".join(f"{float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Eigenvalues -- {joined}"


class G16L716ThermochemistryRotSymNumComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.rotsymnum"
    gaussian_block = "l716.thermochemistry.rotsymnum"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("rotational_symmetry_number"))
        if value is not None:
            return f" Rotational symmetry number {int(value)}."
        return self.raw_text


class G16L716ThermochemistryRotTempComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.rottemp"
    gaussian_block = "l716.thermochemistry.rottemp"

    def _render_fakeg(self, **kwargs) -> str:
        rotational_temperatures = _restore_raw_payload(self.payload.get("rotational_temperatures"))
        if rotational_temperatures is None:
            return self.raw_text
        arr = (
            rotational_temperatures.m
            if hasattr(rotational_temperatures, "m")
            else rotational_temperatures
        )
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Rotational temperatures (Kelvin){joined}"


class G16L716ThermochemistryRotConstsComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.rotconsts"
    gaussian_block = "l716.thermochemistry.rotconsts"

    def _render_fakeg(self, **kwargs) -> str:
        rotational_constants = _restore_raw_payload(self.payload.get("rotational_constants"))
        if rotational_constants is None:
            return self.raw_text
        arr = rotational_constants.m if hasattr(rotational_constants, "m") else rotational_constants
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Rotational constants (GHZ):{joined}"


class G16L716ThermochemistryVibTempComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.vibtemp"
    gaussian_block = "l716.thermochemistry.vibtemp"

    def _render_fakeg(self, **kwargs) -> str:
        vibrational_temperatures = _restore_raw_payload(
            self.payload.get("vibrational_temperatures")
        )
        if vibrational_temperatures is None:
            return self.raw_text
        arr = (
            vibrational_temperatures.m
            if hasattr(vibrational_temperatures, "m")
            else vibrational_temperatures
        )
        joined = "".join(f" {float(v):10.4f}" for v in np.asarray(arr).reshape(-1))
        return f" Vibrational temperatures:{joined}"


class G16L716ThermochemistryZPEComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.zpe"
    gaussian_block = "l716.thermochemistry.zpe"

    def _render_fakeg(self, **kwargs) -> str:
        zpe = _restore_raw_payload(self.payload.get("ZPVE"))
        if zpe is not None:
            return f"Zero-point correction= {_format_float(zpe, 6)} (Hartree/Particle)"
        return self.raw_text


class G16L716ThermoPropsComponent(G16BaseComponent):
    component_name = "l716.thermoprops"
    gaussian_block = "l716.thermoprops"


class G16L716ThermochemistryEnergyComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.energy"
    gaussian_block = "l716.thermochemistry.energy"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCE"))
        if value is not None:
            lines.append(f"Thermal correction to Energy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("U_0"))
        if value is not None:
            lines.append(f"Sum of electronic and zero-point Energies= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("U_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Energies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16L716ThermochemistryEnthalpyComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.enthalpy"
    gaussian_block = "l716.thermochemistry.enthalpy"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCH"))
        if value is not None:
            lines.append(f"Thermal correction to Enthalpy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("H_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Enthalpies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16L716ThermochemistryGibbsComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.gibbs"
    gaussian_block = "l716.thermochemistry.gibbs"

    def _render_fakeg(self, **kwargs) -> str:
        lines = []
        value = _restore_raw_payload(self.payload.get("TCG"))
        if value is not None:
            lines.append(f"Thermal correction to Gibbs Free Energy= {_format_float(value, 6)}")
        value = _restore_raw_payload(self.payload.get("G_T"))
        if value is not None:
            lines.append(f"Sum of electronic and thermal Free Energies= {_format_float(value, 6)}")
        return "\n".join(lines) if lines else self.raw_text


class G16L716ThermochemistryEntropyComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.entropy"
    gaussian_block = "l716.thermochemistry.entropy"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("S"))
        if value is not None:
            return f"Total Entropy= {_format_float(value, 4)} cal/mol-K"
        return self.raw_text


class G16L716ThermochemistryHeatCapacityComponent(G16BaseComponent):
    component_name = "l716.thermochemistry.heatcapacity"
    gaussian_block = "l716.thermochemistry.heatcapacity"

    def _render_fakeg(self, **kwargs) -> str:
        value = _restore_raw_payload(self.payload.get("C_V"))
        if value is not None:
            return f"Total Heat Capacity= {_format_float(value, 4)} cal/mol-K"
        return self.raw_text


class G16L716PolarizabilityComponent(G16BaseComponent):
    component_name = "l716.polarizability"
    gaussian_block = "l716.polarizability"
    allowed_child_component_names = (
        "l716.dipole",
        "l716.polarizability.detail",
    )
    required_frame_fields = ("polarizability",)

    def build_synthetic_children(self) -> tuple[G16SyntheticChildSpec, ...]:
        specs: list[G16SyntheticChildSpec] = []
        if self.payload.get("dipole") is not None:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716DipoleComponent,
                    payload={"dipole": self.payload.get("dipole")},
                )
            )
        detail_payload = {
            key: value
            for key, value in self.payload.items()
            if key != "dipole" and value is not None
        }
        if detail_payload:
            specs.append(
                G16SyntheticChildSpec(
                    component_cls=G16L716PolarizabilityDetailComponent,
                    payload=detail_payload,
                )
            )
        return tuple(specs)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (polarizability := getattr(frame, "polarizability", None)) is not None:
            return dict(_payload_mapping(polarizability))
        return {}


class G16L716DipoleComponent(G16BaseComponent):
    component_name = "l716.dipole"
    gaussian_block = "l716.dipole"


class G16L716PolarizabilityDetailComponent(G16BaseComponent):
    component_name = "l716.polarizability.detail"
    gaussian_block = "l716.polarizability.detail"


class G16L716ForcesComponent(G16BaseComponent):
    component_name = "l716.forces"
    gaussian_block = "l716.forces"
    required_frame_fields = ("forces",)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (forces := getattr(frame, "forces", None)) is not None:
            return {"forces": forces}
        return {}


class G16L716SecondDerivComponent(G16BaseComponent):
    component_name = "l716.secondderiv"
    gaussian_block = "l716.secondderiv"
    required_frame_fields = ("hessian",)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (hessian := getattr(frame, "hessian", None)) is not None:
            return {"hessian": hessian}
        return {}


class G16L103OptimizationComponent(G16BaseComponent):
    component_name = "l103.optimization"
    gaussian_block = "l103.optimization"
    required_frame_fields = ("geometry_optimization_status",)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (opt := getattr(frame, "geometry_optimization_status", None)) is not None:
            return {"geometry_optimization_status": opt}
        return {}


class G16L9999ArchiveComponent(G16BaseComponent):
    component_name = "l9999.archive"
    gaussian_block = "l9999.archive"
    allowed_child_component_names = (
        "l9999.archive.metadata",
        "l9999.archive.energies",
        "l9999.archive.thermochemistry",
        "l9999.archive.polarizability",
        "l9999.archive.hessian",
    )
    required_frame_fields = ("job_type",)
    optional_frame_fields = (
        "functional",
        "basis_set",
        "keywords",
        "title_card",
        "charge",
        "multiplicity",
        "qm_software_version",
        "atoms",
        "coords",
        "energies",
        "thermal_informations",
        "polarizability",
        "hessian",
    )

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        for key in (
            "job_type",
            "functional",
            "basis_set",
            "keywords",
            "title_card",
            "charge",
            "multiplicity",
            "qm_software_version",
            "atoms",
            "coords",
            "energies",
            "thermal_informations",
            "polarizability",
            "hessian",
        ):
            if (value := getattr(frame, key, None)) is not None:
                payload[key] = value
        return payload


class G16L9999FinalComponent(G16BaseComponent):
    component_name = "l9999.final"
    gaussian_block = "l9999.final"
    allowed_child_component_names = ()
    required_frame_fields = ("status",)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (status := getattr(frame, "status", None)) is not None:
            return {"status": status}
        return {}


class G16JobCPUComponent(G16BaseComponent):
    component_name = "jobcpu"
    gaussian_block = "jobcpu"
    allowed_child_component_names = ()
    required_frame_fields = ("running_time",)

    def _render_fakeg(self, **kwargs) -> str:
        running_time = _restore_raw_payload(
            self.payload.get("job_cpu_time") or self.payload.get("running_time")
        )
        if running_time is None:
            return self.raw_text
        seconds = float(
            running_time.to("second").m if hasattr(running_time, "to") else running_time
        )
        return _format_job_time_seconds(seconds)

    @classmethod
    def build_payload_from_frame(cls, frame: Any) -> dict[str, Any]:
        if (running_time := getattr(frame, "running_time", None)) is not None:
            return {"running_time": running_time, "job_cpu_time": running_time}
        return {}


class G16LinkEnterComponent(G16BaseComponent):
    component_name = "link.enter"
    gaussian_block = "link.enter"
    allowed_child_component_names = ()


class G16LeaveComponent(G16BaseComponent):
    component_name = "link.leave"
    gaussian_block = "link.leave"
    allowed_child_component_names = ()


G16_COMPONENT_REGISTRY: tuple[type[G16BaseComponent], ...] = (
    G16L1HeaderComponent,
    G16L101TitleComponent,
    G16L1KeywordsComponent,
    G16L202OrientComponent,
    G16L202RotConstComponent,
    G16L502CycleComponent,
    G16L601PopAnalComponent,
    G16L716FreqComponent,
    G16L716ThermochemistryComponent,
    G16L716PolarizabilityComponent,
    G16L716ForcesComponent,
    G16L716SecondDerivComponent,
    G16L103OptimizationComponent,
    G16L9999ArchiveComponent,
    G16L9999FinalComponent,
    G16JobCPUComponent,
    G16LinkEnterComponent,
    G16LeaveComponent,
)


def get_g16log_component_classes() -> tuple[type[G16BaseComponent], ...]:
    return G16_COMPONENT_REGISTRY


class G16ComponentNode(BaseDataClassWithUnit):
    node_name: str = Field(default="g16.frame")
    component: G16BaseComponent | None = Field(default=None)
    children: list[G16ComponentNode] = Field(default_factory=list)
    include_in_aggregation: bool = Field(default=True, exclude=True, repr=False)
    include_in_render: bool = Field(default=True, exclude=True, repr=False)

    def iter_nodes(self) -> Iterator[G16ComponentNode]:
        yield self
        for child in self.children:
            yield from child.iter_nodes()

    def iter_components(self, *, include_synthetic: bool = False) -> Iterator[G16BaseComponent]:
        if self.component is not None and (self.include_in_aggregation or include_synthetic):
            yield self.component
        for child in self.children:
            yield from child.iter_components(include_synthetic=include_synthetic)

    def render_fakeg(self, **kwargs: Any) -> str:
        if self.component is not None:
            if not self.include_in_render:
                return ""
            return self.component.render_fakeg(**kwargs)
        return "\n".join(
            rendered
            for rendered in (child.render_fakeg(**kwargs).strip("\n") for child in self.children)
            if rendered
        )

    def render_subtree(self, **kwargs: Any) -> str:
        if self.component is not None:
            own_render = self.component.render_fakeg(**kwargs).strip("\n")
            child_renders = [
                rendered
                for rendered in (
                    child.render_subtree(**kwargs).strip("\n") for child in self.children
                )
                if rendered
            ]
            parts = [part for part in [own_render, *child_renders] if part]
            return "\n".join(parts)
        return "\n".join(
            rendered
            for rendered in (child.render_subtree(**kwargs).strip("\n") for child in self.children)
            if rendered
        )


class G16ComponentTree(BaseDataClassWithUnit):
    root: G16ComponentNode = Field(default_factory=G16ComponentNode)
    source_text: str = Field(default="", exclude=True, repr=False)
    only_extract_structure: bool = Field(default=False, exclude=True, repr=False)
    synthetic_children_expanded: bool = Field(default=False, exclude=True, repr=False)
    payloads_rawified: bool = Field(default=False, exclude=True, repr=False)

    def iter_nodes(self) -> Iterator[G16ComponentNode]:
        return self.root.iter_nodes()

    def iter_components(self, *, include_synthetic: bool = False) -> Iterator[G16BaseComponent]:
        return self.root.iter_components(include_synthetic=include_synthetic)

    def component_names(self, *, include_synthetic: bool = False) -> list[str]:
        return [
            component.component_name
            for component in self.iter_components(include_synthetic=include_synthetic)
        ]

    def find_nodes(self, node_name: str) -> list[G16ComponentNode]:
        return [node for node in self.iter_nodes() if node.node_name == node_name]

    def find_nodes_by_prefix(self, node_name_prefix: str) -> list[G16ComponentNode]:
        return [node for node in self.iter_nodes() if node.node_name.startswith(node_name_prefix)]

    def render_node(self, node_name: str, **kwargs: Any) -> str:
        nodes = self.find_nodes(node_name)
        return "\n".join(
            rendered
            for rendered in (node.render_subtree(**kwargs).strip("\n") for node in nodes)
            if rendered
        )

    def render_fakeg(self, **kwargs: Any) -> str:
        return self.root.render_fakeg(**kwargs)

    def validate_contracts(self) -> list[str]:
        issues: list[str] = []
        for node in self.iter_nodes():
            component = node.component
            if component is None:
                continue
            child_names = [
                child.component.component_name
                for child in node.children
                if child.component is not None
            ]
            allowed = set(component.allowed_child_component_names)
            required = set(component.required_child_component_names)
            repeatable = set(component.repeatable_child_component_names)
            if allowed:
                for child_name in child_names:
                    if child_name not in allowed:
                        issues.append(
                            f"{component.component_name}: child {child_name} not allowed by contract"
                        )
            for required_name in required:
                if required_name not in child_names:
                    issues.append(
                        f"{component.component_name}: required child {required_name} missing"
                    )
            seen_counts: dict[str, int] = {}
            for child_name in child_names:
                seen_counts[child_name] = seen_counts.get(child_name, 0) + 1
            for child_name, count in seen_counts.items():
                if count > 1 and child_name not in repeatable:
                    issues.append(
                        f"{component.component_name}: child {child_name} repeated {count} times but not repeatable"
                    )
        return issues


class G16ComponentTreeBuilder:
    @staticmethod
    def _append_synthetic_child(
        parent_node: G16ComponentNode,
        child_spec: G16SyntheticChildSpec,
    ) -> G16ComponentNode | None:
        normalized_payload = {
            key: value for key, value in child_spec.payload.items() if _has_meaningful_value(value)
        }
        if not normalized_payload:
            return None
        parent_component = parent_node.component
        child_component = child_spec.component_cls(
            span_start=parent_component.span_start if parent_component is not None else 0,
            span_end=parent_component.span_end if parent_component is not None else 0,
            raw_text="",
            payload=dict(normalized_payload),
        )
        child_node = G16ComponentNode(
            node_name=child_spec.node_name or child_component.component_name,
            component=child_component,
            include_in_aggregation=child_spec.include_in_aggregation,
            include_in_render=child_spec.include_in_render,
        )
        parent_node.children.append(child_node)
        return child_node

    @classmethod
    def expand_synthetic_children(cls, tree: G16ComponentTree) -> G16ComponentTree:
        if tree.synthetic_children_expanded:
            return tree
        for node in tree.root.children:
            component = node.component
            if component is None:
                continue
            for child_spec in component.build_synthetic_children():
                cls._append_synthetic_child(node, child_spec)
        tree.synthetic_children_expanded = True
        return tree

    @classmethod
    def from_frame_data(
        cls, frame: Any, *, expand_synthetic_children: bool = True
    ) -> G16ComponentTree:
        root = G16ComponentNode(node_name="g16.frame")
        for component_cls in G16_COMPONENT_REGISTRY:
            if not component_cls.can_build_from_frame(frame):
                continue
            payload = component_cls.build_payload_from_frame(frame)
            if not payload:
                continue
            component = component_cls(raw_text="", payload=payload)
            root.children.append(
                G16ComponentNode(node_name=component.component_name, component=component)
            )
        tree = G16ComponentTree(root=root, source_text="", only_extract_structure=False)
        if expand_synthetic_children:
            cls.expand_synthetic_children(tree)
        return tree
