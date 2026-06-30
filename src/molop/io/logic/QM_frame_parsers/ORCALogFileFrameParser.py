from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any, cast

import numpy as np
from rdkit import Chem

from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    ElectronicState,
    ElectronicStates,
    Energies,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    Polarizability,
    Status,
    Vibrations,
)
from molop.io.base_models.FrameParser import BaseFrameParser
from molop.io.logic.QM_frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.unit import atom_ureg


_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_COORD_HEADER_RE = re.compile(
    r"^-+\s*\nCARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\s*$",
    re.MULTILINE,
)
_COORD_ROW_RE = re.compile(
    rf"^\s*(?P<symbol>[A-Z][a-z]?)\s+"
    rf"(?P<x>{_FLOAT_RE})\s+(?P<y>{_FLOAT_RE})\s+(?P<z>{_FLOAT_RE})\s*$"
)
_FINAL_ENERGY_RE = re.compile(rf"FINAL SINGLE POINT ENERGY\s+(?P<energy>{_FLOAT_RE})")
_SCF_ENERGY_RE = re.compile(rf"Total Energy\s*:\s*(?P<energy>{_FLOAT_RE})\s+Eh")
_MP2_ENERGY_RE = re.compile(
    rf"(?:MP2 TOTAL ENERGY:\s*(?P<mp2_total>{_FLOAT_RE})\s*Eh|E\(MP2\)\s*=\s*(?P<mp2_corr>{_FLOAT_RE}))",
    re.I,
)
_MP3_ENERGY_RE = re.compile(rf"E\(MP3\)\s*=\s*(?P<energy>{_FLOAT_RE})", re.I)
_CCSD_ENERGY_RE = re.compile(rf"E\(CCSD(?:\(T\))?\)\s*=\s*(?P<energy>{_FLOAT_RE})", re.I)
_GRADIENT_HEADER_RE = re.compile(
    r"^-+\s*\nCARTESIAN GRADIENT(?: \(NUMERICAL\))?\s*\n-+\s*$",
    re.MULTILINE,
)
_GRADIENT_ROW_RE = re.compile(
    rf"^\s*\d+\s+[A-Z][a-z]?\s*:\s+(?P<x>{_FLOAT_RE})\s+"
    rf"(?P<y>{_FLOAT_RE})\s+(?P<z>{_FLOAT_RE})\s*$"
)
_MULLIKEN_CHARGE_RE = re.compile(rf"^\s*\d+\s+[A-Z][a-z]?\s*:\s*(?P<value>{_FLOAT_RE})\s*$")
_HIRSHFELD_ROW_RE = re.compile(
    rf"^\s*\d+\s+[A-Z][a-z]?\s+(?P<charge>{_FLOAT_RE})\s+(?P<spin>{_FLOAT_RE})\s*$"
)
_DIPOLE_RE = re.compile(
    rf"Total Dipole Moment\s*:\s*(?P<x>{_FLOAT_RE})\s+(?P<y>{_FLOAT_RE})\s+(?P<z>{_FLOAT_RE})"
)
_POLAR_ISO_RE = re.compile(rf"Isotropic polarizability\s*:\s*(?P<value>{_FLOAT_RE})", re.I)
_RUN_TIME_RE = re.compile(
    r"TOTAL RUN TIME:\s*"
    r"(?P<days>\d+)\s+days\s+"
    r"(?P<hours>\d+)\s+hours\s+"
    r"(?P<minutes>\d+)\s+minutes?\s+"
    r"(?P<seconds>\d+)\s+seconds?\s+"
    r"(?P<msec>\d+)\s+msec",
    re.I,
)
_FREQUENCY_RE = re.compile(
    rf"^\s*(?P<idx>\d+):\s*(?P<frequency>{_FLOAT_RE})\s+cm\*\*-1(?:\s+(?P<label>\S+))?\s*$"
)
_STATE_RE = re.compile(
    rf"STATE\s+(?P<root>\d+):\s+E=\s*(?P<au>{_FLOAT_RE})\s+au\s+"
    rf"(?P<ev>{_FLOAT_RE})\s+eV(?P<tail>.*)$",
    re.I,
)
_ABS_TRANSITION_RE = re.compile(
    rf"^\s*\S+\s+->\s+(?P<root>\d+)-(?P<label>\S+)\s+"
    rf"(?P<ev>{_FLOAT_RE})\s+(?P<cm>{_FLOAT_RE})\s+(?P<nm>{_FLOAT_RE})\s+"
    rf"(?P<fosc>{_FLOAT_RE})"
)
_ROOT_TRANSITION_RE = re.compile(
    rf"^\s*(?P<root>\d+)\s+(?P<ev>{_FLOAT_RE})\s+(?P<cm>{_FLOAT_RE})\s+"
    rf"(?P<nm>{_FLOAT_RE})\s+(?P<fosc>{_FLOAT_RE})"
)


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _extract_until_blank_or_rule(text: str, start: int) -> str:
    lines = text[start:].splitlines()
    collected: list[str] = []
    seen_payload = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            if seen_payload:
                break
            collected.append(line)
            continue
        if seen_payload and stripped.startswith("----"):
            break
        collected.append(line)
        if not stripped.startswith("----"):
            seen_payload = True
    return "\n".join(collected)


def _last_section_after_header(text: str, header: str) -> str | None:
    idx = text.rfind(header)
    if idx < 0:
        return None
    return _extract_until_blank_or_rule(text, idx + len(header))


def _parse_labeled_charge_block(text: str, header: str) -> list[float]:
    start = text.rfind(header)
    if start < 0:
        return []
    block = _extract_until_blank_or_rule(text, start + len(header))
    values: list[float] = []
    for line in block.splitlines():
        matched = _MULLIKEN_CHARGE_RE.match(line)
        if matched is not None:
            values.append(_as_float(matched.group("value")))
    return values


def _parse_hirshfeld(text: str) -> tuple[list[float], list[float]]:
    start = text.rfind("HIRSHFELD ANALYSIS")
    if start < 0:
        return [], []
    block = _extract_until_blank_or_rule(text, start + len("HIRSHFELD ANALYSIS"))
    charges: list[float] = []
    spins: list[float] = []
    for line in block.splitlines():
        matched = _HIRSHFELD_ROW_RE.match(line)
        if matched is None:
            continue
        charges.append(_as_float(matched.group("charge")))
        spins.append(_as_float(matched.group("spin")))
    return charges, spins


def _parse_coords(text: str) -> tuple[list[int], Any] | tuple[None, None]:
    matches = list(_COORD_HEADER_RE.finditer(text))
    if not matches:
        return None, None
    start = matches[-1].end()
    rows: list[tuple[int, list[float]]] = []
    pt = Chem.GetPeriodicTable()
    for line in text[start:].splitlines():
        if not line.strip():
            if rows:
                break
            continue
        if line.strip().startswith("-"):
            if rows:
                break
            continue
        matched = _COORD_ROW_RE.match(line)
        if matched is None:
            if rows:
                break
            continue
        symbol = matched.group("symbol")
        atomic_number = pt.GetAtomicNumber(symbol)
        if atomic_number <= 0:
            continue
        rows.append(
            (
                atomic_number,
                [
                    _as_float(matched.group("x")),
                    _as_float(matched.group("y")),
                    _as_float(matched.group("z")),
                ],
            )
        )
    if not rows:
        return None, None
    atoms = [row[0] for row in rows]
    coords = np.asarray([row[1] for row in rows], dtype=float) * atom_ureg.angstrom
    return atoms, coords


def _parse_energies(text: str) -> Energies | None:
    energy_dict: dict[str, Any] = {}
    if matches := list(_SCF_ENERGY_RE.finditer(text)):
        energy_dict["reference_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    for matched in _MP2_ENERGY_RE.finditer(text):
        value = matched.group("mp2_total") or matched.group("mp2_corr")
        if value is not None:
            energy_dict["mp2_energy"] = _as_float(value) * atom_ureg.hartree
    if matches := list(_MP3_ENERGY_RE.finditer(text)):
        energy_dict["mp3_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    if matches := list(_CCSD_ENERGY_RE.finditer(text)):
        energy_dict["ccsd_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    if matches := list(_FINAL_ENERGY_RE.finditer(text)):
        final_energy = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
        if "ccsd_energy" in energy_dict:
            energy_dict["ccsd_energy"] = final_energy
        elif "mp3_energy" in energy_dict:
            energy_dict["mp3_energy"] = final_energy
        elif "mp2_energy" in energy_dict:
            energy_dict["mp2_energy"] = final_energy
        elif "reference_energy" in energy_dict:
            energy_dict["reference_energy"] = final_energy
        else:
            energy_dict["electronic_energy"] = final_energy
    if not energy_dict:
        return None
    return Energies.model_validate(energy_dict)


def _parse_forces(text: str, num_atoms: int | None) -> Any | None:
    matches = list(_GRADIENT_HEADER_RE.finditer(text))
    if not matches:
        if "NORM OF THE MP2 GRADIENT" in text and num_atoms:
            return np.zeros((num_atoms, 3), dtype=float) * atom_ureg.Unit("hartree / bohr")
        return None
    start = matches[-1].end()
    rows: list[list[float]] = []
    for line in text[start:].splitlines():
        if not line.strip():
            if rows:
                break
            continue
        matched = _GRADIENT_ROW_RE.match(line)
        if matched is None:
            if rows:
                break
            continue
        rows.append(
            [
                _as_float(matched.group("x")),
                _as_float(matched.group("y")),
                _as_float(matched.group("z")),
            ]
        )
    if not rows:
        return None
    return np.asarray(rows, dtype=float) * atom_ureg.Unit("hartree / bohr")


def _parse_vibrations(text: str, num_atoms: int | None) -> Vibrations | None:
    start = text.rfind("VIBRATIONAL FREQUENCIES")
    if start < 0:
        return None
    block = text[start:]
    frequencies: list[float] = []
    for line in block.splitlines():
        if frequencies and line.strip().startswith("NORMAL MODES"):
            break
        matched = _FREQUENCY_RE.match(line)
        if matched is None:
            continue
        frequencies.append(_as_float(matched.group("frequency")))
    if not frequencies:
        return None
    if num_atoms is not None and num_atoms > 1:
        expected_counts = [num_atoms * 3 - 6, num_atoms * 3 - 5]
        for expected_count in expected_counts:
            if expected_count > 0 and len(frequencies) > expected_count:
                leading = frequencies[: len(frequencies) - expected_count]
                if all(abs(value) < 1.0e-6 for value in leading):
                    frequencies = frequencies[-expected_count:]
                    break
    return Vibrations(frequencies=np.asarray(frequencies, dtype=float) * atom_ureg.cm_1)


def _parse_populations(text: str) -> ChargeSpinPopulations | None:
    pop_dict: dict[str, Any] = {}
    mulliken = _parse_labeled_charge_block(text, "MULLIKEN ATOMIC CHARGES")
    if mulliken:
        pop_dict["mulliken_charges"] = mulliken
    lowdin = _parse_labeled_charge_block(text, "LOEWDIN ATOMIC CHARGES")
    if lowdin:
        pop_dict["lowdin_charges"] = lowdin
    hirshfeld_charges, hirshfeld_spins = _parse_hirshfeld(text)
    if hirshfeld_charges:
        pop_dict["hirshfeld_charges"] = hirshfeld_charges
    if hirshfeld_spins:
        pop_dict["hirshfeld_spins"] = hirshfeld_spins
    if not pop_dict:
        return None
    try:
        return ChargeSpinPopulations.model_validate(pop_dict)
    except Exception:
        # ORCA can print only a subset of analyses for some jobs. Keep the
        # successfully parsed first population if lengths differ.
        for key in ("mulliken_charges", "lowdin_charges", "hirshfeld_charges"):
            if key in pop_dict:
                return ChargeSpinPopulations.model_validate({key: pop_dict[key]})
    return None


def _parse_polarizability(text: str) -> Polarizability | None:
    polar_dict: dict[str, Any] = {}
    if matches := list(_DIPOLE_RE.finditer(text)):
        matched = matches[-1]
        polar_dict["dipole"] = (
            np.asarray(
                [
                    _as_float(matched.group("x")),
                    _as_float(matched.group("y")),
                    _as_float(matched.group("z")),
                ],
                dtype=float,
            )
            * atom_ureg.debye
        )
        polar_dict["electric_dipole_moment"] = polar_dict["dipole"]
    if matches := list(_POLAR_ISO_RE.finditer(text)):
        polar_dict["isotropic_polarizability"] = (
            _as_float(matches[-1].group("value")) * atom_ureg.bohr**3
        )
    tensor = _parse_polarizability_tensor(text)
    if tensor is not None:
        polar_dict["polarizability_tensor"] = tensor
    if not polar_dict:
        return None
    return Polarizability.model_validate(polar_dict)


def _parse_polarizability_tensor(text: str) -> Any | None:
    start = text.rfind("The raw cartesian tensor (atomic units):")
    if start < 0:
        return None
    lines = text[start:].splitlines()[1:4]
    if len(lines) < 3:
        return None
    tensor: list[list[float]] = []
    for line in lines:
        values = re.findall(_FLOAT_RE, line)
        if len(values) < 3:
            return None
        tensor.append([_as_float(value) for value in values[:3]])
    return np.asarray(tensor, dtype=float) * atom_ureg.bohr**3


def _parse_status(text: str) -> Status | None:
    if "****ORCA TERMINATED NORMALLY****" in text:
        return Status(normal_terminated=True, scf_converged=True)
    if "ORCA finished by error termination" in text or "ORCA TERMINATED ABNORMALLY" in text:
        return Status(normal_terminated=False, scf_converged=False)
    return None


def _parse_running_time(text: str) -> Any | None:
    if matches := list(_RUN_TIME_RE.finditer(text)):
        matched = matches[-1]
        seconds = (
            int(matched.group("days")) * 86400
            + int(matched.group("hours")) * 3600
            + int(matched.group("minutes")) * 60
            + int(matched.group("seconds"))
            + int(matched.group("msec")) / 1000
        )
        return seconds * atom_ureg.second
    return None


def _parse_geometry_optimization_status(text: str) -> GeometryOptimizationStatus | None:
    if "THE OPTIMIZATION HAS CONVERGED" in text or "OPTIMIZATION RUN DONE" in text:
        return GeometryOptimizationStatus(geometry_optimized=True)
    if "GEOMETRY OPTIMIZATION CYCLE" not in text:
        return None
    status_dict: dict[str, Any] = {"geometry_optimized": False}
    metric_map = {
        "Energy change": "energy_change",
        "RMS gradient": "rms_force",
        "MAX gradient": "max_force",
        "RMS step": "rms_displacement",
        "MAX step": "max_displacement",
    }
    for label, field in metric_map.items():
        pattern = re.compile(rf"{re.escape(label)}\s+({_FLOAT_RE})", re.I)
        if matches := list(pattern.finditer(text)):
            status_dict[field] = abs(_as_float(matches[-1].group(1)))
    return GeometryOptimizationStatus.model_validate(status_dict)


def _parse_solvent(text: str) -> ImplicitSolvation | None:
    solvent_dict: dict[str, Any] = {}
    if "CPCM SOLVATION MODEL" in text:
        solvent_dict["solvent_model"] = "CPCM"
    if "Your calculation utilizes the SMD solvation module" in text:
        solvent_dict["solvent_model"] = "SMD"
    if match := re.search(r"Solvent:\s*\.\.\.\s*(?P<solvent>\S+)", text):
        solvent_dict["solvent"] = match.group("solvent").strip().lower()
    if match := re.search(rf"Epsilon\s*\.\.\.\s*(?P<epsilon>{_FLOAT_RE})", text):
        solvent_dict["solvent_epsilon"] = _as_float(match.group("epsilon"))
    if not solvent_dict:
        return None
    return ImplicitSolvation.model_validate(solvent_dict)


def _parse_electronic_states(text: str, method: str | None) -> ElectronicStates | None:
    states: dict[int, ElectronicState] = {}
    for matched in _STATE_RE.finditer(text):
        root = int(matched.group("root"))
        tail = matched.group("tail")
        multiplicity = None
        if mult_match := re.search(r"\bMult\s+(\d+)", tail):
            multiplicity = int(mult_match.group(1))
        irrep = None
        if sym_match := re.search(r"\bSym:\s*(\S+)", tail):
            irrep = sym_match.group(1)
        states[root] = ElectronicState(
            state_index=root,
            root=root,
            label=f"STATE {root}",
            multiplicity=multiplicity,
            irrep=irrep,
            method=method,
            excitation_energy=_as_float(matched.group("ev")) * atom_ureg.eV,
            source="ORCA output",
        )

    for matched in _ABS_TRANSITION_RE.finditer(text):
        root = int(matched.group("root"))
        state = states.get(root)
        transition_dipole = _transition_dipole_from_line(matched.string, matched.start())
        if state is None:
            state = ElectronicState(
                state_index=root,
                root=root,
                label=matched.group("label"),
                method=method,
                excitation_energy=_as_float(matched.group("ev")) * atom_ureg.eV,
                source="ORCA absorption spectrum",
            )
        state.oscillator_strength = _as_float(matched.group("fosc"))
        if transition_dipole is not None:
            state.transition_dipole = transition_dipole
        states[root] = state

    for matched in _ROOT_TRANSITION_RE.finditer(text):
        root = int(matched.group("root"))
        if root in states:
            continue
        states[root] = ElectronicState(
            state_index=root,
            root=root,
            label=f"ROOT {root}",
            method=method,
            excitation_energy=_as_float(matched.group("ev")) * atom_ureg.eV,
            oscillator_strength=_as_float(matched.group("fosc")),
            source="ORCA absorption spectrum",
        )

    if not states:
        if any(
            marker in text
            for marker in (
                "TD-DFT/TDA EXCITED STATES",
                "ADC(2) RESULTS",
                "EOM-CCSD RESULTS",
                "STEOM-CCSD RESULTS",
                "ROCIS-EXCITATION SPECTRA",
            )
        ):
            return ElectronicStates(
                states=[
                    ElectronicState(
                        state_index=1,
                        root=1,
                        label="ORCA electronic state",
                        method=method,
                        source="ORCA output",
                    )
                ]
            )
        return None
    return ElectronicStates(states=[states[root] for root in sorted(states)])


def _transition_dipole_from_line(text: str, start: int) -> Any | None:
    line_end = text.find("\n", start)
    if line_end < 0:
        line_end = len(text)
    values = re.findall(_FLOAT_RE, text[start:line_end])
    if len(values) < 8:
        return None
    try:
        return np.asarray([_as_float(value) for value in values[-3:]], dtype=float) * atom_ureg.debye
    except ValueError:
        return None


class ORCALogFileFrameParserMixin:
    def _parse_frame(self) -> Mapping[str, Any]:
        text = cast(Any, self)._block
        infos: dict[str, Any] = {"qm_software": "ORCA"}

        atoms, coords = _parse_coords(text)
        if atoms is not None and coords is not None:
            infos["atoms"] = atoms
            infos["coords"] = coords
        num_atoms = len(atoms) if atoms is not None else None

        energies = _parse_energies(text)
        if energies is not None:
            infos["energies"] = energies
        if cast(Any, self).only_extract_structure:
            return infos
        forces = _parse_forces(text, num_atoms)
        if forces is not None:
            infos["forces"] = forces
        vibrations = _parse_vibrations(text, num_atoms)
        if vibrations is not None:
            infos["vibrations"] = vibrations
        populations = _parse_populations(text)
        if populations is not None:
            infos["charge_spin_populations"] = populations
        polarizability = _parse_polarizability(text)
        if polarizability is not None:
            infos["polarizability"] = polarizability
        status = _parse_status(text)
        if status is not None:
            infos["status"] = status
        running_time = _parse_running_time(text)
        if running_time is not None:
            infos["running_time"] = running_time
        opt_status = _parse_geometry_optimization_status(text)
        if opt_status is not None:
            infos["geometry_optimization_status"] = opt_status
        solvent = _parse_solvent(text)
        if solvent is not None:
            infos["solvent"] = solvent
        method = getattr(self, "_orca_method", None)
        electronic_states = _parse_electronic_states(text, method)
        if electronic_states is not None:
            infos["electronic_states"] = electronic_states
        return infos


class ORCALogFileFrameParserMemory(
    ORCALogFileFrameParserMixin, BaseFrameParser[ORCALogFileFrameMemory]
):
    _file_frame_class_ = ORCALogFileFrameMemory


class ORCALogFileFrameParserDisk(ORCALogFileFrameParserMixin, BaseFrameParser[ORCALogFileFrameDisk]):
    _file_frame_class_ = ORCALogFileFrameDisk
