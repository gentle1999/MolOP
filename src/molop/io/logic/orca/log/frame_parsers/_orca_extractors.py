from __future__ import annotations

from typing import Any

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
    Vibrations,
)
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns
from molop.unit import atom_ureg


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


def _parse_labeled_charge_block(text: str, header: str) -> list[float]:
    start = text.rfind(header)
    if start < 0:
        return []
    block = _extract_until_blank_or_rule(text, start + len(header))
    values: list[float] = []
    for line in block.splitlines():
        matched = orca_log_patterns.MULLIKEN_CHARGE_ROW.match(line)
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
        matched = orca_log_patterns.HIRSHFELD_ROW.match(line)
        if matched is None:
            continue
        charges.append(_as_float(matched.group("charge")))
        spins.append(_as_float(matched.group("spin")))
    return charges, spins


def extract_orca_coords(text: str) -> tuple[list[int], Any] | tuple[None, None]:
    matches = orca_log_patterns.COORD_HEADER.find_matches(text)
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
        matched = orca_log_patterns.COORD_ROW.match(line)
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


def extract_orca_energies(text: str) -> Energies | None:
    energy_dict: dict[str, Any] = {}
    if matches := orca_log_patterns.SCF_ENERGY.find_matches(text):
        energy_dict["reference_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    for matched in orca_log_patterns.MP2_ENERGY.find_matches(text):
        value = matched.group("mp2_total") or matched.group("mp2_corr")
        if value is not None:
            energy_dict["mp2_energy"] = _as_float(value) * atom_ureg.hartree
    if matches := orca_log_patterns.MP3_ENERGY.find_matches(text):
        energy_dict["mp3_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    if matches := orca_log_patterns.CCSD_ENERGY.find_matches(text):
        energy_dict["ccsd_energy"] = _as_float(matches[-1].group("energy")) * atom_ureg.hartree
    if matches := orca_log_patterns.FINAL_ENERGY.find_matches(text):
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


def extract_orca_forces(text: str, num_atoms: int | None) -> Any | None:
    matches = orca_log_patterns.GRADIENT_HEADER.find_matches(text)
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
        matched = orca_log_patterns.GRADIENT_ROW.match(line)
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


def extract_orca_vibrations(text: str, num_atoms: int | None) -> Vibrations | None:
    start = text.rfind("VIBRATIONAL FREQUENCIES")
    if start < 0:
        return None
    block = text[start:]
    frequencies: list[float] = []
    for line in block.splitlines():
        if frequencies and line.strip().startswith("NORMAL MODES"):
            break
        matched = orca_log_patterns.FREQUENCY.match(line)
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


def extract_orca_populations(text: str) -> ChargeSpinPopulations | None:
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


def extract_orca_polarizability(text: str) -> Polarizability | None:
    polar_dict: dict[str, Any] = {}
    if matches := orca_log_patterns.DIPOLE.find_matches(text):
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
    if matches := orca_log_patterns.POLAR_ISOTROPIC.find_matches(text):
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
        values = [matched.group("value") for matched in orca_log_patterns.FLOAT.find_matches(line)]
        if len(values) < 3:
            return None
        tensor.append([_as_float(value) for value in values[:3]])
    return np.asarray(tensor, dtype=float) * atom_ureg.bohr**3


def extract_orca_geometry_optimization_status(text: str) -> GeometryOptimizationStatus | None:
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
        metric_pattern = orca_log_patterns.optimization_metric(label)
        if matches := metric_pattern.find_matches(text):
            status_dict[field] = abs(_as_float(matches[-1].group("value")))
    return GeometryOptimizationStatus.model_validate(status_dict)


def extract_orca_solvent(text: str) -> ImplicitSolvation | None:
    solvent_dict: dict[str, Any] = {}
    if "CPCM SOLVATION MODEL" in text:
        solvent_dict["solvent_model"] = "CPCM"
    if "Your calculation utilizes the SMD solvation module" in text:
        solvent_dict["solvent_model"] = "SMD"
    if matches := orca_log_patterns.SOLVENT_NAME.find_matches(text):
        solvent_dict["solvent"] = matches[-1].group("solvent").strip().lower()
    if matches := orca_log_patterns.SOLVENT_EPSILON.find_matches(text):
        solvent_dict["solvent_epsilon"] = _as_float(matches[-1].group("epsilon"))
    if not solvent_dict:
        return None
    return ImplicitSolvation.model_validate(solvent_dict)


def extract_orca_electronic_states(text: str, method: str | None) -> ElectronicStates | None:
    states: dict[int, ElectronicState] = {}
    for matched in orca_log_patterns.STATE.find_matches(text):
        root = int(matched.group("root"))
        tail = matched.group("tail")
        multiplicity = None
        if mult_matches := orca_log_patterns.STATE_MULTIPLICITY.find_matches(tail):
            multiplicity = int(mult_matches[-1].group("multiplicity"))
        irrep = None
        if sym_matches := orca_log_patterns.STATE_IRREP.find_matches(tail):
            irrep = sym_matches[-1].group("irrep")
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

    for matched in orca_log_patterns.ABS_TRANSITION.find_matches(text):
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

    for matched in orca_log_patterns.ROOT_TRANSITION.find_matches(text):
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
    values = [
        matched.group("value")
        for matched in orca_log_patterns.FLOAT.find_matches(text[start:line_end])
    ]
    if len(values) < 8:
        return None
    try:
        return (
            np.asarray([_as_float(value) for value in values[-3:]], dtype=float) * atom_ureg.debye
        )
    except ValueError:
        return None
