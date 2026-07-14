from __future__ import annotations

from typing import Any, Literal

import numpy as np
from rdkit import Chem

from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    ElectronicState,
    ElectronicStates,
    Energies,
    EnergyObservation,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    Polarizability,
    Vibrations,
)
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns
from molop.unit import atom_ureg


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _decimal_places(value: str) -> int:
    exponent_positions = [
        position for marker in ("E", "e", "D", "d") if (position := value.find(marker)) >= 0
    ]
    mantissa = value[: min(exponent_positions)] if exponent_positions else value
    return len(mantissa.partition(".")[2])


def _normalized_optional_text(value: Any, *, upper: bool = False) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    return normalized.upper() if upper else normalized


def _orca_concrete_method(model_chemistry: Any) -> str:
    method_family = _normalized_optional_text(
        getattr(model_chemistry, "method_family", None), upper=True
    )
    method = _normalized_optional_text(getattr(model_chemistry, "method", None), upper=True)
    functional = _normalized_optional_text(getattr(model_chemistry, "functional", None))
    dispersion = _normalized_optional_text(
        getattr(model_chemistry, "dispersion_correction", None), upper=True
    )
    if method_family == "DFT" and functional is not None:
        if dispersion is not None and functional.upper().endswith(f"-{dispersion}"):
            functional = functional[: -(len(dispersion) + 1)]
        return functional
    return method or method_family or "electronic"


def _orca_coupled_cluster_methods(model_chemistry: Any) -> tuple[str, str]:
    concrete_method = _orca_concrete_method(model_chemistry)
    upper_method = concrete_method.upper()
    if upper_method == "CCSD(T)" or upper_method.startswith("DLPNO-CCSD(T)"):
        triples_position = upper_method.index("(T)")
        ccsd_method = (
            concrete_method[:triples_position] + concrete_method[triples_position + len("(T)") :]
        )
        return ccsd_method, concrete_method
    if upper_method == "CCSD" or upper_method.startswith("DLPNO-CCSD"):
        return concrete_method, f"{concrete_method}(T)"
    return "CCSD", "CCSD(T)"


def _energy_observation(
    method: str,
    quantity_semantics: Literal["total_energy", "correlation_correction", "component"],
    value: str,
    source_label: str,
) -> EnergyObservation:
    return EnergyObservation(
        method=method,
        quantity_semantics=quantity_semantics,
        value=_as_float(value) * atom_ureg.hartree,
        source_label=source_label,
    )


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


def extract_orca_coords(
    text: str,
    *,
    capture_source_evidence: bool = False,
) -> tuple[list[int] | None, Any | None, int | None]:
    matches = orca_log_patterns.COORD_HEADER.find_matches(text)
    if not matches:
        return None, None, None
    start = matches[-1].end()
    rows: list[tuple[int, list[float]]] = []
    decimal_places: list[int] | None = [] if capture_source_evidence else None
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
        if decimal_places is not None:
            decimal_places.extend(_decimal_places(matched.group(axis)) for axis in ("x", "y", "z"))
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
        return None, None, None
    atoms = [row[0] for row in rows]
    coords = np.asarray([row[1] for row in rows], dtype=float) * atom_ureg.angstrom
    coordinate_decimal_places = min(decimal_places) if decimal_places else None
    return atoms, coords, coordinate_decimal_places


def extract_orca_energies(
    text: str,
    *,
    capture_source_evidence: bool = False,
    model_chemistry: Any = None,
) -> Energies | None:
    energy_dict: dict[str, Any] = {}
    observations: list[EnergyObservation] = []

    def observe(
        method: str,
        quantity_semantics: Literal["total_energy", "correlation_correction", "component"],
        value: str,
        source_label: str,
    ) -> None:
        if capture_source_evidence:
            observations.append(
                _energy_observation(method, quantity_semantics, value, source_label)
            )

    if scf_matches := orca_log_patterns.SCF_ENERGY.find_matches(text):
        value = scf_matches[-1].group("energy")
        energy_dict["reference_energy"] = _as_float(value) * atom_ureg.hartree
        observe("reference", "total_energy", value, "Total Energy")

    mp2_value: str | None = None
    mp2_source_label: str | None = None
    for matched in orca_log_patterns.MP2_ENERGY.find_matches(text):
        value = matched.group("mp2_total") or matched.group("mp2_equation_total")
        if value is not None:
            energy_dict["mp2_energy"] = _as_float(value) * atom_ureg.hartree
            mp2_value = value
            mp2_source_label = (
                "MP2 TOTAL ENERGY" if matched.group("mp2_total") is not None else "E(MP2)"
            )
    if mp2_value is not None and mp2_source_label is not None:
        observe("MP2", "total_energy", mp2_value, mp2_source_label)

    if capture_source_evidence and (
        matches := orca_log_patterns.MP2_CORRELATION_ENERGY.find_matches(text)
    ):
        matched = matches[-1]
        for group_name, source_label in (
            ("labeled", "MP2 CORRELATION ENERGY"),
            ("equation", "EC(MP2)"),
            ("component", "E(MP2)"),
        ):
            if value := matched.group(group_name):
                observe("MP2", "correlation_correction", value, source_label)
                break

    if mp3_matches := orca_log_patterns.MP3_ENERGY.find_matches(text):
        value = mp3_matches[-1].group("energy")
        energy_dict["mp3_energy"] = _as_float(value) * atom_ureg.hartree
        observe("MP3", "total_energy", value, "E(MP3)")
    if capture_source_evidence:
        if matches := orca_log_patterns.MP3_CORRELATION_ENERGY.find_matches(text):
            observe(
                "MP3",
                "correlation_correction",
                matches[-1].group("energy"),
                "EC(MP3)",
            )
        if matches := orca_log_patterns.MP3_COMPONENT_ENERGY.find_matches(text):
            observe("MP3", "component", matches[-1].group("energy"), "E3")

    ccsd_method, ccsd_t_method = (
        _orca_coupled_cluster_methods(model_chemistry)
        if capture_source_evidence
        else ("CCSD", "CCSD(T)")
    )
    ccsd_matches = orca_log_patterns.CCSD_ENERGY.find_matches(text)
    ccsd_total_matches = (
        orca_log_patterns.CCSD_TOTAL_ENERGY.find_matches(text)
        if capture_source_evidence or not ccsd_matches
        else []
    )
    selected_ccsd_matches = ccsd_matches or ccsd_total_matches
    if selected_ccsd_matches:
        energy_dict["ccsd_energy"] = (
            _as_float(selected_ccsd_matches[-1].group("energy")) * atom_ureg.hartree
        )
    if ccsd_matches:
        observe(ccsd_method, "total_energy", ccsd_matches[-1].group("energy"), "E(CCSD)")
    if ccsd_total_matches:
        observe(
            ccsd_method,
            "total_energy",
            ccsd_total_matches[-1].group("energy"),
            "E(TOT)",
        )

    if ccsd_t_matches := orca_log_patterns.CCSD_T_ENERGY.find_matches(text):
        value = ccsd_t_matches[-1].group("energy")
        energy_dict["ccsd_t_energy"] = _as_float(value) * atom_ureg.hartree
        observe(ccsd_t_method, "total_energy", value, "E(CCSD(T))")

    if capture_source_evidence:
        if matches := orca_log_patterns.CCSD_CORRELATION_ENERGY.find_matches(text):
            observe(
                ccsd_method,
                "correlation_correction",
                matches[-1].group("energy"),
                "E(CORR)",
            )
        if matches := orca_log_patterns.TRIPLES_CORRECTION.find_matches(text):
            observe(
                ccsd_t_method,
                "component",
                matches[-1].group("energy"),
                "Triples Correction (T)",
            )
        if matches := orca_log_patterns.SCALED_TRIPLES_CORRECTION.find_matches(text):
            observe(
                ccsd_t_method,
                "component",
                matches[-1].group("energy"),
                "Scaled triples correction (T)",
            )
        if matches := orca_log_patterns.FINAL_CORRELATION_ENERGY.find_matches(text):
            concrete_method = _orca_concrete_method(model_chemistry).upper()
            observe(
                ccsd_t_method if "CCSD(T)" in concrete_method else ccsd_method,
                "correlation_correction",
                matches[-1].group("energy"),
                "Final correlation energy",
            )

    if final_matches := orca_log_patterns.FINAL_ENERGY.find_matches(text):
        value = final_matches[-1].group("energy")
        energy_dict["electronic_energy"] = _as_float(value) * atom_ureg.hartree
        if capture_source_evidence:
            observe(
                _orca_concrete_method(model_chemistry),
                "total_energy",
                value,
                "FINAL SINGLE POINT ENERGY",
            )
    if observations:
        energy_dict["observations"] = observations
    if not energy_dict:
        return None
    return Energies.model_validate(energy_dict)


def extract_orca_forces(text: str, num_atoms: int | None) -> Any | None:
    candidates = [
        (matched.start(), matched.end(), orca_log_patterns.GRADIENT_ROW)
        for matched in orca_log_patterns.GRADIENT_HEADER.find_matches(text)
    ]
    candidates.extend(
        (matched.start(), matched.end(), orca_log_patterns.MP2_GRADIENT_ROW)
        for matched in orca_log_patterns.MP2_GRADIENT_HEADER.find_matches(text)
    )
    if not candidates:
        return None
    _, start, row_pattern = max(candidates, key=lambda candidate: candidate[0])
    rows: list[list[float]] = []
    for line in text[start:].splitlines():
        if not line.strip():
            if rows:
                break
            continue
        matched = row_pattern.match(line)
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
    if num_atoms is not None and len(rows) != num_atoms:
        return None
    return -np.asarray(rows, dtype=float) * atom_ureg.Unit("hartree / bohr")


def extract_orca_vibrations(text: str, num_atoms: int | None) -> Vibrations | None:
    start = text.rfind("VIBRATIONAL FREQUENCIES")
    if start < 0:
        return None
    block = text[start:]
    frequencies: list[float] = []
    mode_indices: list[int] = []
    for line in block.splitlines():
        if frequencies and line.strip().startswith("NORMAL MODES"):
            break
        matched = orca_log_patterns.FREQUENCY.match(line)
        if matched is None:
            continue
        frequencies.append(_as_float(matched.group("frequency")))
        mode_indices.append(int(matched.group("idx")))
    if not frequencies:
        return None
    if num_atoms is not None and num_atoms > 1:
        expected_counts = [num_atoms * 3 - 6, num_atoms * 3 - 5]
        for expected_count in expected_counts:
            if expected_count > 0 and len(frequencies) > expected_count:
                leading = frequencies[: len(frequencies) - expected_count]
                if all(abs(value) < 1.0e-6 for value in leading):
                    frequencies = frequencies[-expected_count:]
                    mode_indices = mode_indices[-expected_count:]
                    break
    return Vibrations(
        frequencies=np.asarray(frequencies, dtype=float) * atom_ureg.cm_1,
        mode_indices=mode_indices,
    )


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


def extract_orca_geometry_optimization_status(
    text: str,
    *,
    capture_source_evidence: bool = False,
) -> GeometryOptimizationStatus | None:
    geometry_optimized = "THE OPTIMIZATION HAS CONVERGED" in text or "OPTIMIZATION RUN DONE" in text
    metric_matches = orca_log_patterns.OPTIMIZATION_CONVERGENCE_METRIC.find_matches(text)
    if not geometry_optimized and "GEOMETRY OPTIMIZATION CYCLE" not in text and not metric_matches:
        return None
    status_dict: dict[str, Any] = {"geometry_optimized": geometry_optimized}
    source_converged: dict[str, bool | None] = {}
    source_labels: dict[str, str] = {}
    metric_map = {
        "Energy change": "energy_change",
        "RMS gradient": "rms_force",
        "MAX gradient": "max_force",
        "RMS step": "rms_displacement",
        "MAX step": "max_displacement",
    }
    metric_units = {
        "energy_change": atom_ureg.hartree,
        "rms_force": atom_ureg.hartree / atom_ureg.bohr,
        "max_force": atom_ureg.hartree / atom_ureg.bohr,
        "rms_displacement": atom_ureg.bohr,
        "max_displacement": atom_ureg.bohr,
    }
    for matched in metric_matches:
        label = matched.group("label")
        field = metric_map[label]
        unit = metric_units[field]
        status_dict[field] = abs(_as_float(matched.group("value"))) * unit
        status_dict[f"{field}_threshold"] = abs(_as_float(matched.group("threshold"))) * unit
        if capture_source_evidence:
            source_converged[field] = matched.group("converged") == "YES"
            source_labels[field] = label
    if source_converged:
        status_dict["source_converged"] = source_converged
    if source_labels:
        status_dict["source_labels"] = source_labels
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
