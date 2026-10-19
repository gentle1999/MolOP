from __future__ import annotations

import shlex
from collections.abc import Iterator
from typing import Any

import numpy as np
from rdkit import Chem

from molop.io.base_models.DataClasses import (
    BondOrderMatrix,
    BondOrders,
    ChargeSpinPopulations,
    Energies,
    EnergyObservation,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    MolecularOrbitals,
    Polarizability,
    QMModelChemistry,
    QMResourceRequest,
    QMTaskRequest,
    SinglePointProperties,
    Status,
    ThermalInformations,
    Vibrations,
)
from molop.io.base_models.SearchPattern import MolOPMatch, MolOPPattern
from molop.io.codec_exceptions import FormatMismatchError, UnsupportedFormatError
from molop.io.logic.xtb.output.parsers._xtb_output_patterns import xtb_output_patterns
from molop.unit import atom_ureg


def _to_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _float_tokens(value: str) -> list[str]:
    return xtb_output_patterns.FLOAT_TOKEN.find_named_group_values(value, "value")


def _is_float_token(value: str) -> bool:
    try:
        _to_float(value)
    except ValueError:
        return False
    return True


def iter_xtb_version_matches(file_content: str) -> Iterator[MolOPMatch]:
    matches = [
        *xtb_output_patterns.MODERN_VERSION.find_matches(file_content),
        *xtb_output_patterns.LEGACY_VERSION.find_matches(file_content),
    ]
    yield from sorted(matches, key=lambda matched: matched.start())


def extract_xtb_version(file_content: str) -> str | None:
    matches = tuple(iter_xtb_version_matches(file_content))
    if not matches:
        return None
    return " ".join(matches[-1].group("version").split())


def extract_xtb_major_version(file_content: str) -> int | None:
    version = extract_xtb_version(file_content)
    if version is None:
        return None
    major = version.split(".", 1)[0]
    return int(major) if major.isdigit() else None


def ensure_xtb_output_content(file_content: str) -> None:
    prefix = file_content[:20_000]
    if not xtb_output_patterns.BANNER.find_matches(prefix):
        raise FormatMismatchError("Not an xTB output file: missing xTB banner.")
    major_version = extract_xtb_major_version(prefix)
    if major_version is None:
        raise FormatMismatchError("Not an xTB output file: missing version declaration.")
    if major_version < 5:
        raise UnsupportedFormatError(
            f"xTB output major version {major_version} is below the supported minimum 5."
        )


def _last_group(pattern: MolOPPattern, text: str, group: str = "value") -> str | None:
    values = pattern.find_named_group_values(text, group)
    if not values:
        return None
    return values[-1].strip()


def extract_xtb_program_call(file_content: str) -> str:
    return _last_group(xtb_output_patterns.PROGRAM_CALL, file_content) or ""


def extract_xtb_input_file_name(file_content: str) -> str:
    return _last_group(xtb_output_patterns.COORDINATE_FILE, file_content) or ""


def _command_tokens(program_call: str) -> list[str]:
    if not program_call:
        return []
    try:
        return shlex.split(program_call)
    except ValueError:
        return program_call.split()


def _option_value(tokens: list[str], *names: str) -> str | None:
    lowered_names = tuple(name.lower() for name in names)
    for index, token in enumerate(tokens):
        lowered = token.lower()
        for name in lowered_names:
            if lowered == name and index + 1 < len(tokens):
                return tokens[index + 1]
            if lowered.startswith(f"{name}="):
                return token.split("=", 1)[1]
    return None


def extract_xtb_method(file_content: str) -> str | None:
    if values := xtb_output_patterns.HAMILTONIAN.find_named_group_values(file_content, "value"):
        return values[-1]
    spaced_matches = xtb_output_patterns.SPACED_HAMILTONIAN.find_matches(file_content)
    if spaced_matches:
        level = "".join((spaced_matches[-1].group("level") or "").split())
        return f"GFN{level}-xTB" if level else "GFN-xTB"
    return None


def extract_xtb_charge(file_content: str, program_call: str) -> int:
    tokens = _command_tokens(program_call)
    if value := _option_value(tokens, "--charge", "--chrg", "-c"):
        return round(_to_float(value))
    if values := xtb_output_patterns.SETUP_CHARGE.find_named_group_values(file_content, "value"):
        return round(_to_float(values[-1]))
    if values := xtb_output_patterns.TOTAL_CHARGE.find_named_group_values(file_content, "value"):
        return round(_to_float(values[-1]))
    return 0


def extract_xtb_multiplicity(file_content: str, program_call: str) -> int:
    tokens = _command_tokens(program_call)
    if value := _option_value(tokens, "--uhf", "-u"):
        return max(1, round(_to_float(value)) + 1)
    if values := xtb_output_patterns.SETUP_SPIN.find_named_group_values(file_content, "value"):
        return max(1, round(_to_float(values[-1])) + 1)
    return 1


def extract_xtb_tasks(program_call: str) -> list[QMTaskRequest]:
    tokens = [token.lower() for token in _command_tokens(program_call)]

    def has(*options: str) -> bool:
        return any(
            token == option or token.startswith(f"{option}=")
            for token in tokens
            for option in options
        )

    tasks: list[QMTaskRequest] = []
    if has("--opt", "--ohess"):
        source = [token for token in tokens if token in {"--opt", "--ohess"}]
        tasks.append(QMTaskRequest(task_type="opt", derivative_order=1, source_keywords=source))
    if has("--hess", "--ohess"):
        source = [token for token in tokens if token in {"--hess", "--ohess"}]
        tasks.append(QMTaskRequest(task_type="freq", derivative_order=2, source_keywords=source))
    if has("--grad"):
        tasks.append(
            QMTaskRequest(task_type="gradient", derivative_order=1, source_keywords=["--grad"])
        )
    properties = [
        option.removeprefix("--")
        for option in ("--vip", "--vea", "--vipea", "--vfukui", "--fukui")
        if has(option)
    ]
    if not tasks or properties:
        tasks.append(
            QMTaskRequest(
                task_type="sp",
                derivative_order=0,
                properties=properties,
                source_keywords=[f"--{name}" for name in properties],
            )
        )
    return tasks


def extract_xtb_solvent(file_content: str) -> ImplicitSolvation | None:
    model = _last_group(xtb_output_patterns.SOLVATION_MODEL, file_content)
    solvent = _last_group(xtb_output_patterns.SOLVENT, file_content)
    if model is None and solvent is None:
        return None
    return ImplicitSolvation(solvent_model=model, solvent=solvent)


def extract_xtb_temperature(file_content: str):
    if values := xtb_output_patterns.VIB_TEMPERATURE.find_named_group_values(file_content, "value"):
        return _to_float(values[-1]) * atom_ureg.kelvin
    if values := xtb_output_patterns.SOLVENT_TEMPERATURE.find_named_group_values(
        file_content, "value"
    ):
        return _to_float(values[-1]) * atom_ureg.kelvin
    return None


def extract_xtb_electron_temperature(file_content: str):
    if values := xtb_output_patterns.ELECTRON_TEMPERATURE.find_named_group_values(
        file_content, "value"
    ):
        return _to_float(values[-1]) * atom_ureg.kelvin
    return None


def extract_xtb_metadata(file_content: str) -> dict[str, Any]:
    program_call = extract_xtb_program_call(file_content)
    input_file_name = extract_xtb_input_file_name(file_content)
    method = extract_xtb_method(file_content)
    solvent = extract_xtb_solvent(file_content)
    model_chemistry = QMModelChemistry(
        method_family="SEMIEMPIRICAL",
        method=method,
        solvation_model=solvent.solvent_model if solvent else None,
        solvent=solvent.solvent if solvent else None,
        raw_keywords=program_call,
    )
    metadata: dict[str, Any] = {
        "qm_software": "xTB",
        "qm_software_version": extract_xtb_version(file_content) or "",
        "input_file_name": input_file_name,
        "keywords": program_call,
        "method": method or "",
        "model_chemistry": model_chemistry,
        "task_requests": extract_xtb_tasks(program_call),
        "charge": extract_xtb_charge(file_content, program_call),
        "multiplicity": extract_xtb_multiplicity(file_content, program_call),
    }
    if solvent is not None:
        metadata["solvent"] = solvent
    if temperature := extract_xtb_temperature(file_content):
        metadata["temperature"] = temperature
    if electron_temperature := extract_xtb_electron_temperature(file_content):
        metadata["electron_temperature"] = electron_temperature
    if threads := _last_group(xtb_output_patterns.OMP_THREADS, file_content):
        num_cpu = int(threads)
        metadata["request_num_cpu"] = num_cpu
        metadata["resource_request"] = QMResourceRequest(
            num_cpu=num_cpu,
            raw=f"omp threads: {num_cpu}",
            options={"omp_threads": num_cpu},
        )
    return metadata


def _decimal_places(token: str) -> int:
    mantissa = token.replace("d", "e").replace("D", "E").split("E", 1)[0]
    return len(mantissa.partition(".")[2])


def _atoms_and_coords(
    rows: list[tuple[str, str, str, str]],
    *,
    unit,
) -> tuple[list[int], Any, int] | tuple[None, None, None]:
    if not rows:
        return None, None, None
    periodic_table = Chem.GetPeriodicTable()
    atoms: list[int] = []
    values: list[tuple[float, float, float]] = []
    precision = 0
    for symbol, x, y, z in rows:
        atomic_number = periodic_table.GetAtomicNumber(symbol.capitalize())
        if atomic_number <= 0:
            return None, None, None
        atoms.append(atomic_number)
        values.append((_to_float(x), _to_float(y), _to_float(z)))
        precision = max(precision, *(_decimal_places(token) for token in (x, y, z)))
    coords = (np.asarray(values, dtype=float) * unit).to(atom_ureg.angstrom)
    return atoms, coords, precision


def _extract_final_structure(file_content: str):
    matches = xtb_output_patterns.FINAL_STRUCTURE.find_matches(file_content)
    if not matches:
        return None, None, None
    block = file_content[matches[-1].end() :]
    end_markers = [
        position
        for marker in ("Bond Distances", "FINAL SINGLEPOINT CALCULATION")
        if (position := block.find(marker)) >= 0
    ]
    if end_markers:
        block = block[: min(end_markers)]

    if "$coord" in block:
        coord_body = block.split("$coord", 1)[1].split("$end", 1)[0]
        rows = []
        for line in coord_body.splitlines():
            tokens = line.split()
            if len(tokens) >= 4 and _is_float_token(tokens[0]):
                rows.append((tokens[3], tokens[0], tokens[1], tokens[2]))
        return _atoms_and_coords(rows, unit=atom_ureg.bohr)

    lines = block.splitlines()
    for index, line in enumerate(lines):
        if "V2000" not in line:
            continue
        tokens = line.split()
        if not tokens or not tokens[0].isdigit():
            continue
        atom_count = int(tokens[0])
        rows = []
        for atom_line in lines[index + 1 : index + 1 + atom_count]:
            atom_tokens = atom_line.split()
            if len(atom_tokens) < 4:
                break
            rows.append((atom_tokens[3], atom_tokens[0], atom_tokens[1], atom_tokens[2]))
        if len(rows) == atom_count:
            return _atoms_and_coords(rows, unit=atom_ureg.angstrom)

    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped.isdigit():
            continue
        atom_count = int(stripped)
        rows = []
        for atom_line in lines[index + 2 : index + 2 + atom_count]:
            tokens = atom_line.split()
            if len(tokens) < 4 or not tokens[0].isalpha() or len(tokens[0]) > 3:
                break
            rows.append((tokens[0], tokens[1], tokens[2], tokens[3]))
        if len(rows) == atom_count:
            return _atoms_and_coords(rows, unit=atom_ureg.angstrom)
    return None, None, None


def _extract_legacy_setup_structure(file_content: str):
    headers = xtb_output_patterns.LEGACY_COORDINATE_HEADER.find_matches(file_content)
    if not headers:
        return None, None, None
    rows: list[tuple[str, str, str, str]] = []
    for line in file_content[headers[-1].end() :].splitlines():
        tokens = line.split()
        if not tokens:
            if rows:
                break
            continue
        if len(tokens) < 8 or not tokens[0].isdigit():
            if rows:
                break
            continue
        rows.append((tokens[1], tokens[-3], tokens[-2], tokens[-1]))
    return _atoms_and_coords(rows, unit=atom_ureg.bohr)


def extract_xtb_coords(file_content: str):
    atoms, coords, precision = _extract_final_structure(file_content)
    if atoms is not None:
        return atoms, coords, precision, "final structure"
    atoms, coords, precision = _extract_legacy_setup_structure(file_content)
    if atoms is not None:
        return atoms, coords, precision, "legacy calculation setup"
    return None, None, None, None


def extract_xtb_final_energy(
    file_content: str,
    *,
    method: str | None,
    capture_source_evidence: bool,
) -> Energies | None:
    candidates: list[tuple[int, float, str]] = []
    for pattern in xtb_output_patterns.FINAL_ENERGIES:
        candidates.extend(
            (matched.start(), _to_float(matched.group("value")), matched.group().strip())
            for matched in pattern.find_matches(file_content)
        )
    if not candidates:
        return None
    _, value, source_label = max(candidates, key=lambda candidate: candidate[0])
    energy = value * atom_ureg.hartree
    observations = []
    if capture_source_evidence:
        observations.append(
            EnergyObservation(
                method=method or "xTB",
                quantity_semantics="total_energy",
                value=energy,
                source_label=source_label,
            )
        )
    return Energies(electronic_energy=energy, observations=observations)


def extract_xtb_running_time(file_content: str):
    matches = xtb_output_patterns.TOTAL_WALL_TIME.find_matches(file_content)
    matched = matches[-1] if matches else None
    if matched is None:
        matches = xtb_output_patterns.ANY_WALL_TIME.find_matches(file_content)
        matched = matches[-1] if matches else None
    if matched is None:
        return None
    seconds = (
        int(matched.group("days")) * 86_400
        + int(matched.group("hours")) * 3_600
        + int(matched.group("minutes")) * 60
        + _to_float(matched.group("seconds"))
    )
    return seconds * atom_ureg.second


def extract_xtb_status(file_content: str, *, include_termination: bool = True) -> Status:
    lower = file_content.lower()
    failure_markers = (
        "scc is not converged",
        "scc failed",
        "scf not converged",
        "convergence criteria cannot be satisfied",
    )
    if any(marker in lower for marker in failure_markers):
        scf_converged = False
    elif (
        "convergence criteria satisfied" in lower
        or extract_xtb_final_energy(file_content, method=None, capture_source_evidence=False)
        is not None
    ):
        scf_converged = True
    else:
        scf_converged = None
    normal_terminated = None
    if include_termination:
        normal_terminated = "* finished run on" in lower
    return Status(scf_converged=scf_converged, normal_terminated=normal_terminated)


def extract_xtb_optimization_status(file_content: str):
    lower = file_content.lower()
    requested = "--opt" in lower or "--ohess" in lower or "geometry optimization" in lower
    converged = "geometry optimization converged" in lower
    failed = "geometry optimization failed" in lower or "ancopt failed" in lower
    if not requested and not converged and not failed:
        return None, None, None

    values: dict[str, Any] = {
        "geometry_optimized": True if converged else False if failed else None
    }
    if matches := xtb_output_patterns.ENERGY_CONVERGENCE.find_named_group_values(
        file_content, "value"
    ):
        values["energy_change_threshold"] = abs(_to_float(matches[-1])) * atom_ureg.hartree
    if matches := xtb_output_patterns.ENERGY_CHANGE.find_named_group_values(file_content, "value"):
        values["energy_change"] = abs(_to_float(matches[-1])) * atom_ureg.hartree
    gradient_threshold = None
    if matches := xtb_output_patterns.GRADIENT_CONVERGENCE.find_named_group_values(
        file_content, "value"
    ):
        gradient_threshold = abs(_to_float(matches[-1])) * atom_ureg.Unit("hartree / bohr")
    gradient_norm = None
    if matches := xtb_output_patterns.GRADIENT_NORM.find_named_group_values(file_content, "value"):
        gradient_norm = abs(_to_float(matches[-1])) * atom_ureg.Unit("hartree / bohr")
    return GeometryOptimizationStatus.model_validate(values), gradient_norm, gradient_threshold


def _extract_indexed_values(block: str) -> list[float]:
    indexed = [
        (int(matched.group("index")), _to_float(matched.group("value")))
        for matched in xtb_output_patterns.INDEXED_VALUE.find_matches(block)
    ]
    if not indexed or any(index < 1 for index, _ in indexed):
        return []
    if len({index for index, _ in indexed}) != len(indexed):
        return []
    result = [np.nan] * max(index for index, _ in indexed)
    for index, value in indexed:
        result[index - 1] = value
    return result


def extract_xtb_vibrations(file_content: str) -> Vibrations | None:
    """Read projected modes and unit-qualified Raman activities from stdout.

    Basis: xTB v6.7.1 hessian.F90 and the locally calculated PTB water fixture;
    see docs/en/reference/formats/xtbout.md for units and companion-file limits.
    """
    starts = xtb_output_patterns.FREQUENCY_HEADER.find_matches(file_content)
    if not starts:
        return None
    block = file_content[starts[-1].end() :]
    reduced_matches = xtb_output_patterns.REDUCED_MASS_HEADER.find_matches(block)
    reduced_start = reduced_matches[0] if reduced_matches else None
    ir_matches = xtb_output_patterns.IR_HEADER.find_matches(block)
    ir_start = ir_matches[0] if ir_matches else None
    raman_matches = xtb_output_patterns.RAMAN_HEADER.find_matches(block)
    raman_start = raman_matches[0] if raman_matches else None

    frequency_block = block[: reduced_start.start() if reduced_start else len(block)]
    frequencies = [
        _to_float(token)
        for values in xtb_output_patterns.FREQUENCY_ROW.find_named_group_values(
            frequency_block, "values"
        )
        for token in _float_tokens(values)
    ]
    if not frequencies:
        return None

    leading_external_modes = 0
    for frequency in frequencies:
        if abs(frequency) > 1.0e-2:
            break
        leading_external_modes += 1
    if leading_external_modes in {5, 6}:
        frequencies = frequencies[leading_external_modes:]

    reduced_masses: list[float] = []
    if reduced_start is not None:
        end = ir_start.start() if ir_start else len(block)
        reduced_masses = _extract_indexed_values(block[reduced_start.end() : end])
        if leading_external_modes in {5, 6}:
            reduced_masses = reduced_masses[leading_external_modes:]
    intensities: list[float] = []
    if ir_start is not None:
        end = raman_start.start() if raman_start else len(block)
        intensities = _extract_indexed_values(block[ir_start.end() : end])
        if leading_external_modes in {5, 6}:
            intensities = intensities[leading_external_modes:]

    payload: dict[str, Any] = {
        "frequencies": np.asarray(frequencies) * atom_ureg.cm_1,
        "mode_indices": list(range(len(frequencies))),
        "source_program": "xTB",
    }
    if ir_start is not None:
        end = raman_start.start() if raman_start else len(block)
        printed_indices = sorted(
            {
                int(matched.group("index"))
                for matched in xtb_output_patterns.INDEXED_VALUE.find_matches(
                    block[ir_start.end() : end]
                )
            }
        )
        if leading_external_modes in {5, 6}:
            printed_indices = printed_indices[leading_external_modes:]
        if len(printed_indices) == len(frequencies):
            payload["source_mode_indices"] = printed_indices
    if matches := xtb_output_patterns.HESSIAN_METHOD.find_matches(
        file_content[: starts[-1].start()]
    ):
        payload["frequency_method"] = matches[-1].group("method")
    if len(reduced_masses) == len(frequencies) and np.isfinite(reduced_masses).all():
        payload["reduced_masses"] = np.asarray(reduced_masses) * atom_ureg.amu
    if len(intensities) == len(frequencies) and np.isfinite(intensities).all():
        payload["IR_intensities"] = np.asarray(intensities) * atom_ureg.Unit("km/mol")
    # v6.7.1 calls this "intensities", but hessian.F90 computes activities.
    # Older "(amu)" zero-filled tables have no established activity unit.
    if raman_start is not None and raman_start.group("unit") in {
        "Ä⁴*amu⁻¹",
        "Å⁴*amu⁻¹",
        "A**4/amu",
        "Å^4/amu",
    }:
        values = _extract_indexed_values(
            block[raman_start.end() :].split("output can be read", 1)[0]
        )
        if leading_external_modes in {5, 6}:
            values = values[leading_external_modes:]
        if len(values) == len(frequencies) and np.isfinite(values).all():
            payload["raman_activities"] = np.asarray(values) * atom_ureg.Unit("angstrom**4/amu")
            if matches := xtb_output_patterns.VIB_INTENSITY_METHOD.find_matches(
                file_content[: starts[-1].start()]
            ):
                payload["raman_activity_method"] = matches[-1].group("method")
    return Vibrations.model_validate(payload)


def extract_xtb_bond_orders(file_content: str, num_atoms: int) -> BondOrders | None:
    """Read thresholded WBO neighbors in legacy and modern stdout layouts.

    Basis: xTB v6.7.1 printout; see docs/en/reference/formats/xtbout.md.
    Reciprocal entries can differ in the last printed decimal. Keep the first
    observed value symmetrically, without averaging or filling unprinted zeros.
    """
    headers = xtb_output_patterns.WBO_HEADER.find_matches(file_content)
    if not headers:
        return None
    header = headers[-1]
    lines = file_content[header.end() :].splitlines()
    for index, line in enumerate(lines):
        if line.strip().startswith(("molecular ", "Topologies ", "Writing topology")):
            lines = lines[:index]
            break
    atom_count_source = "frame atoms"
    if num_atoms <= 0:
        # Single-point/hessian stdout can omit coordinates. Complete WBO row
        # numbering still defines the matrix axes; it does not supply geometry.
        indices = {
            int(row.group("atom"))
            for line in lines
            if (
                row := xtb_output_patterns.WBO_MODERN_ROW.match(line)
                or xtb_output_patterns.WBO_LEGACY_ROW.match(line)
            )
            is not None
        }
        if not indices or indices != set(range(1, max(indices) + 1)):
            return None
        num_atoms = max(indices)
        atom_count_source = "source WBO row indices"
    matrix = np.full((num_atoms, num_atoms), np.nan)
    atom: int | None = None
    modern = False
    for line in lines:
        if line.strip().startswith(("molecular ", "Topologies ", "Writing topology")):
            break
        if row := xtb_output_patterns.WBO_MODERN_ROW.match(line):
            modern = True
        elif not modern:
            row = xtb_output_patterns.WBO_LEGACY_ROW.match(line)
        if row is not None:
            atom = int(row.group("atom")) - 1
            pairs = row.group("pairs")
        elif modern and line[:1].isspace():
            pairs = line
        else:
            continue
        if atom is None:
            continue
        pattern = (
            xtb_output_patterns.WBO_MODERN_PAIR if modern else xtb_output_patterns.WBO_LEGACY_PAIR
        )
        for pair in pattern.find_matches(pairs):
            neighbor = int(pair.group("atom")) - 1
            if not (0 <= atom < num_atoms and 0 <= neighbor < num_atoms) or atom == neighbor:
                return None
            if np.isnan(matrix[atom, neighbor]):
                matrix[atom, neighbor] = matrix[neighbor, atom] = _to_float(pair.group("value"))
    if not np.isfinite(matrix).any():
        return None
    return BondOrders(
        bond_orders={
            "wiberg_bond_order": BondOrderMatrix(
                scheme="wiberg",
                values=matrix,
                source_label="Wiberg/Mayer (AO) data: Wiberg bond orders",
                print_threshold=_to_float(header.group("threshold")),
                missing_reason="pair not printed in thresholded source table",
                metadata={
                    "source_program": "xTB",
                    "basis": "AO",
                    "reciprocal_value_policy": "first_printed",
                    "atom_count_source": atom_count_source,
                },
            )
        }
    )


def extract_xtb_thermal_information(file_content: str) -> ThermalInformations | None:
    payload: dict[str, Any] = {}
    patterns = {
        "ZPVE": xtb_output_patterns.ZPVE,
        "H_T": xtb_output_patterns.TOTAL_ENTHALPY,
        "G_T": xtb_output_patterns.TOTAL_FREE_ENERGY,
    }
    for field, pattern in patterns.items():
        if value := _last_group(pattern, file_content):
            payload[field] = _to_float(value) * atom_ureg.Unit("hartree / particle")
    tot_lines = xtb_output_patterns.THERMO_TOTAL_ROW.find_named_group_values(file_content, "values")
    if tot_lines:
        values = _float_tokens(tot_lines[-1])
        if len(values) >= 3:
            payload["C_V"] = _to_float(values[1]) * atom_ureg.Unit("cal/mol/K")
            payload["S"] = _to_float(values[2]) * atom_ureg.Unit("cal/mol/K")
    if payload and (value := _last_group(xtb_output_patterns.MOLECULAR_MASS, file_content)):
        payload["molecular_mass"] = _to_float(value) * atom_ureg.amu
    if payload:
        payload.update(source_program="xTB", temperature=extract_xtb_temperature(file_content))
        for matched in xtb_output_patterns.THERMO_SETTING.find_matches(file_content):
            setting_value = _to_float(matched.group("value"))
            label = matched.group("label")
            if label == "scaling factor":
                payload["frequency_scale_factor"] = setting_value
            elif label == "rotor cutoff":
                payload["rotor_reference_frequency"] = setting_value * atom_ureg.cm_1
            elif label == "imag. cutoff":
                payload["imaginary_frequency_cutoff"] = setting_value * atom_ureg.cm_1
        if "T·S(HO)" in file_content and "T·S(FR)" in file_content:
            payload["entropy_model"] = "QRRHO"
    return ThermalInformations.model_validate(payload) if payload else None


def extract_xtb_rotation_constants(file_content: str):
    matches = xtb_output_patterns.ROTATION_CONSTANTS.find_matches(file_content)
    if not matches:
        return None
    values = np.asarray([_to_float(matches[-1].group(axis)) for axis in ("a", "b", "c")])
    return values * 29.9792458 * atom_ureg.gigahertz


def extract_xtb_populations(
    file_content: str, atoms: list[int] | None = None
) -> ChargeSpinPopulations | None:
    """Read indexed charge tables, allowing symmetry notices interleaved in stdout.

    Basis: xTB charge tables and the locally calculated PTB water fixture;
    see docs/en/reference/formats/xtbout.md. Incomplete or misidentified atom
    rows cannot be assigned to the frame's source atom order.
    """
    populations: dict[str, Any] = {}
    table_indices: list[int] = []
    table_atoms: list[int] = []
    gfn1_headers = xtb_output_patterns.GFN1_CHARGE_HEADER.find_matches(file_content)
    if gfn1_headers:
        mulliken: list[float] = []
        cm5: list[float] = []
        for line in file_content[gfn1_headers[-1].end() :].splitlines():
            row_matches = xtb_output_patterns.GFN1_CHARGE_ROW.find_matches(line)
            if not row_matches:
                if mulliken:
                    break
                continue
            matched = row_matches[0]
            table_indices.append(int(matched.group("index")))
            table_atoms.append(Chem.GetPeriodicTable().GetAtomicNumber(matched.group("symbol")))
            mulliken.append(_to_float(matched.group("mulliken")))
            cm5.append(_to_float(matched.group("cm5")))
        if mulliken:
            populations = {
                "mulliken_charges": {
                    "scheme": "mulliken",
                    "quantity": "charge",
                    "values": mulliken,
                    "source_label": "Mulliken/CM5 charges",
                },
                "cm5_charges": {
                    "scheme": "cm5",
                    "quantity": "charge",
                    "values": cm5,
                    "source_label": "Mulliken/CM5 charges",
                },
            }

    gfn2_headers = xtb_output_patterns.GFN2_CHARGE_HEADER.find_matches(file_content)
    if gfn2_headers:
        mulliken = []
        table_indices = []
        table_atoms = []
        for line in file_content[gfn2_headers[-1].end() :].splitlines():
            row_matches = xtb_output_patterns.GFN2_CHARGE_ROW.find_matches(line)
            if not row_matches:
                if line.strip().startswith(
                    ("Molecule has the following symmetry elements:", "It seems to be the")
                ):
                    continue
                if mulliken:
                    break
                continue
            matched = row_matches[0]
            table_indices.append(int(matched.group("index")))
            table_atoms.append(int(matched.group("atomic_number")))
            if Chem.GetPeriodicTable().GetElementSymbol(table_atoms[-1]) != matched.group("symbol"):
                return None
            mulliken.append(_to_float(matched.group("charge")))
        if mulliken:
            populations = {
                "mulliken_charges": {
                    "scheme": "mulliken",
                    "quantity": "charge",
                    "values": mulliken,
                    "source_label": "GFN2 atom-resolved q column",
                }
            }
    if table_indices != list(range(1, len(table_indices) + 1)):
        return None
    if atoms and table_atoms != atoms:
        return None
    return (
        ChargeSpinPopulations.model_validate({"populations": populations}) if populations else None
    )


def _shared_xtb_orbitals(
    file_content: str, energies: Any, occupancies: list[float | None], **frontiers: Any
) -> MolecularOrbitals:
    multiplicity = extract_xtb_multiplicity(file_content, extract_xtb_program_call(file_content))
    high_spin = (
        multiplicity > 1
        and all(value in (0, 1, 2) for value in occupancies)
        and sum(value == 1 for value in occupancies) == multiplicity - 1
    )
    return MolecularOrbitals.from_spatial_orbitals(
        energies, occupancies, closed_shell=multiplicity == 1, high_spin=high_spin, **frontiers
    )


def _extract_modern_orbitals(file_content: str) -> MolecularOrbitals | None:
    starts = xtb_output_patterns.MODERN_ORBITAL_HEADER.find_matches(file_content)
    if not starts:
        return None
    block = file_content[starts[-1].end() :]
    ends = xtb_output_patterns.HL_GAP.find_matches(block)
    if ends:
        block = block[: ends[0].start()]

    energies: list[float] = []
    occupancies: list[float | None] = []
    frontiers: dict[str, int] = {}
    for line in block.splitlines():
        if "..." in line:
            continue
        tokens = line.replace("(HOMO)", "").replace("(LUMO)", "").split()
        if not tokens or not tokens[0].isdigit():
            continue
        numeric = [token for token in tokens[1:] if _is_float_token(token)]
        if len(numeric) < 2:
            continue
        index = int(tokens[0])
        occupancy = _to_float(numeric[0]) if len(numeric) >= 3 else 0.0
        energy = _to_float(numeric[-2])
        while len(energies) < index - 1:
            energies.append(np.nan)
            occupancies.append(None)
        if len(energies) != index - 1:
            return None
        energies.append(energy)
        occupancies.append(occupancy)
        for marker, field in (("(HOMO)", "reported_homo_index"), ("(LUMO)", "reported_lumo_index")):
            if marker in line:
                frontiers[field] = index - 1
    if not energies:
        return None
    energy_array = np.asarray(energies) * atom_ureg.hartree
    return _shared_xtb_orbitals(file_content, energy_array, occupancies, **frontiers)


def _extract_legacy_orbitals(file_content: str) -> MolecularOrbitals | None:
    starts = xtb_output_patterns.LEGACY_ORBITAL_HEADER.find_matches(file_content)
    if not starts:
        return None
    block = file_content[starts[-1].end() :]
    ends = xtb_output_patterns.SCC_ENERGY_HEADER.find_matches(block)
    if ends:
        block = block[: ends[0].start()]
    pairs = xtb_output_patterns.LEGACY_ORBITAL_ROWS.find_matches(block)
    if not pairs:
        return None
    occupancies: list[float | None] = []
    energies_ev: list[float] = []
    for matched in pairs:
        occupancies.extend(
            _to_float(token) for token in _float_tokens(matched.group("occupancies"))
        )
        energies_ev.extend(_to_float(token) for token in _float_tokens(matched.group("energies")))
    if not energies_ev or len(occupancies) != len(energies_ev):
        return None
    energy_array = (np.asarray(energies_ev) * atom_ureg.electron_volt).to(atom_ureg.hartree)
    return _shared_xtb_orbitals(file_content, energy_array, occupancies)


def extract_xtb_orbitals(file_content: str) -> MolecularOrbitals | None:
    return _extract_modern_orbitals(file_content) or _extract_legacy_orbitals(file_content)


def extract_xtb_dipole(file_content: str) -> Polarizability | None:
    matches = xtb_output_patterns.DIPOLE.find_matches(file_content)
    if not matches:
        return None
    return Polarizability(
        dipole=np.asarray([_to_float(matches[-1].group(axis)) for axis in ("x", "y", "z")])
        * atom_ureg.debye
    )


def extract_xtb_single_point_properties(file_content: str) -> SinglePointProperties | None:
    payload: dict[str, Any] = {}
    scalar_patterns = {
        "vip": xtb_output_patterns.VIP,
        "vea": xtb_output_patterns.VEA,
        "gei": xtb_output_patterns.GEI,
    }
    for field, pattern in scalar_patterns.items():
        if value := _last_group(pattern, file_content):
            payload[field] = _to_float(value) * atom_ureg.Unit("eV / particle")

    fukui_headers = xtb_output_patterns.FUKUI_HEADER.find_matches(file_content)
    if fukui_headers:
        positive: list[float] = []
        negative: list[float] = []
        zero: list[float] = []
        for line in file_content[fukui_headers[-1].end() :].splitlines():
            row_matches = xtb_output_patterns.FUKUI_ROW.find_matches(line)
            if not row_matches:
                if positive:
                    break
                continue
            matched = row_matches[0]
            positive.append(_to_float(matched.group("positive")))
            negative.append(_to_float(matched.group("negative")))
            zero.append(_to_float(matched.group("zero")))
        if positive:
            payload.update(
                fukui_positive=positive,
                fukui_negative=negative,
                fukui_zero=zero,
            )
    return SinglePointProperties.model_validate(payload) if payload else None


__all__ = [
    "ensure_xtb_output_content",
    "extract_xtb_bond_orders",
    "extract_xtb_coords",
    "extract_xtb_dipole",
    "extract_xtb_final_energy",
    "extract_xtb_major_version",
    "extract_xtb_metadata",
    "extract_xtb_optimization_status",
    "extract_xtb_orbitals",
    "extract_xtb_populations",
    "extract_xtb_rotation_constants",
    "extract_xtb_running_time",
    "extract_xtb_single_point_properties",
    "extract_xtb_status",
    "extract_xtb_thermal_information",
    "extract_xtb_version",
    "extract_xtb_vibrations",
    "iter_xtb_version_matches",
]
