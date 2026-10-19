"""ORCA result extraction grounded in manual semantics and output fixtures.

Source IDs in docstrings resolve in
``docs/en/reference/formats/orcaout.md`` (Chinese peer under
``docs/zh``). That format page includes official manuals, property mappings and derivations. Fixture paths are relative to
``tests/test_files/orca/output_files``.
"""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
from rdkit import Chem

from molop.io.base_models.DataClasses import (
    NMR,
    BondOrderMatrix,
    BondOrders,
    ChargeSpinPopulations,
    ElectronicState,
    ElectronicStates,
    Energies,
    EnergyObservation,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    MolecularOrbitals,
    Polarizability,
    ShieldingTensor,
    ThermalInformations,
    TotalSpin,
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


def _parse_labeled_charge_block(
    text: str, header: str, atoms: list[int] | None = None
) -> tuple[list[float], list[float]]:
    start = text.rfind(header)
    if start < 0:
        return [], []
    charges: list[float] = []
    spins: list[float] = []
    for line in text[start + len(header) :].splitlines():
        matched = orca_log_patterns.MULLIKEN_CHARGE_ROW.match(line)
        if matched is not None:
            index = int(matched.group("index"))
            if index != len(charges) or (
                atoms is not None
                and (
                    index >= len(atoms)
                    or Chem.GetPeriodicTable().GetElementSymbol(atoms[index])
                    != matched.group("symbol")
                )
            ):
                return [], []
            charges.append(_as_float(matched.group("value")))
            if spin := matched.group("spin"):
                spins.append(_as_float(spin))
        elif charges:
            break
    if atoms is not None and len(charges) != len(atoms):
        return [], []
    return charges, spins


def _parse_hirshfeld(text: str, atoms: list[int] | None = None) -> tuple[list[float], list[float]]:
    start = text.rfind("HIRSHFELD ANALYSIS")
    if start < 0:
        return [], []
    charges: list[float] = []
    spins: list[float] = []
    for line in text[start + len("HIRSHFELD ANALYSIS") :].splitlines():
        matched = orca_log_patterns.HIRSHFELD_ROW.match(line)
        if matched is not None:
            index = int(matched.group("index"))
            if index != len(charges) or (
                atoms is not None
                and (
                    index >= len(atoms)
                    or Chem.GetPeriodicTable().GetElementSymbol(atoms[index])
                    != matched.group("symbol")
                )
            ):
                return [], []
            charges.append(_as_float(matched.group("charge")))
            spins.append(_as_float(matched.group("spin")))
        elif charges:
            break
    if atoms is not None and len(charges) != len(atoms):
        return [], []
    return charges, spins


def extract_orca_coords(
    text: str,
    *,
    capture_source_evidence: bool = False,
) -> tuple[list[int] | None, Any | None, int | None]:
    """Read printed Angstrom coordinates in source atom order.

    Basis: ORCA-COORD; ``local/h2o_orca_v5_charges.out`` coordinate headers.
    Selecting the last table and retaining print precision are MolOP policies.
    """
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


def extract_orca_atomic_masses(
    text: str,
    *,
    expected_atom_count: int | None = None,
    atoms: list[int] | None = None,
) -> Any | None:
    """Extract the MASS column, not coordinate magnitudes, from the A.U. table.

    Basis: ORCA-COORD / ORCA-UNITS and the column labels in
    ``local/h2o_orca_v5_charges.out``. Printed masses are in amu, not bohr.
    """

    matches = orca_log_patterns.COORD_AU_HEADER.find_matches(text)
    if not matches:
        return None

    rows: list[tuple[int, float]] = []
    symbols: list[str] = []
    for line in text[matches[-1].end() :].splitlines():
        matched = orca_log_patterns.COORD_AU_MASS_ROW.match(line)
        if matched is None:
            if rows:
                break
            continue
        rows.append(
            (
                int(matched.group("index")),
                _as_float(matched.group("mass")),
            )
        )
        symbols.append(matched.group("symbol"))
    if not rows:
        return None
    indices = [index for index, _mass in rows]
    if indices != list(range(len(indices))):
        return None
    if atoms is not None and symbols != [
        Chem.GetPeriodicTable().GetElementSymbol(atom) for atom in atoms
    ]:
        return None
    if expected_atom_count is not None and len(rows) != expected_atom_count:
        return None
    return np.asarray([mass for _index, mass in rows], dtype=float) * atom_ureg.amu


def extract_orca_energies(
    text: str,
    *,
    capture_source_evidence: bool = False,
    model_chemistry: Any = None,
) -> Energies | None:
    """Map observed total/reference/correlation labels without conflating them.

    Basis: ORCA-ENERGY; ``cclib/basicORCA6.0/water_{mp2,mp3,ccsd,ccsd_t}.out``.
    The spelling and separator distinguish E(MP2) totals from corrections.
    """
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


def extract_orca_forces(
    text: str, num_atoms: int | None, atoms: list[int] | None = None
) -> Any | None:
    """Negate the Cartesian gradient: force = -dE/dR, in hartree/bohr.

    Basis: ORCA-ENERGY / ORCA-OPT; ``local/h2_grad_orca.out``,
    ``local/numerical_orca.out`` and the printed MP2 gradient grammar.
    """
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
    symbols: list[str] = []
    index_base = 0 if row_pattern is orca_log_patterns.MP2_GRADIENT_ROW else 1
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
        if int(matched.group("index")) != len(rows) + index_base:
            return None
        if row_pattern is orca_log_patterns.GRADIENT_ROW:
            symbols.append(matched.group("symbol"))
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
    if (
        symbols
        and atoms is not None
        and symbols != [Chem.GetPeriodicTable().GetElementSymbol(atom) for atom in atoms]
    ):
        return None
    return -np.asarray(rows, dtype=float) * atom_ureg.Unit("hartree / bohr")


def extract_orca_vibrations(text: str, num_atoms: int | None) -> Vibrations | None:
    """Align frequency, displacement and spectrum rows by source mode number.

    Basis: ORCA-VIB; ``local/H2O_hess_orca.out`` (legacy) and
    ``cclib/basicORCA6.0/dvb_raman.out`` (modern IR and Raman).
    Leading-zero trimming and complete-array requirements are MolOP policies.
    """
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
    payload: dict[str, Any] = {
        "frequencies": np.asarray(frequencies, dtype=float) * atom_ureg.cm_1,
        "mode_indices": list(range(len(frequencies))),
        "source_mode_indices": mode_indices,
        "source_program": "ORCA",
        "orientation": "source",
    }
    modes = _parse_orca_normal_modes(block, mode_indices, num_atoms)
    if modes:
        payload.update(
            vibration_modes=modes,
            normalization="source_program",
            mass_weighting="not_mass_weighted",
        )
    ir_header, ir_rows = _parse_orca_spectrum(block, "IR SPECTRUM")
    # ORCA-VIB documents the change after 4.2.1: modern Int is km/mol,
    # modern T**2 is a.u. Legacy H2O_hess_orca.out uses the T**2-labeled
    # third column. Detect the actual header, rather than assuming a version.
    intensity_column = 1 if "km/mol" in ir_header else 0
    if all(index in ir_rows and len(ir_rows[index]) > intensity_column for index in mode_indices):
        payload["IR_intensities"] = (
            np.asarray([ir_rows[index][intensity_column] for index in mode_indices])
            * atom_ureg.km
            / atom_ureg.mol
        )
    _, raman_rows = _parse_orca_spectrum(block, "RAMAN SPECTRUM")
    if all(index in raman_rows and len(raman_rows[index]) == 2 for index in mode_indices):
        payload["raman_activities"] = (
            np.asarray([raman_rows[index][0] for index in mode_indices])
            * atom_ureg.angstrom**4
            / atom_ureg.amu
        )
        payload["raman_depolarization_ratios"] = (
            np.asarray([raman_rows[index][1] for index in mode_indices]) * atom_ureg.dimensionless
        )
    return Vibrations.model_validate(payload)


def _parse_orca_normal_modes(
    text: str, mode_indices: list[int], num_atoms: int | None
) -> list[Any]:
    """Preserve printed Cartesian vectors without another mass transform.

    Basis: NORMAL MODES text explicitly states M(i,i)=1/sqrt(m[i]) in both
    ORCA-VIB fixtures. The common model's displacement unit is a convention
    for these normalized vectors, not an observed absolute amplitude.
    """
    start = text.find("NORMAL MODES")
    if start < 0 or num_atoms is None:
        return []
    columns: list[int] = []
    pending_rows = False
    displacements: dict[int, dict[int, float]] = {}
    for line in text[start:].splitlines()[1:]:
        if matched := orca_log_patterns.MODE_COLUMNS.match(line):
            columns = [int(value) for value in matched.group("indices").split()]
            if len(columns) != len(set(columns)):
                return []
            pending_rows = True
            continue
        tokens = line.split()
        if (
            pending_rows
            and len(tokens) == len(columns)
            and all(token[0].isalpha() for token in tokens)
        ):
            # Symmetry-enabled ORCA 6 tables put irreps below the indices.
            continue
        if columns and (matched := orca_log_patterns.MODE_ROW.match(line)):
            pending_rows = False
            values = [_as_float(value) for value in matched.group("values").split()]
            if len(values) != len(columns):
                return []
            row = int(matched.group("index"))
            for index, value in zip(columns, values, strict=True):
                displacements.setdefault(index, {})[row] = value
        elif displacements and line.strip() and set(line.strip()) != {"-"}:
            break
    expected_rows = set(range(3 * num_atoms))
    if any(set(displacements.get(index, {})) != expected_rows for index in mode_indices):
        return []
    # Preserve ORCA's normalized Cartesian displacement vectors; do not
    # apply another mass factor or independently normalize each atom.
    return [
        np.asarray([displacements[index][row] for row in range(3 * num_atoms)]).reshape(
            num_atoms, 3
        )
        * atom_ureg.angstrom
        for index in mode_indices
    ]


def _parse_orca_spectrum(text: str, header: str) -> tuple[str, dict[int, list[float]]]:
    start = text.rfind(header)
    if start < 0:
        return "", {}
    header_lines: list[str] = []
    rows: dict[int, list[float]] = {}
    block = text[start + len(header) :]
    if sections := orca_log_patterns.SECTION_HEADER.find_matches(block):
        block = block[: sections[0].start()]
    for line in block.splitlines():
        if matched := orca_log_patterns.SPECTRUM_ROW.match(line):
            rows[int(matched.group("index"))] = [
                _as_float(value) for value in matched.group("values").split()
            ]
        elif rows:
            break
        else:
            header_lines.append(line)
    return "\n".join(header_lines), rows


def extract_orca_orbitals(text: str) -> MolecularOrbitals | None:
    """Read the last printed orbital table, retaining fractional occupations.

    Basis: ORCA-OUTPUT / ORCA-METHOD; ORBITAL ENERGIES in the water, ROCIS
    and ORCA 6.1 coupling fixtures. ``dvb_rocis.out`` explicitly limits virtual
    printing. Spin splitting follows the printed reference type and high-spin
    evidence; NaN/None padding is a MolOP representation policy.
    """
    start = text.rfind("ORBITAL ENERGIES")
    if start < 0:
        return None
    channels: dict[str, list[tuple[float, float | None, str | None]]] = {"alpha": [], "beta": []}
    channel = "alpha"
    for line in text[start:].splitlines()[1:]:
        if "SPIN DOWN ORBITALS" in line:
            channel = "beta"
        elif matched := orca_log_patterns.ORBITAL_ROW.match(line):
            if int(matched.group("index")) != len(channels[channel]):
                return None
            channels[channel].append(
                (
                    _as_float(matched.group("energy")),
                    _as_float(matched.group("occupancy")),
                    matched.group("symmetry"),
                )
            )
        elif channels[channel] and line.strip():
            break
    if not channels["alpha"]:
        return None
    # A single restricted table contains total, rather than per-spin, OCC.
    if not channels["beta"] and "SPIN UP ORBITALS" not in text[start:]:
        return _split_orca_restricted_orbitals(channels["alpha"], text[:start])
    payload: dict[str, Any] = {"orbital_basis": "spin_resolved", "occupancy_origin": "observed"}
    orbital_count = max(len(rows) for rows in channels.values())
    has_symmetries = any(row[2] is not None for rows in channels.values() for row in rows)
    for spin, rows in channels.items():
        missing_count = orbital_count - len(rows)
        # A limited virtual-orbital printout can leave one spin shorter.
        # Align source indices without inventing energies or occupations.
        payload[f"{spin}_energies"] = (
            np.asarray([row[0] for row in rows] + [np.nan] * missing_count) * atom_ureg.hartree
        )
        payload[f"{spin}_occupancies"] = [row[1] for row in rows] + [None] * missing_count
        if has_symmetries:
            payload[f"{spin}_symmetries"] = [row[2] for row in rows] + [None] * missing_count
    return MolecularOrbitals.model_validate(payload)


def _split_orca_restricted_orbitals(
    rows: list[tuple[float, float | None, str | None]], text: str
) -> MolecularOrbitals:
    """ORCA-METHOD: distinguish half-electron RHF from high-spin ROHF/ROKS.

    General ROHF can also be configuration/spin averaged or a coupled CSF.
    Without high-spin evidence, open-shell spin occupations remain unknown.
    """
    reference = ""
    if matches := orca_log_patterns.HF_TYPE.find_matches(text):
        reference = matches[-1].group("kind").upper()
    case = ""
    if matches := orca_log_patterns.ROHF_CASE.find_matches(text):
        case = matches[-1].group("case").upper()
    multiplicity = None
    if matches := orca_log_patterns.REFERENCE_MULTIPLICITY.find_matches(text):
        multiplicity = int(matches[-1].group("value"))
    singly_occupied = sum(occupancy == 1 for _, occupancy, _ in rows)
    integer_occupations = all(occupancy in (0, 1, 2) for _, occupancy, _ in rows)
    open_shell = reference in {"ROHF", "ROKS"}
    high_spin = (
        open_shell
        and integer_occupations
        and (
            case == "HIGHSPIN"
            or (not case and (reference == "ROKS" or multiplicity == singly_occupied + 1))
        )
    )
    return MolecularOrbitals.from_spatial_orbitals(
        np.asarray([energy for energy, _, _ in rows]) * atom_ureg.hartree,
        [occupancy for _, occupancy, _ in rows],
        closed_shell=reference in {"RHF", "RKS"},
        high_spin=high_spin,
        symmetries=[symmetry for _, _, symmetry in rows] if any(row[2] for row in rows) else None,
    )


def extract_orca_bond_orders(text: str, num_atoms: int | None) -> BondOrders | None:
    """Read printed Mayer atom-pair entries into a symmetric matrix.

    Basis: ORCA-POP, MAYER_BONDORDERTHRESH; ``local/H2O_hess_orca.out``.
    Unprinted pairs remain NaN; the source print threshold is retained.
    """
    start = text.rfind("Mayer bond orders larger than")
    if start < 0 or num_atoms is None:
        return None
    matrix = np.full((num_atoms, num_atoms), np.nan)
    found = False
    for line in text[start:].splitlines()[1:]:
        matches = orca_log_patterns.MAYER_BOND.find_matches(line)
        if not matches:
            if found or line.strip():
                break
            continue
        for matched in matches:
            first, second = int(matched.group("first")), int(matched.group("second"))
            if first >= num_atoms or second >= num_atoms or first == second:
                return None
            matrix[first, second] = matrix[second, first] = _as_float(matched.group("value"))
            found = True
    if not found:
        return None
    thresholds = orca_log_patterns.MAYER_THRESHOLD.find_matches(text[start:])
    return BondOrders(
        bond_orders={
            "mayer_bond_order": BondOrderMatrix(
                scheme="mayer",
                values=matrix,
                source_label=text[start:].splitlines()[0].strip(),
                print_threshold=_as_float(thresholds[0].group("value")) if thresholds else None,
                missing_reason="not printed in thresholded source table",
            ),
        },
    )


def extract_orca_total_spin(text: str) -> TotalSpin | None:
    """Retain observed <S**2> and derive an effective spin quantum number.

    Basis: ORCA-SPIN; ``cclib/basicORCA6.0/dvb_rocis.out`` contamination block.
    sqrt(<S**2> + 1/4) - 1/2 is derived, not the separately printed ideal spin.
    """
    matches = orca_log_patterns.SPIN_SQUARE.find_matches(text)
    if not matches:
        return None
    spin_square = _as_float(matches[-1].group("value"))
    return TotalSpin(
        spin_square=spin_square,
        spin_quantum_number=np.sqrt(max(spin_square, 0.0) + 0.25) - 0.5,
    )


def extract_orca_thermochemistry(text: str) -> dict[str, Any]:
    """Project the last temperature block into the shared thermal model.

    Basis: ORCA-THERMO / ORCA-UNITS; ``local/H2O_hess_orca.out``.
    Entropy rows hold T*S, so divide by T. TCE/TCH/TCG are total corrections
    relative to E(el), not ORCA's individual incremental correction rows.
    """
    start = text.rfind("THERMOCHEMISTRY AT")
    if start < 0:
        return {}
    block = text[start:]
    units = {"Eh": "hartree/particle", "AMU": "amu", "K": "K", "atm": "atm"}
    rows = {
        matched.group("label").strip(): atom_ureg.Quantity(
            _as_float(matched.group("value")), units[matched.group("unit")]
        )
        for matched in orca_log_patterns.THERMAL_ROW.find_matches(block)
    }
    payload: dict[str, Any] = {}
    for label, field in (
        ("Zero point energy", "ZPVE"),
        ("Non-thermal (ZPE) correction", "ZPVE"),
        ("Total correction", "TCE"),
        ("Total thermal energy", "U_T"),
        ("Total Enthalpy", "H_T"),
        ("Final Gibbs free energy", "G_T"),
        ("Total Mass", "molecular_mass"),
    ):
        if label in rows:
            payload[field] = rows[label]
    if "Electronic energy" in rows:
        electronic = rows["Electronic energy"]
        if "ZPVE" in payload:
            payload["U_0"] = electronic + payload["ZPVE"]
        for total, correction in (("U_T", "TCE"), ("H_T", "TCH"), ("G_T", "TCG")):
            if total in payload:
                payload[correction] = payload[total] - electronic
    if matches := orca_log_patterns.GIBBS_CORRECTION.find_matches(block):
        payload["TCG"] = _as_float(matches[-1].group("value")) * atom_ureg.Unit("hartree/particle")
    temperature = rows.get("Temperature")
    if temperature is not None and temperature.magnitude > 0 and "Final entropy term" in rows:
        payload["S"] = rows["Final entropy term"] / temperature
    if matches := orca_log_patterns.ROTATIONAL_CONSTANTS.find_matches(block):
        payload["rotational_constants"] = (
            np.asarray([_as_float(matches[-1].group(axis)) for axis in ("x", "y", "z")])
            * atom_ureg.cm_1
            * atom_ureg.speed_of_light
        ).to("GHz")
    if matches := orca_log_patterns.SYMMETRY_NUMBER.find_matches(block):
        payload["rotational_symmetry_number"] = int(matches[-1].group("value"))
    result: dict[str, Any] = {}
    if payload:
        payload.update(
            source_program="ORCA", temperature=temperature, pressure=rows.get("Pressure")
        )
        if "Vibrational entropy computed according to the QRRHO" in block:
            payload["entropy_model"] = "QRRHO"
        elif "Vibrational entropy computed according to the RRHO" in block:
            payload["entropy_model"] = "RRHO"
        if matches := orca_log_patterns.FREQUENCY_SCALE.find_matches(text[:start]):
            payload["frequency_scale_factor"] = _as_float(matches[-1].group("value"))
        if matches := orca_log_patterns.THERMAL_CUTOFF.find_matches(block):
            payload["low_frequency_cutoff"] = _as_float(matches[-1].group("value")) * atom_ureg.cm_1
        result["thermal_informations"] = ThermalInformations.model_validate(payload)
    for label, field in (("Temperature", "temperature"), ("Pressure", "pressure")):
        if label in rows:
            result[field] = rows[label]
    return result


def extract_orca_nmr(text: str, num_atoms: int | None) -> NMR | None:
    """Read source-order shielding tensors and isotropic J couplings.

    Basis: ORCA-NMR; ``cclib/basicORCA6.0/dvb_nmr.out`` and
    ``cclib/basicORCA6.1/dvb_coupling_nmr.out``. Use zero-based output indices;
    do not import indexing conventions from unrelated eprnmr input options.
    """
    nuclei = orca_log_patterns.NMR_NUCLEUS.find_matches(text)
    shieldings: dict[int, ShieldingTensor] = {}
    summary_start = text.rfind("CHEMICAL SHIELDING SUMMARY (ppm)")
    summary = {}
    for line in text[summary_start:].splitlines() if summary_start >= 0 else []:
        if matched := orca_log_patterns.SHIELDING_SUMMARY_ROW.match(line):
            summary[int(matched.group("index"))] = matched
        elif summary:
            break
    for position, nucleus in enumerate(nuclei):
        index = int(nucleus.group("index"))
        if num_atoms is not None and index >= num_atoms:
            continue
        end = nuclei[position + 1].start() if position + 1 < len(nuclei) else len(text)
        block = text[nucleus.end() : end]
        start = block.find("Total shielding tensor (ppm):")
        if start < 0:
            continue
        tensor = []
        for line in block[start:].splitlines()[1:]:
            matched = orca_log_patterns.TENSOR_ROW.match(line)
            if matched is None:
                break
            tensor.append([_as_float(matched.group(axis)) for axis in ("x", "y", "z")])
            if len(tensor) == 3:
                break
        if len(tensor) != 3:
            continue
        payload: dict[str, Any] = {
            "atom_index": index,
            "atom_symbol": nucleus.group("symbol"),
            "shielding_tensor": np.asarray(tensor) * atom_ureg.ppm,
            "orientation": "source",
        }
        if values := orca_log_patterns.SHIELDING_PRINCIPAL_VALUES.find_matches(block):
            payload["principal_values"] = (
                np.asarray([_as_float(values[0].group(axis)) for axis in ("x", "y", "z")])
                * atom_ureg.ppm
            )
            payload["isotropic"] = _as_float(values[0].group("isotropic")) * atom_ureg.ppm
        if index in summary:
            for field in ("isotropic", "anisotropy"):
                payload[field] = _as_float(summary[index].group(field)) * atom_ureg.ppm
            payload["anisotropy_convention"] = "ORCA printed shielding anisotropy"
        shieldings[index] = ShieldingTensor.model_validate(payload)
    payload = (
        {"shielding_tensors": [shieldings[index] for index in sorted(shieldings)]}
        if shieldings
        else {}
    )
    payload.update(_parse_orca_couplings(text, num_atoms))
    if not payload:
        return None
    if matches := orca_log_patterns.NMR_GAUGE.find_matches(text):
        payload["gauge"] = matches[-1].group("gauge")
    return NMR.model_validate(payload)


def _parse_orca_couplings(text: str, num_atoms: int | None) -> dict[str, Any]:
    """Read J[...] iso= rows and the blocked isotropic J summary in Hz.

    Basis: ORCA-NMR and the ORCA 6.1 coupling fixture. FC/SD/PSO/DSO map to
    shared fields; full tensors and SD/FC cross terms remain unrepresented.
    """
    observed: dict[str, dict[tuple[int, int], float]] = {}
    for matched in orca_log_patterns.COUPLING_ISOTROPIC.find_matches(text):
        first, second = int(matched.group("first")), int(matched.group("second"))
        if num_atoms is not None and max(first, second) >= num_atoms:
            continue
        observed.setdefault(matched.group("component"), {})[first, second] = _as_float(
            matched.group("value")
        )
    start = text.rfind("SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz)")
    columns: list[int] = []
    for line in text[start:].splitlines()[1:] if start >= 0 else []:
        if matched := orca_log_patterns.COUPLING_ROW.match(line):
            values = [_as_float(value) for value in matched.group("values").split()]
            row = int(matched.group("index"))
            if len(values) == len(columns):
                for column, value in zip(columns, values, strict=True):
                    if num_atoms is None or max(row, column) < num_atoms:
                        observed.setdefault("Total", {})[row, column] = value
        elif matches := orca_log_patterns.COUPLING_COLUMNS.find_matches(line):
            columns = [int(matched.group("index")) for matched in matches]
        elif columns and line.strip() and set(line.strip()) != {"-"}:
            break
    indices = sorted({index for rows in observed.values() for pair in rows for index in pair})
    if not indices:
        return {}
    positions = {index: position for position, index in enumerate(indices)}
    matrices = {}
    for component, rows in observed.items():
        matrix = np.full((len(indices), len(indices)), np.nan)
        for (first, second), value in rows.items():
            matrix[positions[first], positions[second]] = value
            matrix[positions[second], positions[first]] = value
        matrices[component] = matrix * atom_ureg.Hz
    payload: dict[str, Any] = {"coupling_atom_indices": indices}
    if "Total" in matrices:
        payload["spin_spin_coupling_j"] = matrices.pop("Total")
    if matrices:
        payload["spin_spin_coupling_j_components"] = matrices
    return payload


def extract_orca_populations(
    text: str, atoms: list[int] | None = None
) -> ChargeSpinPopulations | None:
    """Read atomic charge and spin columns, keeping their scheme and labels.

    Basis: ORCA-POP; ``local/h2o_orca_v5_charges.out`` (closed shell) and
    ``cclib/basicORCA6.0/dvb_rocis.out`` (open shell).
    """
    populations: dict[str, Any] = {}

    def add_series(
        key: str,
        values: list[float],
        *,
        scheme: str,
        quantity: str,
        source_label: str,
        spin_channel: str | None = None,
    ) -> None:
        if values:
            populations[key] = {
                "scheme": scheme,
                "quantity": quantity,
                "values": values,
                "spin_channel": spin_channel,
                "source_label": source_label,
            }

    mulliken_charges, mulliken_spins = _parse_labeled_charge_block(
        text, "MULLIKEN ATOMIC CHARGES", atoms
    )
    add_series(
        "mulliken_charges",
        mulliken_charges,
        scheme="mulliken",
        quantity="charge",
        source_label="MULLIKEN ATOMIC CHARGES",
    )
    add_series(
        "mulliken_spins",
        mulliken_spins,
        scheme="mulliken",
        quantity="spin_density",
        spin_channel="total",
        source_label="MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS",
    )
    lowdin_charges, lowdin_spins = _parse_labeled_charge_block(
        text, "LOEWDIN ATOMIC CHARGES", atoms
    )
    add_series(
        "lowdin_charges",
        lowdin_charges,
        scheme="lowdin",
        quantity="charge",
        source_label="LOEWDIN ATOMIC CHARGES",
    )
    add_series(
        "lowdin_spins",
        lowdin_spins,
        scheme="lowdin",
        quantity="spin_density",
        spin_channel="total",
        source_label="LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS",
    )
    hirshfeld_charges, hirshfeld_spins = _parse_hirshfeld(text, atoms)
    add_series(
        "hirshfeld_charges",
        hirshfeld_charges,
        scheme="hirshfeld",
        quantity="charge",
        source_label="HIRSHFELD ANALYSIS",
    )
    add_series(
        "hirshfeld_spins",
        hirshfeld_spins,
        scheme="hirshfeld",
        quantity="spin_density",
        spin_channel="total",
        source_label="HIRSHFELD ANALYSIS",
    )
    return (
        ChargeSpinPopulations.model_validate({"populations": populations}) if populations else None
    )


def extract_orca_polarizability(text: str) -> Polarizability | None:
    """Read dipoles in e*bohr and polarizabilities in bohr cubed.

    Basis: ORCA-ELECTRIC / ORCA-UNITS; ORCA 5.0/6.0 ``Trp_polar.out`` fixtures.
    """
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
            * atom_ureg.elementary_charge
            * atom_ureg.bohr
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
    """Map printed energy/gradient/step criteria and convergence markers.

    Basis: ORCA-OPT; ``local/opt_orca.out`` and ORCA 6.0 ``dvb_gopt.out``.
    Absolute magnitudes and retaining source YES/NO are model conventions.
    """
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
    """Read explicit CPCM/SMD headings, solvent names and dielectric constants.

    Basis: ORCA-SOLVENT; ORCA 6.0 ``water_hf_solvent_{cpcm,smd}.log`` fixtures.
    """
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
    """Read STATE rows and a bounded ground-state electric absorption table.

    Basis: ORCA-EXCITED; TD/ADC/EOM/STEOM/ROCIS fixtures. CD, velocity,
    combined multipole and SOC-corrected spectra describe different observables.
    Left-right strengths take precedence over separate left/right CC strengths.
    """
    states: list[ElectronicState] = []
    tolerances: list[float] = []
    for matched in orca_log_patterns.STATE.find_matches(text):
        root = int(matched.group("root"))
        tail = matched.group("tail")
        multiplicity = None
        if mult_matches := orca_log_patterns.STATE_MULTIPLICITY.find_matches(tail):
            multiplicity = int(mult_matches[-1].group("multiplicity"))
        irrep = None
        if sym_matches := orca_log_patterns.STATE_IRREP.find_matches(tail):
            irrep = sym_matches[-1].group("irrep")
        states.append(
            ElectronicState(
                state_index=len(states),
                is_ground_state=False,
                root=root,
                label=f"STATE {root}",
                multiplicity=multiplicity,
                irrep=irrep,
                method=method,
                excitation_energy=_as_float(matched.group("ev")) * atom_ureg.eV,
                source="ORCA output",
            )
        )
        # STATE energies have fewer printed decimals than spectrum energies.
        tolerances.append(0.5 * 10 ** -_decimal_places(matched.group("ev")) + 1e-6)

    for row in _parse_orca_electric_absorption(text):
        energy = row["energy"]
        candidates = [
            state
            for state, tolerance in zip(states, tolerances, strict=True)
            if state.excitation_energy is not None
            and abs(state.excitation_energy.m_as("eV") - energy) <= tolerance
            and (row["multiplicity"] is None or state.multiplicity in (None, row["multiplicity"]))
            and (row["irrep"] is None or state.irrep in (None, row["irrep"]))
        ]
        # ORCA spectrum root numbers restart in each multiplicity/irrep block.
        # Use the root only to disambiguate already compatible states.
        same_root = [state for state in candidates if state.root == row["root"]]
        if len(same_root) == 1:
            state = same_root[0]
        elif len(candidates) == 1:
            state = candidates[0]
        elif candidates:
            continue
        else:
            state = ElectronicState(
                state_index=len(states),
                is_ground_state=False,
                root=row["root"],
                label=row["label"],
                multiplicity=row["multiplicity"],
                irrep=row["irrep"],
                method=method,
                excitation_energy=energy * atom_ureg.eV,
                source="ORCA absorption spectrum",
            )
            states.append(state)
            tolerances.append(1e-6)
        state.oscillator_strength = row["oscillator_strength"]
        state.transition_dipole = row["transition_dipole"]
        if state.multiplicity is None:
            state.multiplicity = row["multiplicity"]
        if state.irrep is None:
            state.irrep = row["irrep"]
        state.properties["absorption_root"] = row["root"]
        state.properties["absorption_label"] = row["label"]
        state.properties["absorption_convention"] = row["convention"]
        state.properties["initial_root"] = row["initial_root"]
        state.properties["initial_label"] = row["initial_label"]
    return ElectronicStates(states=states) if states else None


def _parse_orca_electric_absorption(text: str) -> list[dict[str, Any]]:
    """Read contiguous rows under an exact electric-dipole absorption heading.

    ORCA 6 transition labels include multiplicity/irrep; legacy state tables use
    wavenumbers without an eV column. Column units come from the table header.
    """
    tables: dict[str, list[dict[str, Any]]] = {}
    convention = "electric_dipole"
    lines = text.splitlines()
    for index, line in enumerate(lines):
        title = line.strip()
        if title.startswith("SPECTRUM FOR ") and title.endswith(" TRANSITION MOMENTS"):
            convention = (
                title.removeprefix("SPECTRUM FOR ").removesuffix(" TRANSITION MOMENTS").lower()
            )
        if title != "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS":
            continue
        rows: list[dict[str, Any]] = []
        header: list[str] = []
        for data_line in lines[index + 1 :]:
            header_text = "\n".join(header)
            modern = orca_log_patterns.ABS_TRANSITION.match(data_line)
            numeric_pattern = (
                orca_log_patterns.ROOT_TRANSITION
                if "(eV)" in header_text
                else orca_log_patterns.LEGACY_ROOT_TRANSITION
            )
            matched = modern or numeric_pattern.match(data_line)
            # Require units and oscillator-strength column, not just numbers.
            if matched is not None and "fosc" in header_text and "(nm)" in header_text:
                ev = matched.group("ev") if modern or "(eV)" in header_text else None
                energy = (
                    _as_float(ev)
                    if ev is not None
                    else atom_ureg.Quantity(_as_float(matched.group("cm")), "cm_1")
                    .to("eV", "spectroscopy")
                    .magnitude
                )
                multiplicity = int(matched.group("multiplicity")) if modern else None
                irrep = matched.group("irrep") if modern else None
                root = int(matched.group("root"))
                dipole = matched.group("dipole")
                vector = None
                if dipole and "(au)" in header_text:
                    vector = (
                        np.asarray([_as_float(value) for value in dipole.split()[1:]])
                        * atom_ureg.elementary_charge
                        * atom_ureg.bohr
                    )
                rows.append(
                    {
                        "root": root,
                        "initial_root": 0 if modern else None,
                        "initial_label": modern.group("initial_label") if modern else None,
                        "label": matched.group("label") if modern else f"ROOT {root}",
                        "multiplicity": multiplicity,
                        "irrep": irrep,
                        "energy": energy,
                        "oscillator_strength": _as_float(matched.group("fosc")),
                        "transition_dipole": vector,
                        "convention": convention,
                    }
                )
            elif (
                rows
                or (title in data_line and header)
                or data_line.strip()
                and not (
                    set(data_line.strip()) == {"-"}
                    or any(
                        token in data_line for token in ("Transition", "State", "Energy", "(nm)")
                    )
                )
            ):
                break
            else:
                header.append(data_line)
        if rows:
            tables[convention] = rows
    for preferred in ("left-right", "electric_dipole", "right", "left"):
        if preferred in tables:
            return tables[preferred]
    return []
