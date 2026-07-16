from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
from rdkit import Chem
from scipy.constants import alpha

from molop.io.base_models.DataClasses import (
    NMR,
    ChargeSpinPopulations,
    Energies,
    EnergyObservation,
    GeometryOptimizationStatus,
    ImplicitSolvation,
    MolecularOrbitals,
    Polarizability,
    QMTaskRequest,
    Status,
    ThermalInformations,
    TotalSpin,
    Vibrations,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_patterns import fchk_patterns
from molop.io.logic.gaussian.fchk.output.parsers._fchk_records import (
    FCHKRecord,
    parse_fchk_header,
    parse_fchk_records,
)
from molop.io.logic.gaussian.input.GaussianRoute import (
    build_gaussian_excited_state_requests,
    build_gaussian_model_chemistry,
    build_gaussian_task_requests,
)
from molop.io.logic.gaussian.input.GaussianRouteParsing import (
    parse_gaussian_route_semantic,
)
from molop.unit import atom_ureg
from molop.utils.functions import fill_symmetric_matrix


METADATA_LABELS = frozenset(
    {
        "Route",
        "Gaussian Version",
        "Charge",
        "Multiplicity",
        "Job Status",
        "Total Energy",
    }
)

FRAME_LABELS = frozenset(
    {
        "Atomic numbers",
        "Current cartesian coordinates",
        "Number of alpha electrons",
        "Number of beta electrons",
        "SCF Energy",
        "Total Energy",
        "MP2 Energy",
        "MP3 Energy",
        "MP4 Energy",
        "MP4D Energy",
        "MP4DQ Energy",
        "MP4SDQ Energy",
        "Cluster Energy",
        "Cluster Energy with triples",
        "Cartesian Gradient",
        "Cartesian Force Constants",
        "Alpha Orbital Energies",
        "Beta Orbital Energies",
        "Mulliken Charges",
        "Mulliken Spin Densities",
        "APT Charges",
        "Hirshfeld Charges",
        "Hirshfeld Spin Densities",
        "CM5 Charges",
        "ESP Charges",
        "NPA Charges",
        "NPA Spins",
        "NPA Alpha Spin Densities",
        "NPA Beta Spin Densities",
        "S**2",
        "Dipole Moment",
        "Polarizability",
        "Quadrupole Moment",
        "NMR shielding",
        "Isotropic Spin-Spin Couplings",
        "Number of Normal Modes",
        "Vib-E2",
        "Vib-Modes",
        "Thermal Energy",
        "Thermal Enthalpy",
        "Thermal Free Energy",
        "Job Status",
    }
)

_FCHK_SHIELDING_AU_TO_PPM = alpha**2 * 1_000_000.0


def _value(records: Mapping[str, FCHKRecord], label: str) -> Any | None:
    record = records.get(label)
    return record.value if record is not None else None


def _int_value(records: Mapping[str, FCHKRecord], label: str) -> int | None:
    value = _value(records, label)
    return int(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def _float_value(records: Mapping[str, FCHKRecord], label: str) -> float | None:
    value = _value(records, label)
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def _string_value(records: Mapping[str, FCHKRecord], label: str) -> str:
    value = _value(records, label)
    return value.strip() if isinstance(value, str) else ""


def _int_array(records: Mapping[str, FCHKRecord], label: str) -> list[int]:
    value = _value(records, label)
    if not isinstance(value, list):
        return []
    return [int(item) for item in value]


def _float_array(records: Mapping[str, FCHKRecord], label: str) -> list[float]:
    value = _value(records, label)
    if not isinstance(value, list):
        return []
    return [float(item) for item in value]


def _fallback_tasks(job_type: str) -> list[QMTaskRequest]:
    normalized = job_type.lower()
    tasks: list[QMTaskRequest] = []
    if "opt" in normalized:
        tasks.append(QMTaskRequest(task_type="opt", derivative_order=1))
    if "freq" in normalized:
        tasks.append(QMTaskRequest(task_type="freq", derivative_order=2))
    if not tasks:
        tasks.append(QMTaskRequest(task_type="sp", derivative_order=0))
    return tasks


def _enrich_route_model_chemistry_from_header(semantic_route: Any, method: str, basis: str) -> None:
    header_method = method.removesuffix("-FC")
    header_semantic = parse_gaussian_route_semantic(f"# {header_method}")
    target = semantic_route.model_chemistry
    source = header_semantic.model_chemistry
    if target.method_token is None:
        target.method_token = source.method_token
        target.method_family = source.method_family
        target.functional = source.functional
    if target.spin_qualifier is None and source.spin_qualifier is not None:
        target.method_token = source.method_token
        target.method_family = source.method_family
        target.functional = source.functional
        target.spin_qualifier = source.spin_qualifier
    if target.basis_set is None and basis:
        target.basis_set = basis


def extract_fchk_metadata(file_content: str) -> dict[str, Any]:
    header = parse_fchk_header(file_content)
    records = parse_fchk_records(file_content, wanted_labels=METADATA_LABELS)
    route = _string_value(records, "Route")
    semantic_route = parse_gaussian_route_semantic(route)
    _enrich_route_model_chemistry_from_header(semantic_route, header.method, header.basis_set)
    model_chemistry = build_gaussian_model_chemistry(
        semantic_route,
        keywords=route,
        legacy_method=header.method,
        legacy_basis_set=header.basis_set,
    )
    tasks = build_gaussian_task_requests(semantic_route) or _fallback_tasks(header.job_type)
    status = extract_fchk_status(records)
    metadata: dict[str, Any] = {
        "qm_software": "Gaussian",
        "qm_software_version": _string_value(records, "Gaussian Version"),
        "title_card": header.title,
        "keywords": route,
        "semantic_route": semantic_route,
        "model_chemistry": model_chemistry,
        "task_requests": tasks,
        "excited_state_requests": build_gaussian_excited_state_requests(semantic_route),
        "method": model_chemistry.method or model_chemistry.method_family or "",
        "functional": model_chemistry.functional or "",
        "basis_set": model_chemistry.basis_set or header.basis_set,
        "charge": _int_value(records, "Charge") or 0,
        "multiplicity": _int_value(records, "Multiplicity") or 1,
        "status": status,
    }
    if model_chemistry.solvation_model or model_chemistry.solvent:
        metadata["solvent"] = ImplicitSolvation(
            solvent_model=model_chemistry.solvation_model,
            solvent=model_chemistry.solvent,
        )
    return metadata


def extract_fchk_structure(
    records: Mapping[str, FCHKRecord],
) -> tuple[list[int], Any] | tuple[None, None]:
    atoms = _int_array(records, "Atomic numbers")
    coordinate_values = _float_array(records, "Current cartesian coordinates")
    if not atoms or len(coordinate_values) != len(atoms) * 3:
        return None, None
    coords = (np.asarray(coordinate_values).reshape(-1, 3) * atom_ureg.bohr).to(atom_ureg.angstrom)
    return atoms, coords


def extract_fchk_energies(
    records: Mapping[str, FCHKRecord],
    *,
    method: str | None,
    capture_source_evidence: bool,
) -> Energies | None:
    label_fields = (
        ("SCF Energy", "reference_energy"),
        ("MP2 Energy", "mp2_energy"),
        ("MP3 Energy", "mp3_energy"),
        ("MP4 Energy", "mp4_energy"),
        ("MP4D Energy", "mp4_energy"),
        ("MP4DQ Energy", "mp4_energy"),
        ("MP4SDQ Energy", "mp4_energy"),
        ("Cluster Energy", "ccsd_energy"),
        ("Cluster Energy with triples", "ccsd_t_energy"),
        ("Total Energy", "electronic_energy"),
    )
    payload: dict[str, Any] = {}
    observations: list[EnergyObservation] = []
    for label, field in label_fields:
        value = _float_value(records, label)
        if value is None:
            continue
        energy = value * atom_ureg.hartree
        payload[field] = energy
        if capture_source_evidence:
            observations.append(
                EnergyObservation(
                    method=method or "Gaussian",
                    quantity_semantics="total_energy" if label == "Total Energy" else "component",
                    value=energy,
                    source_label=label,
                )
            )
    if observations:
        payload["observations"] = observations
    return Energies.model_validate(payload) if payload else None


def extract_fchk_forces(records: Mapping[str, FCHKRecord], num_atoms: int) -> Any | None:
    values = _float_array(records, "Cartesian Gradient")
    if len(values) != num_atoms * 3:
        return None
    return -np.asarray(values).reshape(num_atoms, 3) * atom_ureg.hartree / atom_ureg.bohr


def extract_fchk_hessian(records: Mapping[str, FCHKRecord], num_atoms: int) -> Any | None:
    values = _float_array(records, "Cartesian Force Constants")
    size = num_atoms * 3
    if len(values) != size * (size + 1) // 2:
        return None
    return (
        fill_symmetric_matrix(np.asarray(values, dtype=float))
        * atom_ureg.hartree
        / atom_ureg.bohr**2
    )


def extract_fchk_orbitals(records: Mapping[str, FCHKRecord]) -> MolecularOrbitals | None:
    alpha = _float_array(records, "Alpha Orbital Energies")
    if not alpha:
        return None
    beta = _float_array(records, "Beta Orbital Energies") or list(alpha)
    alpha_electrons = _int_value(records, "Number of alpha electrons") or 0
    beta_electrons = _int_value(records, "Number of beta electrons") or 0
    return MolecularOrbitals(
        alpha_energies=np.asarray(alpha) * atom_ureg.hartree,
        beta_energies=np.asarray(beta) * atom_ureg.hartree,
        alpha_occupancies=[1.0 if index < alpha_electrons else 0.0 for index in range(len(alpha))],
        beta_occupancies=[1.0 if index < beta_electrons else 0.0 for index in range(len(beta))],
    )


def extract_fchk_populations(
    records: Mapping[str, FCHKRecord], num_atoms: int
) -> ChargeSpinPopulations | None:
    record_specs = {
        "Mulliken Charges": ("mulliken_charges", "mulliken", "charge", None),
        "Mulliken Spin Densities": (
            "mulliken_spins",
            "mulliken",
            "spin_density",
            "total",
        ),
        "APT Charges": ("apt_charges", "apt", "charge", None),
        "Hirshfeld Charges": ("hirshfeld_charges", "hirshfeld", "charge", None),
        "Hirshfeld Spin Densities": (
            "hirshfeld_spins",
            "hirshfeld",
            "spin_density",
            "total",
        ),
        "CM5 Charges": ("cm5_charges", "cm5", "charge", None),
        "ESP Charges": ("esp_charges", "esp", "charge", None),
        "NPA Charges": ("npa_charges", "npa", "charge", None),
        "NPA Spins": ("npa_spins", "npa", "spin_density", "total"),
        "NPA Alpha Spin Densities": (
            "npa_alpha_spin_densities",
            "npa",
            "spin_density",
            "alpha",
        ),
        "NPA Beta Spin Densities": (
            "npa_beta_spin_densities",
            "npa",
            "spin_density",
            "beta",
        ),
    }
    populations: dict[str, Any] = {}
    for source_label, (key, scheme, quantity, spin_channel) in record_specs.items():
        values = _float_array(records, source_label)
        if len(values) == num_atoms:
            populations[key] = {
                "scheme": scheme,
                "quantity": quantity,
                "values": values,
                "spin_channel": spin_channel,
                "source_label": source_label,
            }
    return (
        ChargeSpinPopulations.model_validate({"populations": populations}) if populations else None
    )


def extract_fchk_total_spin(records: Mapping[str, FCHKRecord]) -> TotalSpin | None:
    spin_square = _float_value(records, "S**2")
    if spin_square is None:
        return None
    return TotalSpin(
        spin_square=spin_square,
        spin_quantum_number=math.sqrt(max(spin_square, 0.0) + 0.25) - 0.5,
    )


def extract_fchk_polarizability(records: Mapping[str, FCHKRecord]) -> Polarizability | None:
    payload: dict[str, Any] = {}
    atomic_charge = atom_ureg.atomic_unit_of_current * atom_ureg.atomic_unit_of_time
    dipole = _float_array(records, "Dipole Moment")
    if len(dipole) == 3:
        payload["dipole"] = np.asarray(dipole) * atomic_charge * atom_ureg.bohr
    tensor = _float_array(records, "Polarizability")
    if len(tensor) == 6:
        payload["polarizability_tensor"] = np.asarray(tensor) * atom_ureg.bohr**3
    quadrupole = _float_array(records, "Quadrupole Moment")
    if len(quadrupole) == 6:
        payload["quadrupole"] = np.asarray(quadrupole) * atomic_charge * atom_ureg.bohr**2
    return Polarizability.model_validate(payload) if payload else None


def extract_fchk_nmr(
    records: Mapping[str, FCHKRecord],
    atoms: list[int],
    *,
    route: str,
) -> NMR | None:
    num_atoms = len(atoms)
    shielding_values = _float_array(records, "NMR shielding")
    if num_atoms <= 0 or len(shielding_values) != num_atoms * 9:
        return None

    route_matches = fchk_patterns.NMR_ROUTE.find_matches(route)
    route_options = route_matches[0].group("options").upper() if route_matches else ""
    gauge = next(
        (
            candidate
            for candidate in ("GIAO", "CSGT", "IGAIM", "SINGLEORIGIN")
            if candidate in route_options
        ),
        None,
    )
    periodic_table = Chem.GetPeriodicTable()
    tensors = (
        np.asarray(shielding_values, dtype=float).reshape(num_atoms, 3, 3).transpose(0, 2, 1)
        * _FCHK_SHIELDING_AU_TO_PPM
    )
    shielding_tensors: list[dict[str, Any]] = []
    for atom_index, (atomic_number, tensor) in enumerate(zip(atoms, tensors, strict=True)):
        principal_values = np.linalg.eigvalsh((tensor + tensor.T) / 2.0)
        shielding_tensors.append(
            {
                "atom_index": atom_index,
                "atom_symbol": periodic_table.GetElementSymbol(atomic_number),
                "shielding_tensor": tensor * atom_ureg.ppm,
                "isotropic": float(np.trace(tensor) / 3.0) * atom_ureg.ppm,
                "anisotropy": float(
                    principal_values[2] - (principal_values[0] + principal_values[1]) / 2.0
                )
                * atom_ureg.ppm,
                "principal_values": principal_values * atom_ureg.ppm,
                "anisotropy_convention": "Gaussian",
                "orientation": "unknown",
            }
        )
    payload: dict[str, Any] = {"gauge": gauge, "shielding_tensors": shielding_tensors}

    coupling_values = _float_array(records, "Isotropic Spin-Spin Couplings")
    packed_size = num_atoms * (num_atoms + 1) // 2
    if len(coupling_values) == packed_size * 4:
        contributions = np.asarray(coupling_values, dtype=float).reshape(4, packed_size)
        component_matrices = {
            component_name: fill_symmetric_matrix(contribution) * atom_ureg.Hz
            for component_name, contribution in zip(
                ("FC", "SD", "PSO", "DSO"), contributions, strict=True
            )
        }
        payload["coupling_atom_indices"] = list(range(num_atoms))
        payload["spin_spin_coupling_k_components"] = component_matrices
        payload["spin_spin_coupling_k"] = sum(
            component_matrices.values(), start=np.zeros((num_atoms, num_atoms)) * atom_ureg.Hz
        )
    return NMR.model_validate(payload)


def extract_fchk_vibrations(records: Mapping[str, FCHKRecord], num_atoms: int) -> Vibrations | None:
    num_modes = _int_value(records, "Number of Normal Modes") or 0
    vib_e2 = _float_array(records, "Vib-E2")
    if num_modes <= 0 or len(vib_e2) < num_modes:
        return None
    payload: dict[str, Any] = {
        "frequencies": np.asarray(vib_e2[:num_modes]) * atom_ureg.cm_1,
        "mode_indices": list(range(num_modes)),
    }
    if len(vib_e2) >= num_modes * 2:
        payload["reduced_masses"] = np.asarray(vib_e2[num_modes : num_modes * 2]) * atom_ureg.amu
    if len(vib_e2) >= num_modes * 3:
        payload["force_constants"] = (
            np.asarray(vib_e2[num_modes * 2 : num_modes * 3]) * atom_ureg.mdyne / atom_ureg.angstrom
        )
    if len(vib_e2) >= num_modes * 4:
        payload["IR_intensities"] = (
            np.asarray(vib_e2[num_modes * 3 : num_modes * 4]) * atom_ureg.km / atom_ureg.mol
        )
    mode_values = _float_array(records, "Vib-Modes")
    mode_size = num_atoms * 3
    if mode_size and len(mode_values) == num_modes * mode_size:
        payload["vibration_modes"] = [
            np.asarray(mode_values[index * mode_size : (index + 1) * mode_size]).reshape(
                num_atoms, 3
            )
            * atom_ureg.angstrom
            for index in range(num_modes)
        ]
        payload["axis_order"] = ("mode", "atom", "cartesian")
        payload["atom_order"] = "source"
        payload["normalization"] = "source_program"
        payload["mass_weighting"] = "unknown"
    return Vibrations.model_validate(payload)


def extract_fchk_thermal_information(
    records: Mapping[str, FCHKRecord],
) -> ThermalInformations | None:
    payload: dict[str, Any] = {}
    for label, field in (
        ("Thermal Energy", "U_T"),
        ("Thermal Enthalpy", "H_T"),
        ("Thermal Free Energy", "G_T"),
    ):
        value = _float_value(records, label)
        if value is not None:
            payload[field] = value * atom_ureg.Unit("hartree / particle")
    return ThermalInformations.model_validate(payload) if payload else None


def extract_fchk_status(records: Mapping[str, FCHKRecord]) -> Status:
    job_status = _int_value(records, "Job Status")
    has_energy = _float_value(records, "Total Energy") is not None
    return Status(
        scf_converged=True if has_energy else None,
        normal_terminated=(job_status == 1) if job_status is not None else None,
    )


def extract_fchk_optimization_status(
    records: Mapping[str, FCHKRecord], task_requests: list[QMTaskRequest]
) -> GeometryOptimizationStatus | None:
    if not any(task.enabled and task.task_type == "opt" for task in task_requests):
        return None
    status = extract_fchk_status(records)
    return GeometryOptimizationStatus(geometry_optimized=status.normal_terminated)


def parse_fchk_frame_records(file_content: str) -> dict[str, FCHKRecord]:
    return parse_fchk_records(file_content, wanted_labels=FRAME_LABELS)


__all__ = [
    "extract_fchk_energies",
    "extract_fchk_forces",
    "extract_fchk_hessian",
    "extract_fchk_metadata",
    "extract_fchk_nmr",
    "extract_fchk_optimization_status",
    "extract_fchk_orbitals",
    "extract_fchk_polarizability",
    "extract_fchk_populations",
    "extract_fchk_status",
    "extract_fchk_structure",
    "extract_fchk_thermal_information",
    "extract_fchk_total_spin",
    "extract_fchk_vibrations",
    "parse_fchk_frame_records",
]
