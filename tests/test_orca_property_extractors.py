"""Numerical contracts for ORCA electronic, thermal and spectroscopic results."""

import json
from pathlib import Path

import numpy as np
import pytest

from molop.io import AutoParser
from molop.io.base_models.DataClasses import Vibrations
from molop.io.logic.orca.log.frame_parsers._orca_extractors import (
    extract_orca_atomic_masses,
    extract_orca_bond_orders,
    extract_orca_electronic_states,
    extract_orca_forces,
    extract_orca_nmr,
    extract_orca_orbitals,
    extract_orca_polarizability,
    extract_orca_populations,
    extract_orca_thermochemistry,
    extract_orca_total_spin,
    extract_orca_vibrations,
)
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory
from molop.unit import atom_ureg


CORPUS = Path(__file__).parent / "test_files" / "orca" / "output_files"


@pytest.fixture(params=["memory", "auto"])
def parse_orca(request):
    def parse(relative_path):
        path = CORPUS / relative_path
        if request.param == "memory":
            return ORCALogFileParserMemory().parse(path.read_text())[-1]
        batch = AutoParser(str(path), parser_detection="orcaout", n_jobs=1)
        assert len(batch) == 1
        return batch[0][-1]

    return parse


@pytest.mark.orca_output_corpus
def test_orca_water_orbitals_mayer_and_thermochemistry_are_numeric(parse_orca):
    frame = parse_orca("local/H2O_hess_orca.out")
    orbitals = frame.molecular_orbitals
    assert orbitals is not None
    assert len(orbitals) == 24
    assert orbitals.HOMO_id == 4
    assert orbitals.LUMO_id == 5
    assert sum(orbitals.alpha_occupancies) == 5
    assert sum(orbitals.beta_occupancies) == 5
    assert (
        sum(orbitals.alpha_occupancies) + sum(orbitals.beta_occupancies)
        == sum(frame.atoms) - frame.charge
    )
    assert orbitals.HOMO_energy.m_as("hartree") == pytest.approx(-0.299506)
    assert orbitals.LUMO_energy.m_as("hartree") == pytest.approx(0.058915)

    bonds = frame.bond_orders
    assert bonds is not None
    assert bonds["mayer_bond_order"].scheme == "mayer"
    assert bonds["mayer_bond_order"].print_threshold == pytest.approx(0.1)
    assert bonds["mayer_bond_order"].source_label.startswith("Mayer bond orders larger than")
    np.testing.assert_allclose(
        bonds["mayer_bond_order"].values,
        [[np.nan, 0.9893, 0.9885], [0.9893, np.nan, np.nan], [0.9885, np.nan, np.nan]],
    )
    assert frame.temperature.m_as("K") == pytest.approx(298.15)
    assert frame.pressure.m_as("atm") == pytest.approx(1)
    thermal = frame.thermal_informations
    assert thermal is not None
    for field, value in (
        ("ZPVE", 0.01978992),
        ("TCE", 0.02262617),
        ("TCH", 0.02357038),
        ("TCG", 0.00145673),
        ("U_0", -76.25286716),
        ("U_T", -76.25003091),
        ("H_T", -76.24908670),
        ("G_T", -76.27120035),
    ):
        assert getattr(thermal, field).m_as("hartree/particle") == pytest.approx(value)
    assert thermal.S.m_as("hartree/particle/K") == pytest.approx(0.02211365 / 298.15)
    assert thermal.molecular_mass.m_as("amu") == pytest.approx(18.02)
    expected_constants = (
        np.array([31.740488, 12.269240, 8.848763]) * atom_ureg.cm_1 * atom_ureg.speed_of_light
    ).m_as("GHz")
    np.testing.assert_allclose(thermal.rotational_constants.m_as("GHz"), expected_constants)
    assert thermal.rotational_symmetry_number == 1
    # Public exports must work after populating additional quantities.
    json.dumps(frame.to_unitless_dump_with_unit_keys(), allow_nan=False)
    mol = frame.qm_embedded_rdmol()
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumBonds() == 2


@pytest.mark.orca_output_corpus
def test_orca_water_ir_and_normal_modes_align_with_source_indices(parse_orca):
    frame = parse_orca("local/H2O_hess_orca.out")
    vibrations = frame.vibrations
    assert vibrations is not None
    assert vibrations.mode_indices == [0, 1, 2]
    assert vibrations.source_mode_indices == [6, 7, 8]
    np.testing.assert_allclose(vibrations.frequencies.m_as("cm^-1"), [1567.61, 3467.7, 3651.46])
    np.testing.assert_allclose(
        vibrations.IR_intensities.m_as("km/mol"), [79.882633, 7.256436, 69.094092]
    )
    assert len(vibrations.vibration_modes) == 3
    assert all(mode.shape == (3, 3) for mode in vibrations.vibration_modes)
    np.testing.assert_allclose(
        vibrations.vibration_modes[0].magnitude,
        [[0.000188, -0.073304, 0], [-0.400098, 0.581752, 0], [0.397107, 0.581738, 0]],
    )
    assert vibrations.vibration_modes[2].magnitude[0, 0] == pytest.approx(0.073903)
    assert vibrations.axis_order == ("mode", "atom", "cartesian")
    assert vibrations.atom_order == "source"
    assert vibrations.normalization == "source_program"
    assert vibrations.mass_weighting == "not_mass_weighted"


@pytest.mark.orca_output_corpus
def test_orca_raman_activities_depolarization_and_modern_ir_are_structured(parse_orca):
    vibrations = parse_orca("cclib/basicORCA6.0/dvb_raman.out").vibrations
    assert vibrations is not None
    assert len(vibrations) == len(vibrations.vibration_modes) == 54
    assert len(vibrations.raman_activities) == len(vibrations.raman_depolarization_ratios) == 54
    np.testing.assert_allclose(vibrations.IR_intensities.m_as("km/mol")[:4], [0.02, 0, 0.37, 0.13])
    np.testing.assert_allclose(
        vibrations.raman_activities.m_as("angstrom**4/amu")[:3], [0.000035, 2.788888, 0.000008]
    )
    np.testing.assert_allclose(
        vibrations.raman_depolarization_ratios.magnitude[:3], [0.749920, 0.750000, 0.443520]
    )
    assert vibrations[1].raman_activity.m_as("angstrom**4/amu") == pytest.approx(2.788888)
    spectrum = vibrations.to_raman_spectrum()
    assert len(spectrum.bands) == 54
    assert spectrum.bands[0].metadata["mode_index"] == 0
    assert spectrum.bands[0].metadata["source_mode_index"] == 6
    assert spectrum.bands[0].metadata["depolarization_ratio"] == pytest.approx(0.749920)


@pytest.mark.orca_output_corpus
def test_orca_rocis_spin_and_truncated_spin_channels_are_structured(parse_orca):
    frame = parse_orca("cclib/basicORCA6.0/dvb_rocis.out")
    assert frame.total_spin is not None
    assert frame.total_spin.spin_square == pytest.approx(1.346022)
    orbitals = frame.molecular_orbitals
    assert orbitals is not None
    assert sum(orbitals.alpha_occupancies) == 35
    assert sum(value for value in orbitals.beta_occupancies if value is not None) == 34
    assert len(orbitals.alpha_energies) == 46
    assert len(orbitals.beta_energies) == 46
    assert np.isnan(orbitals.beta_energies.m_as("hartree")[-1])
    assert orbitals.beta_occupancies[-1] is None
    assert orbitals[45].beta_energy is None
    assert orbitals[45].beta_occupancy is None
    assert orbitals[45].alpha_energy.m_as("hartree") == pytest.approx(0.467608)
    assert orbitals.SOMO_ids == [34]
    dumped = orbitals.to_unitless_dump_with_unit_keys()
    assert dumped["beta_energies (hartree)"][-1] is None
    assert dumped["beta_occupancies"][-1] is None
    assert json.loads(json.dumps(dumped, allow_nan=False))["beta_energies (hartree)"][-1] is None
    table = orbitals.to_orbital_table("beta")
    assert table.row(45)["energy"] is None
    assert table.row(45)["occupancy"] is None
    assert table.row(44)["energy"].m_as("hartree") == pytest.approx(0.468058)
    json.dumps(orbitals.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)


@pytest.mark.orca_output_corpus
@pytest.mark.parametrize(
    "relative_path, first_isotropic",
    [
        ("cclib/basicORCA6.0/dvb_nmr.out", 114.398),
        ("cclib/basicORCA6.1/dvb_coupling_nmr.out", 114.114),
    ],
)
def test_orca_nmr_shielding_and_coupling_tables_are_numeric(
    parse_orca, relative_path, first_isotropic
):
    frame = parse_orca(relative_path)
    nmr = frame.nmr
    assert nmr is not None
    assert nmr.gauge == "GIAO"
    assert [tensor.atom_index for tensor in nmr.shielding_tensors] == list(range(20))
    first = nmr.shielding_tensors[0]
    assert first.atom_symbol == "C"
    assert first.orientation == "source"
    assert first.isotropic.m_as("ppm") == pytest.approx(first_isotropic)
    for shielding in nmr.shielding_tensors:
        assert np.trace(shielding.shielding_tensor.m_as("ppm")) / 3 == pytest.approx(
            shielding.isotropic.m_as("ppm"), abs=1e-3
        )
    if "coupling" in relative_path:
        assert nmr.coupling_atom_indices == list(range(20))
        assert nmr.spin_spin_coupling_j.shape == (20, 20)
        assert nmr.spin_spin_coupling_j.m_as("Hz")[0, 1] == pytest.approx(9.485)
        assert nmr.spin_spin_coupling_j.m_as("Hz")[2, 6] == pytest.approx(133.980)
        np.testing.assert_allclose(
            nmr.spin_spin_coupling_j.magnitude, nmr.spin_spin_coupling_j.magnitude.T
        )
        assert nmr.spin_spin_coupling_j_components["FC"].m_as("Hz")[0, 1] == pytest.approx(7.440)
        orbitals = frame.molecular_orbitals
        assert orbitals is not None
        assert orbitals.alpha_symmetries[0] == "1-Bu"
        assert orbitals.beta_symmetries[0] == "1-Bu"
    else:
        np.testing.assert_allclose(
            first.shielding_tensor.m_as("ppm"),
            [[50.006, 4.530, 0], [6.563, 80.454, 0], [0, 0, 212.734]],
        )
        assert first.anisotropy.m_as("ppm") == pytest.approx(147.504)
        np.testing.assert_allclose(first.principal_values.m_as("ppm"), [49.029, 81.431, 212.734])
    json.dumps(nmr.to_unitless_dump_with_unit_keys(), allow_nan=False)


def test_orca_missing_property_sections_do_not_fabricate_results():
    assert extract_orca_orbitals("") is None
    assert extract_orca_bond_orders("", 3) is None
    assert extract_orca_total_spin("") is None
    assert extract_orca_thermochemistry("") == {}
    assert extract_orca_nmr("", 3) is None
    assert extract_orca_polarizability("") is None
    assert extract_orca_electronic_states("", None) is None


@pytest.mark.parametrize(
    "header,rows",
    [
        ("MULLIKEN ATOMIC CHARGES", ["0 O : -0.8", "1 H : 0.4", "2 H : 0.4"]),
        ("LOEWDIN ATOMIC CHARGES", ["0 O : -0.8", "1 H : 0.4", "2 H : 0.4"]),
        ("HIRSHFELD ANALYSIS", ["0 O -0.8 0.0", "1 H 0.4 0.0", "2 H 0.4 0.0"]),
    ],
)
@pytest.mark.parametrize("mutation", ["missing_atom", "wrong_index", "wrong_element"])
def test_orca_population_rows_must_match_source_atoms(header, rows, mutation):
    valid = extract_orca_populations(header + "\n" + "\n".join(rows), [8, 1, 1])
    assert valid is not None
    assert valid.population_items()[0][1].values == [-0.8, 0.4, 0.4]
    rows = list(rows)
    if mutation == "missing_atom":
        rows.pop(1)
    elif mutation == "wrong_index":
        rows[1] = rows[1].replace("1 H", "2 H")
    else:
        rows[1] = rows[1].replace("1 H", "1 O")
    assert extract_orca_populations(header + "\n" + "\n".join(rows), [8, 1, 1]) is None


@pytest.mark.parametrize("mp2", [False, True])
def test_orca_gradient_rows_do_not_ignore_source_atom_indices(mp2):
    if mp2:
        text = "The final MP2 gradient\n0: 0.1 0.2 0.3\n1: -0.1 -0.2 -0.3\n"
        invalid = text.replace("1:", "0:")
    else:
        text = "------------------\nCARTESIAN GRADIENT\n------------------\n1 H : 0.1 0.2 0.3\n2 O : -0.1 -0.2 -0.3\n"
        invalid = text.replace("2 O", "1 O")
        assert extract_orca_forces(text.replace("2 O", "2 H"), 2, [1, 8]) is None
    forces = extract_orca_forces(text, 2, [1, 8])
    assert forces is not None
    np.testing.assert_allclose(forces.magnitude.sum(axis=0), 0)
    assert extract_orca_forces(invalid, 2, [1, 8]) is None


def test_orca_mass_rows_require_matching_source_atom_labels_and_zero_based_indices():
    text = "----------------------------\nCARTESIAN COORDINATES (A.U.)\n----------------------------\n0 H 1.0 0 1.008 0.0 0.0 0.0\n1 O 8.0 0 15.999 1.0 0.0 0.0\n"
    masses = extract_orca_atomic_masses(text, atoms=[1, 8])
    assert masses is not None
    np.testing.assert_allclose(masses.magnitude, [1.008, 15.999])
    assert extract_orca_atomic_masses(text, atoms=[8, 1]) is None
    assert extract_orca_atomic_masses(text.replace("0 H", "1 H").replace("1 O", "2 O")) is None


def test_orca_orbitals_preserve_fractional_occupations_and_latest_table():
    table = "Hartree-Fock type HFTyp .... RHF\nORBITAL ENERGIES\n NO OCC E(Eh) E(eV)\n 0 1.5 -1D+0 -27.2\n 1 0.5 0.2 5.44\n"
    orbitals = extract_orca_orbitals(table + table.replace("-1D+0", "-2D+0"))
    assert orbitals is not None
    np.testing.assert_allclose(orbitals.alpha_energies.m_as("hartree"), [-2, 0.2])
    assert orbitals.alpha_occupancies == orbitals.beta_occupancies == [0.75, 0.25]
    assert extract_orca_orbitals(table.replace(" 1 0.5", " 3 0.5")) is None


@pytest.mark.parametrize(
    "reference, evidence, expected_alpha, expected_beta",
    [
        ("ROHF", "ROHF_case HIGHSPIN", [1, 1, 0], [1, 0, 0]),
        ("ROHF", "Multiplicity Mult .... 2", [1, 1, 0], [1, 0, 0]),
        ("ROKS", "", [1, 1, 0], [1, 0, 0]),
        ("ROHF", "", [1, None, 0], [1, None, 0]),
        ("ROHF", "ROHF_case SAHF", [1, None, 0], [1, None, 0]),
        ("ROHF", "ROHF_case CAHF", [1, None, 0], [1, None, 0]),
        ("ROHF", "ROHF_case AF_CSF", [1, None, 0], [1, None, 0]),
        ("RHF", "Multiplicity Mult .... 2", [1, 0.5, 0], [1, 0.5, 0]),
    ],
)
def test_orca_restricted_spin_occupations_follow_reference_and_evidence(
    reference, evidence, expected_alpha, expected_beta
):
    orbitals = extract_orca_orbitals(f"""Hartree-Fock type HFTyp .... {reference}
{evidence}
ORBITAL ENERGIES
NO OCC E(Eh) E(eV)
0 2.0 -1.0 -27.2
1 1.0 -0.2 -5.4
2 0.0 0.1 2.7
""")
    assert orbitals is not None
    assert orbitals.alpha_occupancies == expected_alpha
    assert orbitals.beta_occupancies == expected_beta
    assert orbitals.SOMO_ids == ([1] if expected_beta[1] == 0 else [])
    np.testing.assert_allclose(orbitals.alpha_energies.magnitude, [-1, -0.2, 0.1])
    np.testing.assert_allclose(orbitals.beta_energies.magnitude, [-1, -0.2, 0.1])


@pytest.mark.orca_output_corpus
@pytest.mark.parametrize(
    "relative_path, vector_au, magnitude_debye",
    [
        ("cclib/basicORCA5.0/Trp_polar.out", [-0.47508, 0.37562, -0.87797], 2.71106),
        (
            "cclib/basicORCA6.0/Trp_polar.out",
            [-0.475081306, 0.375626728, -0.877962624],
            2.711058506,
        ),
        ("local/H2O_hess_orca.out", [0.00165, -0.74767, 0], 1.90043),
    ],
)
def test_orca_dipole_vector_converts_atomic_units_to_printed_debye(
    parse_orca, relative_path, vector_au, magnitude_debye
):
    polar = parse_orca(relative_path).polarizability
    assert polar is not None
    np.testing.assert_allclose(polar.dipole.m_as("elementary_charge * bohr"), vector_au)
    # ORCA's displayed conversion differs slightly from Pint's constants.
    assert np.linalg.norm(polar.dipole.m_as("debye")) == pytest.approx(magnitude_debye, rel=3e-5)
    np.testing.assert_allclose(polar.electric_dipole_moment.magnitude, polar.dipole.magnitude)


@pytest.mark.orca_output_corpus
def test_orca_td_states_are_not_matrix_rows_or_cd_velocity_observables(parse_orca):
    states = parse_orca("cclib/basicORCA6.0/dvb_td.out").electronic_states
    assert states is not None
    assert len(states) == 10
    assert [state.root for state in states] == list(range(1, 11))
    assert [state.multiplicity for state in states] == [1] * 5 + [3] * 5
    assert states[0].excitation_energy.m_as("eV") == pytest.approx(5.352)
    assert states[0].oscillator_strength == pytest.approx(0.005038390)
    assert states[1].oscillator_strength == pytest.approx(1.170940055)
    np.testing.assert_allclose(
        states[0].transition_dipole.m_as("elementary_charge * bohr"),
        [0.16288, 0.10906, 0],
    )
    assert states[5].excitation_energy.m_as("eV") == pytest.approx(3.129)
    assert states[5].oscillator_strength == 0
    assert states[5].properties["absorption_root"] == 1
    assert states[0].properties["absorption_root"] == 1


@pytest.mark.orca_output_corpus
@pytest.mark.parametrize(
    "filename, count, energy, strength, convention",
    [
        ("dvb_rocis.out", 4, 291.659506, 0.030923917, "electric_dipole"),
        ("dvb_adc2.log", 2, 5.554128, 0.001260182, "left-right"),
        ("dvb_eom_ccsd.log", 2, 5.471765, 0.000528270, "left-right"),
        ("dvb_steom_ccsd.log", 2, 4.690735, 0.000317410, "left-right"),
    ],
)
def test_orca_absorption_table_keeps_its_observable_and_convention(
    parse_orca, filename, count, energy, strength, convention
):
    states = parse_orca(f"cclib/basicORCA6.0/{filename}").electronic_states
    assert states is not None
    assert len(states) == count
    assert states[0].excitation_energy.m_as("eV") == pytest.approx(energy)
    assert states[0].oscillator_strength == pytest.approx(strength)
    assert states[0].properties["absorption_convention"] == convention


def test_orca_spectrum_requires_section_units_and_does_not_fabricate_roots():
    rows = "1 43168.0 231.7 0.005038390 0.03842 0.16288 0.10906 0.0\n"
    assert extract_orca_electronic_states("OVERLAP MATRIX\n" + rows, None) is None
    assert extract_orca_electronic_states("TD-DFT/TDA EXCITED STATES\n", None) is None
    heading = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
    assert extract_orca_electronic_states(heading + rows, None) is None
    legacy = (
        heading
        + "State Energy Wavelength fosc T2 TX TY TZ\n(cm-1) (nm) (au**2) (au) (au) (au)\n"
        + rows
    )
    states = extract_orca_electronic_states(legacy + "\n" + rows.replace("1 ", "99 ", 1), None)
    assert states is not None
    assert len(states) == 1
    assert states[0].root == 1
    assert states[0].excitation_energy.m_as("eV") == pytest.approx(5.352149878)
    np.testing.assert_allclose(
        states[0].transition_dipole.m_as("elementary_charge * bohr"), [0.16288, 0.10906, 0]
    )


def test_orca_mayer_validates_atom_indices_and_uses_latest_table():
    table = "Mayer bond orders larger than 0.1\nB( 0-O , 1-H ) : 0.9 B( 0-O , 2-H ) : 8D-1\n\n"
    bonds = extract_orca_bond_orders(table + table.replace("0.9", "0.7"), 3)
    assert bonds is not None
    assert bonds["mayer_bond_order"].values[0, 1] == pytest.approx(0.7)
    assert bonds["mayer_bond_order"].values[2, 0] == pytest.approx(0.8)
    assert extract_orca_bond_orders(table, 2) is None
    assert extract_orca_bond_orders(table, None) is None


def test_orca_incomplete_unrestricted_table_does_not_become_restricted():
    orbitals = extract_orca_orbitals(
        "ORBITAL ENERGIES\nSPIN UP ORBITALS\nNO OCC E(Eh) E(eV)\n0 1.0 -0.5 -13.6\n"
    )
    assert orbitals is not None
    assert orbitals.alpha_occupancies == [1]
    assert orbitals.beta_occupancies == [None]
    assert np.isnan(orbitals.beta_energies.m_as("hartree")).all()
    assert orbitals.beta_HOMO_id is None
    assert orbitals.beta_LUMO_id is None
    assert orbitals.beta_HOMO_LUMO_gap is None
    assert orbitals.SOMO_ids == []
    assert orbitals.to_unitless_dump_with_unit_keys()["beta_energies (hartree)"] == [None]


def test_orca_pads_shorter_alpha_channel_and_its_symmetries():
    orbitals = extract_orca_orbitals("""ORBITAL ENERGIES
SPIN UP ORBITALS
NO OCC E(Eh) E(eV) Irrep
0 1.0 -0.5 -13.6 1-Ag
SPIN DOWN ORBITALS
NO OCC E(Eh) E(eV) Irrep
0 1.0 -0.6 -16.3 1-Ag
1 0.0 0.1 2.7 1-Bu
""")
    assert orbitals is not None
    assert len(orbitals.alpha_energies) == len(orbitals.beta_energies) == 2
    assert np.isnan(orbitals.alpha_energies.m_as("hartree")[1])
    assert orbitals.alpha_occupancies == [1, None]
    assert orbitals.alpha_symmetries == ["1-Ag", None]
    assert orbitals.beta_symmetries == ["1-Ag", "1-Bu"]
    assert orbitals.SOMO_ids == []


def test_orca_symmetry_labels_in_normal_mode_blocks():
    text = """VIBRATIONAL FREQUENCIES
0: 10.0 cm**-1
1: 20.0 cm**-1
NORMAL MODES
------------
0
Ag
0 0.1
1 0.2
2 0.3
1
Bu
0 0.4
1 0.5
2 0.6
"""
    vibrations = extract_orca_vibrations(text, 1)
    assert vibrations is not None
    assert len(vibrations.vibration_modes) == 2
    np.testing.assert_allclose(vibrations.vibration_modes[1].magnitude, [[0.4, 0.5, 0.6]])


def test_orca_last_temperature_block_and_entropy_energy_semantics():
    block = """THERMOCHEMISTRY AT 300K
Temperature ... 300 K
Pressure ... 1 atm
Electronic energy ... -10 Eh
Zero point energy ... 0.1 Eh
Total thermal energy -9.8 Eh
Total Enthalpy ... -9.7 Eh
Final entropy term ... 0.03 Eh
Final Gibbs free energy ... -9.73 Eh
"""
    result = extract_orca_thermochemistry(block + block.replace("300", "600"))
    assert result["temperature"].m_as("K") == 600
    thermal = result["thermal_informations"]
    assert thermal.TCE.m_as("hartree/particle") == pytest.approx(0.2)
    assert thermal.S.m_as("hartree/particle/K") == pytest.approx(0.03 / 600)


def test_orca_partial_ir_does_not_consume_following_raman_table():
    text = """VIBRATIONAL FREQUENCIES
 0: 100.0 cm**-1
 1: 200.0 cm**-1
IR SPECTRUM
-----------
Mode freq Int
RAMAN SPECTRUM
--------------
Mode freq Activity Depolarization
 0: 100.0 1.2 0.5
 1: 200.0 2.3 0.7
"""
    vibrations = extract_orca_vibrations(text, None)
    assert vibrations is not None
    assert len(vibrations.IR_intensities) == 0
    np.testing.assert_allclose(vibrations.raman_activities.magnitude, [1.2, 2.3])


def test_orca_partial_normal_modes_and_incomplete_spectra_do_not_misalign():
    text = """VIBRATIONAL FREQUENCIES
 6: 100.0 cm**-1
 7: 200.0 cm**-1
NORMAL MODES
------------
 6 7
 0 0.1 0.2
IR SPECTRUM
-----------
 Mode freq Int
 7: 200 5.0
"""
    vibrations = extract_orca_vibrations(text, 3)
    assert vibrations is not None
    assert vibrations.vibration_modes == []
    assert len(vibrations.IR_intensities) == 0


def test_orca_coupling_only_output_preserves_sparse_source_atom_indices():
    nmr = extract_orca_nmr("J[2,5](FC) 1 1 1 iso= 2D+1\nJ[2,5](Total) 1 1 1 iso= 21.5", 6)
    assert nmr is not None
    assert nmr.coupling_atom_indices == [2, 5]
    np.testing.assert_allclose(
        nmr.spin_spin_coupling_j.m_as("Hz"), [[np.nan, 21.5], [21.5, np.nan]]
    )
    assert nmr.spin_spin_coupling_j_components["FC"].m_as("Hz")[1, 0] == 20
    assert extract_orca_nmr("Nucleus 0C :\nTotal shielding tensor (ppm):\n1 0 0\n", 1) is None


def test_raman_fields_survive_iteration_filtering_and_public_export():
    vibrations = Vibrations(
        frequencies=np.array([-10, 20]) * atom_ureg.cm_1,
        mode_indices=[6, 7],
        raman_activities=np.array([2, 3]) * atom_ureg.angstrom**4 / atom_ureg.amu,
        raman_depolarization_ratios=np.array([0.5, 0.7]) * atom_ureg.dimensionless,
    )
    imaginary = vibrations.imaginary_vibrations
    assert imaginary.mode_indices == [6]
    assert imaginary[0].raman_activity.magnitude == 2
    assert imaginary[0].raman_depolarization_ratio.magnitude == 0.5
    dumped = vibrations.to_unitless_dump_with_unit_keys()
    assert any(key.startswith("raman_activities (") for key in dumped)
    json.dumps(dumped, allow_nan=False)
    with pytest.raises(ValueError, match="raman_activities must match frequencies"):
        vibrations.model_copy(
            update={"raman_activities": np.array([1]) * atom_ureg.angstrom**4 / atom_ureg.amu}
        ).validate_mode_metadata()


@pytest.mark.orca_output_corpus
def test_orca_structure_only_skips_additional_properties():
    path = CORPUS / "local/H2O_hess_orca.out"
    frame = ORCALogFileParserMemory(only_extract_structure=True).parse(path.read_text())[-1]
    assert len(frame.atoms) == 3
    for name in (
        "molecular_orbitals",
        "bond_orders",
        "thermal_informations",
        "vibrations",
        "nmr",
        "total_spin",
    ):
        assert getattr(frame, name) is None
