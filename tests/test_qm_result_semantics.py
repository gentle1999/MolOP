"""Cross-format contracts for result identity, missing data and provenance."""

import json
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

from molop.io.base_models.DataClasses import (
    NMR,
    BondOrderMatrix,
    BondOrders,
    ElectronicState,
    ElectronicStates,
    MolecularOrbitals,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_extractors import extract_fchk_orbitals
from molop.io.logic.gaussian.fchk.output.parsers._fchk_records import FCHKRecord
from molop.io.logic.gaussian.fchk.output.parsers.G16FchkFileParser import G16FchkFileParserMemory
from molop.io.logic.gaussian.log.frame_parsers._g16_extractors import _extract_nmr_coupling_matrix
from molop.io.logic.gaussian.log.frame_parsers._g16_shared import (
    _extract_molecular_orbital_payload_from_text,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.log.frame_parsers._orca_extractors import extract_orca_nmr
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory
from molop.io.logic.xtb.output.parsers.XTBOutputFileParser import XTBOutputFileParserMemory
from molop.unit import atom_ureg


FILES = Path(__file__).parent / "test_files"


def test_state_identity_does_not_depend_on_position_or_program_root():
    excited = ElectronicState(
        state_index=0, root=1, is_ground_state=False, excitation_energy=2 * atom_ureg.eV
    )
    ground = ElectronicState(
        state_index=4, root=0, is_ground_state=True, energy=-1 * atom_ureg.hartree
    )
    states = ElectronicStates(states=[excited, ground])
    assert states.ground_state is ground
    assert states.excited_states == [excited]
    transition = states.to_transitions()[0]
    assert (transition.initial_state, transition.final_state) == (4, 0)
    json.dumps(states.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)
    with pytest.raises(ValidationError, match="Only one"):
        ElectronicStates(states=[ground, ground.model_copy(update={"state_index": 5})])


@pytest.mark.orca_output_corpus
@pytest.mark.parametrize("filename, count", [("dvb_td.out", 10), ("dvb_adc2.log", 2)])
def test_orca_excited_only_results_export_every_observed_transition(filename, count):
    path = FILES / "orca/output_files/cclib/basicORCA6.0" / filename
    states = ORCALogFileParserMemory().parse(path.read_text(encoding="utf-8"))[-1].electronic_states
    assert states.ground_state is None
    assert len(states.excited_states) == count
    assert len(states.to_spectrum().transitions) == count
    assert [state.state_index for state in states] == list(range(count))
    assert all(transition.initial_state is None for transition in states.to_transitions())
    assert states.to_transitions()[0].properties["absorption_root"] == 1
    json.dumps(states.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)


@pytest.mark.parametrize("format_id", ["gaussian", "orca"])
def test_sparse_couplings_keep_unobserved_pairs_distinct_from_printed_zero(format_id):
    if format_id == "orca":
        nmr = extract_orca_nmr("J[0,1](Total) 1 1 1 iso=0.0\nJ[1,2](Total) 1 1 1 iso=2.5", 3)
    else:
        matrix = _extract_nmr_coupling_matrix(" 1 2 3\n 2 0.000000D+00\n 3 2.500000D+00", 3)
        # Both values in this synthetic lower-triangle table belong to column 1.
        nmr = NMR(spin_spin_coupling_j=matrix * atom_ureg.Hz, coupling_atom_indices=[0, 1, 2])
    matrix = nmr.spin_spin_coupling_j.m_as("Hz")
    assert matrix[0, 1] == matrix[1, 0] == 0
    unknown = (0, 2) if format_id == "orca" else (1, 2)
    assert np.isnan(matrix[unknown])
    dumped = nmr.to_unitless_dump_with_unit_keys()
    assert dumped["spin_spin_coupling_j (hertz)"][unknown[0]][unknown[1]] is None
    assert nmr.to_unitless_dump()["spin_spin_coupling_j"][0][1] == 0
    binary = nmr.to_unitless_dump_with_unit_keys(array_mode="ndarray")
    assert np.isnan(binary["spin_spin_coupling_j (hertz)"][unknown])
    json.dumps(dumped, allow_nan=False)


def test_missing_matrix_serializers_do_not_accept_infinities_or_mutate_data():
    matrix = np.array([[np.nan, 0.0], [0.0, np.nan]])
    bonds = BondOrders(
        bond_orders={
            "mayer_bond_order": BondOrderMatrix(scheme="mayer", values=matrix, print_threshold=0.1)
        }
    )
    assert bonds.to_unitless_dump()["bond_orders"]["mayer_bond_order"]["values"] == [
        [None, 0],
        [0, None],
    ]
    assert (
        bonds.to_property_bundle().tensor_properties["mayer_bond_order"].metadata["print_threshold"]
        == 0.1
    )
    json.dumps(bonds.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)
    assert np.isnan(matrix[0, 0])
    invalid = NMR(spin_spin_coupling_j=np.array([[np.inf]]) * atom_ureg.Hz)
    with pytest.raises(ValueError, match="non-finite"):
        invalid.to_unitless_dump_with_unit_keys()


def test_gaussian_restricted_sharing_requires_reference_evidence():
    text = "Alpha occ. eigenvalues -- -0.5\nAlpha virt. eigenvalues -- 0.1\n"
    closed = MolecularOrbitals.model_validate(
        _extract_molecular_orbital_payload_from_text(text, reference_method="RHF")
    )
    partial = MolecularOrbitals.model_validate(
        _extract_molecular_orbital_payload_from_text(text, reference_method="UHF")
    )
    assert closed.beta_occupancies == [1, 0]
    assert partial.beta_occupancies == [None, None]
    assert np.isnan(partial.beta_energies.magnitude).all()
    records = {
        "Alpha Orbital Energies": FCHKRecord("Alpha Orbital Energies", "R", [-0.5, 0.1], 2),
        "Number of alpha electrons": FCHKRecord("Number of alpha electrons", "I", 1),
        "Number of beta electrons": FCHKRecord("Number of beta electrons", "I", 1),
    }
    assert extract_fchk_orbitals(records, reference_method="RHF").beta_occupancies == [1, 0]
    assert extract_fchk_orbitals(records, reference_method="UHF").beta_occupancies == [None, None]


def test_real_gaussian_log_and_fchk_closed_shell_channels_agree():
    log = G16LogFileParserMemory().parse(
        (FILES / "g16log/test_nmr.log").read_text(encoding="utf-8")
    )[-1]
    fchk = G16FchkFileParserMemory().parse(
        (FILES / "g16fchk/molecule.fchk").read_text(encoding="utf-8")
    )[-1]
    a, b = log.molecular_orbitals, fchk.molecular_orbitals
    assert a.orbital_basis == b.orbital_basis == "shared_spatial"
    assert a.alpha_occupancies == a.beta_occupancies == b.alpha_occupancies == b.beta_occupancies
    np.testing.assert_allclose(a.beta_energies.magnitude, b.beta_energies.magnitude, atol=1e-5)


def test_xtb_frontiers_survive_omitted_orbitals_and_fractional_smearing():
    frame = XTBOutputFileParserMemory().parse(
        (FILES / "xtbout/6-6-1-hess.out").read_text(encoding="utf-8")
    )[-1]
    orbitals = frame.molecular_orbitals
    assert (orbitals.HOMO_id, orbitals.LUMO_id) == (22, 23)
    assert orbitals.HOMO_LUMO_gap.m_as("eV") == pytest.approx(0.28295295, abs=5e-6)
    assert orbitals.SOMO_ids == []
    assert orbitals.total_occupancies[22:24] == [1.9917, 0.0084]
    assert orbitals.alpha_occupancies[22] == orbitals.beta_occupancies[22] == 1.9917 / 2
    assert orbitals.alpha_occupancies[1] is orbitals.beta_occupancies[1] is None
    assert orbitals.to_orbital_table().row(1)["energy"] is None
    json.dumps(orbitals.to_unitless_dump_with_unit_keys(), allow_nan=False)


@pytest.mark.parametrize(
    "parser, relative, source_indices, entropy_model",
    [
        (G16LogFileParserMemory, "g16log/2-modfreq.log", [1, 2, 3], "RRHO"),
        (G16FchkFileParserMemory, "g16fchk/dsgdb9nsd_000001-3.fchk", [None] * 3, None),
        (XTBOutputFileParserMemory, "xtbout/6-6-1-hess.out", [7, 8, 9], "QRRHO"),
        pytest.param(
            ORCALogFileParserMemory,
            "orca/output_files/local/H2O_hess_orca.out",
            [6, 7, 8],
            "QRRHO",
            marks=pytest.mark.orca_output_corpus,
        ),
    ],
)
def test_mode_identity_and_thermal_conditions_survive_container_interfaces(
    parser, relative, source_indices, entropy_model
):
    frame = parser().parse((FILES / relative).read_text(encoding="utf-8"))[-1]
    modes = frame.vibrations
    assert modes.mode_indices == list(range(len(modes)))
    assert modes.source_mode_indices[:3] == source_indices
    first = modes[0]
    assert first.mode_index == 0
    assert first.source_mode_index == source_indices[0]
    assert first.normalization == modes.normalization
    assert first.mass_weighting == modes.mass_weighting
    assert modes[:1][0].source_mode_index == first.source_mode_index
    assert modes.to_ir_spectrum().bands[0].metadata["source_mode_index"] == first.source_mode_index
    for item in modes.imaginary_vibrations:
        assert item.source_mode_index == modes.source_mode_indices[item.mode_index]
    thermal = frame.thermal_informations
    assert thermal.entropy_model == entropy_model
    bundle = thermal.to_property_bundle()
    assert bundle.metadata["entropy_model"] == entropy_model
    if entropy_model is not None:
        assert thermal.temperature.m_as("K") == pytest.approx(298.15)
        assert bundle.scalar_properties["temperature"].m_as("K") == pytest.approx(298.15)
    if frame.qm_software == "xTB":
        assert thermal.rotor_reference_frequency.m_as("cm^-1") == 50
        assert thermal.imaginary_frequency_cutoff.m_as("cm^-1") == -20
        assert thermal.frequency_scale_factor == 1
    json.dumps(bundle.to_unitless_dump_with_unit_keys(), allow_nan=False)
