"""Cross-format numerical contracts backed by locally calculated Raman fixtures."""

import json
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

from molop.io import AutoParser
from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.DataClasses import Vibrations
from molop.io.logic.gaussian.fchk.output.parsers._fchk_extractors import extract_fchk_vibrations
from molop.io.logic.gaussian.fchk.output.parsers._fchk_records import parse_fchk_records
from molop.io.logic.gaussian.fchk.output.parsers.G16FchkFileParser import G16FchkFileParserMemory
from molop.io.logic.gaussian.log.frame_parsers._g16_extractors import (
    ParseState,
    extract_bond_orders_from_state,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.xtb.output.parsers._xtb_output_extractors import extract_xtb_bond_orders
from molop.io.logic.xtb.output.parsers.XTBOutputFileParser import XTBOutputFileParserMemory
from molop.unit import atom_ureg


FIXTURES = Path(__file__).parent / "test_files"


@pytest.fixture(params=["memory", "auto"])
def parse_fixture(request):
    parsers = {
        "g16log": G16LogFileParserMemory,
        "g16fchk": G16FchkFileParserMemory,
        "xtbout": XTBOutputFileParserMemory,
    }

    def parse(relative_path):
        path = FIXTURES / relative_path
        if request.param == "auto":
            return AutoParser(str(path), parser_detection=relative_path.split("/", 1)[0])[0][-1]
        return parsers[relative_path.split("/", 1)[0]]().parse(path.read_text(encoding="utf-8"))[-1]

    return parse


def test_gaussian_closed_shell_nbo_matrices_reach_common_container(parse_fixture):
    frame = parse_fixture("g16log/anion_0167_opt_g16_nbo_sp.log")
    orders = frame.bond_orders
    assert orders is not None
    assert set(orders.bond_order_names) == {
        "wiberg_bond_order",
        "atom_atom_overlap_bond_order",
        "mo_bond_order",
    }
    for key, expected in (
        ("wiberg_bond_order", 1.1790),
        ("atom_atom_overlap_bond_order", 0.9822),
        ("mo_bond_order", 0.4525),
    ):
        record = orders[key]
        assert record.values.shape == (10, 10)
        np.testing.assert_allclose(record.values, record.values.T)
        assert record.values[0, 1] == expected
        assert record.values[0, 0] == 0
        assert record.spin_channel == "total"
        assert record.print_threshold is None
        assert record.source_label
    assert orders["mo_bond_order"].values[0, 2] == -0.1565
    assert orders["wiberg_bond_order"].values[0, 9] == 0.0008
    assert orders["wiberg_bond_order"].metadata["basis"] == "NAO"
    json.dumps(orders.to_property_bundle().to_unitless_dump_with_unit_keys(), allow_nan=False)


def test_gaussian_nbo_spin_channels_are_independent_observations(parse_fixture):
    orders = parse_fixture("g16log/radical_0554_opt_g16_nbo_sp.log").bond_orders
    assert orders is not None and len(orders) == 9
    for key, channel, expected in (
        ("wiberg_bond_order", "total", 1.0200),
        ("wiberg_bond_order_for_alpha_spin", "alpha", 0.2505),
        ("wiberg_bond_order_for_beta_spin", "beta", 0.2610),
    ):
        record = orders[key]
        assert record.spin_channel == channel
        assert record.values.shape == (27, 27)
        assert record.values[0, 1] == expected
    assert orders["mo_bond_order_for_alpha_spin"].values[0, 1] == 0.4580
    assert orders["mo_bond_order_for_beta_spin"].values[0, 1] == 0.4744


def test_truncated_nbo_matrix_preserves_unknown_cells_and_printed_zero():
    content = """Wiberg bond index matrix in the NAO basis:
     Atom    1       2
     ---- ------  ------
   1.  H  0.0000  0.9000
"""
    orders = extract_bond_orders_from_state(ParseState(content), 2)
    assert orders is not None
    record = orders["wiberg_bond_order"]
    assert record.values[0, 0] == 0
    assert np.isnan(record.values[1]).all()
    assert record.missing_reason
    assert orders.to_unitless_dump()["bond_orders"]["wiberg_bond_order"]["values"][1] == [
        None,
        None,
    ]
    assert (
        extract_bond_orders_from_state(
            ParseState(content.replace("Atom    1       2", "Atom    1       3")), 2
        )
        is None
    )


def test_xtb_legacy_wiberg_neighbors_are_thresholded_not_mayer(parse_fixture):
    frame = parse_fixture("xtbout/opt_point_charge_xtb.out")
    orders = frame.bond_orders
    assert orders is not None and orders.bond_order_names == ["wiberg_bond_order"]
    record = orders["wiberg_bond_order"]
    assert record.scheme == "wiberg"
    assert record.print_threshold == 0.10
    assert record.values.shape == (5, 5)
    assert record.values[0, 1] == 0.998
    assert record.values[0, 4] == 0.999
    assert np.isnan(record.values[1, 2])
    assert np.isnan(record.values.diagonal()).all()
    assert record.metadata["basis"] == "AO"
    assert record.metadata["reciprocal_value_policy"] == "first_printed"
    json.dumps(orders.to_unitless_dump(), allow_nan=False)


def test_xtb_modern_wrapped_neighbors_parse_without_inventing_geometry(parse_fixture):
    frame = parse_fixture("xtbout/6-6-1-hess.out")
    assert frame.atoms == []
    record = frame.bond_orders["wiberg_bond_order"]
    assert record.values.shape == (14, 14)
    assert record.values[0, 9] == 0.947
    assert record.values[2, 4] == 0.176
    assert record.values[5, 13] == 0.935
    assert record.metadata["atom_count_source"] == "source WBO row indices"
    assert np.isnan(record.values[0, 2])
    assert len(frame.vibrations.raman_activities) == 0  # old, undefined "(amu)" output


def test_xtb_wbo_invalid_atom_reference_is_not_silently_clipped():
    content = (FIXTURES / "xtbout/water_raman_ptb_6_7_1.out").read_text(encoding="utf-8")
    assert extract_xtb_bond_orders(content, 2) is None


def test_gaussian_raman_activity_and_both_polarizations_match_fchk(parse_fixture):
    log = parse_fixture("g16log/water_raman_hf_sto3g.log")
    fchk = parse_fixture("g16fchk/water_raman_hf_sto3g.fchk")
    for frame in (log, fchk):
        vibrations = frame.vibrations
        assert vibrations is not None and len(vibrations) == 3
        np.testing.assert_allclose(
            vibrations.raman_activities.m_as("angstrom**4/amu"),
            [9.2638, 47.7964, 21.5447],
            atol=5e-5,
            rtol=0,
        )
        np.testing.assert_allclose(
            vibrations.raman_depolarization_ratios_polarized.magnitude,
            [0.7246, 0.1791, 0.7500],
            atol=5e-5,
            rtol=0,
        )
        np.testing.assert_allclose(
            vibrations.raman_depolarization_ratios_unpolarized.magnitude,
            [0.8403, 0.3038, 0.8571],
            atol=5e-5,
            rtol=0,
        )
        assert len(vibrations.raman_depolarization_ratios) == 0
        mode = vibrations[0]
        assert mode.raman_depolarization_ratio is None
        assert mode.raman_depolarization_ratio_polarized is not None
        assert mode.raman_depolarization_ratio_unpolarized is not None
        metadata = vibrations.to_raman_spectrum().bands[0].metadata
        assert metadata["depolarization_ratio_polarized"] == pytest.approx(
            mode.raman_depolarization_ratio_polarized.magnitude
        )
        assert metadata["depolarization_ratio_unpolarized"] == pytest.approx(
            mode.raman_depolarization_ratio_unpolarized.magnitude
        )
        json.dumps(vibrations.to_unitless_dump(), allow_nan=False)


def test_gaussian_linear_one_mode_keeps_displacements_and_raman(parse_fixture):
    frame = parse_fixture("g16log/hydrogen_raman_hf_sto3g.log")
    vibrations = frame.vibrations
    assert vibrations is not None and len(vibrations) == 1
    assert vibrations.vibration_modes[0].shape == (2, 3)
    assert np.linalg.norm(vibrations.vibration_modes[0].magnitude) > 0
    assert len(vibrations.raman_activities) == 1
    assert vibrations.raman_activities[0].magnitude > 0
    assert len(vibrations.raman_depolarization_ratios_polarized) == 1


def test_xtb_ptb_raman_uses_stdout_units_and_preserves_method_identity(parse_fixture):
    frame = parse_fixture("xtbout/water_raman_ptb_6_7_1.out")
    vibrations = frame.vibrations
    assert vibrations is not None
    np.testing.assert_allclose(vibrations.frequencies.m_as("cm_1"), [1538.34, 3645.62, 3659.41])
    np.testing.assert_allclose(
        vibrations.raman_activities.m_as("angstrom**4/amu"), [3.15, 32.12, 18.37]
    )
    np.testing.assert_allclose(vibrations.IR_intensities.magnitude, [78.71, 67.05, 126.10])
    assert vibrations.source_mode_indices == [7, 8, 9]
    assert vibrations.frequency_method == "GFN2-xTB"
    assert vibrations.raman_activity_method == "PTB"
    np.testing.assert_allclose(
        frame.charge_spin_populations["mulliken_charges"].values, [-0.565, 0.282, 0.282]
    )
    assert vibrations[0].raman_activity_method == "PTB"
    assert vibrations.to_raman_spectrum().bands[0].metadata["raman_activity_method"] == "PTB"
    assert len(vibrations.raman_depolarization_ratios) == 0
    assert len(vibrations.raman_depolarization_ratios_polarized) == 0
    assert len(vibrations.raman_depolarization_ratios_unpolarized) == 0
    assert frame.bond_orders["wiberg_bond_order"].values[0, 1] == 0.920


@pytest.mark.parametrize("mutation", ["missing_atom", "wrong_index", "wrong_element"])
def test_xtb_charge_table_cannot_shift_values_between_atoms(mutation):
    content = (FIXTURES / "xtbout/water_raman_ptb_6_7_1.out").read_text(encoding="utf-8")
    row = "     2   1 H        0.806     0.282     0.776     1.384"
    assert row in content
    replacement = {
        "missing_atom": "",
        "wrong_index": row.replace("     2   1 H", "     3   1 H"),
        "wrong_element": row.replace("     2   1 H", "     2   8 O"),
    }[mutation]
    frame = XTBOutputFileParserMemory().parse(content.replace(row, replacement))[-1]
    assert frame.charge_spin_populations is None
    assert len(frame.vibrations) == 3


def test_fchk_raman_requires_computation_evidence_even_with_allocated_slots():
    records = parse_fchk_records(
        (FIXTURES / "g16fchk/water_raman_hf_sto3g.fchk").read_text(encoding="utf-8")
    )
    assert len(extract_fchk_vibrations(records, 3).raman_activities) == 3
    # A computed, physically zero activity remains observed when derivatives exist.
    records["Vib-E2"].value[12:21] = [0.0] * 9
    observed = extract_fchk_vibrations(records, 3)
    np.testing.assert_equal(observed.raman_activities.magnitude, np.zeros(3))
    del records["Polarizability Derivatives"]
    assert len(extract_fchk_vibrations(records, 3).raman_activities) == 0


@pytest.mark.parametrize(
    "label,field",
    [
        ("Raman Activ --", "raman_activities"),
        ("Depolar (P) --", "raman_depolarization_ratios_polarized"),
        ("Red. masses --", "reduced_masses"),
        ("Frc consts  --", "force_constants"),
        ("IR Inten    --", "IR_intensities"),
    ],
)
def test_incomplete_gaussian_raman_columns_do_not_discard_frequencies(label, field):
    path = FIXTURES / "g16log/water_raman_hf_sto3g.log"
    lines = path.read_text(encoding="utf-8").splitlines()
    lines = [line.rsplit(None, 1)[0] if line.strip().startswith(label) else line for line in lines]
    vibrations = G16LogFileParserMemory().parse("\n".join(lines))[-1].vibrations
    assert len(vibrations.frequencies) == 3
    assert len(getattr(vibrations, field)) == 0
    assert len(vibrations.raman_depolarization_ratios_unpolarized) == 3
    json.dumps(vibrations.to_unitless_dump(), allow_nan=False)


@pytest.mark.parametrize(
    "token,field",
    [
        ("8: 32.12", "raman_activities"),
        ("8: 67.05", "IR_intensities"),
        ("8:  1.67", "reduced_masses"),
    ],
)
def test_xtb_missing_index_does_not_shift_or_export_unmasked_nan(token, field):
    content = (FIXTURES / "xtbout/water_raman_ptb_6_7_1.out").read_text(encoding="utf-8")
    assert token in content
    content = content.replace(token, "")
    vibrations = XTBOutputFileParserMemory().parse(content)[-1].vibrations
    assert len(vibrations.frequencies) == 3
    assert len(getattr(vibrations, field)) == 0
    json.dumps(vibrations.to_unitless_dump(), allow_nan=False)


@pytest.mark.parametrize("mutation", ["truncated", "wrong_index", "wrong_element"])
def test_gaussian_displacements_require_complete_matching_atom_identities(mutation):
    content = (FIXTURES / "g16log/water_raman_hf_sto3g.log").read_text(encoding="utf-8")
    lines = content.splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().startswith("Atom  AN"))
    row = start + 2
    tokens = lines[row].split()
    assert tokens[:2] == ["2", "1"]
    if mutation == "truncated":
        del lines[row]
    else:
        tokens[0 if mutation == "wrong_index" else 1] = "3" if mutation == "wrong_index" else "8"
        lines[row] = " ".join(tokens)
    frame = G16LogFileParserMemory().parse("\n".join(lines))[-1]
    assert len(frame.vibrations) == 3
    assert frame.vibrations.vibration_modes == []
    assert len(frame.vibrations.IR_intensities) == 3


def test_vibration_arrays_normalize_units_before_collection_access_and_export():
    vibrations = Vibrations(
        frequencies=np.array([-10000, 20000, 30000]) * atom_ureg.Unit("1/m"),
        reduced_masses=np.array([1, 2, 3]) * atom_ureg.kg,
        force_constants=np.array([1, 2, 3]) * atom_ureg.Unit("N/m"),
        IR_intensities=np.array([1, 2, 3]) * atom_ureg.Unit("m/mol"),
        vibration_modes=[np.ones((3, 3)) * atom_ureg.bohr for _ in range(3)],
    )
    assert vibrations.frequencies.units == atom_ureg.cm_1
    np.testing.assert_allclose(vibrations.frequencies.magnitude, [-100, 200, 300])
    assert vibrations.reduced_masses.units == atom_ureg.amu
    assert vibrations.force_constants.units == atom_ureg.Unit("mdyne/angstrom")
    np.testing.assert_allclose(vibrations.IR_intensities.magnitude, [0.001, 0.002, 0.003])
    np.testing.assert_allclose(vibrations.vibration_modes[0].magnitude, 0.529177210544)
    assert vibrations.vibration_modes[0].units == atom_ureg.angstrom
    assert vibrations[0].frequency.magnitude == pytest.approx(-100)
    assert vibrations.imaginary_vibrations.frequencies.magnitude.tolist() == [-100]
    assert "frequencies (reciprocal_centimeter)" in vibrations.to_unitless_dump_with_unit_keys()


@pytest.mark.parametrize("field", ["reduced_masses", "force_constants", "IR_intensities"])
def test_vibration_scalar_arrays_cannot_silently_shift_mode_identity(field):
    with pytest.raises(ValidationError, match="must match frequencies"):
        Vibrations(
            frequencies=np.array([100, 200]) * atom_ureg.cm_1,
            **{field: getattr(Vibrations(), field).units * np.ones(1)},
        )


def test_vibration_displacements_must_match_the_chemical_frame_atom_count():
    with pytest.raises(ValidationError, match="frame atom count"):
        BaseCalcFrame(
            atoms=[8, 1, 1],
            coords=np.zeros((3, 3)) * atom_ureg.angstrom,
            vibrations=Vibrations(
                frequencies=np.array([100, 200, 300]) * atom_ureg.cm_1,
                vibration_modes=[np.ones((2, 3)) * atom_ureg.angstrom for _ in range(3)],
            ),
        )


def test_fchk_thermal_corrections_are_derived_and_export_their_basis(parse_fixture):
    log = parse_fixture("g16log/water_raman_hf_sto3g.log")
    fchk = parse_fixture("g16fchk/water_raman_hf_sto3g.fchk")
    thermal = fchk.thermal_informations
    assert thermal is not None
    for correction, total in (("TCE", "U_T"), ("TCH", "H_T"), ("TCG", "G_T")):
        actual = getattr(thermal, correction).m_as("hartree/particle")
        expected = getattr(log.thermal_informations, correction).m_as("hartree/particle")
        assert actual == pytest.approx(expected, abs=1e-6)
        assert thermal.derivations[correction] == f"{total} - Total Energy"
    assert thermal.entropy_model is None
    assert thermal.to_property_bundle().metadata["derivations"] == thermal.derivations


def test_polarized_raman_fields_survive_filtering_and_validate_mode_alignment():
    vibrations = Vibrations(
        frequencies=np.array([-100, 200]) * atom_ureg.cm_1,
        raman_depolarization_ratios_polarized=np.array([0.1, 0.2]) * atom_ureg.dimensionless,
        raman_depolarization_ratios_unpolarized=np.array([0.3, 0.4]) * atom_ureg.dimensionless,
        raman_activity_method="PTB",
    )
    imaginary = vibrations.imaginary_vibrations
    assert imaginary[0].raman_depolarization_ratio_polarized.magnitude == 0.1
    assert imaginary[0].raman_depolarization_ratio_unpolarized.magnitude == 0.3
    assert imaginary[0].raman_activity_method == "PTB"
    with pytest.raises(ValidationError, match="must match frequencies"):
        Vibrations(
            frequencies=np.array([100]) * atom_ureg.cm_1,
            raman_depolarization_ratios_polarized=np.array([0.1, 0.2]) * atom_ureg.dimensionless,
        )
