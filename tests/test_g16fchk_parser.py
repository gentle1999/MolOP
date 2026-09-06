from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.codec_exceptions import FormatMismatchError
from molop.io.logic.gaussian.fchk.output.models.G16FchkFile import G16FchkFileDisk
from molop.io.logic.gaussian.fchk.output.parsers._fchk_extractors import (
    extract_fchk_atomic_masses,
    extract_fchk_populations,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_records import (
    FCHKRecord,
    parse_fchk_header,
    parse_fchk_records,
)
from molop.io.logic.gaussian.fchk.output.parsers.G16FchkFileParser import (
    G16FchkFileParserDisk,
    G16FchkFileParserMemory,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserDisk
from molop.unit import atom_ureg


FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "g16fchk"
LOG_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "g16log"


def test_g16fchk_all_maintained_fixtures_parse() -> None:
    fixture_paths = sorted(FIXTURE_DIR.glob("*.fchk"))

    assert len(fixture_paths) >= 13
    for fixture_path in fixture_paths:
        parsed = G16FchkFileParserDisk().parse(str(fixture_path), release_file_content=True)

        assert len(parsed.frames) == 1, fixture_path.name
        assert parsed.qm_software == "Gaussian", fixture_path.name
        assert parsed.frames[0].atoms, fixture_path.name
        assert parsed.frames[0].energies is not None, fixture_path.name
        assert parsed.frames[0].energies.total_energy is not None, fixture_path.name


def test_g16fchk_frequency_fixture_exposes_structured_results() -> None:
    parsed = G16FchkFileParserDisk(capture_source_evidence=True).parse(
        str(FIXTURE_DIR / "dsgdb9nsd_000007-6.fchk"),
        release_file_content=False,
    )
    frame = parsed.frames[0]

    assert parsed.qm_software_version == "ES64L-G16RevC.01"
    assert parsed.title_card == "dsgdb9nsd_000007-6"
    assert parsed.method == "DFT"
    assert parsed.functional == "B3LYP-GD3BJ"
    assert parsed.basis_set == "6-311+G(d,p)"
    assert parsed.charge == 0
    assert parsed.multiplicity == 2
    assert [task.task_type for task in parsed.task_requests] == ["opt", "freq"]
    assert frame.atoms == [6, 6, 1, 1, 1, 1, 1]
    assert frame.coords.shape == (7, 3)
    assert frame.atomic_masses_source == "gaussian_fchk_real_atomic_weights"
    assert frame.atomic_masses is not None
    np.testing.assert_allclose(
        frame.atomic_masses.m_as("amu"),
        [12.0, 12.0, 1.00782504, 1.00782504, 1.00782504, 1.00782504, 1.00782504],
    )
    assert frame.coordinate_source == "observed"
    assert frame.energies is not None
    assert frame.energies.total_energy.m_as("hartree") == pytest.approx(-79.1896514864202)
    assert {item.source_label for item in frame.energies.observations} >= {
        "SCF Energy",
        "Total Energy",
    }
    assert frame.forces is not None
    assert frame.forces.shape == (7, 3)
    assert frame.force_source_field == "Cartesian Gradient"
    assert frame.hessian is not None
    assert frame.hessian.shape == (21, 21)
    assert np.allclose(frame.hessian.magnitude, frame.hessian.magnitude.T)
    assert frame.molecular_orbitals is not None
    assert len(frame.molecular_orbitals.alpha_energies) == 74
    assert len(frame.molecular_orbitals.beta_energies) == 74
    assert frame.charge_spin_populations is not None
    assert len(frame.charge_spin_populations["mulliken_charges"].values) == 7
    assert frame.total_spin is not None
    assert frame.total_spin.spin_square == pytest.approx(0.7538265141283933)
    assert frame.polarizability is not None
    assert frame.polarizability.dipole.shape == (3,)
    assert frame.polarizability.polarizability_tensor.shape == (6,)
    assert frame.polarizability.quadrupole.shape == (6,)
    assert frame.vibrations is not None
    assert len(frame.vibrations) == 15
    assert frame.vibrations.frequencies[0].m_as("cm^-1") == pytest.approx(106.305208)
    assert frame.vibrations.reduced_masses[0].m_as("amu") == pytest.approx(1.00815858)
    assert len(frame.vibrations.vibration_modes) == 15
    assert frame.vibrations.vibration_modes[0].shape == (7, 3)
    assert frame.thermal_informations is not None
    assert frame.thermal_informations.U_T.m_as("hartree/particle") == pytest.approx(
        -79.12666220220767
    )
    assert frame.status is not None
    assert frame.status.scf_converged is True
    assert frame.status.normal_terminated is True
    assert frame.geometry_optimization_status is not None
    assert frame.geometry_optimization_status.geometry_optimized is True


def test_g16fchk_correlated_energy_fields_are_preserved() -> None:
    parsed = G16FchkFileParserDisk().parse(str(FIXTURE_DIR / "H_ccsd_t.fchk"))
    energies = parsed.frames[0].energies

    assert parsed.method == "CCSD(T)"
    assert parsed.model_chemistry.spin_treatment == "U"
    assert energies is not None
    assert energies.reference_energy.m_as("hartree") == pytest.approx(-0.499821176023958)
    assert energies.mp2_energy.m_as("hartree") == pytest.approx(-0.499821176023958)
    assert energies.mp3_energy.m_as("hartree") == pytest.approx(-0.499821176023958)
    assert energies.mp4_energy.m_as("hartree") == pytest.approx(-0.499821176023958)
    assert energies.ccsd_energy.m_as("hartree") == pytest.approx(-0.499821176023958)
    assert energies.ccsd_t_energy.m_as("hartree") == pytest.approx(-0.499821176023958)


def test_g16fchk_npa_population_is_structured_when_present() -> None:
    parsed = G16FchkFileParserDisk().parse(str(FIXTURE_DIR / "radical_0554_opt_g16_nbo_sp.fchk"))
    populations = parsed.frames[0].charge_spin_populations

    assert populations is not None
    assert len(populations["mulliken_charges"].values) == 27
    assert len(populations["npa_charges"].values) == 27


def test_g16fchk_population_records_support_spin_and_extensible_schemes() -> None:
    records = {
        label: FCHKRecord(label=label, data_type="R", value=values, count=2)
        for label, values in {
            "Mulliken Charges": [-0.2, 0.2],
            "Mulliken Spin Densities": [0.6, 0.4],
            "ESP Charges": [-0.3, 0.3],
            "NPA Charges": [-0.1, 0.1],
            "NPA Spins": [0.7, 0.3],
        }.items()
    }

    populations = extract_fchk_populations(records, num_atoms=2)

    assert populations is not None
    assert populations["mulliken_charges"].values == [-0.2, 0.2]
    assert populations["mulliken_spins"].values == [0.6, 0.4]
    assert populations["npa_charges"].values == [-0.1, 0.1]
    assert populations["esp_charges"].values == [-0.3, 0.3]
    assert populations["npa_spins"].values == [0.7, 0.3]


def test_g16fchk_atomic_mass_extractor_requires_one_mass_per_atom() -> None:
    records = {
        "Real atomic weights": FCHKRecord(
            label="Real atomic weights",
            data_type="R",
            value=[12.0, 1.00782504],
            count=2,
        )
    }

    masses = extract_fchk_atomic_masses(records, num_atoms=2)

    assert masses is not None
    assert masses.units == atom_ureg.amu
    np.testing.assert_allclose(masses.magnitude, [12.0, 1.00782504])
    assert extract_fchk_atomic_masses(records, num_atoms=1) is None
    assert extract_fchk_atomic_masses({}, num_atoms=2) is None


@pytest.mark.parametrize(
    ("fchk_name", "log_name", "atom_count"),
    [
        ("molecule.fchk", "test_nmr.log", 5),
        ("ethanol.fchk", "test_nmr_coupling.log", 9),
    ],
)
def test_g16fchk_nmr_shielding_matches_corresponding_log(
    fchk_name: str, log_name: str, atom_count: int
) -> None:
    fchk_frame = G16FchkFileParserDisk().parse(str(FIXTURE_DIR / fchk_name)).frames[0]
    log_frame = G16LogFileParserDisk().parse(str(LOG_FIXTURE_DIR / log_name)).frames[-1]

    assert fchk_frame.nmr is not None
    assert log_frame.nmr is not None
    assert fchk_frame.nmr.gauge == log_frame.nmr.gauge == "GIAO"
    assert len(fchk_frame.nmr.shielding_tensors) == atom_count
    assert len(log_frame.nmr.shielding_tensors) == atom_count
    for fchk_shielding, log_shielding in zip(
        fchk_frame.nmr.shielding_tensors,
        log_frame.nmr.shielding_tensors,
        strict=True,
    ):
        assert fchk_shielding.atom_index == log_shielding.atom_index
        assert fchk_shielding.atom_symbol == log_shielding.atom_symbol
        np.testing.assert_allclose(
            fchk_shielding.shielding_tensor.m_as("ppm"),
            log_shielding.shielding_tensor.m_as("ppm"),
            atol=6e-5,
        )
        assert fchk_shielding.isotropic is not None
        assert log_shielding.isotropic is not None
        assert fchk_shielding.isotropic.m_as("ppm") == pytest.approx(
            log_shielding.isotropic.m_as("ppm"), abs=6e-5
        )
        assert fchk_shielding.anisotropy is not None
        assert log_shielding.anisotropy is not None
        assert fchk_shielding.anisotropy.m_as("ppm") == pytest.approx(
            log_shielding.anisotropy.m_as("ppm"), abs=6e-5
        )
        assert fchk_shielding.principal_values is not None
        assert log_shielding.principal_values is not None
        np.testing.assert_allclose(
            fchk_shielding.principal_values.m_as("ppm"),
            log_shielding.principal_values.m_as("ppm"),
            atol=6e-5,
        )


def test_g16fchk_nmr_spin_spin_components_reconstruct_total_k() -> None:
    frame = G16FchkFileParserDisk().parse(str(FIXTURE_DIR / "ethanol.fchk")).frames[0]
    log_frame = (
        G16LogFileParserDisk().parse(str(LOG_FIXTURE_DIR / "test_nmr_coupling.log")).frames[-1]
    )

    assert frame.nmr is not None
    assert log_frame.nmr is not None
    assert frame.nmr.coupling_atom_indices == list(range(9))
    assert frame.nmr.spin_spin_coupling_k is not None
    assert log_frame.nmr.spin_spin_coupling_k is not None
    assert set(frame.nmr.spin_spin_coupling_k_components) == {"FC", "SD", "PSO", "DSO"}
    assert set(log_frame.nmr.spin_spin_coupling_k_components) == {"FC", "SD", "PSO", "DSO"}
    for component_name in ("FC", "SD", "PSO", "DSO"):
        np.testing.assert_allclose(
            frame.nmr.spin_spin_coupling_k_components[component_name].m_as("Hz"),
            log_frame.nmr.spin_spin_coupling_k_components[component_name].m_as("Hz"),
            atol=5e-5,
        )
    np.testing.assert_allclose(
        sum(
            component.m_as("Hz") for component in frame.nmr.spin_spin_coupling_k_components.values()
        ),
        frame.nmr.spin_spin_coupling_k.m_as("Hz"),
    )
    np.testing.assert_allclose(
        frame.nmr.spin_spin_coupling_k.m_as("Hz"),
        log_frame.nmr.spin_spin_coupling_k.m_as("Hz"),
        atol=5e-5,
    )
    assert frame.nmr.spin_spin_coupling_j is None
    assert frame.nmr.spin_spin_coupling_j_components == {}


def test_g16fchk_only_extract_structure_skips_property_arrays() -> None:
    parsed = G16FchkFileParserDisk(only_extract_structure=True).parse(
        str(FIXTURE_DIR / "dsgdb9nsd_000007-6.fchk")
    )
    frame = parsed.frames[0]

    assert frame.atoms == [6, 6, 1, 1, 1, 1, 1]
    assert frame.coords.shape == (7, 3)
    assert frame.energies is None
    assert frame.molecular_orbitals is None
    assert frame.vibrations is None


def test_g16fchk_autoparser_preserves_single_source_span() -> None:
    fixture_path = FIXTURE_DIR / "dsgdb9nsd_000001-3.fchk"
    source = fixture_path.read_bytes().decode("utf-8")
    batch = AutoParser(str(fixture_path), capture_source_evidence=True)
    parsed = cast(G16FchkFileDisk, batch[0])

    assert parsed.detected_format_id == "g16fchk"
    assert parsed.source_format == "g16fchk"
    assert parsed.source_segments[0].source_span.start_char == 0
    assert parsed.source_segments[0].source_span.end_char == len(source)
    assert parsed.frames[0].source_span == parsed.source_segments[0].source_span


def test_g16fchk_record_decoder_handles_fixed_width_character_arrays() -> None:
    source = (FIXTURE_DIR / "H_anion.fchk").read_text(encoding="utf-8")
    header = parse_fchk_header(source)
    records = parse_fchk_records(source, wanted_labels={"Full Title", "Route", "Gaussian Version"})

    assert header.job_type == "Freq"
    assert header.method == "RMP2-FC"
    assert header.basis_set == "def2SVP"
    assert records["Full Title"].value == "Title: H_anion"
    assert "scrf(SMD,solvent=ethanol)" in str(records["Route"].value)
    assert records["Gaussian Version"].value == "ES64L-G16RevC.01"


def test_g16fchk_probe_does_not_claim_gaussian_log() -> None:
    log_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"

    assert G16FchkFileParserDisk.probe_file_format(log_path) is False
    with pytest.raises(FormatMismatchError):
        G16FchkFileParserMemory().parse(log_path.read_text(encoding="utf-8"))
