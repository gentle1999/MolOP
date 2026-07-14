from pathlib import Path
from typing import cast

import pytest

from molop.io import AutoParser


FIXTURE = Path("tests/test_files/g16log/2-TS1-Opt.log")


def _parse_fixture():
    return AutoParser(FIXTURE.as_posix(), parser_detection="g16log", n_jobs=1)[0]


def test_g16log_autoparser_exposes_stable_file_and_frame_fields() -> None:
    parsed_file = _parse_fixture()
    frame = parsed_file[-1]

    assert parsed_file.qm_software == "Gaussian"
    assert parsed_file.qm_software_version == "ES64L-G16RevC.01"
    assert parsed_file.method == frame.method == "DFT"
    assert parsed_file.basis_set == frame.basis_set == "pseudopotential"
    assert parsed_file.functional == frame.functional == "B3LYP"
    assert parsed_file.charge == frame.charge == 0
    assert parsed_file.multiplicity == frame.multiplicity == 1

    assert frame.frame_id == len(parsed_file) - 1
    assert len(frame.atoms) == 101
    assert frame.coords.shape == (101, 3)
    assert frame.standard_coords is not None
    assert frame.standard_coords.shape == (101, 3)
    assert frame.standard_orientation_transformation_matrix is not None

    assert parsed_file.status.normal_terminated is True
    assert parsed_file.is_error is False
    assert frame.status is not None
    assert frame.status.scf_converged is True
    assert frame.status.normal_terminated is None
    assert frame.is_error is False
    assert frame.is_normal is True


def test_g16log_public_qm_result_fields_are_available() -> None:
    frame = _parse_fixture()[-1]

    assert frame.energies is not None
    assert frame.energies.reference_energy is not None
    assert frame.energies.total_energy is not None
    assert frame.energies.total_energy.to("hartree").m == pytest.approx(
        frame.energies.reference_energy.to("hartree").m
    )

    assert frame.thermal_informations is not None
    assert frame.thermal_informations.ZPVE is not None
    assert frame.thermal_informations.G_T is not None
    assert frame.thermal_informations.rotational_constants is not None
    assert frame.thermal_informations.vibrational_temperatures is not None

    assert frame.vibrations is not None
    assert len(frame.vibrations.frequencies) == 297
    assert frame.vibrations.num_imaginary == 1
    assert frame.is_TS is True

    assert frame.molecular_orbitals is not None
    assert len(frame.molecular_orbitals.alpha_energies) == 1050
    assert frame.forces is not None
    assert frame.forces.shape == (101, 3)
    assert frame.forces_axis_order == ("atom", "cartesian")
    assert frame.forces_atom_order == "source"
    assert frame.forces_orientation == "unknown"
    assert frame.hessian is not None
    assert frame.hessian.shape == (303, 303)
    assert frame.hessian_axis_order == ("atom_cartesian", "atom_cartesian")
    assert frame.hessian_atom_order == "source"
    assert frame.hessian_orientation == "unknown"
    assert frame.polarizability is not None
    assert frame.polarizability.dipole is not None


def test_g16log_force_and_hessian_conventions_are_machine_readable() -> None:
    frame = _parse_fixture()[-1]

    assert frame.forces is not None
    assert frame.forces_axis_order == ("atom", "cartesian")
    assert frame.forces_atom_order == "source"
    assert frame.forces_orientation == "unknown"
    assert frame.hessian is not None
    assert frame.hessian_axis_order == ("atom_cartesian", "atom_cartesian")
    assert frame.hessian_atom_order == "source"
    assert frame.hessian_orientation == "unknown"


def test_g16log_internal_component_state_is_not_serialized_as_public_model_data() -> None:
    frame = _parse_fixture()[-1]
    dumped = frame.model_dump()

    assert "component_tree" not in dumped
    assert "_component_tree" not in dumped
    assert "component_tree_input" not in type(frame).model_fields


def test_g16log_summary_dataframe_contract_for_user_entrypoints() -> None:
    parsed_file = _parse_fixture()
    default_df = parsed_file.to_summary_df()
    brief_df = parsed_file.to_summary_df(frame="all")
    full_df = parsed_file.to_summary_df(frame="all", brief=False)

    assert len(default_df) == 1
    assert default_df[("General", "FrameID", "")].tolist() == [len(parsed_file) - 1]
    assert len(brief_df) == len(parsed_file)
    assert ("General", "FrameID", "") in brief_df.columns
    assert ("Calc Parameter", "Software", "") in brief_df.columns
    assert ("Status", "IsError", "") in brief_df.columns
    assert ("Energy", "total_energy", "hartree") not in brief_df.columns

    assert ("Energy", "reference_energy", "hartree") in full_df.columns
    assert ("Energy", "total_energy", "hartree") in full_df.columns
    assert ("Thermal", "G_T", "kilocalorie / mole") in full_df.columns
    assert ("Vibration", "num_imaginary", "") in full_df.columns
    assert brief_df.to_json(orient="records")

    batch = AutoParser(FIXTURE.as_posix(), parser_detection="g16log", n_jobs=1)
    batch_df = batch.to_summary_df(frame="all", n_jobs=1, flatten_columns=True)
    assert len(batch_df) == len(parsed_file)
    assert "General.FrameID" in batch_df.columns
    assert "Status.IsError" in batch_df.columns


def test_g16log_summary_keeps_ts_smiles_when_rdkit_canon_smiles_cannot_reparse() -> None:
    parsed_file = AutoParser(
        "tests/test_files/g16log/2-modfreq.log",
        parser_detection="g16log",
        n_jobs=1,
        only_last_frame=True,
    )[0]

    summary = parsed_file[-1].to_summary_dict(brief=False)

    assert parsed_file[-1].is_TS is True
    assert summary[("General", "PreCanonicalSMILES", "")]
    assert summary[("General", "PostCanonicalSMILES", "")]


def test_g16log_fakeg_file_transform_is_file_level_and_reparseable() -> None:
    parsed_file = _parse_fixture()

    rendered = cast(str, parsed_file.format_transform("fakeg"))
    rendered_all = cast(str, parsed_file.format_transform("fakeg", frame="all"))

    assert isinstance(rendered, str)
    assert "Standard orientation:" in rendered
    assert "SCF Done:" in rendered
    assert "Frequencies --" in rendered
    assert "Zero-point correction=" in rendered
    assert "- Thermochemistry -" in rendered
    assert "Normal termination of Gaussian" not in rendered

    assert rendered_all.count("Standard orientation:") > rendered.count("Standard orientation:")
    assert rendered_all.count("Normal termination of Gaussian") == 2
    assert "Job cpu time:" not in rendered_all
