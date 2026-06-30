from __future__ import annotations

from pathlib import Path

import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import ORCAInpFileFrameDisk


ORCA_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "orca"
ORCA_SINGLE_POINT_INPUT_DIR = ORCA_FIXTURE_DIR / "single_point_inputs"
ORCA_SCF_STABILITY_INPUT_DIR = ORCA_FIXTURE_DIR / "scf_stability_inputs"
ORCA_OPTIMIZATION_INPUT_DIR = ORCA_FIXTURE_DIR / "optimization_inputs"
ORCA_FREQUENCY_INPUT_DIR = ORCA_FIXTURE_DIR / "frequency_inputs"
ORCA_EXCITED_STATE_INPUT_DIR = ORCA_FIXTURE_DIR / "excited_state_inputs"
ORCA_MRCI_INPUT_DIR = ORCA_FIXTURE_DIR / "mrci_inputs"
ORCA_SOLVATOR_INPUT_DIR = ORCA_FIXTURE_DIR / "solvator_inputs"
ORCA_INPUT_FIXTURE_DIRS = (
    ORCA_SINGLE_POINT_INPUT_DIR,
    ORCA_SCF_STABILITY_INPUT_DIR,
    ORCA_OPTIMIZATION_INPUT_DIR,
    ORCA_FREQUENCY_INPUT_DIR,
    ORCA_EXCITED_STATE_INPUT_DIR,
    ORCA_MRCI_INPUT_DIR,
    ORCA_SOLVATOR_INPUT_DIR,
)


@pytest.mark.parametrize(
    "fixture_path",
    [path for fixture_dir in ORCA_INPUT_FIXTURE_DIRS for path in sorted(fixture_dir.glob("*.inp"))],
)
def test_orcainp_fixture_parse_structured_frame(fixture_path: Path) -> None:
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]
    assert len(file_model) > 0

    frame = file_model[-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.qm_software == "ORCA"
    assert frame.qm_software_version == "Any"
    assert frame.keywords or frame.blocks
    assert frame.method or frame.basis_set or frame.functional or frame.blocks
    assert frame.task_requests
    assert frame.model_chemistry.raw_keywords == frame.keywords
    if frame.resources_raw:
        assert frame.resource_request.raw == frame.resources_raw
    if "D3BJ" in frame.keywords:
        assert frame.dispersion_correction == "D3BJ"
        assert frame.model_chemistry.dispersion_correction == "D3BJ"
    assert frame.keyword_lines or frame.blocks
    assert frame.geometry is not None
    if "Print[" in frame.resources_raw:
        assert frame.output_print_settings
    if frame.geometry.ctype in {"xyz", "cart", "cartesian"}:
        if "..." not in fixture_path.read_text(encoding="utf-8"):
            assert len(frame.atoms) > 0
            assert tuple(frame.coords.shape) == (len(frame.atoms), 3)
    elif frame.geometry.ctype in {"int", "internal", "gzmt"}:
        assert frame.geometry.atoms
        assert frame.atoms
        if frame.geometry.internal_coords is not None:
            assert tuple(frame.coords.shape) == (len(frame.atoms), 3)
    else:
        assert frame.geometry.external_path


def test_orcainp_single_point_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(
        ORCA_SINGLE_POINT_INPUT_DIR.glob("orca6_energygradients_*.inp")
    )
    assert len(manual_fixtures) == 66


def test_orcainp_scf_stability_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(
        ORCA_SCF_STABILITY_INPUT_DIR.glob("orca6_scfstability_*.inp")
    )
    assert len(manual_fixtures) == 1


def test_orcainp_optimization_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(
        ORCA_OPTIMIZATION_INPUT_DIR.glob("orca6_optimizations_*.inp")
    )
    assert len(manual_fixtures) == 24


def test_orcainp_frequency_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(ORCA_FREQUENCY_INPUT_DIR.glob("orca6_frequencies_*.inp"))
    assert len(manual_fixtures) == 2


def test_orcainp_excited_state_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(ORCA_EXCITED_STATE_INPUT_DIR.glob("*.inp"))
    assert len(manual_fixtures) == 13


def test_orcainp_mrci_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(ORCA_MRCI_INPUT_DIR.glob("*.inp"))
    assert len(manual_fixtures) == 7


def test_orcainp_solvator_inputs_include_orca6_manual_fixtures() -> None:
    manual_fixtures = sorted(ORCA_SOLVATOR_INPUT_DIR.glob("*.inp"))
    assert len(manual_fixtures) == 1


@pytest.mark.parametrize(
    ("fixture_name", "expected"),
    [
        (
            "tddft_triplets_int.inp",
            {
                "method": "DFT",
                "functional": "B3LYP",
                "family": "TDDFT",
                "nroots": 10,
                "triplets": True,
            },
        ),
        (
            "tddft_followiroot_xyz.inp",
            {
                "method": "DFT",
                "functional": "wB97X",
                "family": "TDDFT",
                "nroots": 5,
                "iroot": 3,
                "followiroot": True,
            },
        ),
        (
            "cis_iroot1_int.inp",
            {
                "method": "CIS",
                "basis_set": "DEF2-SVP",
                "family": "CIS",
                "nroots": 1,
                "iroot": 1,
                "reference_method": "HF",
            },
        ),
        (
            "rocis_xyz.inp",
            {
                "method": "ROCIS",
                "basis_set": "SVP",
                "family": "ROCIS",
                "nroots": 2,
            },
        ),
        (
            "mcrpa_int.inp",
            {
                "method": "CASSCF",
                "basis_set": "DEF2-SVP",
                "family": "MCRPA",
                "reference_method": "CASSCF",
                "nroots": 8,
            },
        ),
        (
            "mdci_eomccsd_xyz.inp",
            {
                "method": "EOM-CCSD",
                "basis_set": "cc-pVDZ",
                "family": "EOM-CCSD",
                "reference_method": "RHF",
                "nroots": 9,
            },
        ),
        (
            "mdci_ipeomccsd_xyz.inp",
            {
                "method": "IP-EOM-CCSD",
                "basis_set": "cc-pVDZ",
                "family": "IP-EOM-CCSD",
                "doalpha": True,
                "nroots": 7,
            },
        ),
        (
            "mdci_adc2_xyz.inp",
            {
                "method": "ADC2",
                "basis_set": "cc-pVDZ",
                "family": "ADC2",
                "nroots": 9,
            },
        ),
        (
            "mdci_steomccsd_xyz.inp",
            {
                "method": "STEOM-CCSD",
                "basis_set": "cc-pVDZ",
                "family": "STEOM-CCSD",
                "do_dbfilter": True,
                "nroots": 9,
            },
        ),
        (
            "mdci_ihfsmrccsd_xyz.inp",
            {
                "method": "IH-FSMR-CCSD",
                "basis_set": "cc-pVDZ",
                "family": "IH-FSMR-CCSD",
                "nroots": 6,
            },
        ),
        (
            "mdci_btpno_eomccsd_xyz.inp",
            {
                "method": "BT-PNO-EOM-CCSD",
                "basis_set": "def2-TZVP",
                "family": "BT-PNO-EOM-CCSD",
                "nroots": 9,
            },
        ),
        (
            "mdci_steom_dlpno_xyz.inp",
            {
                "method": "STEOM-DLPNO-CCSD",
                "basis_set": "def2-TZVP",
                "family": "STEOM-DLPNO-CCSD",
                "rootwise": True,
                "do_store_steom": True,
                "add_l2_term": True,
                "nroots": 6,
            },
        ),
    ],
)
def test_orcainp_excited_state_manual_fixtures_are_structured(
    fixture_name: str, expected: dict[str, object]
) -> None:
    fixture_path = ORCA_EXCITED_STATE_INPUT_DIR / fixture_name
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)

    dump = frame.excited_state_semantic.model_dump()
    assert dump["enabled"] is True
    assert dump["family"] == expected["family"]
    assert len(frame.excited_state_requests) == 1
    request = frame.excited_state_requests[0]
    assert request.enabled is True
    assert request.family == expected["family"]
    assert any(task.task_type == "excited_state" for task in frame.task_requests)
    assert frame.method == expected["method"]
    if "functional" in expected:
        assert frame.functional == expected["functional"]
    if "basis_set" in expected:
        assert frame.basis_set == expected["basis_set"]
    if "reference_method" in expected:
        assert dump["reference_method"] == expected["reference_method"]
        assert request.reference_method == expected["reference_method"]
    if "nroots" in expected:
        assert dump["nroots"] == expected["nroots"]
        assert request.nroots == expected["nroots"]
    for key in (
        "iroot",
        "followiroot",
        "triplets",
        "doalpha",
        "do_dbfilter",
        "rootwise",
        "do_store_steom",
        "add_l2_term",
    ):
        if key in expected:
            assert dump[key] == expected[key]
    if "iroot" in expected:
        assert request.root == expected["iroot"]
    if "followiroot" in expected:
        assert request.follow_root == expected["followiroot"]
    if "triplets" in expected:
        assert request.triplets == expected["triplets"]


@pytest.mark.parametrize(
    ("fixture_name", "expected"),
    [
        (
            "mrci_simple_int.inp",
            {
                "method": "MRCI",
                "basis_set": "def2-SVP",
                "reference_method": None,
                "ci_type": "MRCI",
                "ewin": (-3.0, 1000.0),
                "tsel": 1e-6,
                "tpre": 1e-2,
                "solver": "Diag",
                "int_mode": "FullTrafo",
                "use_ivos": True,
                "all_singles": True,
                "new_blocks": 3,
            },
        ),
        (
            "sorci_h2co_moread.inp",
            {
                "method": "SORCI",
                "reference_method": "CASSCF",
                "ci_type": "SORCI",
                "tsel": 1e-6,
                "tpre": 1e-4,
                "tnat": 1e-5,
                "int_mode": "FullTrafo",
                "all_singles": True,
                "do_nat_orbs": True,
                "new_blocks": 3,
            },
        ),
        (
            "mracpf2a_doddcimp2_xyz.inp",
            {
                "method": "MRACPF2a",
                "basis_set": "def2-SVP",
                "reference_method": "ROHF",
                "ci_type": "MRACPF2a",
                "ewin": (-3.0, 1000.0),
                "tsel": 0.0,
                "solver": "DIIS",
                "int_mode": "FullTrafo",
                "use_ivos": True,
                "do_ddcimp2": True,
                "new_blocks": 2,
            },
        ),
        (
            "mracpf_correlation_energy_multijob.inp",
            {
                "method": "MRACPF",
                "basis_set": "aug-SVP",
                "reference_method": "CASSCF",
                "ci_type": None,
                "tsel": 1e-8,
                "tpre": 1e-6,
                "new_blocks": 0,
            },
        ),
        (
            "mrciq_hf_scan_xyz.inp",
            {
                "method": "MRCI+Q",
                "reference_method": "CASSCF",
                "ci_type": None,
                "tsel": 1e-8,
                "tpre": 1e-5,
                "new_blocks": 0,
            },
        ),
        (
            "ri_mrmp2_stilbene_scan_int.inp",
            {
                "method": "MRMP2",
                "basis_set": "def2-TZVP",
                "reference_method": "CASSCF",
                "ci_type": None,
                "tsel": 1e-8,
                "new_blocks": 0,
            },
        ),
    ],
)
def test_orcainp_mrci_manual_fixtures_are_structured(
    fixture_name: str, expected: dict[str, object]
) -> None:
    fixture_path = ORCA_MRCI_INPUT_DIR / fixture_name
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)

    dump = frame.multi_reference_semantic.model_dump()
    assert dump["enabled"] is True
    assert frame.method == expected["method"]
    assert len(frame.multireference_requests) == 1
    request = frame.multireference_requests[0]
    assert request.enabled is True
    assert request.method == expected["method"]
    assert any(task.task_type == "multi_reference" for task in frame.task_requests)
    if "basis_set" in expected:
        assert frame.basis_set == expected["basis_set"]
    for key in (
        "reference_method",
        "ci_type",
        "ewin",
        "tsel",
        "tpre",
        "tnat",
        "solver",
        "int_mode",
        "use_ivos",
        "all_singles",
        "do_nat_orbs",
        "do_ddcimp2",
    ):
        if key in expected:
            assert dump[key] == expected[key]
    assert request.reference_method == expected["reference_method"]
    assert request.ci_type == expected["ci_type"]
    for key in ("tsel", "tpre", "tnat"):
        if key in expected:
            assert request.thresholds[key] == expected[key]
    assert len(dump["new_blocks"]) == expected["new_blocks"]
    assert len(request.state_blocks) == expected["new_blocks"]
    assert frame.excited_state_semantic.enabled is False


def test_orcainp_mrci_multijob_manual_fixture_splits_frames_and_structures_last_job() -> None:
    fixture_path = ORCA_MRCI_INPUT_DIR / "o2_mrci_multijob.inp"
    batch = AutoParser(str(fixture_path), parser_detection="orcainp")
    file_model = batch[0]
    assert len(file_model) == 2

    first_frame = file_model[0]
    second_frame = file_model[1]
    assert isinstance(first_frame, ORCAInpFileFrameDisk)
    assert isinstance(second_frame, ORCAInpFileFrameDisk)

    assert first_frame.method == "RI-MP2"
    assert second_frame.method == "MRCI"
    assert second_frame.excited_state_semantic.enabled is False
    assert len(second_frame.multireference_requests) == 1
    assert second_frame.multireference_requests[0].method == "MRCI"
    assert len(second_frame.multireference_requests[0].state_blocks) == 3
    dump = second_frame.multi_reference_semantic.model_dump()
    assert dump["enabled"] is True
    assert dump["ci_type"] == "MRCI"
    assert dump["reference_method"] == "CASSCF"
    assert dump["tsel"] == 1e-7
    assert dump["tpre"] == 1e-5
    assert len(dump["new_blocks"]) == 3
    assert dump["new_blocks"][0]["multiplicity"] == 3
    assert dump["new_blocks"][0]["irrep"] == "1"
    assert dump["new_blocks"][0]["refs"] == "cas(8,6)"


def test_orcainp_mrci_manual_scan_fixture_keeps_coordinate_parameters() -> None:
    fixture_path = ORCA_MRCI_INPUT_DIR / "ri_mrmp2_stilbene_scan_int.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.coordinate_parameters is not None
    dump = frame.geometry.coordinate_parameters.model_dump()
    assert dump["items"] == [
        {
            "name": "DIHED",
            "raw_value": "90,270, 19",
            "start": 90.0,
            "stop": 270.0,
            "steps": 19,
        }
    ]


def test_orcainp_mrci_manual_xyz_scan_fixture_uses_first_parameter_value_for_frame_coords() -> None:
    fixture_path = ORCA_MRCI_INPUT_DIR / "mrciq_hf_scan_xyz.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.coordinate_parameters is not None
    assert frame.geometry.coordinate_parameters[0].name == "R"
    assert frame.geometry.coordinate_parameters[0].start == 0.85
    assert tuple(frame.coords.shape) == (2, 3)
    assert frame.coords.magnitude[1, 2] == pytest.approx(0.85)


def test_orcainp_orca6_manual_internal_coords_fill_frame_atoms_and_coords() -> None:
    fixture_path = (
        ORCA_SINGLE_POINT_INPUT_DIR
        / "orca6_energygradients_017_mp2_def2_tzvp_tightscf.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "int"
    assert frame.geometry.internal_coords is not None
    assert frame.atoms == [6, 8, 1, 1]
    assert tuple(frame.coords.shape) == (4, 3)
    assert frame.method == "MP2"
    assert frame.basis_set == "def2-TZVP"


def test_orcainp_orca6_manual_output_print_without_equals_is_structured() -> None:
    fixture_path = ORCA_SINGLE_POINT_INPUT_DIR / "orca6_energygradients_090_rhf_sv_p.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.output_print_settings[0].target == "P_ReducedOrbPopMO_L"
    assert frame.output_print_settings[0].value == "1"
    assert frame.method == "RHF"
    assert frame.basis_set == "SV(P)"


def test_orcainp_orca6_manual_dispersion_is_functional_suffix() -> None:
    fixture_path = ORCA_SINGLE_POINT_INPUT_DIR / "orca6_energygradients_078_blyp_d3_def2_qzvpp_opt.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.method == "DFT"
    assert frame.functional == "BLYP-D3"
    assert frame.dispersion_correction == "D3"
    assert frame.basis_set == "def2-QZVPP"
    assert frame.model_chemistry.functional == "BLYP-D3"
    assert frame.model_chemistry.dispersion_correction == "D3"
    assert frame.model_chemistry.basis_set == "def2-QZVPP"


def test_orcainp_orca6_manual_dlpno_double_hybrid_dispersion_is_functional_suffix() -> None:
    fixture_path = (
        ORCA_SINGLE_POINT_INPUT_DIR
        / "orca6_energygradients_034_dlpno_b2plyp_d3_normalpno_def2_tzvp_def2_tzvp_c_.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.method == "DFT"
    assert frame.functional == "DLPNO-B2PLYP-D3"
    assert frame.dispersion_correction == "D3"
    assert frame.basis_set == "def2-TZVP"
    assert frame.model_chemistry.functional == "DLPNO-B2PLYP-D3"
    assert frame.model_chemistry.dispersion_correction == "D3"


def test_orcainp_orca6_scf_stability_block_is_structured() -> None:
    fixture_path = (
        ORCA_SCF_STABILITY_INPUT_DIR
        / "orca6_scfstability_000_bhlyp_def2_svp_nori.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.method == "DFT"
    assert frame.functional == "BHLYP"
    assert frame.basis_set == "def2-SVP"
    assert frame.geometry is not None
    assert frame.geometry.ctype == "xyz"
    assert frame.atoms == [1, 1]
    assert tuple(frame.coords.shape) == (2, 3)
    assert [block.name for block in frame.blocks] == ["scf"]
    assert [line.text.strip() for line in frame.blocks[0].lines] == [
        "guess hcore",
        "HFTyp UHF",
        "STABPerform true",
    ]
    assert "STABPerform true" in frame.resources_raw


def test_orcainp_orca6_optimization_constraints_block_is_not_truncated() -> None:
    fixture_path = (
        ORCA_OPTIMIZATION_INPUT_DIR
        / "orca6_optimizations_010_rks_b3lyp_g_sv_p_opt.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "int"
    assert [block.name for block in frame.blocks] == ["geom"]
    assert "{ B 0 1 1.25 C }" in frame.resources_raw
    assert "{ A 2 0 3 120.0 C }" in frame.resources_raw
    assert frame.blocks[0].raw_text.strip().endswith("end")


def test_orcainp_orca6_optimization_hess_internal_block_is_not_truncated() -> None:
    fixture_path = (
        ORCA_OPTIMIZATION_INPUT_DIR
        / "orca6_optimizations_020_geomopt_hess_internal_scanname.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "xyz"
    assert [block.name for block in frame.blocks] == ["geom"]
    assert 'XYZ1 "scanName.003.xyz"' in frame.resources_raw
    assert 'GBW2 "ScanName.005.xyz"' in frame.resources_raw
    assert "Update BFGS" in frame.resources_raw
    assert "ProjectTR false" in frame.resources_raw


def test_orcainp_orca6_optimization_compound_block_keeps_nested_steps() -> None:
    fixture_path = (
        ORCA_OPTIMIZATION_INPUT_DIR
        / "orca6_optimizations_007_compound_two_step_opt.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "int"
    assert frame.keyword_lines == []
    assert [block.name for block in frame.blocks] == ["compound"]
    assert "! AM1 NumFreq" in frame.resources_raw
    assert "!B3LYP def2-svp def2/J Opt" in frame.resources_raw
    assert "%geom" in frame.resources_raw
    assert 'InHessName "two_step_opt_Compound_1.hess"' in frame.resources_raw
    assert frame.blocks[0].raw_text.strip().endswith("End")


def test_orcainp_orca6_optimization_scan_block_keeps_all_scan_lines() -> None:
    fixture_path = (
        ORCA_OPTIMIZATION_INPUT_DIR / "orca6_optimizations_027_b3lyp_g_sv_p_opt.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "int"
    assert frame.resources_raw.count("will be scanned") == 3
    assert "B 0 1 = 1.35, 1.10, 12" in frame.resources_raw
    assert "B 0 2 = 1.20, 1.00, 5" in frame.resources_raw
    assert "A 2 0 1 = 140, 100, 5" in frame.resources_raw


def test_orcainp_orca6_optimization_pdbfile_geometry_is_structured() -> None:
    fixture_path = ORCA_OPTIMIZATION_INPUT_DIR / "orca6_optimizations_055_l_opt.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "pdbfile"
    assert frame.geometry.external_path == "CHMH.pdb"
    assert frame.charge == 0
    assert frame.multiplicity == 1
    assert frame.atoms == []


def test_orcainp_orca6_optimization_fragment_mixed_basis_is_structured() -> None:
    fixture_path = (
        ORCA_OPTIMIZATION_INPUT_DIR
        / "orca6_optimizations_016_fragment_constraints_bp86.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "xyz"
    assert frame.charge == 1
    assert frame.multiplicity == 2
    assert len(frame.geometry.atoms) == 80
    assert len(frame.atoms) == 80
    assert frame.has_mixed_basis
    assert frame.geometry.atoms[0].symbol == "Fe"
    assert frame.geometry.atoms[0].fragment_id == 1
    assert frame.geometry.atoms[0].basis_set == "TZVP"
    atom_basis_sets = [
        basis
        for basis in frame.model_chemistry.basis_sets
        if basis.scope == "atom" and basis.role == "orbital"
    ]
    assert atom_basis_sets
    assert atom_basis_sets[0].name == "TZVP"
    assert atom_basis_sets[0].atom_indices == [0]
    assert frame.model_chemistry.options["has_mixed_basis"] is True
    assert "ConnectFragments" in frame.resources_raw
    assert "optimizeHydrogens true" in frame.resources_raw


def test_orcainp_orca6_optimization_neb_block_and_xtb_geometry_are_structured() -> None:
    fixture_path = ORCA_OPTIMIZATION_INPUT_DIR / "orca6_optimizations_061_xtb_neb_ts.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "xyz"
    assert [block.name for block in frame.blocks] == ["neb"]
    assert 'neb_end_xyzfile "final.xyz"' in frame.resources_raw
    assert len(frame.atoms) == 8
    assert tuple(frame.coords.shape) == (8, 3)


def test_orcainp_orca6_frequency_opt_then_freq_example_is_structured() -> None:
    fixture_path = (
        ORCA_FREQUENCY_INPUT_DIR
        / "orca6_frequencies_000_bp_def2_tzvp_opt_anfreq_numfreq.inp"
    )
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "xyz"
    assert frame.method == "DFT"
    assert frame.functional == "BP"
    assert frame.basis_set == "def2-TZVP"
    assert "AnFreq" in frame.keywords
    assert "NumFreq" in frame.keywords
    assert [block.name for block in frame.blocks] == ["freq"]
    assert "CentralDiff true" in frame.resources_raw
    assert "Increment   0.005" in frame.resources_raw
    assert len(frame.atoms) == 4
    assert tuple(frame.coords.shape) == (4, 3)


def test_orcainp_orca6_frequency_restart_example_is_structured() -> None:
    fixture_path = ORCA_FREQUENCY_INPUT_DIR / "orca6_frequencies_002_sto_3g_numfreq_restart.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.geometry is not None
    assert frame.geometry.ctype == "int"
    assert frame.method == ""
    assert frame.functional == ""
    assert frame.basis_set == "STO-3G"
    assert "NumFreq" in frame.keywords
    assert [block.name for block in frame.blocks] == ["freq"]
    assert "Restart true" in frame.resources_raw
    assert len(frame.atoms) == 4
    assert tuple(frame.coords.shape) == (4, 3)


def test_orcainp_orca6_solvator_fixture_is_structured() -> None:
    fixture_path = ORCA_SOLVATOR_INPUT_DIR / "orca6_solvator_basic.inp"
    frame = AutoParser(str(fixture_path), parser_detection="orcainp")[0][-1]
    assert isinstance(frame, ORCAInpFileFrameDisk)
    assert frame.model_chemistry.solvation_model == "CPCM"
    assert frame.model_chemistry.solvent == "water"
    assert frame.explicit_solvent_requests
    request = frame.explicit_solvent_requests[0]
    assert request.enabled is True
    assert request.nsolv == 12
    assert request.cluster_mode == "sphere"
    assert request.droplet is True
