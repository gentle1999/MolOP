from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from molop.io import AutoParser
from molop.io.base_models.FrameParser import FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.gaussian.input.frame_models.GJFFileFrame import (
    GJFFileFrameMemory,
    GJFRouteSection,
)
from molop.io.logic.gaussian.input.frame_parsers.GJFFileFrameParser import GJFFileFrameParserMemory
from molop.io.logic.gaussian.input.GaussianInputPatterns import g16_input_patterns
from molop.io.logic.gaussian.input.GaussianLink0 import (
    GaussianLink0Commands,
    render_gaussian_link0_shared_memory_line,
)
from molop.io.logic.gaussian.input.GaussianRoute import (
    GaussianRouteSemantic,
    build_gaussian_model_chemistry,
)
from molop.io.logic.gaussian.input.GaussianRouteParsing import parse_gaussian_route_semantic
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import G16LogFileFrameMemory
from molop.io.logic.gaussian.log.models.G16LogFile import G16LogFileMemory
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.unit import atom_ureg


def test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities() -> None:
    semantic = parse_gaussian_route_semantic(
        "#p b3lyp/6-31g(d) opt scrf=(smd,solvent=ethanol) em=gd3bj"
    )
    assert semantic.dieze_tag == "#P"
    assert semantic.model_chemistry.method_token == "b3lyp"
    assert semantic.model_chemistry.method_family == "DFT"
    assert semantic.model_chemistry.functional == "b3lyp"
    assert semantic.model_chemistry.basis_set == "6-31g(d)"
    assert semantic.model_chemistry.basis_family == "pople"
    assert semantic.model_chemistry.basis_has_diffuse is False
    assert semantic.model_chemistry.basis_polarization == ["d"]
    assert "opt" in semantic.job_types
    assert "GeometryOptimization" in semantic.capabilities
    assert semantic.solvation_model == "scrf=(smd,solvent=ethanol)"
    assert "Dispersion" in semantic.capabilities
    assert semantic.empirical_dispersion == "gd3bj"
    assert semantic.option_maps["scrf"].params == {"smd": None, "solvent": "ethanol"}
    assert semantic.option_maps["em"].scalar_value == "gd3bj"
    assert semantic.opt_options.enabled is True
    assert semantic.scrf_options.enabled is True
    assert semantic.scrf_options.model == "smd"
    assert semantic.scrf_options.model_family == "smd"
    assert semantic.scrf_options.smd is True
    assert semantic.scrf_options.solvent == "ethanol"
    assert semantic.diagnostics.confidence >= 0.8
    assert semantic.to_route_dict()["em"] == "gd3bj"


def test_dieze_tag_allows_spacing_but_only_as_first_token() -> None:
    spaced = parse_gaussian_route_semantic("# p b3lyp/6-31g(d) opt")
    assert spaced.dieze_tag == "#P"

    non_leading = parse_gaussian_route_semantic("b3lyp/6-31g(d) #p opt")
    assert non_leading.dieze_tag is None
    assert "#p" in non_leading.unknown_tokens


def test_gjf_route_section_uses_shared_semantic_model() -> None:
    route = GJFRouteSection.from_str("#p b3lyp def2svp opt")
    semantic = route.semantic_route
    assert semantic.model_chemistry.method_token == "b3lyp"
    assert semantic.model_chemistry.method_family == "DFT"
    assert semantic.model_chemistry.basis_set == "def2svp"
    assert semantic.model_chemistry.basis_family == "def2"
    assert semantic.job_types == ["opt"]
    assert route.semantic_route.raw_route == route.route


def test_gjf_patterns_expose_named_groups_for_parser_fields() -> None:
    route_match = g16_input_patterns.ROUTE.search("#p hf/sto-3g\n\n")
    title_match = g16_input_patterns.TITLE.search("water\n")
    charge_match = g16_input_patterns.CHARGE_MULTIPLICITY.search("0 1\n")
    atom_match = g16_input_patterns.ATOMS.match("H 0.0 0.0 0.7")
    modredundant_ref_match = g16_input_patterns.MODREDUNDANT_ATOM_REF.match("*")

    assert route_match is not None
    assert route_match.group("route") == "#p hf/sto-3g"
    assert title_match is not None
    assert title_match.group("title") == "water"
    assert charge_match is not None
    assert charge_match.group("charge") == "0"
    assert charge_match.group("multiplicity") == "1"
    assert atom_match is not None
    assert atom_match.group("symbol") == "H"
    assert float(atom_match.group("x")) == 0.0
    assert float(atom_match.group("z")) == 0.7
    assert modredundant_ref_match is not None
    assert modredundant_ref_match.group("atom_ref") == "*"


def test_gjf_route_section_stores_explicit_semantic_field() -> None:
    semantic = parse_gaussian_route_semantic("#p hf/3-21g sp")
    route = GJFRouteSection.model_validate({"route": "#p hf/3-21g sp", "semantic_route": semantic})
    assert route.semantic_route.model_chemistry.method_token == "hf"
    assert route.semantic_route.model_chemistry.method_family == "HF"
    assert route.semantic_route.model_chemistry.functional is None
    assert route.semantic_route.job_types == ["sp"]


def test_gjf_route_section_direct_model_construction_does_not_parse_semantics() -> None:
    route = GJFRouteSection.model_validate({"route": "b3lyp/def2svp opt"})

    assert route.route == "# b3lyp/def2svp opt"
    assert route.semantic_route.raw_route == ""
    assert route.semantic_route.model_chemistry.method_token is None
    assert route.semantic_route.job_types == []


def test_gjf_frame_parser_returns_model_ready_sections() -> None:
    parser = GJFFileFrameParserMemory()
    block = """%nprocshared=4
%nosave
#p b3lyp/def2svp opt

water

0 1
O 0.0 0.0 0.0
H 0.0 0.0 0.9
H 0.8 0.0 0.0
"""
    result = parser._parse_block_to_result(block)
    assert isinstance(result, ModelParseResult)
    assert result.has_value("molecule_specifications") is True
    payload = parser._parse_frame(block, context=FrameParseContext(additional_data={}))

    assert payload["link0_commands"].cpu_request() == 4
    assert [(link0.key, link0.value) for link0 in payload["link0_commands"].link0_keywords] == [
        ("nprocshared", "4"),
        ("nosave", None),
    ]
    assert payload["route_section"].semantic_route.model_chemistry.method_family == "DFT"
    assert payload["route_section"].semantic_route.model_chemistry.basis_set == "def2svp"
    assert payload["title_card"].title_card == "water"
    assert payload["molecule_specifications"].total_charge == 0
    assert payload["molecule_specifications"].atomic_numbers() == [8, 1, 1]


def test_gaussian_link0_container_supports_gjf_and_g16_fakeg_projection() -> None:
    link0 = GaussianLink0Commands.from_str("%nprocshared=4\n%mem=1GB\n")

    assert link0.cpu_request() == 4
    assert link0.shared_memory_cpu_value() == "4"
    assert link0.shared_memory_cpu_render_line() == (
        " Will use up to    4 processors via shared memory."
    )
    assert render_gaussian_link0_shared_memory_line("%nprocshared=4\n%mem=1GB\n") == (
        " Will use up to    4 processors via shared memory."
    )


def test_gjf_frame_parser_returns_model_ready_additional_sections() -> None:
    parser = GJFFileFrameParserMemory()
    block = """#p b3lyp/def2svp opt=modredundant

water

0 1
O 0.0 0.0 0.0
H 0.0 0.0 0.9
H 0.8 0.0 0.0

B 1 2 F
"""

    payload = parser._parse_frame(block, context=FrameParseContext(additional_data={}))

    assert payload["additional_sections"].strip() == "B 1 2 F"
    assert payload["additional_section_diagnostics"] == []
    parsed_sections = payload["parsed_additional_sections"]
    assert len(parsed_sections) == 1
    assert parsed_sections[0].section_type == "modredundant"
    assert parsed_sections[0].lines[0].coordinate_type == "B"
    assert parsed_sections[0].lines[0].atom_refs == ["1", "2"]
    assert parsed_sections[0].lines[0].action == "F"


def test_hf_route_does_not_populate_functional() -> None:
    semantic = parse_gaussian_route_semantic("#p hf/3-21g sp")

    assert semantic.model_chemistry.method_token == "hf"
    assert semantic.model_chemistry.method_family == "HF"
    assert semantic.model_chemistry.functional is None
    assert semantic.model_chemistry.basis_set == "3-21g"


def test_spin_prefixed_dft_route_keeps_raw_token_and_normalizes_model_chemistry() -> None:
    semantic = parse_gaussian_route_semantic("#p RB3LYP/6-31G(d) freq")

    assert semantic.model_chemistry.method_token == "rb3lyp"
    assert semantic.model_chemistry.method_family == "DFT"
    assert semantic.model_chemistry.functional == "b3lyp"
    assert semantic.model_chemistry.spin_qualifier == "R"

    model = build_gaussian_model_chemistry(
        semantic,
        keywords=semantic.raw_route,
        legacy_method="DFT",
        legacy_functional="RB3LYP",
    )
    assert model.method == "B3LYP"
    assert model.functional == "B3LYP"
    assert model.spin_treatment == "R"


@pytest.mark.parametrize(
    ("route_method", "expected_method", "expected_spin", "expected_functional"),
    [
        ("UHF", "HF", "U", None),
        ("ROHF", "HF", "RO", None),
        ("R2SCAN", "R2SCAN", None, "R2SCAN"),
    ],
)
def test_gaussian_method_normalization_does_not_confuse_spin_prefixes(
    route_method: str,
    expected_method: str,
    expected_spin: str | None,
    expected_functional: str | None,
) -> None:
    semantic = parse_gaussian_route_semantic(f"#p {route_method}/def2svp sp")
    model = build_gaussian_model_chemistry(semantic, keywords=semantic.raw_route)

    assert model.method == expected_method
    assert model.spin_treatment == expected_spin
    assert model.functional == expected_functional


def test_gjf_frame_populates_qm_metadata_from_shared_semantic_route() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16gjf" / "test_solvent.gjf"
    batch = AutoParser(str(fixture_path))
    frame = cast(Any, batch[0][0])
    assert frame.method == "DFT"
    assert frame.basis_set.lower() == "def2svp"
    assert frame.functional == "B3LYP-GD3BJ"
    assert frame.route_section.semantic_route.model_chemistry.basis_set == "def2svp"
    assert frame.model_chemistry.method_family == "DFT"
    assert frame.model_chemistry.method == "B3LYP"
    assert frame.model_chemistry.functional == "B3LYP-GD3BJ"
    assert frame.model_chemistry.dispersion_correction == "GD3BJ"
    assert frame.model_chemistry.basis_set == "def2svp"
    assert {"opt", "freq", "population_analysis"} <= {
        task.task_type for task in frame.task_requests
    }


def test_g16log_frame_exposes_shared_semantic_route() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))
    frame = cast(Any, batch[0][0])
    semantic = frame.semantic_route
    assert semantic.model_chemistry.method_token == "ccsd"
    assert semantic.model_chemistry.method_family == "CCSD"
    assert semantic.model_chemistry.basis_set == "aug-cc-pvtz"
    assert "opt" in semantic.job_types
    assert frame.dieze_tag == semantic.dieze_tag


def test_g16log_segment_metadata_result_populates_semantic_route() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"

    result = G16LogFileParserMemory()._parse_segment_metadata_result(fixture_path.read_text())
    semantic = result.model_data()["semantic_route"]

    assert isinstance(semantic, GaussianRouteSemantic)
    assert semantic.model_chemistry.method_family == "CCSD"
    assert semantic.model_chemistry.basis_set == "aug-cc-pvtz"


def test_g16log_link1_frames_use_stable_concrete_method_names() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "H2O.log"
    parsed_file = AutoParser(str(fixture_path))[0]

    assert {frame.model_chemistry.method for frame in parsed_file} == {"B3LYP"}
    assert parsed_file[-1].semantic_route.model_chemistry.method_token == "rb3lyp"
    assert parsed_file[-1].model_chemistry.functional == "B3LYP"
    assert parsed_file[-1].model_chemistry.spin_treatment == "R"


def test_g16log_models_do_not_parse_semantic_route_from_raw_keywords() -> None:
    file_model = G16LogFileMemory.model_validate({"keywords": "#p b3lyp/def2svp opt"})
    frame_model = G16LogFileFrameMemory.model_validate(
        {
            "keywords": "#p b3lyp/def2svp opt",
            "atoms": [8],
            "coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
            "standard_coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
        }
    )

    assert file_model.semantic_route.raw_route == ""
    assert file_model.semantic_route.model_chemistry.method_token is None
    assert frame_model.semantic_route.raw_route == ""
    assert frame_model.semantic_route.model_chemistry.method_token is None


def test_g16log_models_project_hf_method_from_semantic_route() -> None:
    semantic = parse_gaussian_route_semantic("#p hf/3-21g sp")
    file_model = G16LogFileMemory.model_validate(
        {
            "keywords": semantic.raw_route,
            "semantic_route": semantic,
            "method": "stale-method",
            "functional": "stale-functional",
            "basis_set": "stale-basis",
        }
    )
    frame_model = G16LogFileFrameMemory.model_validate(
        {
            "keywords": semantic.raw_route,
            "semantic_route": semantic,
            "method": "stale-method",
            "functional": "stale-functional",
            "basis_set": "stale-basis",
            "atoms": [1],
            "coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
            "standard_coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
        }
    )

    assert file_model.method == "HF"
    assert file_model.functional == ""
    assert file_model.basis_set == "3-21g"
    assert file_model.model_chemistry.method_family == "HF"
    assert frame_model.method == "HF"
    assert frame_model.functional == ""
    assert frame_model.basis_set == "3-21g"
    assert frame_model.model_chemistry.method_family == "HF"


def test_gjf_frame_projects_hf_method_from_shared_semantic_route() -> None:
    frame = GJFFileFrameMemory.model_validate(
        {
            "route_section": GJFRouteSection.from_str("#p hf/3-21g sp"),
            "atoms": [1],
            "coords": np.array([[0.0, 0.0, 0.0]]) * atom_ureg.angstrom,
        }
    )

    assert frame.method == "HF"
    assert frame.functional == ""
    assert frame.basis_set == "3-21g"
    assert frame.model_chemistry.method_family == "HF"
    assert frame.model_chemistry.functional is None
    assert frame.model_chemistry.basis_set == "3-21g"


def test_gjf_route_section_to_dict_projects_from_semantic_route() -> None:
    route = GJFRouteSection.from_str("#p b3lyp/6-31g(d) opt scrf=(smd,solvent=ethanol) em=gd3bj")
    projected = route.to_dict()
    assert projected["b3lyp"] is None
    assert projected["6-31g(d)"] is None
    assert projected["opt"] is None
    assert projected["scrf"] == {"smd": None, "solvent": "ethanol"}
    assert projected["em"] == "gd3bj"


def test_shared_semantic_route_exposes_structured_dispersion_and_solvation() -> None:
    fixture_path = (
        Path(__file__).resolve().parent
        / "test_files"
        / "g16log"
        / "MnCO3C6H6PMe3-mod2-sp-smd-revDSDPBEP86d3.log"
    )
    batch = AutoParser(str(fixture_path))
    frame = cast(Any, batch[0][0])
    semantic = frame.semantic_route
    assert semantic.solvation_model is not None
    assert semantic.solvation_model.lower().startswith("scrf=")
    assert semantic.empirical_dispersion == "gd3bj"
    assert semantic.option_maps["scrf"].params["solvent"] == "nitromethane"
    assert frame.functional.endswith("-GD3BJ")
    assert frame.model_chemistry.dispersion_correction == "GD3BJ"
    assert frame.model_chemistry.functional.endswith("-GD3BJ")
    assert frame.model_chemistry.solvation_model == "smd"
    assert frame.model_chemistry.solvent == "nitromethane"


def test_semantic_route_captures_geom_and_external_options() -> None:
    semantic = parse_gaussian_route_semantic("#p freq geom=allcheck external='./xtb.sh'")
    assert semantic.checkpoint_geometry_mode == "allcheck"
    assert semantic.option_maps["geom"].scalar_value == "allcheck"
    assert semantic.geom_options.enabled is True
    assert semantic.geom_options.mode == "allcheck"
    assert semantic.geom_options.allcheck is True
    assert semantic.external_program == "'./xtb.sh'"
    assert semantic.option_maps["external"].scalar_value == "'./xtb.sh'"


def test_geom_options_are_strictly_extracted_from_param_map() -> None:
    semantic = parse_gaussian_route_semantic(
        "#p geom=(checkpoint,huge,modify,newdefinition,newredundant,notest,addgic,readallgic,connectivity,modconnectivity,genconnectivity,zmconnectivity,distance,nodistance,cangle,angle,noangle,cdihedral,dihedral,nodihedral,printinputorient,print,step=3,nogic)"
    )
    assert semantic.checkpoint_geometry_mode == "checkpoint"
    assert semantic.geom_options.enabled is True
    assert semantic.geom_options.mode == "checkpoint"
    assert semantic.geom_options.checkpoint is True
    assert semantic.geom_options.huge is True
    assert semantic.geom_options.modify is True
    assert semantic.geom_options.new_definition is True
    assert semantic.geom_options.new_redundant is True
    assert semantic.geom_options.no_test is True
    assert semantic.geom_options.add_gic is True
    assert semantic.geom_options.read_all_gic is True
    assert semantic.geom_options.connectivity is True
    assert semantic.geom_options.mod_connectivity is True
    assert semantic.geom_options.gen_connectivity is True
    assert semantic.geom_options.zm_connectivity is True
    assert semantic.geom_options.distance is True
    assert semantic.geom_options.no_distance is True
    assert semantic.geom_options.cangle is True
    assert semantic.geom_options.angle is True
    assert semantic.geom_options.no_angle is True
    assert semantic.geom_options.cdihedral is True
    assert semantic.geom_options.dihedral is True
    assert semantic.geom_options.no_dihedral is True
    assert semantic.geom_options.print_input_orient is True
    assert semantic.geom_options.print is True
    assert semantic.geom_options.no_gic is True
    assert semantic.geom_options.step == 3
    assert semantic.geom_options.ngeom == 4


def test_double_slash_route_implies_opt_plus_sp_with_layered_model_chemistry() -> None:
    semantic = parse_gaussian_route_semantic("# CCSD/6-31G(d)//B3LYP/6-31G(d)")
    assert semantic.model_chemistry.method_token == "CCSD"
    assert semantic.model_chemistry.method_family == "CCSD"
    assert semantic.model_chemistry.basis_set == "6-31G(d)"
    assert semantic.model_chemistry.low_level is not None
    assert semantic.model_chemistry.low_level.method_token == "B3LYP"
    assert semantic.model_chemistry.low_level.method_family == "DFT"
    assert semantic.model_chemistry.low_level.basis_set == "6-31G(d)"
    assert "opt" in semantic.job_types
    assert "sp" in semantic.job_types
    assert "GeometryOptimization" in semantic.capabilities
    assert "SinglePointEnergy" in semantic.capabilities


def test_shared_route_parser_supports_more_gaussian_job_types_and_params() -> None:
    semantic = parse_gaussian_route_semantic(
        "#p stable volume force td(nstates=10) pop=(full,nbo) admp ircmax oniom(b3lyp/6-31g(d):pm3)"
    )
    for job in ["stable", "volume", "force", "td", "pop", "admp", "ircmax", "oniom"]:
        assert job in semantic.job_types
    for capability in [
        "WavefunctionStability",
        "MolecularVolume",
        "ForceConstants",
        "ExcitedState",
        "PopulationAnalysis",
        "DirectDynamics",
        "ReactionPathMaximum",
        "ONIOM",
    ]:
        assert capability in semantic.capabilities
    assert semantic.option_maps["td"].params == {"nstates": "10"}
    assert semantic.option_maps["pop"].params == {"full": None, "nbo": None}
    assert semantic.option_maps["oniom"].params == {"b3lyp/6-31g(d):pm3": None}
    assert semantic.td_options.enabled is True
    assert semantic.td_options.nstates == 10
    assert semantic.pop_options.enabled is True
    assert semantic.pop_options.full is True
    assert semantic.pop_options.nbo is True


def test_population_options_are_strictly_extracted() -> None:
    semantic = parse_gaussian_route_semantic(
        "#p pop=(none,full,nbo,nboread,nbo6read,nbo7read,hirshfeld,cm5,mk,chelpg,orbitals=3,readradii,readatradii)"
    )
    assert semantic.pop_options.enabled is True
    assert semantic.pop_options.none is True
    assert semantic.pop_options.full is True
    assert semantic.pop_options.nbo is True
    assert semantic.pop_options.nbo_read is True
    assert semantic.pop_options.nbo6_read is True
    assert semantic.pop_options.nbo7_read is True
    assert semantic.pop_options.hirshfeld is True
    assert semantic.pop_options.cm5 is True
    assert semantic.pop_options.mk is True
    assert semantic.pop_options.chelpg is True
    assert semantic.pop_options.orbitals == 3
    assert semantic.pop_options.read_radii is True
    assert semantic.pop_options.read_at_radii is True


def test_single_point_keyword_is_explicitly_supported() -> None:
    semantic = parse_gaussian_route_semantic("#p sp b3lyp/def2svp")
    assert semantic.job_types == ["sp"]
    assert semantic.capabilities == ["SinglePointEnergy"]


def test_typed_common_job_submodels_extract_core_semantics_and_preserve_extra_options() -> None:
    semantic = parse_gaussian_route_semantic(
        "#p opt=(restart,ts,saddle=2,verytight,calcfc,calcall,calchffc,readfc,rcfc,maxcycles=120,maxstep=8,recalcfc=5,tight,noexpert,noeigentest,cartesian) modredundant freq=(anharmonic,readanharm,projected,tprojected,hinderedrotor,vibrot,polar,hpmodes,readisotopes,selectnormalmodes,savenormalmodes,vcd,raman,noraman,cphf=rdfreq,layer=real,atoms=1-3,notatoms=H,temperature=350,pressure=2.0) td(nstates=12,root=2,triplets,tda) scrf=(read,iefpcm,solvent=water)"
    )
    assert semantic.opt_options.enabled is True
    assert semantic.opt_options.restart is True
    assert semantic.opt_options.transition_state is True
    assert semantic.opt_options.saddle_order == 2
    assert semantic.opt_options.very_tight is True
    assert semantic.opt_options.calc_fc is True
    assert semantic.opt_options.calc_all is True
    assert semantic.opt_options.calc_hf_fc is True
    assert semantic.opt_options.tight is True
    assert semantic.opt_options.max_cycles == 120
    assert semantic.opt_options.max_step == 8
    assert semantic.opt_options.recalc_fc == 5
    assert semantic.opt_options.read_fc is True
    assert semantic.opt_options.read_cartesian_fc is True
    assert semantic.opt_options.has_modredundant is True
    assert semantic.opt_options.expert is False
    assert semantic.opt_options.eigen_test is False
    assert semantic.opt_options.coordinate_system == "cartesian"
    assert semantic.opt_options.extra_options == {}

    assert semantic.freq_options.enabled is True
    assert semantic.freq_options.anharmonic is True
    assert semantic.freq_options.read_anharm is True
    assert semantic.freq_options.projected is True
    assert semantic.freq_options.tprojected is True
    assert semantic.freq_options.hindered_rotor is True
    assert semantic.freq_options.vibrot is True
    assert semantic.freq_options.polar is True
    assert semantic.freq_options.hpmodes is True
    assert semantic.freq_options.read_isotopes is True
    assert semantic.freq_options.select_normal_modes is True
    assert semantic.freq_options.save_normal_modes is True
    assert semantic.freq_options.vcd is True
    assert semantic.freq_options.raman is True
    assert semantic.freq_options.no_raman is True
    assert semantic.freq_options.cphf_rd_freq is True
    assert semantic.freq_options.layer == "real"
    assert semantic.freq_options.atoms == "1-3"
    assert semantic.freq_options.not_atoms == "h"
    assert semantic.freq_options.temperature == 350.0
    assert semantic.freq_options.pressure == 2.0

    assert semantic.td_options.enabled is True
    assert semantic.td_options.nstates == 12
    assert semantic.td_options.root == 2
    assert semantic.td_options.triplets is True
    assert semantic.td_options.tda is True

    assert semantic.scrf_options.enabled is True
    assert semantic.scrf_options.read is True
    assert semantic.scrf_options.iefpcm is True
    assert semantic.scrf_options.model_family == "iefpcm"
    assert semantic.scrf_options.solvent == "water"
