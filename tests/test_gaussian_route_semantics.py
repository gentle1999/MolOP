from pathlib import Path
from typing import Any, cast

from molop.io import AutoParser
from molop.io.logic.gaussian_route_models import parse_gaussian_route_semantic
from molop.io.logic.qminput_frame_models.GJFFileFrame import GJFRouteSection


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


def test_gjf_route_section_stores_explicit_semantic_field() -> None:
    semantic = parse_gaussian_route_semantic("#p hf/3-21g sp")
    route = GJFRouteSection.model_validate({"route": "#p hf/3-21g sp", "semantic_route": semantic})
    assert route.semantic_route.model_chemistry.method_token == "hf"
    assert route.semantic_route.job_types == ["sp"]


def test_gjf_frame_populates_qm_metadata_from_shared_semantic_route() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16gjf" / "test_solvent.gjf"
    batch = AutoParser(str(fixture_path))
    frame = cast(Any, batch[0][0])
    assert frame.method == "DFT"
    assert frame.basis_set.lower() == "def2svp"
    assert frame.functional.lower() == "b3lyp"
    assert frame.route_section.semantic_route.model_chemistry.basis_set == "def2svp"


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
