from pathlib import Path

from molop.io.logic.QM_frame_models.G16V3Components import (
    G16V3ComponentTreeBuilder,
    G16V3JobCPUComponent,
    G16V3L1HeaderComponent,
    G16V3L202OrientComponent,
    G16V3L502CycleComponent,
    G16V3L601PopAnalComponent,
    G16V3L716FreqComponent,
    G16V3L716ThermochemistryComponent,
    G16V3L9999ArchiveComponent,
    G16V3L9999FinalComponent,
    G16V3LeaveComponent,
    G16V3LinkEnterComponent,
)
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParserV3 import G16LogFileFrameParserV3Memory
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"
L502_CYCLE_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "02209_00862_00000_0000.log"
)


def _last_frame_block() -> str:
    file_content = FIXTURE.read_text()
    parser = G16LogFileParserMemory()
    return parser._split_file(file_content)[-1]


def test_g16log_v3_tree_segments_major_nodes_in_source_order():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    tree = parser.build_component_tree(block)

    component_names = tree.component_names()

    assert "l202.orient" in component_names
    assert "l502.cycle" in component_names
    assert "l601.popanal" in component_names
    assert "l716.freq" in component_names
    assert "l716.thermochemistry" in component_names
    assert "l103.optimization" in component_names
    assert "l9999.archive" in component_names
    assert "jobcpu" in component_names
    assert "l9999.final" in component_names

    ordered_indices = [
        component_names.index(name)
        for name in [
            "l202.orient",
            "l502.cycle",
            "l601.popanal",
            "l716.freq",
            "l716.thermochemistry",
            "l103.optimization",
            "l9999.archive",
            "jobcpu",
            "l9999.final",
        ]
    ]
    assert ordered_indices == sorted(ordered_indices)


def test_g16log_v3_first_frame_includes_link1_header_component():
    block = G16LogFileParserMemory()._split_file(FIXTURE.read_text())[0]
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None
    component_names = tree.component_names()
    expanded_names = tree.component_names(include_synthetic=True)

    assert "l1.header" in component_names
    assert "l1.options" in expanded_names
    assert "l1.keywords" in expanded_names
    assert "l101.title" in expanded_names
    assert "l101.charge_multiplicity" in expanded_names


def test_g16log_v3_components_expose_class_local_boundary_specs():
    cycle_spec = G16V3L502CycleComponent.primary_boundary_spec()
    jobcpu_spec = G16V3JobCPUComponent.primary_boundary_spec()
    final_spec = G16V3L9999FinalComponent.primary_boundary_spec()
    link_enter_spec = G16V3LinkEnterComponent.primary_boundary_spec()
    l1_header_spec = G16V3L1HeaderComponent.primary_boundary_spec()

    assert cycle_spec.boundary_source == "iochem"
    assert cycle_spec.match_mode == "section"
    assert cycle_spec.include_end is True
    assert r"^\s*Cycle\s+.*" in cycle_spec.start_patterns
    assert r"^\s*SCF Done:.*" in cycle_spec.end_patterns

    assert jobcpu_spec.boundary_source == "iochem"
    assert jobcpu_spec.match_mode == "section"
    assert jobcpu_spec.overlap_policy == "allow_shared_end"
    assert r"^\s*(Job cpu time|Elapsed time):" in jobcpu_spec.start_patterns

    assert final_spec.boundary_source == "iochem"
    assert final_spec.overlap_policy == "allow_shared_end"
    assert r"^\s*Normal termination of Gaussian" in final_spec.start_patterns

    assert link_enter_spec.boundary_source == "iochem"
    assert link_enter_spec.match_mode == "line"
    assert r"^\s*\(Enter\s+.+\)\s*$" in link_enter_spec.start_patterns

    assert l1_header_spec.boundary_source == "iochem"
    assert l1_header_spec.match_mode == "section"
    assert any("Entering Link 1 =" in marker for marker in l1_header_spec.contains_markers)
    assert any("Input orientation:" in pattern for pattern in l1_header_spec.end_patterns)


def test_g16log_v3_components_own_child_expansion_contracts():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    tree = parser.build_component_tree(block)

    payload_by_name = {}
    for component in tree.iter_components():
        component.parse_payload(block, tree)
        payload_by_name[component.component_name] = component

    freq_children = payload_by_name["l716.freq"].build_synthetic_children()
    thermo_children = payload_by_name["l716.thermochemistry"].build_synthetic_children()

    assert any(
        child.component_cls.component_name == "l716.forceconstants" for child in freq_children
    )
    assert any(
        child.component_cls.component_name == "l716.vibration.mode" for child in freq_children
    )
    assert any(
        child.component_cls.component_name == "l716.thermochemistry.zpe"
        for child in thermo_children
    )
    assert any(
        child.component_cls.component_name == "l716.thermochemistry.mass"
        for child in thermo_children
    )


def test_g16log_v3_cycle_component_captures_cycle_body_and_scf_summary():
    block = G16LogFileParserMemory()._split_file(L502_CYCLE_FIXTURE.read_text())[0]
    parser = G16LogFileFrameParserV3Memory()
    tree = parser.build_component_tree(block)

    cycle_components = [
        component
        for component in tree.iter_components()
        if component.component_name == "l502.cycle"
    ]

    assert cycle_components
    assert any(
        "Cycle" in component.raw_text and "SCF Done:" in component.raw_text
        for component in cycle_components
    )


def test_g16log_v3_parsed_tree_expands_l716_child_nodes():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    node_names = [node.node_name for node in tree.iter_nodes()]

    for expected_name in (
        "l716.forceconstants",
        "l716.diagvib",
        "l716.irspectrum",
        "l716.thermochemistry.mass",
        "l716.thermochemistry.temperature",
        "l716.thermochemistry.moi",
        "l716.thermochemistry.rotsymnum",
        "l716.thermochemistry.rottemp",
        "l716.thermochemistry.rotconsts",
        "l716.thermochemistry.vibtemp",
        "l716.thermochemistry.zpe",
        "l716.thermoprops",
        "l716.thermochemistry.energy",
        "l716.thermochemistry.enthalpy",
        "l716.thermochemistry.gibbs",
        "l716.thermochemistry.entropy",
        "l716.thermochemistry.heatcapacity",
    ):
        assert expected_name in node_names

    assert "l716.dipole" not in node_names
    assert "l716.polarizability.detail" not in node_names
    assert any(name.startswith("l716.vibration.mode[") for name in node_names)

    freq_nodes = [node for node in tree.iter_nodes() if node.node_name == "l716.freq"]
    assert freq_nodes
    assert any(
        any(child.node_name.startswith("l716.vibration.mode[") for child in freq_node.children)
        for freq_node in freq_nodes
    )


def test_g16log_v3_component_queries_can_include_synthetic_children():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    default_names = tree.component_names()
    expanded_names = tree.component_names(include_synthetic=True)

    assert "l716.forceconstants" not in default_names
    assert "l716.forceconstants" in expanded_names
    assert any(name == "l716.vibration.mode" for name in expanded_names)


def test_g16log_v3_protocol_tree_contracts_validate_cleanly():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None
    assert tree.validate_contracts() == []


def test_g16log_v3_protocol_builder_from_frame_data_rebuilds_component_tree():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    frame = parser.parse(block)

    tree = G16V3ComponentTreeBuilder.from_frame_data(frame)

    component_names = tree.component_names(include_synthetic=True)
    assert "l202.orient" in component_names
    assert "l502.cycle" in component_names
    assert "l601.popanal" in component_names
    assert "l716.freq" in component_names
    assert "l716.thermochemistry" in component_names
    assert "l9999.archive" in component_names
    assert any(name == "l716.vibration.mode" for name in component_names)


def test_g16log_v3_protocol_builder_requires_declared_frame_fields():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    frame = parser.parse(block)

    frame.vibrations = None
    frame.thermal_informations = None

    tree = G16V3ComponentTreeBuilder.from_frame_data(frame)
    component_names = tree.component_names(include_synthetic=True)

    assert "l716.freq" not in component_names
    assert "l716.thermochemistry" not in component_names
    assert all(name != "l716.vibration.mode" for name in component_names)


def test_g16log_v3_protocol_builder_allows_optional_frame_fields():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    frame = parser.parse(block)

    frame.pressure = None

    tree = G16V3ComponentTreeBuilder.from_frame_data(frame)
    component_names = tree.component_names(include_synthetic=True)

    assert "l716.thermochemistry" in component_names


def test_g16log_v3_components_expose_child_constraint_contracts():
    assert "l202.distmat" in G16V3L202OrientComponent.allowed_child_component_names
    assert "l601.molecular_orbitals" in G16V3L601PopAnalComponent.allowed_child_component_names
    assert "l716.vibration.mode" in G16V3L716FreqComponent.allowed_child_component_names
    assert (
        "l716.thermochemistry.zpe"
        in G16V3L716ThermochemistryComponent.allowed_child_component_names
    )
    assert "l9999.archive.energies" in G16V3L9999ArchiveComponent.allowed_child_component_names
    assert G16V3JobCPUComponent.allowed_child_component_names == ()
    assert G16V3L9999FinalComponent.allowed_child_component_names == ()
    assert G16V3LeaveComponent.allowed_child_component_names == ()


def test_g16log_v3_protocol_tree_reports_contract_violations():
    block = _last_frame_block()
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    freq_node = tree.find_nodes("l716.freq")[0]
    freq_node.children.append(freq_node.children[0])

    issues = tree.validate_contracts()

    assert any("not repeatable" in issue for issue in issues)
