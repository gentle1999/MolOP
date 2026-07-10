from pathlib import Path

from molop.io.logic.gaussian.log.frame_models.G16Components import (
    G16ComponentTreeBuilder,
    G16JobCPUComponent,
    G16L202OrientComponent,
    G16L601PopAnalComponent,
    G16L716FreqComponent,
    G16L716ThermochemistryComponent,
    G16L9999ArchiveComponent,
    G16L9999FinalComponent,
    G16LeaveComponent,
    get_g16log_component_classes,
)
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"


def _parsed_frames():
    return G16LogFileParserMemory().parse(FIXTURE.read_text()).frames


def _last_frame():
    return _parsed_frames()[-1]


def test_g16log_component_tree_builds_major_nodes_from_frame_data():
    tree = _last_frame().component_tree
    component_names = tree.component_names()

    for expected_name in (
        "l202.orient",
        "l202.rotconst",
        "l502.cycle",
        "l601.popanal",
        "l716.freq",
        "l716.thermochemistry",
        "l716.polarizability",
        "l716.forces",
        "l716.secondderiv",
        "l103.optimization",
        "l9999.archive",
        "jobcpu",
        "l9999.final",
    ):
        assert expected_name in component_names

    assert tree.source_text == ""
    assert all(component.raw_text == "" for component in tree.iter_components())


def test_g16log_component_registry_exposes_declared_component_classes():
    component_class_names = tuple(
        component_cls.__name__ for component_cls in get_g16log_component_classes()
    )

    assert component_class_names == (
        "G16L1HeaderComponent",
        "G16L101TitleComponent",
        "G16L1KeywordsComponent",
        "G16L202OrientComponent",
        "G16L202RotConstComponent",
        "G16L502CycleComponent",
        "G16L601PopAnalComponent",
        "G16L716FreqComponent",
        "G16L716ThermochemistryComponent",
        "G16L716PolarizabilityComponent",
        "G16L716ForcesComponent",
        "G16L716SecondDerivComponent",
        "G16L103OptimizationComponent",
        "G16L9999ArchiveComponent",
        "G16L9999FinalComponent",
        "G16JobCPUComponent",
        "G16LinkEnterComponent",
        "G16LeaveComponent",
    )


def test_g16log_first_frame_includes_link1_header_component_from_frame_data():
    tree = _parsed_frames()[0].component_tree

    component_names = tree.component_names()
    expanded_names = tree.component_names(include_synthetic=True)

    assert "l1.header" in component_names
    assert "l1.options" in expanded_names
    assert "l1.keywords" in expanded_names
    assert "l101.title" in expanded_names
    assert "l101.charge_multiplicity" in expanded_names


def test_g16log_components_own_child_expansion_contracts():
    payload_by_name = {
        component.component_name: component
        for component in _last_frame().component_tree.iter_components(include_synthetic=True)
    }

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


def test_g16log_cycle_component_renders_scf_summary_from_frame_data():
    rendered = _last_frame().component_tree.render_node("l502.cycle")

    assert "SCF Done:" in rendered


def test_g16log_component_tree_expands_l716_child_nodes():
    tree = _last_frame().component_tree
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

    assert "l716.dipole" in node_names
    assert "l716.polarizability.detail" in node_names
    assert any(name.startswith("l716.vibration.mode[") for name in node_names)

    freq_nodes = [node for node in tree.iter_nodes() if node.node_name == "l716.freq"]
    assert freq_nodes
    assert any(
        any(child.node_name.startswith("l716.vibration.mode[") for child in freq_node.children)
        for freq_node in freq_nodes
    )


def test_g16log_component_queries_can_include_synthetic_children():
    tree = _last_frame().component_tree

    default_names = tree.component_names()
    expanded_names = tree.component_names(include_synthetic=True)

    assert "l716.forceconstants" not in default_names
    assert "l716.forceconstants" in expanded_names
    assert any(name == "l716.vibration.mode" for name in expanded_names)


def test_g16log_component_tree_contracts_validate_cleanly():
    assert _last_frame().component_tree.validate_contracts() == []


def test_g16log_component_tree_builder_from_frame_data_rebuilds_component_tree():
    tree = G16ComponentTreeBuilder.from_frame_data(_last_frame())

    component_names = tree.component_names(include_synthetic=True)
    assert "l202.orient" in component_names
    assert "l502.cycle" in component_names
    assert "l601.popanal" in component_names
    assert "l716.freq" in component_names
    assert "l716.thermochemistry" in component_names
    assert "l9999.archive" in component_names
    assert any(name == "l716.vibration.mode" for name in component_names)


def test_g16log_component_tree_builder_requires_declared_frame_fields():
    frame = _last_frame()
    frame.vibrations = None
    frame.thermal_informations = None

    tree = G16ComponentTreeBuilder.from_frame_data(frame)
    component_names = tree.component_names(include_synthetic=True)

    assert "l716.freq" not in component_names
    assert "l716.thermochemistry" not in component_names
    assert all(name != "l716.vibration.mode" for name in component_names)


def test_g16log_component_tree_builder_allows_optional_frame_fields():
    frame = _last_frame()
    frame.pressure = None

    tree = G16ComponentTreeBuilder.from_frame_data(frame)
    component_names = tree.component_names(include_synthetic=True)

    assert "l716.thermochemistry" in component_names


def test_g16log_components_expose_child_constraint_contracts():
    assert "l202.distmat" in G16L202OrientComponent.allowed_child_component_names
    assert "l601.molecular_orbitals" in G16L601PopAnalComponent.allowed_child_component_names
    assert "l716.vibration.mode" in G16L716FreqComponent.allowed_child_component_names
    assert (
        "l716.thermochemistry.zpe" in G16L716ThermochemistryComponent.allowed_child_component_names
    )
    assert "l9999.archive.energies" in G16L9999ArchiveComponent.allowed_child_component_names
    assert G16JobCPUComponent.allowed_child_component_names == ()
    assert G16L9999FinalComponent.allowed_child_component_names == ()
    assert G16LeaveComponent.allowed_child_component_names == ()


def test_g16log_component_tree_reports_contract_violations():
    tree = _last_frame().component_tree
    freq_node = tree.find_nodes("l716.freq")[0]
    freq_node.children.append(freq_node.children[0])

    issues = tree.validate_contracts()

    assert any("not repeatable" in issue for issue in issues)
