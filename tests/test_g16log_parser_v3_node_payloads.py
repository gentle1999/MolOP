from pathlib import Path

from molop.io.logic.QM_frame_parsers.G16LogFileFrameParserV3 import G16LogFileFrameParserV3Memory
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"


def _build_tree():
    file_content = FIXTURE.read_text()
    file_parser = G16LogFileParserMemory()
    block = file_parser._split_file(file_content)[-1]
    parser = G16LogFileFrameParserV3Memory()
    tree = parser.build_component_tree(block)
    for component in tree.iter_components():
        component.parse_payload(block, tree)
    return tree


def test_g16log_v3_major_components_parse_payloads():
    tree = _build_tree()
    components = list(tree.iter_components())
    payloads = {component.component_name: component.payload for component in tree.iter_components()}
    jobcpu_payloads = [
        component.payload for component in components if component.component_name == "jobcpu"
    ]
    final_payloads = [
        component.payload for component in components if component.component_name == "l9999.final"
    ]

    assert payloads["l202.orient"]["atoms"]
    assert (
        payloads["l202.orient"].get("coords") is not None
        or payloads["l202.orient"].get("standard_coords") is not None
    )
    assert payloads["l502.cycle"].get("energies") is not None
    assert payloads["l601.popanal"].get("charge_spin_populations") is not None
    assert payloads["l716.freq"].get("vibrations") is not None
    assert payloads["l716.thermochemistry"].get("thermal_informations") is not None
    assert payloads["l9999.archive"].get("charge") is not None
    assert payloads["l9999.archive"].get("multiplicity") is not None
    assert any(payload.get("job_cpu_time") is not None for payload in jobcpu_payloads)
    assert any(payload.get("elapsed_time") is not None for payload in jobcpu_payloads)
    assert any(payload.get("status") is not None for payload in jobcpu_payloads)
    assert any(payload.get("termination_kind") == "Normal" for payload in jobcpu_payloads)
    assert any(
        "Normal termination of Gaussian" in payload.get("termination_line", "")
        for payload in jobcpu_payloads
    )
    assert any(payload.get("status") is not None for payload in final_payloads)
    assert any(payload.get("termination_kind") == "Normal" for payload in final_payloads)


def test_g16log_v3_first_frame_parses_link1_header_payloads():
    file_content = FIXTURE.read_text()
    block = G16LogFileParserMemory()._split_file(file_content)[0]
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    payload_by_name = {
        component.component_name: component.payload
        for component in tree.iter_components(include_synthetic=True)
    }

    assert payload_by_name["l1.header"]["keywords"].strip().startswith("# opt=maxcycle=150")
    assert "%mem=" in payload_by_name["l1.header"]["options"]
    assert payload_by_name["l1.header"]["title_card"].strip() == "B3D3_Pople_Opt_Gas"
    assert payload_by_name["l1.header"]["charge"] == -1
    assert payload_by_name["l1.header"]["multiplicity"] == 1
    assert payload_by_name["l1.keywords"]["keywords"].strip().startswith("# opt=maxcycle=150")
    assert "%nprocshared=" in payload_by_name["l1.options"]["options"]
    assert payload_by_name["l101.title"]["title_card"].strip() == "B3D3_Pople_Opt_Gas"
    assert payload_by_name["l101.charge_multiplicity"]["charge"] == -1
    assert payload_by_name["l101.charge_multiplicity"]["multiplicity"] == 1


def test_g16log_v3_l716_children_decompose_parent_payloads():
    file_content = FIXTURE.read_text()
    block = G16LogFileParserMemory()._split_file(file_content)[-1]
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    payload_by_node = {
        node.node_name: node.component.payload
        for node in tree.iter_nodes()
        if node.component is not None and not node.include_in_aggregation
    }
    assert payload_by_node["l716.forceconstants"].get("force_constants") is not None
    assert payload_by_node["l716.diagvib"].get("vibration_modes") is not None
    assert payload_by_node["l716.irspectrum"].get("IR_intensities") is not None
    assert payload_by_node["l716.thermochemistry.mass"].get("molecular_mass") is not None
    assert payload_by_node["l716.thermochemistry.temperature"].get("temperature") is not None
    assert payload_by_node["l716.thermochemistry.moi"].get("moments_of_inertia") is not None
    assert (
        payload_by_node["l716.thermochemistry.rotsymnum"].get("rotational_symmetry_number")
        is not None
    )
    assert (
        payload_by_node["l716.thermochemistry.rottemp"].get("rotational_temperatures") is not None
    )
    assert payload_by_node["l716.thermochemistry.rotconsts"].get("rotational_constants") is not None
    assert (
        payload_by_node["l716.thermochemistry.vibtemp"].get("vibrational_temperatures") is not None
    )
    assert payload_by_node["l716.thermochemistry.zpe"].get("ZPVE") is not None
    assert payload_by_node["l716.thermochemistry.energy"]
    assert payload_by_node["l716.thermochemistry.enthalpy"]
    assert payload_by_node["l716.thermochemistry.gibbs"]
    assert payload_by_node["l716.thermochemistry.entropy"].get("S") is not None
    assert payload_by_node["l716.thermochemistry.heatcapacity"].get("C_V") is not None
    assert "l716.dipole" not in payload_by_node
    assert "l716.polarizability.detail" not in payload_by_node


def test_g16log_v3_l716_vibration_mode_children_expose_per_mode_semantics():
    file_content = FIXTURE.read_text()
    block = G16LogFileParserMemory()._split_file(file_content)[-1]
    parser = G16LogFileFrameParserV3Memory()
    parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None

    vibration_mode_nodes = [
        node
        for node in tree.iter_nodes()
        if node.node_name.startswith("l716.vibration.mode[") and node.component is not None
    ]

    assert vibration_mode_nodes
    assert vibration_mode_nodes[0].component is not None

    first_payload = vibration_mode_nodes[0].component.payload
    assert first_payload.get("mode_index") == 0
    assert first_payload.get("frequency") is not None
    assert first_payload.get("reduced_mass") is not None
    assert first_payload.get("force_constant") is not None
    assert first_payload.get("IR_intensity") is not None
    assert first_payload.get("vibration_mode") is not None
    assert "is_imaginary" in first_payload

    freq_nodes = [node for node in tree.iter_nodes() if node.node_name == "l716.freq"]
    assert freq_nodes
    assert "l716.vibration.mode[0]" in [child.node_name for child in freq_nodes[0].children]
