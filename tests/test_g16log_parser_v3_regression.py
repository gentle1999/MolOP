from pathlib import Path

import pytest

from molop.io.logic.QM_frame_models.G16LogFileFrame import G16LogFileFrameMemory
from molop.io.logic.QM_frame_models.G16V3Components import aggregate_g16log_v3_tree
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParser import G16LogFileFrameParserMemory
from molop.io.logic.QM_frame_parsers.G16LogFileFrameParserV3 import G16LogFileFrameParserV3Memory
from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURES = [
    Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log",
    Path(__file__).resolve().parent / "test_files" / "g16log" / "1-INT1-Sp.log",
]

ARCHIVE_ONLY_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
CCSD_FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "CH3-ccsd-sp.log"


def _last_frame_block(fixture: Path) -> str:
    file_content = fixture.read_text()
    parser = G16LogFileParserMemory()
    return parser._split_file(file_content)[-1]


def test_g16log_v3_matches_v1_on_representative_fixtures():
    for fixture in FIXTURES:
        block = _last_frame_block(fixture)
        v1_frame = G16LogFileFrameParserMemory().parse(block)
        v3_frame = G16LogFileFrameParserV3Memory().parse(block)

        assert v3_frame.atoms == v1_frame.atoms
        assert v3_frame.coords is not None
        assert v1_frame.coords is not None
        assert v3_frame.standard_coords is not None or v1_frame.standard_coords is None
        assert v3_frame.energies is not None
        assert v1_frame.energies is not None
        assert v3_frame.energies.scf_energy is not None
        assert v1_frame.energies.scf_energy is not None
        assert v3_frame.energies.scf_energy.to("hartree").m == pytest.approx(
            v1_frame.energies.scf_energy.to("hartree").m
        )
        assert bool(v3_frame.vibrations) == bool(v1_frame.vibrations)
        assert bool(v3_frame.molecular_orbitals) == bool(v1_frame.molecular_orbitals)
        assert bool(v3_frame.charge_spin_populations) == bool(v1_frame.charge_spin_populations)
        assert bool(v3_frame.polarizability) == bool(v1_frame.polarizability)
        assert bool(v3_frame.hessian is not None) == bool(v1_frame.hessian is not None)
        assert bool(v3_frame.forces is not None) == bool(v1_frame.forces is not None)
        assert bool(v3_frame.geometry_optimization_status) == bool(
            v1_frame.geometry_optimization_status
        )


def test_g16log_v3_does_not_invent_archive_energies_on_terminal_frame():
    block = _last_frame_block(ARCHIVE_ONLY_FIXTURE)
    v1_frame = G16LogFileFrameParserMemory().parse(block)
    v3_frame = G16LogFileFrameParserV3Memory().parse(block)

    assert v1_frame.energies is None
    assert v3_frame.energies is None
    assert v1_frame.status is None
    assert v3_frame.status is None
    assert v1_frame.temperature is None
    assert v3_frame.temperature is None
    assert v1_frame.pressure is None
    assert v3_frame.pressure is None


def test_g16log_v3_preserves_live_scf_energy_over_archive_value():
    block = _last_frame_block(CCSD_FIXTURE)
    v1_frame = G16LogFileFrameParserMemory().parse(block)
    v3_frame = G16LogFileFrameParserV3Memory().parse(block)

    assert v1_frame.energies is not None
    assert v3_frame.energies is not None
    assert v1_frame.energies.scf_energy is not None
    assert v3_frame.energies.scf_energy is not None
    assert v3_frame.energies.scf_energy.to("hartree").m == pytest.approx(
        v1_frame.energies.scf_energy.to("hartree").m
    )


def test_g16log_v3_component_aggregation_contract_matches_frame_output():
    block = _last_frame_block(FIXTURES[0])
    parser = G16LogFileFrameParserV3Memory()
    frame = parser.parse(block)
    tree = parser.last_component_tree

    assert tree is not None
    component_names = {component.component_name for component in tree.iter_components()}
    assert "l502.cycle" in component_names
    assert "l601.popanal" in component_names
    assert "l716.freq" in component_names
    assert all(hasattr(component, "aggregate_into") for component in tree.iter_components())
    assert frame.energies is not None
    assert frame.molecular_orbitals is not None
    assert frame.vibrations is not None


def test_g16log_v3_parsed_frame_exposes_component_tree():
    block = _last_frame_block(FIXTURES[0])
    frame = G16LogFileFrameParserV3Memory().parse(block)

    assert frame.component_tree is not None
    assert "l1.header" in frame.component_tree.component_names()
    assert "l1.keywords" in frame.component_tree.component_names(include_synthetic=True)


def test_g16log_file_parser_attaches_component_tree_to_frames():
    parsed = G16LogFileParserMemory().parse(FIXTURES[0].read_text())

    assert parsed[0].component_tree is not None
    assert "l1.header" in parsed[0].component_tree.component_names()


def test_g16log_frame_model_aggregates_from_component_tree_input():
    block = _last_frame_block(FIXTURES[0])
    parser = G16LogFileFrameParserV3Memory()
    tree = parser.build_component_tree(block)
    for component in tree.iter_components():
        component.parse_payload(block, tree)

    assert tree.validate_contracts() == []

    frame = G16LogFileFrameMemory.model_validate(
        {
            "frame_content": block,
            "component_tree_input": tree,
        }
    )
    expected = aggregate_g16log_v3_tree(tree)

    assert frame.component_tree is not None
    assert frame.keywords == expected["keywords"]
    assert frame.charge == expected["charge"]
    assert frame.multiplicity == expected["multiplicity"]
