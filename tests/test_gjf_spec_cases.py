from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from molop.io import AutoParser
from molop.io.logic.qminput_parsers.GJFFileParser import GJFFileParserDisk


_G16GJF_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "g16gjf"


def _parse_first_frame(fixture_name: str):
    fixture_path = _G16GJF_FIXTURE_DIR / fixture_name
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]
    assert len(file_model) > 0
    return cast(Any, file_model[0])


def test_zmat_labeled_variable_blocks_are_applied_to_atom_lines() -> None:
    frame = _parse_first_frame("spec_zmat_variables_labeled.gjf")
    atom_specs = frame.molecule_specifications.molecule_fragments[0].atom_specifications
    assert len(atom_specs) == 3
    assert "r1" not in atom_specs[1].coords_part
    assert "r2" not in atom_specs[2].coords_part
    assert "a1" not in atom_specs[2].coords_part
    assert "0.9600" in atom_specs[1].coords_part
    assert "0.9700" in atom_specs[2].coords_part
    assert "104.5000" in atom_specs[2].coords_part


def test_zmat_unlabeled_variable_blocks_are_applied_by_position() -> None:
    frame = _parse_first_frame("spec_zmat_variables_unlabeled.gjf")
    atom_specs = frame.molecule_specifications.molecule_fragments[0].atom_specifications
    assert len(atom_specs) == 3
    assert "r1" not in atom_specs[1].coords_part
    assert "r2" not in atom_specs[2].coords_part
    assert "a1" not in atom_specs[2].coords_part
    assert "0.9600" in atom_specs[1].coords_part
    assert "0.9700" in atom_specs[2].coords_part
    assert "104.5000" in atom_specs[2].coords_part


def test_zmat_generated_cartesian_geometry_matches_input_variables() -> None:
    frame = _parse_first_frame("spec_zmat_variables_labeled.gjf")
    coords = frame.molecule_specifications.coords().magnitude
    assert coords.shape == (3, 3)

    oh1 = np.linalg.norm(coords[1] - coords[0])
    oh2 = np.linalg.norm(coords[2] - coords[0])
    assert abs(oh1 - 0.96) < 1e-3
    assert abs(oh2 - 0.97) < 1e-3

    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]
    cosang = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    angle = float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))
    assert abs(angle - 104.5) < 1e-2


def test_modredundant_additional_section_detected_as_typed_section() -> None:
    frame = _parse_first_frame("spec_modredundant_basic.gjf")
    assert frame.parsed_additional_sections
    first_section = frame.parsed_additional_sections[0]
    assert first_section.section_type == "modredundant"
    assert len(first_section.lines) == 2


def test_gic_additional_section_detected_as_typed_section() -> None:
    frame = _parse_first_frame("spec_gic_basic.gjf")
    assert frame.parsed_additional_sections
    first_section = frame.parsed_additional_sections[0]
    assert first_section.section_type == "gic"
    assert len(first_section.lines) == 3


def test_gic_lines_are_structurally_parsed() -> None:
    frame = _parse_first_frame("spec_gic_structured.gjf")
    first_section = frame.parsed_additional_sections[0]
    assert first_section.section_type == "gic"
    assert len(first_section.lines) == 4

    line0 = first_section.lines[0]
    assert line0.label == "rOH"
    assert line0.label_options == ["active", "min=0.8", "max=1.2"]
    assert [opt.key for opt in line0.parsed_label_options] == ["active", "min", "max"]
    assert line0.parsed_label_options[0].is_flag is True
    assert line0.parsed_label_options[1].value == "0.8"
    assert line0.parsed_label_options[2].value == "1.2"
    assert line0.expression_kind == "function"
    assert line0.function_name == "R"
    assert line0.function_args == ["1", "2"]

    line2 = first_section.lines[2]
    assert line2.label == "combo"
    assert line2.function_name == "DotDiff"
    assert line2.function_args == ["R(1,2)", "R(1,3)"]

    line3 = first_section.lines[3]
    assert line3.is_standalone_option is True
    assert line3.expression_kind == "standalone"
    assert line3.standalone_keyword == "Atom"
    assert line3.standalone_args == ["2", "Freeze"]
    assert line3.standalone_target == "2"
    assert line3.standalone_action == "freeze"
    assert line3.normalized_state == "frozen"


def test_render_prefers_structured_additional_sections() -> None:
    frame = _parse_first_frame("spec_gic_structured.gjf")
    frame.additional_sections = ""
    rendered = frame._render()
    assert "rOH(active, min=0.8, max=1.2)=R(1, 2)" in rendered
    assert "aHOH=A(2, 1, 3)" in rendered
    assert "combo=DotDiff(R(1,2), R(1,3))" in rendered
    assert "Atom 2 Freeze" in rendered


def test_render_accepts_structured_override_inputs() -> None:
    frame = _parse_first_frame("spec_gic_structured.gjf")
    rendered = frame._render(
        molecule_specifications=frame.molecule_specifications,
        parsed_additional_sections=frame.parsed_additional_sections,
    )
    assert "rOH(active, min=0.8, max=1.2)=R(1, 2)" in rendered
    assert "Atom 2 Freeze" in rendered


def test_zmat_trailing_marker_0_1_is_preserved_and_rendered_without_units() -> None:
    frame = _parse_first_frame("spec_zmat_alt_marker.gjf")
    atom_specs = frame.molecule_specifications.molecule_fragments[0].atom_specifications
    assert atom_specs[2].coords_part.endswith(" 0")
    assert atom_specs[3].coords_part.endswith(" 1")

    rendered = frame._render()
    assert "angstrom" not in rendered
    assert "degree" not in rendered
    assert "104.500000    0" in rendered
    assert "120.000000    1" in rendered


def test_zmat_trailing_marker_1_uses_second_angle_geometry_semantics() -> None:
    frame = _parse_first_frame("spec_zmat_alt_marker.gjf")
    coords = frame.molecule_specifications.coords().magnitude

    bond = np.linalg.norm(coords[3] - coords[0])
    assert abs(bond - 1.5) < 1e-3

    v_ref1 = coords[1] - coords[0]
    v_ref2 = coords[2] - coords[0]
    v_new = coords[3] - coords[0]

    angle1 = float(
        np.degrees(
            np.arccos(
                np.clip(
                    np.dot(v_new, v_ref1) / (np.linalg.norm(v_new) * np.linalg.norm(v_ref1)),
                    -1.0,
                    1.0,
                )
            )
        )
    )
    angle2 = float(
        np.degrees(
            np.arccos(
                np.clip(
                    np.dot(v_new, v_ref2) / (np.linalg.norm(v_new) * np.linalg.norm(v_ref2)),
                    -1.0,
                    1.0,
                )
            )
        )
    )

    assert abs(angle1 - 110.0) < 1e-2
    assert abs(angle2 - 120.0) < 1e-2


def test_free_format_delimiters_are_parsed_and_rendered_canonically() -> None:
    frame = _parse_first_frame("spec_free_format_delimiters.gjf")
    atom_specs = frame.molecule_specifications.molecule_fragments[0].atom_specifications
    assert len(atom_specs) == 3
    assert "0.9600" in atom_specs[1].coords_part
    assert "0.9700" in atom_specs[2].coords_part
    assert "104.5000" in atom_specs[2].coords_part

    molecule_render = frame.molecule_specifications._render()
    assert "," not in molecule_render
    assert "\t" not in molecule_render
    assert "/" not in molecule_render


def test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic() -> None:
    frame = _parse_first_frame("spec_mixed_gic_modredundant_invalid.gjf")
    assert frame.parsed_additional_sections
    first_section = frame.parsed_additional_sections[0]
    assert first_section.section_type == "unknown"
    assert frame.additional_section_diagnostics
    assert any(
        diagnostic.section_type == "mixed-additional-section"
        and "mixes GIC and ModRedundant" in diagnostic.message
        for diagnostic in frame.additional_section_diagnostics
    )


def test_nbo_additional_section_detected_as_typed_section() -> None:
    frame = _parse_first_frame("spec_pop_nbo_section.gjf")
    assert frame.route_section.semantic_route.pop_options.nbo_read is True
    assert frame.route_section.semantic_route.pop_options.hirshfeld is True
    assert frame.route_section.semantic_route.pop_options.orbitals == 3
    assert frame.parsed_additional_sections
    first_section = frame.parsed_additional_sections[0]
    assert first_section.section_type == "nbo"
    assert first_section.header == "$NBO"
    assert first_section.commands == ["FILE=water", "PRINT=3"]
    assert first_section.footer == "$END"
    rendered = frame._render()
    assert "$NBO" in rendered
    assert "FILE=water" in rendered
    assert "$END" in rendered


def test_multifragment_assignments_match_declared_pairs() -> None:
    frame = _parse_first_frame("spec_fragments_valid.gjf")
    assert len(frame.molecule_specifications.molecule_fragments) == 2
    assert [frag.fragment_id for frag in frame.molecule_specifications.molecule_fragments] == [0, 1]
    assert [
        len(frag.atom_specifications) for frag in frame.molecule_specifications.molecule_fragments
    ] == [1, 1]


@pytest.mark.parametrize(
    ("fixture_name", "message"),
    [
        (
            "spec_fragments_invalid_missing_assignment.gjf",
            "require every atom to declare Fragment=n",
        ),
        (
            "spec_fragments_invalid_gap.gjf",
            "Fragment assignments must be contiguous",
        ),
    ],
)
def test_multifragment_invalid_assignments_raise_clear_errors(
    fixture_name: str, message: str
) -> None:
    fixture_path = _G16GJF_FIXTURE_DIR / fixture_name
    with pytest.raises(ValueError, match=message):
        GJFFileParserDisk().parse(str(fixture_path))


@pytest.mark.parametrize(
    ("fixture_name", "message"),
    [
        ("spec_zmat_invalid_bond_length.gjf", "bond length must be positive"),
        ("spec_zmat_invalid_angle.gjf", "bond angle must be between 0 and 180 degrees"),
        (
            "spec_zmat_invalid_alt_second_angle.gjf",
            "Alternate Z-matrix second angle must be between 0 and 180 degrees",
        ),
        (
            "spec_zmat_invalid_forward_reference.gjf",
            "must refer to a previously defined atom",
        ),
        (
            "spec_zmat_invalid_duplicate_reference.gjf",
            "distance and angle references on atom 3 must be different",
        ),
    ],
)
def test_invalid_zmatrix_domains_raise_clear_errors(fixture_name: str, message: str) -> None:
    fixture_path = _G16GJF_FIXTURE_DIR / fixture_name
    with pytest.raises(ValueError, match=message):
        GJFFileParserDisk().parse(str(fixture_path))
