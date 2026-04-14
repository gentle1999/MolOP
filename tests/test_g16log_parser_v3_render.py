from pathlib import Path
from typing import cast

from molop.io.logic.QM_parsers.G16LogFileParser import G16LogFileParserMemory


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "3-m-Py_anion_Opt.log"


def test_g16log_v3_render_contains_expected_gaussian_sections():
    file_content = FIXTURE.read_text()
    frame = G16LogFileParserMemory().parse(file_content)[-1]
    rendered = frame.render_fakeg()

    assert "orientation" in rendered.lower()
    assert "scf done" in rendered.lower()
    assert "frequencies --" in rendered.lower()
    assert "zero-point correction" in rendered.lower()
    assert "- thermochemistry -" in rendered.lower()
    assert "e (thermal)             cv                s" in rendered.lower()
    assert "                      1                      2                      3" in rendered
    assert "                      A                      A                      A" in rendered


def test_g16log_v3_can_render_specific_l716_child_nodes():
    file_content = FIXTURE.read_text()
    frame = G16LogFileParserMemory().parse(file_content)[-1]
    tree = frame.component_tree

    assert tree is not None

    force_render = tree.render_node("l716.forceconstants")
    mode_render = tree.render_node("l716.vibration.mode[0]")
    temp_render = tree.render_node("l716.thermochemistry.temperature")

    assert "frc consts" in force_render.lower()
    assert "frequencies --" in mode_render.lower()
    assert "red. masses --" in mode_render.lower()
    assert "temperature" in temp_render.lower()


def test_g16log_file_model_can_render_fakeg_for_all_frames():
    rendered = cast(str, G16LogFileParserMemory().parse(FIXTURE.read_text()).render_fakeg())

    assert "gaussian 16:" in rendered.lower()
    assert "symbolic z-matrix:" in rendered.lower()
    assert "charge =" in rendered.lower()
    assert "alpha electrons" not in rendered
    assert "beta electrons" not in rendered
    assert "orientation" in rendered.lower()
    assert "scf done" in rendered.lower()
    assert "- thermochemistry -" in rendered.lower()
    assert "e (thermal)             cv                s" in rendered.lower()
    assert "normal termination of gaussian" in rendered.lower()
    assert rendered.count("Standard orientation:") == 5
    assert rendered.count("Step number") == 5
    assert rendered.count("Item               Value     Threshold  Converged?") == 5
    assert rendered.count("Frequencies --") >= 1
    assert rendered.count("Normal termination of Gaussian") == 2
    assert "Job cpu time:" not in rendered


def test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry():
    rendered = cast(str, G16LogFileParserMemory().parse(FIXTURE.read_text()).render_fakeg())
    reparsed = G16LogFileParserMemory().parse(rendered)

    assert len(reparsed.frames) == 5
    assert reparsed.charge == -1
    assert reparsed.multiplicity == 1
    assert reparsed.frames[0].geometry_optimization_status is not None
    assert reparsed.frames[-1].vibrations is not None
    assert reparsed.frames[-1].thermal_informations is not None
