"""CPU usage and elapsed duration remain independent throughout QM I/O."""

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from molop.io import AutoParser
from molop.io.base_models.ChemFile import BaseCalcFile
from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.base_models.ParseContainers import TextParseContext
from molop.io.base_models.timing import accumulate_calculation_timings
from molop.io.logic.gaussian.log.frame_models.G16Components import G16JobCPUComponent
from molop.io.logic.gaussian.log.frame_parsers._g16_shared import _parse_running_time
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import extract_g16_running_time
from molop.io.logic.gaussian.log.parsers._g16_timing import extract_g16_timings
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.log.parsers._orca_log_shared import extract_orca_cpu_timing
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory
from molop.io.logic.xtb.output.parsers._xtb_output_extractors import extract_xtb_timings
from molop.io.logic.xtb.output.parsers.XTBOutputFileParser import XTBOutputFileParserMemory
from molop.unit import atom_ureg


FIXTURES = Path(__file__).parent / "test_files"
JOB_TIMES = "Job cpu time: 0 days 0 hours 0 minutes 100.0 seconds.\nElapsed time: 0 days 0 hours 0 minutes 10.0 seconds.\n"
LINK_TIMES = "Leave Link 502, MaxMem= 100 cpu: 20.0 elap: 2.0\nLeave Link 601, MaxMem= 100 cpu: 30.0 elap: 3.0\n"


@pytest.mark.parametrize("text,cpu,wall", [(JOB_TIMES, 100, 10), (LINK_TIMES, 50, 5)])
def test_gaussian_cpu_and_wall_time_are_never_added(text, cpu, wall):
    timings = extract_g16_timings(text)
    assert timings["cpu_time"].m_as("second") == cpu
    assert timings["wall_time"].m_as("second") == wall
    assert extract_g16_running_time(TextParseContext(text)).m_as("second") == wall
    assert _parse_running_time(text).m_as("second") == wall


def test_gaussian_job_totals_and_link_rows_are_alternative_scopes():
    file_times = extract_g16_timings(LINK_TIMES + JOB_TIMES)
    frame_times = extract_g16_timings(LINK_TIMES + JOB_TIMES, prefer_job_totals=False)
    assert file_times["cpu_time"].m_as("second") == 100
    assert file_times["wall_time"].m_as("second") == 10
    assert frame_times["cpu_time"].m_as("second") == 50
    assert frame_times["wall_time"].m_as("second") == 5


@pytest.mark.parametrize("kind", ["Job cpu time", "Elapsed time"])
@pytest.mark.parametrize("seconds", ["0", "0.0", "12.5"])
def test_gaussian_single_time_kind_and_zero_do_not_invent_the_other(kind, seconds):
    text = f"{kind}: 0 days 0 hours 0 minutes {seconds} seconds.\n"
    observed = "cpu_time" if kind == "Job cpu time" else "wall_time"
    absent = "wall_time" if observed == "cpu_time" else "cpu_time"
    timings = extract_g16_timings(text)
    assert timings[observed].m_as("second") == float(seconds)
    assert timings[absent] is None


def test_gaussian_job_time_converts_day_hour_minute_units_independently():
    timings = extract_g16_timings(
        "Job cpu time: 2 days 3 hours 4 minutes 5.5 seconds.\nElapsed time: 1 days 2 hours 3 minutes 4.5 seconds.\n"
    )
    assert timings["cpu_time"].m_as("second") == 183845.5
    assert timings["wall_time"].m_as("second") == 93784.5


@pytest.mark.parametrize("mode", ["memory", "disk", "last_frame"])
def test_gaussian_real_multistep_totals_follow_selected_segment_scope(mode):
    path = FIXTURES / "g16log/H2O.log"
    if mode == "memory":
        parsed = G16LogFileParserMemory().parse(path.read_text(encoding="utf-8"))
    else:
        parsed = AutoParser(str(path), n_jobs=1, only_last_frame=mode == "last_frame")[0]
    assert parsed.cpu_time.m_as("second") == pytest.approx(46.9 if mode == "last_frame" else 96.7)
    assert parsed.wall_time.m_as("second") == pytest.approx(2.8 if mode == "last_frame" else 5.8)
    assert parsed.running_time == parsed.wall_time
    assert parsed[-1].cpu_time.m_as("second") == pytest.approx(46.5)
    assert parsed[-1].wall_time.m_as("second") == pytest.approx(2.6)
    assert parsed[-1].running_time == parsed[-1].wall_time


def test_gaussian_cpu_only_log_keeps_wall_time_missing():
    content = (FIXTURES / "g16log/H2O.log").read_text(encoding="utf-8")
    content = "\n".join(
        line
        for line in content.splitlines()
        if "Elapsed time:" not in line and "Leave Link" not in line
    )
    parsed = G16LogFileParserMemory().parse(content)
    assert parsed.cpu_time.m_as("second") == pytest.approx(96.7)
    assert parsed.wall_time is None
    assert parsed.running_time is None
    assert parsed[-1].cpu_time.m_as("second") == pytest.approx(46.9)
    assert parsed[-1].wall_time is None
    rendered = parsed[-1].render_fakeg()
    assert "Job cpu time:" in rendered
    assert "Elapsed time:" not in rendered


def test_gaussian_incomplete_segment_time_is_not_reported_as_a_complete_file_total():
    content = (FIXTURES / "g16log/H2O.log").read_text(encoding="utf-8")
    boundary = "Link1:  Proceeding to internal job step number  2."
    first, second = content.split(boundary, 1)
    second = "\n".join(
        line
        for line in second.splitlines()
        if "Elapsed time:" not in line and "Leave Link" not in line
    )
    parsed = G16LogFileParserMemory().parse(first + boundary + second)
    assert parsed.cpu_time.m_as("second") == pytest.approx(96.7)
    assert parsed.wall_time is None
    assert parsed.running_time is None


@pytest.mark.parametrize("model_class", [BaseCalcFrame, BaseCalcFile])
def test_public_time_aliases_validate_units_and_reject_conflicting_wall_values(model_class):
    model = model_class.model_validate(
        {"cpu_time": 3 * atom_ureg.second, "wall_time": 1 * atom_ureg.minute},
        context={"force_unit_transform": True},
    )
    assert model.running_time.m_as("second") == 60
    assert model.wall_time.m_as("second") == 60
    assert model.cpu_time.m_as("second") == 3
    dumped = model.to_unitless_dump_with_unit_keys(exclude_computed_fields=True)
    assert dumped["running_time (second)"] == dumped["wall_time (second)"] == 60
    assert dumped["cpu_time (second)"] == 3
    json.dumps(dumped, allow_nan=False)
    assert model_class(running_time=0 * atom_ureg.second).wall_time.m_as("second") == 0
    assert model_class(cpu_time=1 * atom_ureg.second).running_time is None
    with pytest.raises(ValidationError, match="same wall-clock duration"):
        model_class(running_time=1 * atom_ureg.second, wall_time=2 * atom_ureg.second)


@pytest.mark.parametrize("model_class", [BaseCalcFrame, BaseCalcFile])
def test_time_aliases_respect_disabled_unit_transformation(model_class):
    model = model_class.model_validate(
        {"cpu_time": 3 * atom_ureg.minute, "wall_time": 1 * atom_ureg.minute},
        context={"force_unit_transform": False},
    )
    assert model.running_time == model.wall_time
    assert model.wall_time.units == model.cpu_time.units == atom_ureg.minute


def test_segment_totals_keep_missing_quantities_unknown_and_preserve_zero():
    segments = [
        {"cpu_time": 100 * atom_ureg.second, "wall_time": 10 * atom_ureg.second},
        {"cpu_time": None, "wall_time": 0 * atom_ureg.second},
    ]
    totals = accumulate_calculation_timings(segments)
    assert totals["cpu_time"] is None
    assert totals["running_time"] == totals["wall_time"] == 10 * atom_ureg.second
    segments[1]["cpu_time"] = 20 * atom_ureg.second
    assert accumulate_calculation_timings(segments)["cpu_time"] == 120 * atom_ureg.second
    segments[1]["wall_time"] = None
    totals = accumulate_calculation_timings(segments)
    assert totals["running_time"] is totals["wall_time"] is None


def test_file_summary_preserves_zero_cpu_and_wall_time_as_distinct_columns():
    parsed = BaseCalcFile(cpu_time=0 * atom_ureg.second, wall_time=0 * atom_ureg.second)
    parsed.append(BaseCalcFrame())
    summary = parsed.to_summary_dict()
    assert summary[("Status", "RuningTime", "second")] == 0
    assert summary[("Status", "CPUTime", "second")] == 0
    assert summary[("Status", "WallTime", "second")] == 0


def test_xtb_total_times_do_not_add_cpu_or_component_durations():
    content = "module:\n * wall-time: 0 d, 0 h, 0 min, 3 sec\n * cpu-time: 0 d, 0 h, 0 min, 20 sec\ntotal:\n * wall-time: 0 d, 0 h, 0 min, 10 sec\n * cpu-time: 0 d, 0 h, 0 min, 100 sec\nmodule:\n * wall-time: 0 d, 0 h, 0 min, 1 sec\n * cpu-time: 0 d, 0 h, 0 min, 2 sec\n"
    timings = extract_xtb_timings(content)
    assert timings["cpu_time"].m_as("second") == 100
    assert timings["wall_time"].m_as("second") == 10
    assert timings["running_time"] == timings["wall_time"]


@pytest.mark.parametrize("kind", ["cpu", "wall"])
def test_xtb_missing_total_kind_does_not_reuse_a_component_time(kind):
    other = "cpu" if kind == "wall" else "wall"
    content = f"module:\n * {kind}-time: 0 d, 0 h, 0 min, 123 sec\ntotal:\n * {other}-time: 0 d, 0 h, 0 min, 0 sec\n"
    timings = extract_xtb_timings(content)
    assert timings[f"{kind}_time"] is None
    assert timings[f"{other}_time"].m_as("second") == 0


@pytest.mark.parametrize("copies", [1, 2])
def test_xtb_real_cpu_and_wall_totals_survive_file_and_frame_parsing(copies):
    content = (FIXTURES / "xtbout/water_raman_ptb_6_7_1.out").read_text(encoding="utf-8")
    parsed = XTBOutputFileParserMemory().parse("\n".join([content] * copies))
    assert parsed.cpu_time.m_as("second") == pytest.approx(copies * 0.039)
    assert parsed.wall_time.m_as("second") == pytest.approx(copies * 0.136)
    assert parsed.running_time == parsed.wall_time
    assert parsed[-1].cpu_time.m_as("second") == pytest.approx(0.039)
    assert parsed[-1].wall_time.m_as("second") == pytest.approx(0.136)


def test_xtb_missing_cpu_total_in_one_run_does_not_reuse_the_first_run_total():
    content = (FIXTURES / "xtbout/water_raman_ptb_6_7_1.out").read_text(encoding="utf-8")
    missing_cpu = "\n".join(line for line in content.splitlines() if "0.039 sec" not in line)
    assert missing_cpu != content
    parsed = XTBOutputFileParserMemory().parse(content + "\n" + missing_cpu)
    assert parsed.cpu_time is None
    assert parsed.wall_time.m_as("second") == pytest.approx(0.272)
    assert parsed[-1].cpu_time is None


def test_orca_serial_input_estimates_cpu_time_without_changing_wall_time():
    content = (FIXTURES / "orca/output_files/local/h2_grad_orca.out").read_text(encoding="utf-8")
    parsed = ORCALogFileParserMemory().parse(content)
    assert parsed.running_time == parsed.wall_time
    assert parsed.wall_time.m_as("second") > 0
    assert parsed.cpu_time == parsed.wall_time
    assert parsed.cpu_time_source == "estimated"
    assert parsed.cpu_time_num_processes == 1
    assert parsed[-1].wall_time == parsed[-1].running_time
    assert parsed[-1].cpu_time == parsed[-1].wall_time
    assert parsed[-1].cpu_time_source == "estimated"


@pytest.mark.parametrize(
    "text,kwargs,count",
    [
        ("Program running with 4 parallel MPI-processes", {"requested_num_cpu": 8}, 4),
        ("Program running with 2 parallel MPI-processes", {"inherited_num_cpu": 4}, 2),
        ("", {"requested_num_cpu": 6}, 6),
        ("", {"keywords": "PBE def2-SVP PAL8"}, 8),
        ("", {"serial_input": True}, 1),
        ("", {}, None),
        ("", {"requested_num_cpu": 0, "serial_input": True}, None),
    ],
)
def test_orca_cpu_fallback_uses_known_positive_process_count(text, kwargs, count):
    result = extract_orca_cpu_timing(text, 10 * atom_ureg.second, **kwargs)
    assert result["cpu_time_num_processes"] == count
    if count is None:
        assert result["cpu_time"] is result["cpu_time_source"] is None
    else:
        assert result["cpu_time"].m_as("second") == 10 * count
        assert result["cpu_time_source"] == "estimated"


def test_orca_cpu_fallback_does_not_replace_observed_cpu_time_or_invent_wall_time():
    result = extract_orca_cpu_timing(
        "", 10 * atom_ureg.second, cpu_time=3 * atom_ureg.second, requested_num_cpu=4
    )
    assert result["cpu_time"].m_as("second") == 3
    assert result["cpu_time_source"] == "observed"
    missing_wall = extract_orca_cpu_timing("", None, requested_num_cpu=4)
    assert missing_wall["cpu_time"] is missing_wall["cpu_time_source"] is None
    zero_wall = extract_orca_cpu_timing("", 0 * atom_ureg.second, requested_num_cpu=4)
    assert zero_wall["cpu_time"].m_as("second") == 0
    assert zero_wall["cpu_time_source"] == "estimated"


@pytest.mark.parametrize("mode", ["memory", "disk"])
def test_orca_parallel_cpu_estimate_reaches_files_frames_and_exports(mode):
    path = FIXTURES / "orca/output_files/cclib/basicORCA6.0/dvb_rocis.out"
    parsed = (
        ORCALogFileParserMemory().parse(path.read_text(encoding="utf-8"))
        if mode == "memory"
        else AutoParser(str(path), parser_detection="orcaout", n_jobs=1)[0]
    )
    for obj in (parsed, parsed[-1]):
        assert obj.wall_time.m_as("second") == pytest.approx(30.501)
        assert obj.cpu_time.m_as("second") == pytest.approx(122.004)
        assert obj.running_time == obj.wall_time
        assert obj.cpu_time_source == "estimated"
        assert obj.cpu_time_num_processes == 4
        dumped = obj.to_unitless_dump_with_unit_keys()
        assert dumped["cpu_time_source"] == "estimated"
        assert dumped["cpu_time_num_processes"] == 4
        json.dumps(dumped, allow_nan=False)
    summary = parsed.to_summary_dict()
    assert summary[("Status", "CPUTimeSource", "")] == "estimated"
    assert summary[("Status", "CPUTimeNumProcesses", "")] == 4


def test_orca_file_sums_each_segment_estimate_using_its_own_process_count():
    serial = (FIXTURES / "orca/output_files/local/h2_grad_orca.out").read_text(encoding="utf-8")
    parallel = (FIXTURES / "orca/output_files/cclib/basicORCA6.0/dvb_rocis.out").read_text(
        encoding="utf-8"
    )
    parsed = ORCALogFileParserMemory().parse(
        "******** JOB NUMBER 1 ********\n"
        + serial
        + "\n******** JOB NUMBER 2 ********\n"
        + parallel
    )
    assert parsed.wall_time.m_as("second") == pytest.approx(31.821)
    assert parsed.cpu_time.m_as("second") == pytest.approx(123.324)
    assert parsed.cpu_time_source == "estimated"
    assert parsed.cpu_time_num_processes is None
    assert parsed[-1].cpu_time_num_processes == 4


def test_aggregating_measured_and_estimated_cpu_time_preserves_estimate_marker():
    result = accumulate_calculation_timings(
        [
            {
                "cpu_time": 2 * atom_ureg.second,
                "wall_time": 1 * atom_ureg.second,
                "cpu_time_source": "observed",
            },
            {
                "cpu_time": 8 * atom_ureg.second,
                "wall_time": 2 * atom_ureg.second,
                "cpu_time_source": "estimated",
                "cpu_time_num_processes": 4,
            },
        ]
    )
    assert result["cpu_time"].m_as("second") == 10
    assert result["wall_time"].m_as("second") == 3
    assert result["cpu_time_source"] == "estimated"
    assert result["cpu_time_num_processes"] is None


def test_fakeg_retains_cpu_estimation_provenance():
    frame = BaseCalcFrame(
        cpu_time=40 * atom_ureg.second,
        wall_time=10 * atom_ureg.second,
        cpu_time_source="estimated",
        cpu_time_num_processes=4,
    )
    component = G16JobCPUComponent(payload=G16JobCPUComponent.build_payload_from_frame(frame))
    rendered = component.render_fakeg()
    timings = extract_g16_timings(rendered)
    assert timings["cpu_time"].m_as("second") == 40
    assert timings["wall_time"].m_as("second") == 10
    assert timings["cpu_time_source"] == "estimated"
    assert timings["cpu_time_num_processes"] == 4


def test_fakeg_does_not_infer_uniform_process_count_when_one_marker_lacks_it():
    text = JOB_TIMES + "MolOP CPU time source: estimated; processes=4\n"
    text += JOB_TIMES + "MolOP CPU time source: estimated\n"
    timings = extract_g16_timings(text)
    assert timings["cpu_time_source"] == "estimated"
    assert timings["cpu_time_num_processes"] is None


def test_fakeg_cpu_estimation_provenance_reaches_reparsed_files_and_frames():
    parsed = G16LogFileParserMemory().parse(
        (FIXTURES / "g16log/water_raman_hf_sto3g.log").read_text(encoding="utf-8")
    )
    frame = parsed[-1]
    frame.cpu_time = 40 * atom_ureg.second
    frame.cpu_time_source = "estimated"
    frame.cpu_time_num_processes = 4
    reparsed = G16LogFileParserMemory().parse(parsed.render_fakeg(frame=[len(parsed) - 1]))
    for obj in (reparsed, reparsed[-1]):
        assert obj.cpu_time.m_as("second") == 40
        assert obj.cpu_time_source == "estimated"
        assert obj.cpu_time_num_processes == 4


def test_fakeg_time_labels_preserve_values_and_emit_each_kind_once():
    parsed = G16LogFileParserMemory().parse(
        (FIXTURES / "g16log/water_raman_hf_sto3g.log").read_text(encoding="utf-8")
    )
    frame = parsed[-1]
    rendered = frame.render_fakeg()
    assert rendered.count("Job cpu time:") == rendered.count("Elapsed time:") == 1
    timings = extract_g16_timings(rendered)
    assert timings["cpu_time"].m_as("second") == pytest.approx(
        frame.cpu_time.m_as("second"), abs=0.005
    )
    assert timings["wall_time"].m_as("second") == pytest.approx(
        frame.wall_time.m_as("second"), abs=0.005
    )
    component = G16JobCPUComponent(payload={"wall_time": 0 * atom_ureg.second})
    assert "Elapsed time:" in component.render_fakeg()
    assert "Job cpu time:" not in component.render_fakeg()
