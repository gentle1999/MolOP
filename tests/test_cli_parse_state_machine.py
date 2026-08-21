from __future__ import annotations

import csv
import json
from pathlib import Path

import click
import pytest
from click.shell_completion import _resolve_context, _resolve_incomplete
from click.testing import CliRunner

import molop.cli.state_machine as state_machine
from molop.cli.app import app
from molop.cli.state_machine import (
    OPERATION_REGISTRY,
    BatchInputConfig,
    CliUsageError,
    FilterStateParams,
    OperationCall,
    ToSummaryDfParams,
    build_plan,
    parse_dynamic_options,
)
from molop.config import molopconfig
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.parse_outcomes import BatchParseResult, FileParseOutcome


runner = CliRunner()


def test_top_level_only_exposes_parse_business_command() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "parse" in result.stdout
    for old_command in (
        "summary",
        "transform",
        "visualize",
        "filter-by-codec",
        "sample",
        "stats",
        "groupby",
    ):
        assert old_command not in result.stdout


def test_global_max_jobs_sets_process_worker_ceiling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(molopconfig, "max_jobs", None)
    monkeypatch.setattr(state_machine, "AutoParser", lambda *_args, **_kwargs: FileBatchModelDisk())

    result = runner.invoke(app, ["--max-jobs", "3", "-q", "parse", "dummy.xyz"])

    assert result.exit_code == 0, result.output
    assert molopconfig.max_jobs == 3


def test_global_max_jobs_rejects_non_positive_values() -> None:
    result = runner.invoke(app, ["--max-jobs", "0", "parse", "dummy.xyz"])

    assert result.exit_code != 0
    assert "x>=1" in result.output


def test_parse_operation_is_registered_as_click_command() -> None:
    result = runner.invoke(app, ["parse", "dummy.log", "format-transform", "--help"])

    assert result.exit_code == 0
    assert "--format" in result.stdout
    assert "--output-dir" in result.stdout


def test_parse_help_exposes_public_parse_options() -> None:
    result = runner.invoke(app, ["parse", "--help"])

    assert result.exit_code == 0
    for option in (
        "--input",
        "--total-charge",
        "--total-multiplicity",
        "--only-extract-structure",
        "--only-last-frame",
        "--capture-source-evidence",
        "--source-encoding",
        "--keep-file-content",
        "--force-unit-transform",
        "--graph-reconstruction-backend",
        "--make-dative-bonds",
        "--report",
    ):
        assert option in result.stdout


def test_dynamic_format_options_parse_generic_click_extras() -> None:
    assert parse_dynamic_options(
        ["--graph-policy", "prefer", "--strict", "--no-cache", "--charge", "0"]
    ) == {
        "graph_policy": "prefer",
        "strict": True,
        "cache": False,
        "charge": 0,
    }


def test_format_transform_completes_writer_formats() -> None:
    command = app.commands["parse"].commands["format-transform"]
    format_option = next(param for param in command.params if param.name == "target_format")
    ctx = click.Context(command)

    completions = format_option.shell_complete(ctx, "g")

    assert [item.value for item in completions] == ["gjf"]


def test_format_transform_completes_dynamic_options_for_target_format() -> None:
    command = app.commands["parse"].commands["format-transform"]
    extra_arg = next(param for param in command.params if param.name == "extra_args")
    ctx = click.Context(command)
    ctx.params["target_format"] = "sdf"
    ctx.args = ["--format", "sdf"]

    completions = extra_arg.shell_complete(ctx, "--e")

    assert [item.value for item in completions] == ["--engine"]
    assert (
        completions[0].help
        == "Rendering backend used to generate the target format. Default: 'rdkit'."
    )


def test_format_transform_completes_orcainp_writer_dynamic_options() -> None:
    command = app.commands["parse"].commands["format-transform"]
    extra_arg = next(param for param in command.params if param.name == "extra_args")
    ctx = click.Context(command)
    ctx.params["target_format"] = "orcainp"
    ctx.args = ["--format", "orcainp"]

    completions = extra_arg.shell_complete(ctx, "--max")

    assert [item.value for item in completions] == ["--maxcore"]


def test_format_transform_completes_dynamic_option_values() -> None:
    command = app.commands["parse"].commands["format-transform"]
    extra_arg = next(param for param in command.params if param.name == "extra_args")
    ctx = click.Context(command)
    ctx.params["target_format"] = "sdf"
    ctx.args = ["--format", "sdf", "--engine"]

    completions = extra_arg.shell_complete(ctx, "o")

    assert [item.value for item in completions] == ["openbabel"]
    assert completions[0].help == "Use Open Babel for rendering."


def test_format_transform_dynamic_completion_protocol_includes_help() -> None:
    args = [
        "parse",
        "dummy.log",
        "format-transform",
        "--format",
        "sdf",
    ]
    ctx = _resolve_context(app, {}, "molop", args.copy())
    obj, incomplete = _resolve_incomplete(ctx, args.copy(), "--e")

    completions = obj.shell_complete(ctx, incomplete)
    completions_by_value = {item.value: item for item in completions}

    assert "--embed" in completions_by_value
    assert "--engine" in completions_by_value
    assert (
        completions_by_value["--engine"].help
        == "Rendering backend used to generate the target format. Default: 'rdkit'."
    )


def test_format_transform_completes_dynamic_option_values_via_click_protocol() -> None:
    args = [
        "parse",
        "dummy.log",
        "format-transform",
        "--format",
        "sdf",
        "--engine",
    ]
    ctx = _resolve_context(app, {}, "molop", args.copy())
    obj, incomplete = _resolve_incomplete(ctx, args.copy(), "o")

    completions = obj.shell_complete(ctx, incomplete)

    assert [item.value for item in completions] == ["openbabel"]


def test_completion_show_prints_shell_source() -> None:
    result = runner.invoke(app, ["completion", "show", "--shell", "bash"])

    assert result.exit_code == 0
    assert "_molop_completion" in result.stdout
    assert "_MOLOP_COMPLETE=bash_complete" in result.stdout


def test_completion_install_writes_shell_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))

    result = runner.invoke(app, ["completion", "install", "--shell", "fish"])

    assert result.exit_code == 0
    target = tmp_path / ".config" / "fish" / "completions" / "molop.fish"
    assert target.exists()
    assert "_MOLOP_COMPLETE=fish_complete" in target.read_text(encoding="utf-8")


def test_completion_install_overwrites_existing_shell_registration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    bashrc = tmp_path / ".bashrc"
    bashrc.write_text(
        "\n".join(
            [
                "export KEEP_ME=1",
                "# >>> molop completion >>>",
                "source /tmp/old-molop-completion.bash",
                "# <<< molop completion <<<",
                "alias ll='ls -la'",
                "",
            ]
        ),
        encoding="utf-8",
    )
    target = tmp_path / ".config" / "molop" / "completions" / "bash" / "molop.bash"
    target.parent.mkdir(parents=True)
    target.write_text("old completion source", encoding="utf-8")

    first = runner.invoke(app, ["completion", "install", "--shell", "bash"])
    second = runner.invoke(app, ["completion", "install", "--shell", "bash"])

    assert first.exit_code == 0
    assert second.exit_code == 0
    assert "_MOLOP_COMPLETE=bash_complete" in target.read_text(encoding="utf-8")

    rc_source = bashrc.read_text(encoding="utf-8")
    assert "export KEEP_ME=1" in rc_source
    assert "alias ll='ls -la'" in rc_source
    assert "/tmp/old-molop-completion.bash" not in rc_source
    assert rc_source.count("# >>> molop completion >>>") == 1
    assert str(target) in rc_source


def test_registry_marks_chainable_and_terminal_operations() -> None:
    assert OPERATION_REGISTRY["filter-state"].return_kind == "batch"
    assert OPERATION_REGISTRY["sample"].return_kind == "batch"
    assert OPERATION_REGISTRY["format-transform"].return_kind == "terminal"
    assert OPERATION_REGISTRY["to-summary-df"].return_kind == "terminal"


def test_terminal_operation_cannot_be_followed_by_batch_operation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    called = False

    def _fail_if_called(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("AutoParser must not be called for invalid plans")

    monkeypatch.setattr(state_machine, "AutoParser", _fail_if_called)

    with pytest.raises(CliUsageError, match="must be the last operation"):
        build_plan(
            BatchInputConfig(pattern="tests/test_files/orca/single_point_inputs/h2_grad_orca.inp"),
            [
                OperationCall(
                    spec=OPERATION_REGISTRY["to-summary-df"],
                    params=ToSummaryDfParams(mode="frame"),
                ),
                OperationCall(
                    spec=OPERATION_REGISTRY["filter-state"],
                    params=FilterStateParams(state="normal"),
                ),
            ],
        )

    assert called is False


def test_parse_options_are_forwarded_as_one_immutable_snapshot(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _capture_auto_parser(file_path, **kwargs):
        captured["file_path"] = file_path
        captured.update(kwargs)
        return FileBatchModelDisk()

    monkeypatch.setattr(state_machine, "AutoParser", _capture_auto_parser)
    config = BatchInputConfig(
        pattern="first.log",
        additional_patterns=("second.out", "inputs/*.xyz"),
        parser_detection="g16log",
        n_jobs=3,
        total_charge=-1,
        total_multiplicity=2,
        only_extract_structure=True,
        only_last_frame=True,
        capture_source_evidence=True,
        source_encoding="utf-16",
        release_file_content=False,
        force_unit_transform=True,
        graph_reconstruction_backend="python",
        reconstruction_failure_policy="return_suspicious",
        make_dative_bonds=False,
    )

    state_machine.execute_plan(build_plan(config, []))

    assert captured["file_path"] == ["first.log", "second.out", "inputs/*.xyz"]
    assert captured["n_jobs"] == 3
    assert captured["parser_detection"] == "g16log"
    assert captured["return_report"] is False
    options = captured["parse_options"]
    assert options.total_charge == -1
    assert options.total_multiplicity == 2
    assert options.only_extract_structure is True
    assert options.only_last_frame is True
    assert options.capture_source_evidence is True
    assert options.source_encoding == "utf-16"
    assert options.release_file_content is False
    assert options.force_unit_transform is True
    assert options.graph_reconstruction_backend == "python"
    assert options.reconstruction_failure_policy == "return_suspicious"
    assert options.make_dative_bonds is False


def test_parse_report_cannot_be_chained_before_reading_files(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    called = False

    def _fail_if_called(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("AutoParser must not be called for invalid plans")

    monkeypatch.setattr(state_machine, "AutoParser", _fail_if_called)

    with pytest.raises(CliUsageError, match="--report cannot be combined"):
        build_plan(
            BatchInputConfig(pattern="dummy.log", report=True),
            [
                OperationCall(
                    spec=OPERATION_REGISTRY["sample"],
                    params=state_machine.SampleParams(n=1),
                )
            ],
        )

    assert called is False


def test_parse_report_renders_all_input_outcomes_as_json(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    batch = FileBatchModelDisk()
    report = BatchParseResult(
        batch=batch,
        outcomes=(
            FileParseOutcome(
                file_path="/tmp/first.xyz",
                status="ok",
                detected_format="xyz",
                input_index=0,
            ),
            FileParseOutcome(
                file_path="/tmp/missing.xyz",
                status="missing",
                input_index=1,
            ),
        ),
    )
    captured: dict[str, object] = {}

    def _return_report(file_path, **kwargs):
        captured["file_path"] = file_path
        captured.update(kwargs)
        return report

    monkeypatch.setattr(state_machine, "AutoParser", _return_report)

    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "first.xyz",
            "--input",
            "missing.xyz",
            "--report",
            "--output-format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["file_path"] == ["first.xyz", "missing.xyz"]
    assert captured["return_report"] is True
    payload = json.loads(result.stdout)
    assert payload["summary"] == {"failed": 1, "succeeded": 1, "total": 2}
    assert [outcome["status"] for outcome in payload["outcomes"]] == ["ok", "missing"]


def test_to_summary_df_params_preserve_all_frame_and_options() -> None:
    params = ToSummaryDfParams(
        frame="all",
        n_jobs=2,
        brief=False,
        flatten_columns=True,
        on_missing_frame="error",
    )

    kwargs = params.method_kwargs(BatchInputConfig(pattern="dummy.log", n_jobs=8))

    assert kwargs == {
        "mode": "frame",
        "frame": "all",
        "n_jobs": 2,
        "brief": False,
        "flatten_columns": True,
        "on_missing_frame": "error",
    }


def test_operation_n_jobs_defaults_to_minus_one_independently_of_parse() -> None:
    input_config = BatchInputConfig(pattern="dummy.log", n_jobs=8)

    assert ToSummaryDfParams().method_kwargs(input_config)["n_jobs"] == -1
    assert (
        state_machine.FormatTransformParams(format="xyz").method_kwargs(input_config)["n_jobs"]
        == -1
    )


def test_parse_filter_sample_outputs_batch_paths_json() -> None:
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp",
            "--parser-detection",
            "orcainp",
            "--n-jobs",
            "1",
            "--output-format",
            "json",
            "filter-by-codec",
            "--codec-id",
            "orcainp",
            "sample",
            "--n",
            "1",
            "--seed",
            "1",
        ],
    )

    assert result.exit_code == 0
    paths = json.loads(result.stdout)
    assert len(paths) == 1
    assert paths[0].endswith("h2_grad_orca.inp")


def test_parse_overrides_charge_multiplicity_and_retains_only_last_frame() -> None:
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/g16log/1.log",
            "--parser-detection",
            "g16log",
            "--n-jobs",
            "1",
            "--total-charge",
            "-1",
            "--total-multiplicity",
            "2",
            "--only-last-frame",
            "to-summary-df",
            "--frame",
            "all",
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    rows = json.loads(result.stdout)
    assert len(rows) == 1
    assert rows[0]["General.Charge"] == -1
    assert rows[0]["General.Multiplicity"] == 2


def test_parse_format_transform_writes_output(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp",
            "--parser-detection",
            "orcainp",
            "--n-jobs",
            "1",
            "format-transform",
            "--format",
            "xyz",
            "--output-dir",
            str(out_dir),
            "--graph-policy",
            "prefer",
        ],
    )

    assert result.exit_code == 0
    assert result.stdout == ""
    assert (out_dir / "h2_grad_orca.xyz").exists()


def test_parse_format_transform_no_write_overrides_output_dir(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp",
            "--parser-detection",
            "orcainp",
            "--n-jobs",
            "1",
            "format-transform",
            "--format",
            "xyz",
            "--output-dir",
            str(out_dir),
            "--no-write",
            "--graph-policy",
            "prefer",
        ],
    )

    assert result.exit_code == 0
    assert not out_dir.exists()
    assert "h2_grad_orca.inp" in result.stdout
    assert "H" in result.stdout


def test_parse_format_transform_write_without_output_dir_suppresses_stdout(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "h2.xyz"
    input_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            str(input_path),
            "--parser-detection",
            "xyz",
            "--n-jobs",
            "1",
            "format-transform",
            "--format",
            "gjf",
            "--write",
            "--graph-policy",
            "prefer",
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.stdout == ""
    assert (tmp_path / "h2.gjf").exists()


def test_parse_summary_terminal_writes_csv(tmp_path: Path) -> None:
    out = tmp_path / "summary.csv"
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp",
            "--parser-detection",
            "orcainp",
            "--n-jobs",
            "1",
            "to-summary-df",
            "--out",
            str(out),
        ],
    )

    assert result.exit_code == 0
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert "General.FrameID" in header
    assert "('General', 'FrameID')" not in header


def test_parse_summary_frame_all_flatten_columns(tmp_path: Path) -> None:
    out = tmp_path / "summary.csv"
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/g16log/1.log",
            "--parser-detection",
            "g16log",
            "--n-jobs",
            "1",
            "to-summary-df",
            "--frame",
            "all",
            "--flatten-columns",
            "--out",
            str(out),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(out.read_text(encoding="utf-8").splitlines()))
    assert len(rows) > 1
    assert [row["General.FrameID"] for row in rows] == ["0", "1", "2", "3", "4"]


def test_parse_summary_g16_frame_json_terminal_handles_unit_columns() -> None:
    result = runner.invoke(
        app,
        [
            "-q",
            "parse",
            "tests/test_files/g16log/000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log",
            "--parser-detection",
            "g16log",
            "--n-jobs",
            "1",
            "to-summary-df",
            "--mode",
            "frame",
            "--frame",
            "-1",
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    rows = json.loads(result.output)
    assert rows[0]["Environment.Temperature.kelvin"] == 298.15
    assert rows[0]["Status.IsError"] is None
    assert rows[0]["General.FrameID"] == 33


def test_old_flat_command_is_removed() -> None:
    result = runner.invoke(
        app, ["summary", "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp"]
    )

    assert result.exit_code != 0
    assert "No such command" in result.output
