from __future__ import annotations

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


def test_parse_operation_is_registered_as_click_command() -> None:
    result = runner.invoke(app, ["parse", "dummy.log", "format-transform", "--help"])

    assert result.exit_code == 0
    assert "--format" in result.stdout
    assert "--output-dir" in result.stdout


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
    assert completions[0].help == "Rendering backend used to generate the target format. Default: 'rdkit'."


def test_format_transform_has_no_orcainp_writer_dynamic_options() -> None:
    command = app.commands["parse"].commands["format-transform"]
    extra_arg = next(param for param in command.params if param.name == "extra_args")
    ctx = click.Context(command)
    ctx.params["target_format"] = "orcainp"
    ctx.args = ["--format", "orcainp"]

    completions = extra_arg.shell_complete(ctx, "--no-use")

    assert completions == []


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
    assert out.read_text(encoding="utf-8").strip()


def test_old_flat_command_is_removed() -> None:
    result = runner.invoke(app, ["summary", "tests/test_files/orca/single_point_inputs/h2_grad_orca.inp"])

    assert result.exit_code != 0
    assert "No such command" in result.output
