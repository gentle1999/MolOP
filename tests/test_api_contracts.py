from __future__ import annotations

from inspect import Parameter, signature

import pytest

from molop.cli.app import FORMAT_TRANSFORM_STATIC_OPTIONS
from molop.cli.app import app as cli_app
from molop.cli.state_machine import (
    OPERATION_REGISTRY,
    BatchInputConfig,
    CliUsageError,
    FormatTransformParams,
    OperationCall,
    SampleParams,
    ToSummaryDfParams,
    build_plan,
)
from molop.io import AutoParser
from molop.io._batch_format_transform import BatchFormatTransformMixin
from molop.io.base_models._format_transform import (
    FormatTransformMixin,
    FrameFormatTransformMixin,
)
from molop.io.base_models.ChemFile import BaseCalcFile, BaseChemFile
from molop.io.base_models.ChemFileFrame import BaseCalcFrame
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.FileBatchParserDisk import FileBatchParserDisk
from molop.io.frame_selection import normalize_frame_selector
from molop.utils.progressbar import parallel_map


def test_format_transform_public_signature_contract() -> None:
    frame_sig = signature(FrameFormatTransformMixin.format_transform)
    assert list(frame_sig.parameters) == [
        "self",
        "format",
        "file_path",
        "write_to_disk",
        "kwargs",
    ]
    assert frame_sig.parameters["file_path"].default is None
    assert frame_sig.parameters["write_to_disk"].default is False

    file_sig = signature(FormatTransformMixin.format_transform)
    assert list(file_sig.parameters) == [
        "self",
        "format",
        "frame",
        "file_path",
        "embed_in_one_file",
        "write_to_disk",
        "kwargs",
    ]
    assert file_sig.parameters["frame"].default == -1
    assert "slice" not in str(file_sig.parameters["frame"].annotation)
    assert file_sig.parameters["embed_in_one_file"].default is True
    assert file_sig.parameters["write_to_disk"].default is False

    batch_sig = signature(BatchFormatTransformMixin.format_transform)
    assert list(batch_sig.parameters) == [
        "self",
        "format",
        "output_dir",
        "frame",
        "embed_in_one_file",
        "write_to_disk",
        "n_jobs",
        "kwargs",
    ]
    assert batch_sig.parameters["output_dir"].default is None
    assert batch_sig.parameters["frame"].default == -1
    assert batch_sig.parameters["embed_in_one_file"].default is True
    assert batch_sig.parameters["write_to_disk"].default is False
    assert batch_sig.parameters["n_jobs"].default == -1


def test_to_summary_df_public_signature_contract() -> None:
    file_sig = signature(BaseChemFile.to_summary_df)
    assert file_sig.parameters["frame"].kind is Parameter.KEYWORD_ONLY
    assert file_sig.parameters["frame"].default == -1
    assert "frames" not in file_sig.parameters
    assert file_sig.parameters["flatten_columns"].default is False
    assert file_sig.parameters["on_missing_frame"].default == "skip"

    sig = signature(FileBatchModelDisk.to_summary_df)

    assert sig.parameters["mode"].default == "frame"
    assert sig.parameters["frame"].default == -1
    assert "frames" not in sig.parameters
    assert sig.parameters["n_jobs"].default == -1
    assert sig.parameters["brief"].kind is Parameter.KEYWORD_ONLY
    assert sig.parameters["brief"].default is True
    assert sig.parameters["flatten_columns"].kind is Parameter.KEYWORD_ONLY
    assert sig.parameters["flatten_columns"].default is False
    assert sig.parameters["on_missing_frame"].kind is Parameter.KEYWORD_ONLY
    assert sig.parameters["on_missing_frame"].default == "skip"


def test_ts_endpoint_sampling_signature_contract() -> None:
    callables = (
        BaseCalcFrame.possible_pre_post_ts,
        BaseCalcFrame.save_pre_post_ts,
        BaseCalcFrame.to_diff_rdmol,
        BaseCalcFile.save_pre_post_ts,
        FileBatchModelDisk.save_pre_post_ts,
    )

    for callable_obj in callables:
        parameters = signature(callable_obj).parameters
        assert parameters["min_ratio"].default == 0.6
        assert parameters["max_ratio"].default == 1.4
        assert parameters["steps"].default == 8
        assert parameters["sampling_method"].default == "harmonic_potential"
        assert "ratio" not in parameters
        assert "ratio_attempts" not in parameters


def test_ts_endpoint_additional_sampling_signature_contract() -> None:
    parameters = signature(BaseCalcFrame.additional_pre_post_ts).parameters

    assert parameters["pre_rdmol"].kind is Parameter.POSITIONAL_OR_KEYWORD
    assert parameters["post_rdmol"].kind is Parameter.POSITIONAL_OR_KEYWORD
    assert parameters["min_ratio"].default == 0.6
    assert parameters["max_ratio"].default == 1.4
    assert parameters["steps"].default == 8
    assert parameters["sampling_method"].default == "harmonic_potential"


def test_parallel_execute_public_signature_contract() -> None:
    sig = signature(FileBatchModelDisk.parallel_execute)

    assert sig.parameters["desc"].default == ""
    assert sig.parameters["n_jobs"].default == -1
    assert sig.parameters["return_as"].default == "list"
    assert sig.parameters["return_results"].default is None
    assert sig.parameters["_diskfiles_snapshot"].kind is Parameter.KEYWORD_ONLY
    assert sig.parameters["_diskfiles_snapshot"].default is None


def test_all_public_n_jobs_defaults_are_minus_one() -> None:
    callables = (
        AutoParser,
        FileBatchParserDisk,
        parallel_map,
        FileBatchModelDisk.parallel_execute,
        FileBatchModelDisk.filter_state,
        FileBatchModelDisk.filter_value,
        FileBatchModelDisk.filter_custom,
        FileBatchModelDisk.groupby,
        FileBatchModelDisk.filter_by_codec_id,
        FileBatchModelDisk.to_summary_df,
        FileBatchModelDisk.draw_grid_image,
        FileBatchModelDisk.copy_to,
        FileBatchModelDisk.move_to,
        BatchFormatTransformMixin.format_transform,
    )

    for callable_obj in callables:
        assert signature(callable_obj).parameters["n_jobs"].default == -1

    parse_command = cli_app.commands["parse"]
    commands = (parse_command, *parse_command.commands.values())
    n_jobs_options = [
        parameter
        for command in commands
        for parameter in command.params
        if parameter.name == "n_jobs"
    ]
    assert n_jobs_options
    assert all(option.default == -1 for option in n_jobs_options)

    operation_models = {
        operation.params_model
        for operation in OPERATION_REGISTRY.values()
        if "n_jobs" in operation.params_model.model_fields
    }
    assert operation_models
    assert all(model.model_fields["n_jobs"].default == -1 for model in operation_models)


def test_frame_selector_contract() -> None:
    assert normalize_frame_selector("all", 3) == [0, 1, 2]
    assert normalize_frame_selector(1, 3) == [1]
    assert normalize_frame_selector(-1, 3) == [2]
    assert normalize_frame_selector([0, -1], 3) == [0, 2]

    with pytest.raises(ValueError, match="frame must be an integer"):
        normalize_frame_selector("last", 3)
    with pytest.raises(TypeError, match="sequence must contain only integers"):
        normalize_frame_selector([0, "bad"], 3)  # type: ignore[list-item]
    with pytest.raises(IndexError, match="Frame index 3 is out of range"):
        normalize_frame_selector(3, 3, validate_range=True)


def test_cli_chain_contract_metadata_and_terminal_validation() -> None:
    assert {"--write", "--no-write"}.issubset(FORMAT_TRANSFORM_STATIC_OPTIONS)
    assert OPERATION_REGISTRY["format-transform"].return_kind == "terminal"
    assert OPERATION_REGISTRY["format-transform"].effect == "writes-files"
    assert OPERATION_REGISTRY["to-summary-df"].return_kind == "terminal"

    input_config = BatchInputConfig(pattern="dummy.log")
    terminal = OperationCall(
        spec=OPERATION_REGISTRY["to-summary-df"],
        params=ToSummaryDfParams(),
    )
    chained_after_terminal = OperationCall(
        spec=OPERATION_REGISTRY["sample"],
        params=SampleParams(),
    )
    with pytest.raises(CliUsageError, match="must be the last operation"):
        build_plan(input_config, [terminal, chained_after_terminal])


def test_cli_format_transform_write_contract(tmp_path) -> None:
    input_config = BatchInputConfig(pattern="dummy.xyz", n_jobs=4)
    out_dir = tmp_path / "out"

    render_only = FormatTransformParams(
        format="xyz",
        output_dir=out_dir,
        write_to_disk=False,
    )
    assert render_only.method_kwargs(input_config)["write_to_disk"] is False
    assert not out_dir.exists()

    implicit_write = FormatTransformParams(format="xyz", output_dir=out_dir)
    assert implicit_write.method_kwargs(input_config)["write_to_disk"] is True
    assert out_dir.is_dir()

    source_dir_write = FormatTransformParams(format="xyz", write_to_disk=True)
    kwargs = source_dir_write.method_kwargs(input_config)
    assert kwargs["write_to_disk"] is True
    assert kwargs["output_dir"] is None
