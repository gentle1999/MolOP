from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Literal, Protocol, get_origin, get_type_hints, runtime_checkable

import click
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, field_validator

from molop.cli.shared.frames import parse_frame_selection
from molop.io import AutoParser
from molop.io.codec_types import ParseOptions
from molop.io.FileBatchModelDisk import FileBatchModelDisk
from molop.io.frame_selection import FrameSelector
from molop.io.parse_outcomes import BatchParseResult, FileParseOutcome


ReturnKind = Literal["batch", "terminal"]
EffectKind = Literal["pure", "writes-files", "moves-files", "releases-memory"]
OutputFormat = Literal["text", "json"]


@runtime_checkable
class _ImageWithData(Protocol):
    data: str | bytes


@runtime_checkable
class _ImageWithSave(Protocol):
    def save(self, path: Path) -> object: ...


class CliUsageError(Exception):
    pass


class CliRuntimeError(Exception):
    pass


class CliModel(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


class BatchInputConfig(CliModel):
    pattern: str
    additional_patterns: tuple[str, ...] = ()
    parser_detection: str = "auto"
    n_jobs: int = -1
    output_format: OutputFormat = "text"
    total_charge: int | None = None
    total_multiplicity: int | None = None
    only_extract_structure: bool = False
    only_last_frame: bool = False
    capture_source_evidence: bool = False
    source_encoding: str = "utf-8"
    release_file_content: bool = True
    force_unit_transform: bool | None = None
    graph_reconstruction_backend: Literal["cpp", "python"] | None = None
    make_dative_bonds: bool | None = None
    make_stereochemistry: bool | None = None
    report: bool = False

    def parse_options(self) -> ParseOptions:
        return ParseOptions(
            total_charge=self.total_charge,
            total_multiplicity=self.total_multiplicity,
            only_extract_structure=self.only_extract_structure,
            only_last_frame=self.only_last_frame,
            capture_source_evidence=self.capture_source_evidence,
            source_encoding=self.source_encoding,
            release_file_content=self.release_file_content,
            force_unit_transform=self.force_unit_transform,
            graph_reconstruction_backend=self.graph_reconstruction_backend,
            make_dative_bonds=self.make_dative_bonds,
            make_stereochemistry=self.make_stereochemistry,
        )

    def parser_input(self) -> str | list[str]:
        if not self.additional_patterns:
            return self.pattern
        return [self.pattern, *self.additional_patterns]


class OperationParams(CliModel):
    def method_kwargs(self, input_config: BatchInputConfig) -> dict[str, Any]:
        return self.model_dump(exclude_none=True)


class FilterStateParams(OperationParams):
    state: Literal["ts", "error", "opt", "normal", "thermal", "no-img"]
    negate: bool = False
    n_jobs: int = -1


class FilterValueParams(OperationParams):
    target: Literal["charge", "multiplicity", "format"]
    value: str | float | int
    compare: Literal["==", "!=", ">", "<", ">=", "<="] = "=="
    n_jobs: int = -1

    @field_validator("value", mode="before")
    @classmethod
    def _coerce_value(cls, value: object) -> object:
        if not isinstance(value, str):
            return value
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            return value


class FilterByCodecParams(OperationParams):
    codec_id: str
    negate: bool = False
    n_jobs: int = -1
    on_missing: Literal["keep", "drop", "error"] = "drop"


class SampleParams(OperationParams):
    n: int = 10
    seed: int | None = None


class FormatTransformParams(OperationParams):
    format: str
    output_dir: Path | None = None
    frame: str = "-1"
    embed: bool = True
    write_to_disk: bool | None = None
    n_jobs: int = -1
    format_options: dict[str, Any] = Field(default_factory=dict)

    def method_kwargs(self, input_config: BatchInputConfig) -> dict[str, Any]:
        should_write = (
            self.write_to_disk if self.write_to_disk is not None else self.output_dir is not None
        )
        kwargs: dict[str, Any] = {
            "format": self.format,
            "output_dir": str(self.output_dir) if self.output_dir is not None else None,
            "frame": parse_cli_frame_selection(self.frame),
            "embed_in_one_file": self.embed,
            "write_to_disk": should_write,
            "n_jobs": self.n_jobs,
            **self.format_options,
        }
        if should_write and self.output_dir is not None:
            self.output_dir.mkdir(parents=True, exist_ok=True)
        return kwargs


class ToSummaryDfParams(OperationParams):
    mode: Literal["file", "frame"] = "frame"
    frame: str = "-1"
    n_jobs: int = -1
    out: Path | None = None
    format: Literal["csv", "json"] = "csv"
    brief: bool = True
    flatten_columns: bool = True
    on_missing_frame: Literal["skip", "error"] = "skip"

    def method_kwargs(self, input_config: BatchInputConfig) -> dict[str, Any]:
        frame_selection = parse_cli_frame_selection(self.frame)
        return {
            "mode": self.mode,
            "frame": frame_selection,
            "n_jobs": self.n_jobs,
            "brief": self.brief,
            "flatten_columns": self.flatten_columns,
            "on_missing_frame": self.on_missing_frame,
        }


class DrawGridImageParams(OperationParams):
    out: Path
    mols_per_row: int = 4
    sub_img_width: int = 200
    sub_img_height: int = 200
    max_mols: int = 16
    use_svg: bool | None = None
    n_jobs: int = -1

    def method_kwargs(self, input_config: BatchInputConfig) -> dict[str, Any]:
        use_svg = self.out.suffix.lower() == ".svg" if self.use_svg is None else self.use_svg
        return {
            "molsPerRow": self.mols_per_row,
            "subImgSize": (self.sub_img_width, self.sub_img_height),
            "maxMols": self.max_mols,
            "useSVG": use_svg,
            "n_jobs": self.n_jobs,
        }


class GroupByParams(OperationParams):
    key: Literal["detected_format_id", "file_format", "state"] = "detected_format_id"
    n_jobs: int = -1


class CopyToParams(OperationParams):
    output_dir: Path
    n_jobs: int = -1

    def method_kwargs(self, input_config: BatchInputConfig) -> dict[str, Any]:
        return {
            "output_dir": str(self.output_dir),
            "n_jobs": self.n_jobs,
        }


class MoveToParams(CopyToParams):
    pass


Handler = Callable[[FileBatchModelDisk[Any], OperationParams, BatchInputConfig], object]


class OperationSpec(CliModel):
    cli_name: str
    method_name: str
    params_model: type[OperationParams]
    return_kind: ReturnKind
    effect: EffectKind = "pure"
    handler: Handler | None = Field(default=None, exclude=True)


class OperationCall(CliModel):
    spec: OperationSpec
    params: OperationParams


class BatchPlan(CliModel):
    input: BatchInputConfig
    operations: list[OperationCall]


def parse_cli_frame_selection(frame: str) -> FrameSelector:
    try:
        return parse_frame_selection(frame)
    except ValueError as exc:
        raise CliUsageError(str(exc)) from exc


def infer_return_kind(method: Callable[..., Any]) -> ReturnKind:
    hints = get_type_hints(method)
    return_type = hints.get("return")
    if return_type is None:
        raise CliUsageError(f"Cannot expose {method.__name__}: missing return annotation.")
    if _is_file_batch_return(return_type):
        return "batch"
    return "terminal"


def validate_plan(plan: BatchPlan) -> None:
    if plan.input.report and plan.operations:
        raise CliUsageError("--report cannot be combined with operation commands.")
    for index, operation in enumerate(plan.operations):
        is_last = index == len(plan.operations) - 1
        if operation.spec.return_kind == "terminal" and not is_last:
            raise CliUsageError(
                f"{operation.spec.cli_name} returns non-FileBatchModelDisk "
                "and must be the last operation."
            )


def build_plan(input_config: BatchInputConfig, operations: Sequence[OperationCall]) -> BatchPlan:
    plan = BatchPlan(input=input_config, operations=list(operations))
    validate_plan(plan)
    return plan


def execute_plan(plan: BatchPlan) -> object:
    state: object = AutoParser(
        plan.input.parser_input(),
        n_jobs=plan.input.n_jobs,
        parser_detection=plan.input.parser_detection,
        parse_options=plan.input.parse_options(),
        return_report=plan.input.report,
    )
    for operation in plan.operations:
        if not isinstance(state, FileBatchModelDisk):
            raise CliRuntimeError("Invalid state transition: current state is not a batch.")
        state = _execute_operation(state, operation, plan.input)
    return state


def render_result(result: object, plan: BatchPlan) -> None:
    if _result_is_file_output_only(plan):
        return

    if isinstance(result, BatchParseResult):
        _render_parse_report(result, output_format=plan.input.output_format)
        return

    if isinstance(result, FileBatchModelDisk):
        _render_paths(result.file_paths, output_format=plan.input.output_format)
        return

    last = plan.operations[-1] if plan.operations else None
    params = last.params if last is not None else None

    if isinstance(result, pd.DataFrame):
        if not isinstance(params, ToSummaryDfParams):
            _render_json_or_text(result.to_dict(orient="records"), plan.input.output_format)
            return
        _render_summary_df(result, params, plan.input.output_format)
        return

    if isinstance(result, dict):
        if result and all(isinstance(value, FileBatchModelDisk) for value in result.values()):
            groups = {
                str(key): value.file_paths
                for key, value in result.items()
                if isinstance(value, FileBatchModelDisk)
            }
            _render_json_or_text(groups, plan.input.output_format)
            return
        _render_json_or_text(result, plan.input.output_format)
        return

    if isinstance(params, DrawGridImageParams):
        _write_grid_image(result, params.out)
        click.echo(f"Image written to {params.out}", err=True)
        return

    if result is None:
        return

    if isinstance(result, list):
        _render_json_or_text(result, plan.input.output_format)
        return

    click.echo(str(result))


def _result_is_file_output_only(plan: BatchPlan) -> bool:
    if not plan.operations:
        return False
    params = plan.operations[-1].params
    return isinstance(params, FormatTransformParams) and _format_transform_writes_to_disk(params)


def _format_transform_writes_to_disk(params: FormatTransformParams) -> bool:
    if params.write_to_disk is not None:
        return params.write_to_disk
    return params.output_dir is not None


def _execute_operation(
    batch: FileBatchModelDisk[Any],
    operation: OperationCall,
    input_config: BatchInputConfig,
) -> object:
    if operation.spec.handler is not None:
        return operation.spec.handler(batch, operation.params, input_config)
    method = getattr(batch, operation.spec.method_name)
    return method(**operation.params.method_kwargs(input_config))


def parse_dynamic_options(tokens: Sequence[str]) -> dict[str, Any]:
    options: dict[str, Any] = {}
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if not token.startswith("--"):
            raise CliUsageError(
                f"Unexpected dynamic format argument: {token}. Use --option value or --flag syntax."
            )

        key = token[2:].replace("-", "_")
        if not key:
            raise CliUsageError("Dynamic format option name cannot be empty.")
        index += 1

        if key.startswith("no_"):
            options[key[3:]] = False
            continue

        if index >= len(tokens) or tokens[index].startswith("--"):
            options[key] = True
            continue

        options[key] = _coerce_dynamic_value(tokens[index])
        index += 1

    return options


def _coerce_dynamic_value(value: str) -> object:
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "none" or lowered == "null":
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _is_file_batch_return(return_type: object) -> bool:
    origin = get_origin(return_type)
    if return_type is FileBatchModelDisk or origin is FileBatchModelDisk:
        return True
    return isinstance(return_type, str) and return_type.startswith("FileBatchModelDisk")


def _render_paths(paths: list[str], *, output_format: OutputFormat) -> None:
    if output_format == "json":
        click.echo(json.dumps(paths, indent=2))
        return
    for path in paths:
        click.echo(path)


def _render_json_or_text(value: object, output_format: OutputFormat) -> None:
    if output_format == "json":
        click.echo(json.dumps(value, indent=2, sort_keys=True, default=str))
        return
    if isinstance(value, dict):
        for key, item in value.items():
            click.echo(f"{key}\t{item}")
        return
    if isinstance(value, list):
        for item in value:
            click.echo(item)
        return
    click.echo(str(value))


def _render_parse_report(
    report: BatchParseResult[object, object], *, output_format: OutputFormat
) -> None:
    outcomes = [_parse_outcome_payload(outcome) for outcome in report.outcomes]
    summary = {
        "total": len(outcomes),
        "succeeded": len(report.succeeded),
        "failed": len(report.failures),
    }
    payload = {
        "summary": summary,
        "outcomes": outcomes,
    }
    if output_format == "json":
        click.echo(json.dumps(payload, indent=2, sort_keys=True))
        return

    click.echo(
        f"total={summary['total']}\tsucceeded={summary['succeeded']}\tfailed={summary['failed']}"
    )
    click.echo("status\tformat\tpath\tmessage")
    for outcome in outcomes:
        failure = outcome["failure"]
        message = failure["message"] if isinstance(failure, dict) else ""
        click.echo(
            f"{outcome['status']}\t{outcome['detected_format'] or ''}\t"
            f"{outcome['file_path']}\t{message}"
        )


def _parse_outcome_payload(outcome: FileParseOutcome[object]) -> dict[str, object]:
    failure = outcome.failure
    return {
        "input_index": outcome.input_index,
        "file_path": outcome.file_path,
        "status": outcome.status,
        "detected_format": outcome.detected_format,
        "warnings": [
            {"code": warning.code, "message": warning.message} for warning in outcome.warnings
        ],
        "failure": (
            {
                "kind": failure.kind,
                "message": failure.message,
                "reader_format": failure.reader_format,
                "exception_type": failure.exception_type,
            }
            if failure is not None
            else None
        ),
    }


def _render_summary_df(
    df: pd.DataFrame, params: ToSummaryDfParams, output_format: OutputFormat
) -> None:
    if params.out is not None:
        params.out.parent.mkdir(parents=True, exist_ok=True)
        if params.format == "json":
            df.to_json(params.out, orient="records", force_ascii=True, indent=2)
        else:
            df.to_csv(params.out, index=False)
        click.echo(f"Summary written to {params.out}", err=True)
        return
    if output_format == "json" or params.format == "json":
        click.echo(df.to_json(orient="records", force_ascii=True, indent=2))
    else:
        click.echo(df.to_csv(index=False))


def _write_grid_image(image: object, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".svg":
        content = image.data if isinstance(image, _ImageWithData) else image
        with open(out, "w", encoding="utf-8") as handle:
            handle.write(str(content))
        return
    if isinstance(image, _ImageWithSave):
        image.save(out)
        return
    if isinstance(image, _ImageWithData):
        data = image.data.encode("utf-8") if isinstance(image.data, str) else image.data
        with open(out, "wb") as handle:
            handle.write(data)
        return
    raise CliRuntimeError("Unsupported image object returned by RDKit.")


def _group_key(diskfile: Any, key: str) -> str:
    if key == "detected_format_id":
        fmt = getattr(diskfile, "detected_format_id", "unknown")
        return fmt if fmt is not None else "unknown"
    if key == "file_format":
        return diskfile.file_format
    if key == "state":
        if not diskfile:
            return "unknown"
        last_frame = diskfile[-1]
        if getattr(last_frame, "is_TS", False):
            return "ts"
        if getattr(last_frame, "is_error", False):
            return "error"
        if getattr(last_frame, "is_normal", False):
            return "normal"
        return "unknown"
    return "unknown"


def _groupby_handler(
    batch: FileBatchModelDisk[Any],
    params: OperationParams,
    input_config: BatchInputConfig,
) -> object:
    typed = _ensure_params(params, GroupByParams)
    return batch.groupby(
        key_func=lambda diskfile: _group_key(diskfile, typed.key), n_jobs=typed.n_jobs
    )


def _ensure_params(params: OperationParams, expected_type: type[Any]) -> Any:
    if not isinstance(params, expected_type):
        raise CliRuntimeError(f"Expected {expected_type.__name__}, got {type(params).__name__}")
    return params


def _method(name: str) -> Callable[..., Any]:
    method = getattr(FileBatchModelDisk, name)
    if not callable(method):
        raise TypeError(f"{name} is not callable")
    return method


def _spec(
    cli_name: str,
    method_name: str,
    params_model: type[OperationParams],
    *,
    effect: EffectKind = "pure",
    return_kind: ReturnKind | None = None,
    handler: Handler | None = None,
) -> OperationSpec:
    if return_kind is None:
        return_kind = infer_return_kind(_method(method_name))
    return OperationSpec(
        cli_name=cli_name,
        method_name=method_name,
        params_model=params_model,
        return_kind=return_kind,
        effect=effect,
        handler=handler,
    )


OPERATION_REGISTRY: dict[str, OperationSpec] = {
    spec.cli_name: spec
    for spec in (
        _spec("filter-state", "filter_state", FilterStateParams),
        _spec("filter-value", "filter_value", FilterValueParams),
        _spec("filter-by-codec", "filter_by_codec_id", FilterByCodecParams),
        _spec("sample", "sample", SampleParams),
        _spec(
            "format-transform",
            "format_transform",
            FormatTransformParams,
            effect="writes-files",
        ),
        _spec("to-summary-df", "to_summary_df", ToSummaryDfParams),
        _spec(
            "draw-grid-image",
            "draw_grid_image",
            DrawGridImageParams,
            effect="writes-files",
            return_kind="terminal",
        ),
        _spec(
            "groupby",
            "groupby",
            GroupByParams,
            return_kind="terminal",
            handler=_groupby_handler,
        ),
        _spec(
            "copy-to",
            "copy_to",
            CopyToParams,
            effect="writes-files",
            return_kind="terminal",
        ),
        _spec(
            "move-to",
            "move_to",
            MoveToParams,
            effect="moves-files",
            return_kind="terminal",
        ),
    )
}


def registered_operations_help() -> str:
    return "\n".join(f"  {name}" for name in sorted(OPERATION_REGISTRY))
