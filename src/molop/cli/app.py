from __future__ import annotations

import os
import shlex
from pathlib import Path
from typing import Literal, cast

import click
from click.shell_completion import CompletionItem

from molop.cli.state_machine import (
    OPERATION_REGISTRY,
    BatchInputConfig,
    CliRuntimeError,
    CliUsageError,
    CopyToParams,
    DrawGridImageParams,
    FilterByCodecParams,
    FilterStateParams,
    FilterValueParams,
    FormatTransformParams,
    GroupByParams,
    MoveToParams,
    OperationCall,
    OperationParams,
    SampleParams,
    ToSummaryDfParams,
    build_plan,
    execute_plan,
    parse_dynamic_options,
    render_result,
)
from molop.io.codec_registry import get_supported_writer_formats, get_writer_option_specs


def _version() -> str:
    import importlib.metadata

    try:
        return importlib.metadata.version("molop")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}
PARSE_VALUE_OPTIONS = {"--parser-detection", "--n-jobs", "-j", "--output-format"}
PARSE_FLAG_OPTIONS = {"--help", "-h"}
FORMAT_TRANSFORM_STATIC_OPTIONS = {
    "--format",
    "--output-dir",
    "--frame",
    "--embed",
    "--no-embed",
    "--n-jobs",
    "-h",
    "--help",
}
FORMAT_TRANSFORM_VALUE_OPTIONS = {"--format", "--output-dir", "--frame", "--n-jobs"}
COMPLETION_SHELLS = ("bash", "zsh", "fish")
COMPLETION_MARKER_START = "# >>> molop completion >>>"
COMPLETION_MARKER_END = "# <<< molop completion <<<"
FORMAT_VALUE_COMPLETION_HELP = {
    "graph_policy": {
        "prefer": "Use molecular graph data when available; fall back to coordinates.",
        "strict": "Require molecular graph data before rendering.",
        "coords": "Render from coordinate data only.",
    },
    "engine": {
        "rdkit": "Use RDKit for rendering.",
        "openbabel": "Use Open Babel for rendering.",
    },
    "coords_type": {
        "auto": "Preserve the stored Gaussian coordinate representation.",
        "cartesian": "Render Gaussian coordinates in Cartesian form.",
        "internal": "Render Gaussian coordinates in internal-coordinate form.",
    },
}


class DynamicFormatOptionsArgument(click.Argument):
    def shell_complete(
        self,
        ctx: click.Context,
        incomplete: str,
    ) -> list[CompletionItem]:
        target_format = _current_format_transform_format(ctx)
        if not target_format:
            return []

        option_specs = get_writer_option_specs(target_format)
        tokens = _format_transform_tokens(ctx)
        if _completing_dynamic_option_value(ctx):
            option_name = _last_dynamic_option_name(tokens)
            if option_name is None:
                return []
            return [
                CompletionItem(value, help=_dynamic_value_help(spec.name, value))
                for spec in option_specs
                if spec.option == option_name
                for value in spec.value_candidates
                if value.startswith(incomplete)
            ]

        used_options = _used_dynamic_options(tokens)
        completions: list[CompletionItem] = []
        for spec in option_specs:
            if spec.option not in used_options and spec.option.startswith(incomplete):
                completions.append(CompletionItem(spec.option, help=spec.help))
            no_option = _no_option_name(spec.option)
            if (
                spec.supports_no
                and no_option not in used_options
                and no_option.startswith(incomplete)
            ):
                completions.append(CompletionItem(no_option, help=spec.help))
        return completions


class FormatTransformCommand(click.Command):
    def shell_complete(
        self,
        ctx: click.Context,
        incomplete: str,
    ) -> list[CompletionItem]:
        completions = super().shell_complete(ctx, incomplete)
        if incomplete and not incomplete[0].isalnum():
            dynamic_arg = next(
                (
                    param
                    for param in self.get_params(ctx)
                    if isinstance(param, DynamicFormatOptionsArgument)
                ),
                None,
            )
            if dynamic_arg is not None:
                completions.extend(dynamic_arg.shell_complete(ctx, incomplete))
        return _dedupe_completions(completions)


class ParseChainGroup(click.Group):
    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        return super().parse_args(ctx, self._normalize_parse_options(args))

    def _normalize_parse_options(self, args: list[str]) -> list[str]:
        prefix: list[str] = []
        index = 0
        while index < len(args):
            consumed = _consume_parse_option(args, index)
            if consumed is None:
                break
            values, index = consumed
            prefix.extend(values)

        if index >= len(args) or args[index].startswith("-"):
            return args

        pattern = args[index]
        index += 1
        while index < len(args) and args[index] not in self.commands:
            consumed = _consume_parse_option(args, index)
            if consumed is None:
                break
            values, index = consumed
            prefix.extend(values)

        return [*prefix, pattern, *args[index:]]


def _consume_parse_option(args: list[str], index: int) -> tuple[list[str], int] | None:
    token = args[index]
    if token in PARSE_FLAG_OPTIONS:
        return [token], index + 1
    if token in PARSE_VALUE_OPTIONS:
        if index + 1 >= len(args):
            return [token], index + 1
        return [token, args[index + 1]], index + 2
    if any(
        token.startswith(f"{option}=") for option in PARSE_VALUE_OPTIONS if option.startswith("--")
    ):
        return [token], index + 1
    return None


def _writer_format_complete(
    ctx: click.Context,
    param: click.Parameter,
    incomplete: str,
) -> list[CompletionItem]:
    _ = (ctx, param)
    return [
        CompletionItem(format_id)
        for format_id in get_supported_writer_formats()
        if format_id.startswith(incomplete)
    ]


def _dedupe_completions(completions: list[CompletionItem]) -> list[CompletionItem]:
    seen: set[str] = set()
    deduped: list[CompletionItem] = []
    for item in completions:
        if item.value in seen:
            continue
        seen.add(item.value)
        deduped.append(item)
    return deduped


def _current_format_transform_format(ctx: click.Context) -> str | None:
    value = ctx.params.get("target_format")
    if isinstance(value, str) and value:
        return value

    args = _format_transform_tokens(ctx)
    for index, token in enumerate(args):
        if token == "--format" and index + 1 < len(args):
            return args[index + 1]
        if token.startswith("--format="):
            return token.split("=", 1)[1]
    return None


def _format_transform_tokens(ctx: click.Context) -> list[str]:
    extra_args = ctx.params.get("extra_args")
    if isinstance(extra_args, tuple):
        return list(extra_args)
    if isinstance(extra_args, list):
        return list(extra_args)
    return list(ctx.args)


def _completing_dynamic_option_value(ctx: click.Context) -> bool:
    args = _format_transform_tokens(ctx)
    option_name = _last_dynamic_option_name(args)
    if option_name is None:
        return False
    if not args or args[-1] != option_name:
        return False
    return option_name not in _boolean_dynamic_options(_current_format_from_args(args))


def _last_dynamic_option_name(args: list[str]) -> str | None:
    if not args:
        return None
    token = args[-1]
    if token.startswith("--") and token not in FORMAT_TRANSFORM_STATIC_OPTIONS and "=" not in token:
        return token
    return None


def _used_dynamic_options(args: list[str]) -> set[str]:
    used: set[str] = set()
    for token in args:
        if (
            token.startswith("--")
            and token not in FORMAT_TRANSFORM_STATIC_OPTIONS
            and "=" not in token
        ):
            used.add(token)
    return used


def _current_format_from_args(args: list[str]) -> str | None:
    for index, token in enumerate(args):
        if token == "--format" and index + 1 < len(args):
            return args[index + 1]
        if token.startswith("--format="):
            return token.split("=", 1)[1]
    return None


def _boolean_dynamic_options(format_id: str | None) -> set[str]:
    if not format_id:
        return set()
    return {spec.option for spec in get_writer_option_specs(format_id) if not spec.takes_value}


def _no_option_name(option: str) -> str:
    return f"--no-{option[2:]}" if option.startswith("--") else f"no-{option}"


def _dynamic_value_help(option_name: str, value: str) -> str | None:
    return FORMAT_VALUE_COMPLETION_HELP.get(option_name, {}).get(value)


def _resolve_completion_shell(shell: str) -> str:
    if shell == "auto":
        shell = Path(os.environ.get("SHELL", "")).name
    if shell not in COMPLETION_SHELLS:
        raise click.ClickException(f"Unsupported shell: {shell}")
    return shell


def _completion_source(shell: str) -> str:
    from click.shell_completion import get_completion_class

    comp_cls = get_completion_class(shell)
    if comp_cls is None:
        raise click.ClickException(f"Unsupported shell: {shell}")
    return comp_cls(app, {}, "molop", "_MOLOP_COMPLETE").source()


def _completion_source_path(shell: str) -> Path:
    return Path.home() / ".config" / "molop" / "completions" / shell / f"molop.{shell}"


def _completion_rc_path(shell: str) -> Path:
    if shell == "bash":
        return Path.home() / ".bashrc"
    if shell == "zsh":
        return Path.home() / ".zshrc"
    raise click.ClickException(f"No rc file is managed for shell: {shell}")


def _install_completion_script(shell: str) -> Path:
    source = _completion_source(shell)
    if shell == "fish":
        target = Path.home() / ".config" / "fish" / "completions" / "molop.fish"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(source, encoding="utf-8")
        return target

    target = _completion_source_path(shell)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(source, encoding="utf-8")

    rc_path = _completion_rc_path(shell)
    rc_path.parent.mkdir(parents=True, exist_ok=True)
    source_line = f"source {shlex.quote(str(target))}"
    block = f"{COMPLETION_MARKER_START}\n{source_line}\n{COMPLETION_MARKER_END}\n"
    existing = rc_path.read_text(encoding="utf-8") if rc_path.exists() else ""
    updated = _replace_completion_block(existing, block)
    rc_path.write_text(updated, encoding="utf-8")
    return target


def _replace_completion_block(existing: str, block: str) -> str:
    start = existing.find(COMPLETION_MARKER_START)
    end = existing.find(COMPLETION_MARKER_END)
    if start != -1 and end != -1 and end >= start:
        end += len(COMPLETION_MARKER_END)
        if end < len(existing) and existing[end] == "\n":
            end += 1
        prefix = existing[:start].rstrip("\n")
        suffix = existing[end:].lstrip("\n")
        pieces = [piece for piece in (prefix, block.rstrip("\n"), suffix) if piece]
        return "\n".join(pieces) + "\n"

    if existing and not existing.endswith("\n"):
        existing += "\n"
    return existing + block


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose output.")
@click.option("--quiet", "-q", is_flag=True, help="Enable quiet mode.")
@click.version_option(version=_version(), prog_name="molop")
def _app_callback(verbose: bool = False, quiet: bool = False) -> None:
    """MolOP: Molecule OPerator CLI."""
    if verbose and quiet:
        verbose = False

    if verbose or quiet:
        from molop.config import molopconfig

        if quiet:
            molopconfig.quiet()
        if verbose:
            molopconfig.verbose()


app: click.Group = cast(click.Group, _app_callback)


@app.group("completion", context_settings=CONTEXT_SETTINGS)
def _completion_callback() -> None:
    """Manage shell completion."""


completion: click.Group = cast(click.Group, _completion_callback)


@completion.command("show")
@click.option(
    "--shell",
    type=click.Choice(list(COMPLETION_SHELLS) + ["auto"]),
    default="auto",
    show_default=True,
    help="Shell completion script to print.",
)
def completion_show(shell: Literal["auto", "bash", "zsh", "fish"] = "auto") -> None:
    shell_name = _resolve_completion_shell(shell)
    click.echo(_completion_source(shell_name))


@completion.command("install")
@click.option(
    "--shell",
    type=click.Choice(list(COMPLETION_SHELLS) + ["auto"]),
    default="auto",
    show_default=True,
    help="Shell to install completion for.",
)
def completion_install(shell: Literal["auto", "bash", "zsh", "fish"] = "auto") -> None:
    shell_name = _resolve_completion_shell(shell)
    target = _install_completion_script(shell_name)
    click.echo(str(target), err=True)


@click.group(
    name="parse",
    cls=ParseChainGroup,
    chain=True,
    invoke_without_command=True,
    no_args_is_help=True,
    context_settings=CONTEXT_SETTINGS,
)
@click.argument("pattern")
@click.option(
    "--parser-detection", default="auto", show_default=True, help="Parser detection mode."
)
@click.option("--n-jobs", "-j", default=-1, show_default=True, help="Number of parallel jobs.")
@click.option(
    "--output-format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="Default terminal output format.",
)
def _parse_callback(
    pattern: str,
    parser_detection: str = "auto",
    n_jobs: int = -1,
    output_format: Literal["text", "json"] = "text",
) -> None:
    """Parse files into a FileBatchModelDisk state, then run operation commands."""


parse: ParseChainGroup = cast(ParseChainGroup, _parse_callback)


@parse.result_callback()
def execute_parse_chain(
    operations: list[OperationCall],
    pattern: str,
    parser_detection: str = "auto",
    n_jobs: int = -1,
    output_format: Literal["text", "json"] = "text",
) -> None:
    try:
        plan = build_plan(
            BatchInputConfig(
                pattern=pattern,
                parser_detection=parser_detection,
                n_jobs=n_jobs,
                output_format=output_format,
            ),
            operations,
        )
        render_result(execute_plan(plan), plan)
    except (CliUsageError, CliRuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc


def _operation_call(name: str, params: object) -> OperationCall:
    return OperationCall(spec=OPERATION_REGISTRY[name], params=cast("OperationParams", params))


@parse.command("filter-state")
@click.option(
    "--state",
    required=True,
    type=click.Choice(["ts", "error", "opt", "normal", "thermal", "no-img"]),
    help="Calculation state to keep.",
)
@click.option("--negate", is_flag=True, help="Invert the filter.")
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def filter_state(
    state: Literal["ts", "error", "opt", "normal", "thermal", "no-img"],
    negate: bool = False,
    n_jobs: int | None = None,
) -> OperationCall:
    """Filter the current batch by calculation state."""
    return _operation_call(
        "filter-state", FilterStateParams(state=state, negate=negate, n_jobs=n_jobs)
    )


@parse.command("filter-value")
@click.option(
    "--target",
    required=True,
    type=click.Choice(["charge", "multiplicity", "format"]),
    help="File value target.",
)
@click.option("--value", required=True, help="Value to compare.")
@click.option(
    "--compare",
    type=click.Choice(["==", "!=", ">", "<", ">=", "<="]),
    default="==",
    show_default=True,
    help="Comparison operator.",
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def filter_value(
    target: Literal["charge", "multiplicity", "format"],
    value: str,
    compare: Literal["==", "!=", ">", "<", ">=", "<="] = "==",
    n_jobs: int | None = None,
) -> OperationCall:
    """Filter the current batch by charge, multiplicity, or file format."""
    return _operation_call(
        "filter-value",
        FilterValueParams(target=target, value=value, compare=compare, n_jobs=n_jobs),
    )


@parse.command("filter-by-codec")
@click.option("--codec-id", required=True, help="Detected reader codec id to keep.")
@click.option("--negate", is_flag=True, help="Invert the filter.")
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
@click.option(
    "--on-missing",
    type=click.Choice(["keep", "drop", "error"]),
    default="drop",
    show_default=True,
    help="Behavior when detected codec id is missing.",
)
def filter_by_codec(
    codec_id: str,
    negate: bool = False,
    n_jobs: int | None = None,
    on_missing: Literal["keep", "drop", "error"] = "drop",
) -> OperationCall:
    """Filter the current batch by detected reader codec id."""
    return _operation_call(
        "filter-by-codec",
        FilterByCodecParams(
            codec_id=codec_id,
            negate=negate,
            n_jobs=n_jobs,
            on_missing=on_missing,
        ),
    )


@parse.command("sample")
@click.option("--n", default=10, show_default=True, help="Number of files to sample.")
@click.option("--seed", type=int, default=None, help="Random seed.")
def sample(n: int = 10, seed: int | None = None) -> OperationCall:
    """Randomly sample files from the current batch."""
    return _operation_call("sample", SampleParams(n=n, seed=seed))


@parse.command(
    "format-transform",
    cls=FormatTransformCommand,
    context_settings={
        **CONTEXT_SETTINGS,
        "allow_extra_args": True,
        "ignore_unknown_options": True,
    },
)
@click.option(
    "--format",
    "target_format",
    required=True,
    help="Target writer format id.",
    shell_complete=_writer_format_complete,
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    default=None,
    help="Directory for generated files.",
)
@click.option(
    "--frame", default="-1", show_default=True, help="Frame selection: all, int, or csv ints."
)
@click.option(
    "--embed/--no-embed", default=True, show_default=True, help="Embed selected frames in one file."
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED, cls=DynamicFormatOptionsArgument)
def format_transform(
    target_format: str,
    output_dir: Path | None = None,
    frame: str = "-1",
    embed: bool = True,
    n_jobs: int | None = None,
    extra_args: tuple[str, ...] = (),
) -> OperationCall:
    """Transform the current batch to another file format."""
    return _operation_call(
        "format-transform",
        FormatTransformParams(
            format=target_format,
            output_dir=output_dir,
            frame=frame,
            embed=embed,
            n_jobs=n_jobs,
            format_options=parse_dynamic_options(extra_args),
        ),
    )


@parse.command("to-summary-df")
@click.option(
    "--mode",
    type=click.Choice(["file", "frame"]),
    default="frame",
    show_default=True,
    help="Summary mode.",
)
@click.option(
    "--frame", default="-1", show_default=True, help="Frame selection: all, int, or csv ints."
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
@click.option(
    "--out", type=click.Path(path_type=Path, dir_okay=False), default=None, help="Output file."
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["csv", "json"]),
    default="csv",
    show_default=True,
    help="Output file format.",
)
def to_summary_df(
    mode: Literal["file", "frame"] = "frame",
    frame: str = "-1",
    n_jobs: int | None = None,
    out: Path | None = None,
    output_format: Literal["csv", "json"] = "csv",
) -> OperationCall:
    """Build a summary DataFrame from the current batch."""
    return _operation_call(
        "to-summary-df",
        ToSummaryDfParams(mode=mode, frame=frame, n_jobs=n_jobs, out=out, format=output_format),
    )


@parse.command("draw-grid-image")
@click.option(
    "--out", type=click.Path(path_type=Path, dir_okay=False), required=True, help="Image path."
)
@click.option("--mols-per-row", default=4, show_default=True, help="Molecules per row.")
@click.option("--sub-img-width", default=200, show_default=True, help="Sub-image width.")
@click.option("--sub-img-height", default=200, show_default=True, help="Sub-image height.")
@click.option("--max-mols", default=16, show_default=True, help="Maximum molecules to render.")
@click.option("--use-svg/--no-use-svg", default=None, help="Force SVG output mode.")
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def draw_grid_image(
    out: Path,
    mols_per_row: int = 4,
    sub_img_width: int = 200,
    sub_img_height: int = 200,
    max_mols: int = 16,
    use_svg: bool | None = None,
    n_jobs: int | None = None,
) -> OperationCall:
    """Render a molecule grid image from the current batch."""
    return _operation_call(
        "draw-grid-image",
        DrawGridImageParams(
            out=out,
            mols_per_row=mols_per_row,
            sub_img_width=sub_img_width,
            sub_img_height=sub_img_height,
            max_mols=max_mols,
            use_svg=use_svg,
            n_jobs=n_jobs,
        ),
    )


@parse.command("groupby")
@click.option(
    "--key",
    type=click.Choice(["detected_format_id", "file_format", "state"]),
    default="detected_format_id",
    show_default=True,
    help="Grouping key.",
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def groupby(
    key: Literal["detected_format_id", "file_format", "state"] = "detected_format_id",
    n_jobs: int | None = None,
) -> OperationCall:
    """Group current batch files and print grouped paths."""
    return _operation_call("groupby", GroupByParams(key=key, n_jobs=n_jobs))


@parse.command("copy-to")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    required=True,
    help="Destination directory.",
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def copy_to(output_dir: Path, n_jobs: int | None = None) -> OperationCall:
    """Copy current batch files to a directory."""
    return _operation_call("copy-to", CopyToParams(output_dir=output_dir, n_jobs=n_jobs))


@parse.command("move-to")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    required=True,
    help="Destination directory.",
)
@click.option(
    "--n-jobs", type=int, default=None, help="Number of parallel jobs for this operation."
)
def move_to(output_dir: Path, n_jobs: int | None = None) -> OperationCall:
    """Move current batch files to a directory."""
    return _operation_call("move-to", MoveToParams(output_dir=output_dir, n_jobs=n_jobs))


app.add_command(parse)


if __name__ == "__main__":
    app()
