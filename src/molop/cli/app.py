from pathlib import Path
from typing import Annotated

import typer


app = typer.Typer(
    name="molop",
    help="MolOP: Molecule OPerator CLI",
    add_completion=True,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def version_callback(value: bool):
    if value:
        import importlib.metadata

        try:
            version = importlib.metadata.version("molop")
        except importlib.metadata.PackageNotFoundError:
            version = "unknown"
        typer.echo(f"molop version: {version}")
        raise typer.Exit()


@app.callback()
def main(
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose output.")
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", "-q", help="Enable quiet mode.")] = False,
    _version: Annotated[
        bool | None,
        typer.Option(
            "--version",
            callback=version_callback,
            is_eager=True,
            help="Show the version and exit.",
        ),
    ] = None,
):
    """
    MolOP: Molecule OPerator CLI
    """
    if verbose and quiet:
        # If both are set, quiet takes precedence.
        verbose = False

    if verbose or quiet:
        from molop.config import molopconfig

        if quiet:
            molopconfig.quiet()
        if verbose:
            molopconfig.verbose()


@app.command(help="Generate a summary of the molecules in the given files.")
def summary(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    out: Annotated[Path | None, typer.Option("--out", "-o", help="Output file path.")] = None,
    format: Annotated[
        str, typer.Option("--format", "-f", help="Output format (csv or json).")
    ] = "csv",
    mode: Annotated[
        str, typer.Option("--mode", "-m", help="Summary mode (file or frame).")
    ] = "frame",
    frame: Annotated[
        str, typer.Option("--frame", help="Frame selection (e.g., 'all', '-1', '1,2').")
    ] = "-1",
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    """
    Generate a summary of the molecules in the given files.
    """
    import json
    from collections.abc import Sequence
    from typing import Literal, cast

    import pandas as pd

    from molop.cli.shared.frames import parse_frame_selection
    from molop.io import AutoParser

    if format not in ["csv", "json"]:
        typer.echo(f"Error: Invalid format '{format}'. Supported formats: csv, json", err=True)
        raise typer.Exit(1)

    if mode not in ["file", "frame"]:
        typer.echo(f"Error: Invalid mode '{mode}'. Supported modes: file, frame", err=True)
        raise typer.Exit(1)

    output_path = out if out is not None else Path(f"summary.{format}")

    try:
        frame_selection = parse_frame_selection(frame)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from e

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        typer.echo(f"No files found matching pattern: {pattern}", err=True)
        raise typer.Exit(0)

    frame_ids: int | Sequence[int] | Literal["all"] = -1
    if mode == "frame" and frame_selection == "all":
        max_frames = max(len(f) for f in batch.values())
        frame_ids = list(range(max_frames))
    else:
        frame_ids = -1 if frame_selection == "all" else frame_selection

    df = batch.to_summary_df(
        mode=cast(Literal["file", "frame"], mode),
        frameIDs=cast(int | Sequence[int], frame_ids),
        n_jobs=n_jobs,
    )

    if df.empty:
        typer.echo("Summary is empty.", err=True)
        raise typer.Exit(0)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == "csv":
        df.to_csv(output_path, index=False)
    elif format == "json":
        if isinstance(df.columns, pd.MultiIndex):
            df.columns = [
                "/".join(str(level) for level in col if level).strip("/")
                for col in df.columns.values
            ]

        def serialize_val(val):
            if hasattr(val, "magnitude"):
                val = val.magnitude
            if hasattr(val, "tolist"):
                return val.tolist()
            return val

        # Use applymap for older pandas compatibility if needed,
        # but map is preferred in 2.1+.
        # We use a helper to ensure compatibility.
        def safe_map(df, func):
            if hasattr(df, "map"):
                return df.map(func)
            return df.applymap(func)

        data = safe_map(df, serialize_val).to_dict(orient="records")

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=True, indent=2)

    typer.echo(f"Summary written to {output_path}", err=True)


@app.command(help="Visualize molecules in a grid image.")
def visualize(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    out: Annotated[
        Path, typer.Option("--out", help="Output file path (e.g., grid.png or grid.svg).")
    ],
    max_mols: Annotated[
        int, typer.Option("--max-mols", help="Maximum number of molecules to visualize.")
    ] = 16,
    mols_per_row: Annotated[
        int, typer.Option("--mols-per-row", help="Number of molecules per row.")
    ] = 4,
    sub_img_width: Annotated[
        int, typer.Option("--sub-img-width", help="Width of each sub-image.")
    ] = 200,
    sub_img_height: Annotated[
        int, typer.Option("--sub-img-height", help="Height of each sub-image.")
    ] = 200,
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
):
    """
    Visualize molecules in a grid image.
    """
    from molop.io import AutoParser

    use_svg = out.suffix.lower() == ".svg"

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        typer.echo(f"No files found matching pattern: {pattern}", err=True)
        raise typer.Exit(0)

    img = batch.draw_grid_image(
        molsPerRow=mols_per_row,
        subImgSize=(sub_img_width, sub_img_height),
        maxMols=max_mols,
        useSVG=use_svg,
        n_jobs=n_jobs,
    )

    if img is None:
        raise typer.Exit(0)

    out.parent.mkdir(parents=True, exist_ok=True)

    if use_svg:
        content = img.data if hasattr(img, "data") else img
        with open(out, "w", encoding="utf-8") as f:
            f.write(content)
    else:
        if hasattr(img, "save"):
            img.save(out)
        elif hasattr(img, "data"):
            with open(out, "wb") as f:
                f.write(img.data)
        else:
            typer.echo("Error: Unsupported image object returned by RDKit.", err=True)
            raise typer.Exit(1)

    typer.echo(f"Visualization saved to {out}", err=True)


@app.command(help="")
def transform(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    to: Annotated[
        str, typer.Option("--to", help="Target format (xyz, sdf, cml, gjf, smi, orcainp).")
    ],
    output_dir: Annotated[
        Path, typer.Option("--output-dir", help="Directory to save transformed files.")
    ],
    frame: Annotated[
        str, typer.Option("--frame", help="Frame selection (e.g., 'all', '-1', '1,2').")
    ] = "-1",
    embed: Annotated[
        bool,
        typer.Option("--embed/--no-embed", help="Whether to embed multiple frames in one file."),
    ] = True,
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    from typing import Any, cast

    from molop.cli.shared.frames import parse_frame_selection
    from molop.io import AutoParser, codec_registry

    valid_formats = codec_registry.get_supported_writer_formats()
    if to not in valid_formats:
        typer.echo(
            f"Error: Invalid format '{to}'. Supported formats: {', '.join(valid_formats)}",
            err=True,
        )
        raise typer.Exit(1)

    try:
        frame_selection = parse_frame_selection(frame)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1) from e

    output_dir.mkdir(parents=True, exist_ok=True)

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        typer.echo(f"No files found matching pattern: {pattern}", err=True)
        raise typer.Exit(0)

    batch.format_transform(
        format=cast(Any, to),
        output_dir=str(output_dir),
        frameID=frame_selection,
        embed_in_one_file=embed,
        n_jobs=n_jobs,
    )

    typer.echo(f"Transformation completed. Files saved to {output_dir}", err=True)


@app.command(help="Filter files by their detected codec ID.")
def filter_by_codec(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    codec_id: Annotated[str, typer.Option("--codec-id", help="Codec ID to filter by.")],
    negate: Annotated[bool, typer.Option("--negate", help="Negate the filter.")] = False,
    on_missing: Annotated[
        str, typer.Option("--on-missing", help="Behavior on missing codec ID (keep, drop, error).")
    ] = "drop",
    out_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format (text or json).")
    ] = "text",
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    """
    Filter files by their detected codec ID.
    """
    import json
    from typing import Literal, cast

    from molop.io import AutoParser

    if out_format not in ["text", "json"]:
        typer.echo(f"Error: Invalid format '{out_format}'. Supported formats: text, json", err=True)
        raise typer.Exit(1)

    if on_missing not in ["keep", "drop", "error"]:
        typer.echo(
            f"Error: Invalid on-missing behavior '{on_missing}'. Supported: keep, drop, error",
            err=True,
        )
        raise typer.Exit(1)

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        if out_format == "json":
            typer.echo("[]")
        return

    filtered_batch = batch.filter_by_codec_id(
        codec_id=codec_id,
        negate=negate,
        n_jobs=n_jobs,
        on_missing=cast(Literal["keep", "drop", "error"], on_missing),
    )

    paths = filtered_batch.file_paths

    if out_format == "json":
        typer.echo(json.dumps(paths))
    else:
        for path in paths:
            typer.echo(path)


@app.command(help="Randomly sample N files from the matched files.")
def sample(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    n: Annotated[int, typer.Option("--n", help="Number of files to sample.")],
    seed: Annotated[int | None, typer.Option("--seed", help="Random seed.")] = None,
    out_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format (text or json).")
    ] = "text",
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    """
    Randomly sample N files from the matched files.
    """
    import json

    from molop.io import AutoParser

    if out_format not in ["text", "json"]:
        typer.echo(f"Error: Invalid format '{out_format}'. Supported formats: text, json", err=True)
        raise typer.Exit(1)

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        if out_format == "json":
            typer.echo("[]")
        return

    sampled_batch = batch.sample(n=n, seed=seed)
    paths = sampled_batch.file_paths

    if out_format == "json":
        typer.echo(json.dumps(paths))
    else:
        for path in paths:
            typer.echo(path)


@app.command(help="Show statistics of the molecules in the given files.")
def stats(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    format: Annotated[
        str, typer.Option("--format", "-f", help="Output format (text or json).")
    ] = "text",
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    """
    Show statistics of the molecules in the given files.
    """
    import json
    from collections import Counter

    from molop.io import AutoParser

    if format not in ["text", "json"]:
        typer.echo(f"Error: Invalid format '{format}'. Supported formats: text, json", err=True)
        raise typer.Exit(1)

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        typer.echo(f"No files found matching pattern: {pattern}", err=True)
        raise typer.Exit(0)

    total = len(batch)
    format_counts: Counter[str] = Counter()
    state_counts: Counter[str] = Counter()

    for diskfile in batch.values():
        fmt = getattr(diskfile, "detected_format_id", "unknown")
        if fmt is None:
            fmt = "unknown"
        format_counts[fmt] += 1

        if not diskfile:
            state_counts["unknown"] += 1
            continue

        last_frame = diskfile[-1]
        if getattr(last_frame, "is_TS", False):
            state_counts["ts"] += 1
        elif getattr(last_frame, "is_error", False):
            state_counts["error"] += 1
        elif getattr(last_frame, "is_normal", False):
            state_counts["normal"] += 1
        else:
            state_counts["unknown"] += 1

    if format == "json":
        output = {
            "total": total,
            "by_detected_format_id": dict(sorted(format_counts.items())),
            "by_state": dict(sorted(state_counts.items())),
        }
        typer.echo(json.dumps(output, indent=2, sort_keys=True))
    else:
        typer.echo(f"Total files: {total}")
        typer.echo("\nBy format:")
        for fmt, count in sorted(format_counts.items()):
            typer.echo(f"  {fmt}: {count}")
        typer.echo("\nBy state (last frame):")
        for state, count in sorted(state_counts.items()):
            typer.echo(f"  {state}: {count}")


@app.command(help="Group files by a key and output the groups.")
def groupby(
    pattern: Annotated[str, typer.Argument(help="File pattern to match.")],
    key: Annotated[
        str,
        typer.Option(
            "--key",
            help="Key to group by (detected_format_id, file_format, state).",
        ),
    ] = "detected_format_id",
    format: Annotated[
        str, typer.Option("--format", "-f", help="Output format (text or json).")
    ] = "json",
    parser_detection: Annotated[
        str, typer.Option("--parser-detection", help="Parser detection mode.")
    ] = "auto",
    n_jobs: Annotated[int, typer.Option("--n-jobs", "-j", help="Number of parallel jobs.")] = -1,
):
    """
    Group files by a key and output the groups.
    """

    from molop.io import AutoParser

    if format not in ["text", "json"]:
        typer.echo(f"Error: Invalid format '{format}'. Supported formats: text, json", err=True)
        raise typer.Exit(1)

    if key not in ["detected_format_id", "file_format", "state"]:
        typer.echo(
            f"Error: Invalid key '{key}'. Supported keys: detected_format_id, file_format, state",
            err=True,
        )
        raise typer.Exit(1)

    batch = AutoParser(pattern, n_jobs=n_jobs, parser_detection=parser_detection)

    if not batch:
        if format == "json":
            typer.echo("{}")
        return

    def key_func(diskfile) -> str:
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

    groups = batch.groupby(key_func=key_func, n_jobs=n_jobs)

    import json

    if format == "json":
        output = {k: v.file_paths for k, v in groups.items()}
        typer.echo(json.dumps(output, indent=2, sort_keys=True))
    else:
        for k, v in sorted(groups.items()):
            for path in v.file_paths:
                typer.echo(f"{k}\t{path}")


if __name__ == "__main__":
    app()
