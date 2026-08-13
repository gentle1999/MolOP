from __future__ import annotations

from pathlib import Path
from typing import Literal

import click

from molop.cli.state_machine import OperationCall

app: click.Group
parse: click.Group
completion: click.Group

def completion_show(shell: Literal["auto", "bash", "zsh", "fish"] = "auto") -> None: ...
def completion_install(shell: Literal["auto", "bash", "zsh", "fish"] = "auto") -> None: ...
def execute_parse_chain(
    operations: list[OperationCall],
    pattern: str,
    additional_patterns: tuple[str, ...] = (),
    parser_detection: str = "auto",
    n_jobs: int = -1,
    output_format: Literal["text", "json"] = "text",
    total_charge: int | None = None,
    total_multiplicity: int | None = None,
    only_extract_structure: bool = False,
    only_last_frame: bool = False,
    capture_source_evidence: bool = False,
    source_encoding: str = "utf-8",
    release_file_content: bool = True,
    force_unit_transform: bool | None = None,
    graph_reconstruction_backend: Literal["cpp", "python"] | None = None,
    make_dative_bonds: bool | None = None,
    report: bool = False,
) -> None: ...
def filter_state(
    state: Literal["ts", "error", "opt", "normal", "thermal", "no-img"],
    negate: bool = False,
    n_jobs: int = -1,
) -> OperationCall: ...
def filter_value(
    target: Literal["charge", "multiplicity", "format"],
    value: str,
    compare: Literal["==", "!=", ">", "<", ">=", "<="] = "==",
    n_jobs: int = -1,
) -> OperationCall: ...
def filter_by_codec(
    codec_id: str,
    negate: bool = False,
    n_jobs: int = -1,
    on_missing: Literal["keep", "drop", "error"] = "drop",
) -> OperationCall: ...
def sample(n: int = 10, seed: int | None = None) -> OperationCall: ...
def format_transform(
    target_format: str,
    output_dir: Path | None = None,
    frame: str = "-1",
    embed: bool = True,
    write_to_disk: bool | None = None,
    n_jobs: int = -1,
    extra_args: tuple[str, ...] = (),
) -> OperationCall: ...
def to_summary_df(
    mode: Literal["file", "frame"] = "frame",
    frame: str = "-1",
    n_jobs: int = -1,
    brief: bool = True,
    flatten_columns: bool = True,
    on_missing_frame: Literal["skip", "error"] = "skip",
    out: Path | None = None,
    output_format: Literal["csv", "json"] = "csv",
) -> OperationCall: ...
def draw_grid_image(
    out: Path,
    mols_per_row: int = 4,
    sub_img_width: int = 200,
    sub_img_height: int = 200,
    max_mols: int = 16,
    use_svg: bool | None = None,
    n_jobs: int = -1,
) -> OperationCall: ...
def groupby(
    key: Literal["detected_format_id", "file_format", "state"] = "detected_format_id",
    n_jobs: int = -1,
) -> OperationCall: ...
def copy_to(output_dir: Path, n_jobs: int = -1) -> OperationCall: ...
def move_to(output_dir: Path, n_jobs: int = -1) -> OperationCall: ...
