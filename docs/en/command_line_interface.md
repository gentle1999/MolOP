# Command Line Interface

MolOP provides a modern command-line interface built with [Typer](https://typer.tiangolo.com/). The primary entrypoint is `molop <command>`, which offers a structured and user-friendly way to interact with your molecular data.

For in-repo reproducibility, all examples here use `uv run molop ...`. If MolOP is installed in your environment, you can simply use `molop ...`.

## Global Options

MolOP supports several global options that can be used with any command:

- `-v`, `--verbose`: Enable verbose output for detailed logging.
- `-q`, `--quiet`: Enable quiet mode. This is highly recommended for stable, clean outputs in scripts or documentation.
- `--version`: Show the version and exit.

## Discovering Commands

You can view all available commands and global options by running:

```bash
uv run molop --help
```

```text

 Usage: molop [OPTIONS] COMMAND [ARGS]...                                       

 MolOP: Molecule OPerator CLI                                                   

╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --verbose             -v        Enable verbose output.                       │
│ --quiet               -q        Enable quiet mode.                           │
│ --version                       Show the version and exit.                   │
│ --install-completion            Install completion for the current shell.    │
│ --show-completion               Show completion for the current shell, to    │
│                                 copy it or customize the installation.       │
│ --help                -h        Show this message and exit.                  │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ summary           Generate a summary of the molecules in the given files.    │
│ visualize         Visualize molecules in a grid image.                       │
│ transform         Transform molecular files to another format.               │
│ filter-by-codec   Filter files by their detected codec ID.                   │
│ sample            Randomly sample N files from the matched files.            │
│ stats             Show statistics of the molecules in the given files.       │
│ groupby           Group files by a key and output the groups.                │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Typer Commands

### `summary`

Generate a summary of the molecules in the given files, typically outputting to a CSV or JSON file.

**When to use:** Use this to extract key properties (like energy, SMILES, etc.) from a batch of files into a structured table.

```bash
uv run molop summary --help
```

```text

 Usage: molop summary [OPTIONS] PATTERN                                         

 Generate a summary of the molecules in the given files.                        

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --out               -o      PATH     Output file path.                       │
│ --format            -f      TEXT     Output format (csv or json).            │
│                                      [default: csv]                          │
│ --mode              -m      TEXT     Summary mode (file or frame).           │
│                                      [default: frame]                        │
│ --frame                     TEXT     Frame selection (e.g., 'all', '-1',     │
│                                      '1,2').                                 │
│                                      [default: -1]                           │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q summary "tests/test_files/g16log/2-TS1-Opt.log" --out .sisyphus/tmp/doc_summary.csv --format csv --mode frame --frame -1 --n-jobs 1
```

```text
Summary written to .sisyphus/tmp/doc_summary.csv
```

### `transform`

Transform molecular files from one format to another (e.g., Gaussian log to SDF).

**When to use:** Use this for batch format conversion.

```bash
uv run molop transform --help
```

```text

 Usage: molop transform [OPTIONS] PATTERN                                       

 Transform molecular files to another format.                                   

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --to                                  TEXT     Target format (xyz, sdf,   │
│                                                   cml, gjf, smi, orcainp).   │
│                                                   [required]                 │
│ *  --output-dir                          PATH     Directory to save          │
│                                                   transformed files.         │
│                                                   [required]                 │
│    --frame                               TEXT     Frame selection (e.g.,     │
│                                                   'all', '-1', '1,2').       │
│                                                   [default: -1]              │
│    --embed                 --no-embed             Whether to embed multiple  │
│                                                   frames in one file.        │
│                                                   [default: embed]           │
│    --parser-detection                    TEXT     Parser detection mode.     │
│                                                   [default: auto]            │
│    --n-jobs            -j                INTEGER  Number of parallel jobs.   │
│                                                   [default: -1]              │
│    --help              -h                         Show this message and      │
│                                                   exit.                      │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q transform "tests/test_files/orca/h2_grad_orca.inp" --to sdf --output-dir .sisyphus/tmp/doc_transform_out --frame -1 --embed --parser-detection orcainp --n-jobs 1
```

```text
Transformation completed. Files saved to .sisyphus/tmp/doc_transform_out
```

### `visualize`

Visualize molecules in a grid image (PNG or SVG).

**When to use:** Use this to quickly inspect the structures of a batch of molecules visually.

```bash
uv run molop visualize --help
```

```text

 Usage: molop visualize [OPTIONS] PATTERN                                       

 Visualize molecules in a grid image.                                           

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --out                       PATH     Output file path (e.g., grid.png or  │
│                                         grid.svg).                           │
│                                         [required]                           │
│    --max-mols                  INTEGER  Maximum number of molecules to       │
│                                         visualize.                           │
│                                         [default: 16]                        │
│    --mols-per-row              INTEGER  Number of molecules per row.         │
│                                         [default: 4]                         │
│    --sub-img-width             INTEGER  Width of each sub-image.             │
│                                         [default: 200]                       │
│    --sub-img-height            INTEGER  Height of each sub-image.            │
│                                         [default: 200]                       │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q visualize "tests/test_files/g16log/2-TS1-Opt.log" --out .sisyphus/tmp/doc_grid.svg --max-mols 4 --mols-per-row 2 --sub-img-width 200 --sub-img-height 200 --n-jobs 1 --parser-detection auto
```

```text
Visualization saved to .sisyphus/tmp/doc_grid.svg
```

### `filter-by-codec`

Filter files by their detected codec ID (e.g., `g16log`, `xyz`).

**When to use:** Use this to separate files of different formats or to ensure only specific formats are processed.

```bash
uv run molop filter-by-codec --help
```

```text

 Usage: molop filter-by-codec [OPTIONS] PATTERN                                 

 Filter files by their detected codec ID.                                       

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --codec-id                  TEXT     Codec ID to filter by. [required]    │
│    --negate                             Negate the filter.                   │
│    --on-missing                TEXT     Behavior on missing codec ID (keep,  │
│                                         drop, error).                        │
│                                         [default: drop]                      │
│    --format            -f      TEXT     Output format (text or json).        │
│                                         [default: text]                      │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q filter-by-codec "tests/test_files/g16log/2-TS1-Opt.log" --codec-id g16log --format text --n-jobs 1
```

```text
/Users/tmj/Documents/proj/MolOP/tests/test_files/g16log/2-TS1-Opt.log
```

### `sample`

Randomly sample N files from the matched files.

**When to use:** Use this to pick a representative subset of files for testing or quick inspection.

```bash
uv run molop sample --help
```

```text

 Usage: molop sample [OPTIONS] PATTERN                                          

 Randomly sample N files from the matched files.                                

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ *  --n                         INTEGER  Number of files to sample.           │
│                                         [required]                           │
│    --seed                      INTEGER  Random seed.                         │
│    --format            -f      TEXT     Output format (text or json).        │
│                                         [default: text]                      │
│    --parser-detection          TEXT     Parser detection mode.               │
│                                         [default: auto]                      │
│    --n-jobs            -j      INTEGER  Number of parallel jobs.             │
│                                         [default: -1]                        │
│    --help              -h               Show this message and exit.          │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q sample "tests/test_files/mix_format/*.gjf" --n 2 --seed 1 --format text --n-jobs 1
```

```text
/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/S_Ph_Ni_TS.gjf
/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.gjf
```

### `stats`

Show statistics of the molecules in the given files.

**When to use:** Use this to get a high-level overview of the dataset (e.g., format distribution, calculation status).

```bash
uv run molop stats --help
```

```text

 Usage: molop stats [OPTIONS] PATTERN                                           

 Show statistics of the molecules in the given files.                           

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --format            -f      TEXT     Output format (text or json).           │
│                                      [default: text]                         │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q stats "tests/test_files/g16log/2-TS1-Opt.log"
```

```text
Total files: 1

By format:
  g16log: 1

By state (last frame):
  normal: 1
```

### `groupby`

Group files by a key (e.g., `detected_format_id`, `state`) and output the groups.

**When to use:** Use this to organize files into categories for further processing.

```bash
uv run molop groupby --help
```

```text

 Usage: molop groupby [OPTIONS] PATTERN                                         

 Group files by a key and output the groups.                                    

╭─ Arguments ──────────────────────────────────────────────────────────────────╮
│ *    pattern      TEXT  File pattern to match. [required]                    │
╰──────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────╮
│ --key                       TEXT     Key to group by (detected_format_id,    │
│                                      file_format, state).                    │
│                                      [default: detected_format_id]           │
│ --format            -f      TEXT     Output format (text or json).           │
│                                      [default: json]                         │
│ --parser-detection          TEXT     Parser detection mode. [default: auto]  │
│ --n-jobs            -j      INTEGER  Number of parallel jobs. [default: -1]  │
│ --help              -h               Show this message and exit.             │
╰──────────────────────────────────────────────────────────────────────────────╯
```

**Example:**

```bash
uv run molop -q groupby "tests/test_files/mix_format/*.gjf" --key detected_format_id --format json --n-jobs 1
```

```text
{
  "gjf": [
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/S_Ph_Ni_TS.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/TS_Zy0fwX_ll_ad_14-19_15-16_optts_g16.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/dsgdb9nsd_000007-6.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/irc.gjf",
    "/Users/tmj/Documents/proj/MolOP/tests/test_files/mix_format/opt_point_charge_xtb.gjf"
  ]
}
```

