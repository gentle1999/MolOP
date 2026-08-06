# Transform Behaviors

The `format_transform` method is available on both `FileBatchModelDisk` and
individual file objects. It converts between chemical file formats while
preserving as much structure-level information as the target format allows.

For parameter defaults, return values, and error policy, see
[API Contracts](api_contracts.md).

=== "Render in memory"

    ```python
    from molop import AutoParser

    batch = AutoParser("results/*.out", n_jobs=1)
    rendered = batch.format_transform("xyz", frame=-1, write_to_disk=False)
    ```

    ??? example "Return shape"

        ```text
        dict 1
        ```

    The result is a mapping from each absolute source path to rendered text (or
    a list of texts when `embed_in_one_file=False`). No output file is created.

=== "Write to disk"

    ```python
    from pathlib import Path
    from molop import AutoParser

    Path("structures").mkdir(parents=True, exist_ok=True)
    batch = AutoParser("results/*.out", n_jobs=1)
    batch.format_transform(
        "xyz",
        output_dir="structures",
        write_to_disk=True,
    )
    ```

    ??? example "Created file"

        ```text
        structures/water_mp2.xyz
        ```

    The batch API requires `output_dir` to exist. The CLI creates its
    `--output-dir`; both APIs replace only the source file's final suffix.

!!! warning
    A coordinate writer cannot preserve properties that the target format has
    no representation for. Use a graph-capable target for bonding data and a
    QM-aware target when energies or frequencies must remain available.

## Key Behaviors

- **Frame Selection**: By default, only the last frame (`frame=-1`) is transformed. You can specify `frame="all"` or a sequence of frame IDs to transform more frames.
- **Embedding**: If `embed_in_one_file=True` (default), multiple frames are combined into a single output file if the format supports it (e.g., SDF or multi-frame XYZ).
- **File Output**: Python callers must pass `write_to_disk=True` to write rendered
  content to disk. When writing, `file_path` or batch `output_dir` selects the
  output location; if no path is provided, MolOP writes beside the source file.
  Output names replace only the final suffix, so `name.hash.log` becomes
  `name.hash.xyz` when converting to XYZ.
  When `write_to_disk=False`, `file_path` and `output_dir` are ignored and the
  transform only returns the rendered string or list of strings.
- **Structure Level**:
    - **COORDS (Coordinate Level)**: Formats such as `xyz`, `gjf`, and `orcainp` primarily preserve atomic coordinates and elements.
    - **GRAPH (Graph Level)**: Formats like `sdf`, `smi`, and `cml` preserve bonding information (molecular graph). If the source file only has coordinates (e.g., a `.log` file), MolOP will automatically attempt to reconstruct the molecular graph using its built-in algorithms.
- **Metadata Preservation**:
    - The `gjf` writer preserves structured Gaussian directives and keywords. The `orcainp` reader parses ORCA keywords, blocks, and geometry into frame fields; its registered canonical writer can build a new input from a coordinate-bearing frame plus explicit `keywords`, resources, and blocks, but does not promise source-preserving round trips.
    - Computational properties (energies, frequencies) are generally **NOT** preserved when transforming to simple coordinate formats like XYZ, although some formats like SDF can store them as properties.
