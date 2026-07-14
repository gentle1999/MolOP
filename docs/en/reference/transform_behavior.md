<!--
 * @Author: TMJ
 * @Date: 2026-02-12 12:57:10
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-13 11:56:32
 * @Description: 请填写简介
-->

# Transform Behaviors

The `format_transform` method is available on both `FileBatchModelDisk` and
individual file objects. It converts between chemical file formats while
preserving as much structure-level information as the target format allows.

For parameter defaults, return values, and error policy, see
[API Contracts](api_contracts.md).

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
  - **COORDS (Coordinate Level)**: Formats like `xyz` and `gjf` primarily preserve atomic coordinates and elements. `orcainp` is currently a structured reader only.
  - **GRAPH (Graph Level)**: Formats like `sdf`, `smi`, and `cml` preserve bonding information (molecular graph). If the source file only has coordinates (e.g., a `.log` file), MolOP will automatically attempt to reconstruct the molecular graph using its built-in algorithms.
- **Metadata Preservation**:
  - The `gjf` writer preserves structured Gaussian directives and keywords. The `orcainp` reader parses ORCA keywords, blocks, and geometry into frame fields, but no ORCA writer is currently registered.
  - Computational properties (energies, frequencies) are generally **NOT** preserved when transforming to simple coordinate formats like XYZ, although some formats like SDF can store them as properties.
