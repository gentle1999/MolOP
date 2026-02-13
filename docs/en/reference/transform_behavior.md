<!--
 * @Author: TMJ
 * @Date: 2026-02-12 12:57:10
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-13 11:56:32
 * @Description: 请填写简介
-->

# Transform Behaviors

The `format_transform` method (available on `FileBatchModelDisk` and individual file objects) allows converting between different chemical file formats.

## Key Behaviors

- **Frame Selection**: By default, only the last frame (`frameID=-1`) is transformed. You can specify `frameID="all"` to transform all frames in a file.
- **Embedding**: If `embed_in_one_file=True` (default), multiple frames are combined into a single output file if the format supports it (e.g., SDF or multi-frame XYZ).
- **Structure Level**:
  - **COORDS (Coordinate Level)**: Formats like `xyz`, `gjf`, and `orcainp` primarily preserve atomic coordinates and elements.
  - **GRAPH (Graph Level)**: Formats like `sdf`, `smi`, and `cml` preserve bonding information (molecular graph). If the source file only has coordinates (e.g., a `.log` file), MolOP will automatically attempt to reconstruct the molecular graph using its built-in algorithms.
- **Metadata Preservation**:
  - QM input formats (`gjf`, `orcainp`) attempt to preserve raw directives and keywords.
  - Computational properties (energies, frequencies) are generally **NOT** preserved when transforming to simple coordinate formats like XYZ, although some formats like SDF can store them as properties.
