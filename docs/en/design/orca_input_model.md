# ORCA Input Model and Parser Guide

This page documents the structured reader model for ORCA `.inp` files in MolOP. The ORCA writer is not registered yet; do not rely on `orca_raw_preamble` / `orca_raw_postamble` or raw round-trip compatibility.

## Data Layers

The ORCA file parser splits `$new_job` into multiple `ORCAInpFileFrame` objects, following the same framing idea as Gaussian `--Link1--`:

```text
ORCAInpFile
  frames: list[ORCAInpFileFrame]

ORCAInpFileFrame
  comment_lines
  keyword_lines
  blocks
  geometry
  trailing_lines
  auxiliary_basis_set
  dispersion_correction
  has_mixed_basis
  output_print_settings
  atoms / coords / charge / multiplicity
```

`ORCAInpFileFrame` must directly carry coordinates and primary semantics. `geometry` is the structured source field for the frame, not an extra job wrapper.

## Frame Fields

`comment_lines`

: Top-level `#` comments without the leading marker.

`keyword_lines`

: Top-level `!` keyword lines without the leading marker. The frame validator joins these into the common `keywords` field and derives `method`, `functional`, and `basis_set` when possible. ORCA auxiliary-basis keywords are stored in `auxiliary_basis_set`, such as `def2-SVP/C`. ORCA dispersion keywords are also exposed in `dispersion_correction`, such as `D3BJ` or `D4`; when a DFT functional is recognized, the common `functional` field includes the dispersion suffix, matching the Gaussian behavior (`B3LYP-D3BJ`, `PBE0-D4`).

`blocks`

: Top-level `%...` blocks. Each `ORCABlock` contains:

- `name`: lowercase block name without `%`, such as `pal`, `maxcore`, or `scf`
- `raw_header`: original header line
- `lines`: block body lines without newline terminators

`geometry`

: Structured ORCA geometry input. The parser supports `* xyz`, `* cart`, `* xyzfile`, `* gzmtfile`, and Cartesian coordinates inside `%coords`. Cartesian coordinates are synchronized onto the frame-level `atoms` / `coords` fields.

`trailing_lines`

: Non-empty top-level text that is not classified yet. This is for diagnostics and future extension, not a raw compatibility interface.

`dispersion_correction`

: Dispersion correction recognized from ORCA `!` keyword lines. This is an ORCA input frame field, not a common `BaseQMInputFrame` field. When the frame also has a recognized DFT functional, the dispersion is appended to the common `functional` string as `FUNCTIONAL-DISPERSION`.

`has_mixed_basis`

: `true` when any atom row contains `newgto` or `newauxgto`.

`output_print_settings`

: `Print[ ... ] = ...` settings recognized from `%output` blocks. The original `%output` block remains preserved in `blocks` and `resources_raw`.

## Geometry Model

`ORCAGeometry` fields:

- `ctype`: `xyz`, `cart`, `cartesian`, `int`, `internal`, `gzmt`, `xyzfile`, or `gzmtfile`
- `charge` / `multiplicity`
- `units`, such as `bohr`
- `external_path` for external geometry files
- `atoms` for Cartesian atom rows
- `point_charges` for `Q q x y z` rows
- `source`: `star` or `percent_coords`

`ORCAGeometryAtom` fields:

- `symbol` / `atomic_number`
- `x` / `y` / `z`
- `is_dummy` / `is_ghost`
- `fragment_id`
- `frozen`
- `isotope`
- `nuclear_charge`
- `basis_set` / `auxiliary_basis_set` from atom-row `newgto` / `newauxgto`
- `basis_overrides` for atom-level basis override directives, including directive kind, named basis, and tokens

## Parse Flow

File parser:

1. Read the `.inp` file text.
2. Split top-level `$new_job` sections into frames.
3. Parse each frame with `ORCAInpFileFrameParser`.

Frame parser:

1. Detect a `* ...` geometry section, or a `%coords ... end` geometry block.
2. Scan top-level `%` blocks into `ORCABlock`.
3. Ignore lines occupied by geometry and `%` blocks.
4. Extract `#` comments and `!` keywords from remaining top-level lines.
5. Build top-level `ORCAInpFileFrame` fields.
6. Let the frame validator derive common QM input fields.

## Units and Coordinates

Cartesian input is treated as Å by default. If `units` is `bohr`, `a0`, `au`, or an equivalent token, coordinates and point-charge positions are converted to Å before being stored.

`xyzfile` / `gzmtfile` preserve only `external_path`; the parser does not read external files or fabricate `atoms` / `coords`.

## Development Constraints

- Do not add a `jobs` wrapper layer; `$new_job` maps to multiple frames.
- Do not restore `orca_raw_preamble` / `orca_raw_postamble`.
- Do not register an ORCA writer until a structured renderer exists.
- Write program-specific syntax into ORCA structured fields and common QM structured containers first, then use `project_common_qm_fields()` to project compatibility fields such as `method`, `functional`, and `basis_set`.
- Do not hide semantics in raw strings.
