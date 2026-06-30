# Format Support

<!-- format-support-overview -->

This page is the overview of MolOP's built-in format support. Detailed support
scope is documented on each format page. The machine-checkable source of truth is
`tests/format_feature_coverage/support_matrix.py`, and
`tests/format_feature_coverage/test_support_matrix.py` verifies that the matrix
matches the registered codecs, declares support scope and limits, and points to
real tests.

| Group | Formats | Summary |
| ----- | ------- | ------- |
| Structure formats | [`xyz`](formats/xyz.md), [`sdf`](formats/sdf.md), [`smi`](formats/smi.md), [`cml`](formats/cml.md) | Coordinate/graph structure IO and rendering. |
| QM input formats | [`gjf`](formats/gjf.md), [`orcainp`](formats/orcainp.md) | Gaussian and ORCA input parsing, with Gaussian input rendering currently available. |
| QM output formats | [`g16log`](formats/g16log.md), [`orcaout`](formats/orcaout.md), [`fakeg`](formats/fakeg.md) | Gaussian/ORCA output parsing plus Gaussian-like rendering from parsed Gaussian output. |
| Special readers | [OpenBabel fallback](formats/openbabel-fallback.md) | Unknown-extension fallback through OpenBabel-compatible readers. |

| Format ID | Typical Extensions | Read | Write | Support Page |
| --------- | ------------------ | ---- | ----- | ------------ |
| `xyz` | `.xyz` | Yes | Yes | [XYZ](formats/xyz.md) |
| `sdf` | `.sdf`, `.sd`, `.mol` | Yes | Yes | [SDF/MOL](formats/sdf.md) |
| `smi` | `.smi`, `.txt` | Yes | Yes | [SMILES](formats/smi.md) |
| `gjf` | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` | Yes | Yes | [Gaussian input](formats/gjf.md) |
| `g16log` | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` | Yes | No | [Gaussian output](formats/g16log.md) |
| `orcainp` | `.inp` | Yes | No | [ORCA input](formats/orcainp.md) |
| `orcaout` | `.out`, `.log`, `.orcaout` | Yes | No | [ORCA output](formats/orcaout.md) |
| `cml` | `.cml` | No | Yes | [CML writer](formats/cml.md) |
| `fakeg` | `.fakeg` | No | Yes | [Gaussian-like renderer](formats/fakeg.md) |
| OpenBabel fallback | Unknown or OpenBabel-supported extensions | Yes | No | [OpenBabel fallback](formats/openbabel-fallback.md) |

!!! note
    Suffixes are not always unique. For example `.out` and `.log` can be ORCA
    or Gaussian output, and `.gau` can be Gaussian input or output. Reader
    selection uses extension candidates plus content checks; non-matching readers
    must raise `FormatMismatchError` so that the next candidate can be tried.
