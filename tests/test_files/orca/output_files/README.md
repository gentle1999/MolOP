ORCA output fixtures for parser refactoring.

Purpose:
- Keep ORCA output parsing behavior fenced before replacing the old parser path.
- Cover local historical ORCA outputs plus representative external regression outputs.
- Avoid using ORCA input fixtures as a proxy for output behavior.

Layout:
- `local/`: ORCA output files already present in this repository, copied into an output-specific fixture tree.
- `cclib/`: selected ORCA output files from cclib's regression data.
- `manifest.json`: machine-readable inventory with source, feature tags, and raw text anchors.
- `cclib/LICENSE.cclib`: upstream BSD-3-Clause license text for the cclib fixtures.

cclib source:
- Repository: https://github.com/cclib/cclib
- Commit: `0e5d5dfd7d15b2fbd95e02fdd36bd4fe3336a5d5`
- Upstream directory: `data/ORCA`
- License: BSD-3-Clause, copied in `cclib/LICENSE.cclib`

Selection rule:
- Keep a compact but broad output set instead of mirroring all upstream data.
- Prefer fixtures that exercise different output sections: single point energies,
  post-HF energies, gradients, geometry optimization, frequencies, IR/Raman,
  polarizability, NMR, solvation, empirical dispersion, and excited-state methods.

Parser contract:
- Tests in `tests/test_orca_output_fixtures.py` validate file presence and raw anchors now.
- Structured parsing tests are intentionally marked `xfail` until the ORCA output parser is implemented.
