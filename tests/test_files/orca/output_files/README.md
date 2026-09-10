# ORCA Output Fixtures

## Purpose

- Keep ORCA output parsing behavior fenced before replacing the old parser path.
- Cover local historical ORCA outputs plus representative external regression outputs.
- Avoid using ORCA input fixtures as a proxy for output behavior.

## Layout

- `local/`: ORCA output files already present in this repository, copied into an output-specific fixture tree.
- `cclib/`: selected ORCA output files from cclib's regression data.
- `manifest.json`: machine-readable inventory with source, feature tags, and raw text anchors.
- `cclib/LICENSE.cclib`: upstream BSD-3-Clause license text for the cclib fixtures.

## cclib Source

- Repository: <https://github.com/cclib/cclib>
- Commit: `0e5d5dfd7d15b2fbd95e02fdd36bd4fe3336a5d5`
- Upstream directory: `data/ORCA`
- License: BSD-3-Clause, copied in `cclib/LICENSE.cclib`

## Selection Rule

- Keep a compact but broad output set instead of mirroring all upstream data.
- Prefer fixtures that exercise different output sections: single point energies,
  post-HF energies, gradients, geometry optimization, frequencies, IR/Raman,
  polarizability, NMR, solvation, empirical dispersion, and excited-state methods.

## Parser Contract

- The representative ORCA output smoke test runs with the regular science regression.
- Run `make check-orca` to execute the complete corpus. The tracked `manifest.json` is required;
  a missing manifest is an error rather than a skipped test.
- Tests in `tests/test_orca_output_fixtures.py` validate file presence, raw anchors, and structured
  parser behavior for the covered feature families.
