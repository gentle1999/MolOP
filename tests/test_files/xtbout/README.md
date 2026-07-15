# xTB Output Fixtures

This directory contains xTB command-line standard-output captures used by the
`xtbout` parser tests.

- Existing `6-*`, `xtb_6_*`, and `dsgdb*` files are retained project fixtures
  spanning xTB 6.1 beta through 6.6.1 and several calculation modes.
- `xtb_5_8_1_legacy_contract.out` is a compact parser-contract fixture for the
  xTB 5 legacy print family. It is intentionally limited to stable labels also
  present in the legacy xTB 6.1 output family: version banner, calculation
  setup, Bohr coordinates, SCC convergence, orbital table, total energy,
  optimized structure, and termination timing.

The xTB 5 fixture is not presented as a verbatim vendor regression file. Add a
complete real xTB 5 output here when one is available, preserving its source and
calculation command in this README.
