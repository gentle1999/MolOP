# ORCA Output

<!-- format-support:orcaout -->

| Item | Value |
| ---- | ----- |
| Format ID | `orcaout` |
| Extensions | `.out`, `.log`, `.orcaout` |
| Read | Yes |
| Write | No |
| Registry role | Reader |
| Data level | Coordinates and QM results |

ORCA output parsing for user-facing result attributes such as energies, forces,
vibrations, populations, solvation, excited states, and response properties.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Software metadata, frames, and status -->Software metadata, frames, and status | Fixture-covered | Local and cclib ORCA fixtures expose software name/version, frame counts, and normal termination status across legacy and modern ORCA versions. | Coverage is intentionally fixture-driven; malformed or partial outputs outside this set may remain unsupported. | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_versions_cover_legacy_and_modern_orca`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Energies -->Energies | Partial | Final single-point total energies are exposed in Hartree for fixtures that declare expected final energies. | Only tested ORCA energy fields are advertised; method-specific correction/decomposition tables may remain unstructured. | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Forces and geometry optimization -->Forces and geometry optimization | Partial | Gradient/force arrays and geometry optimization status are exposed for covered ORCA gradient and optimization outputs. | Coverage follows the local/cclib fixture inventory and does not imply every ORCA optimizer diagnostic is structured. | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Vibrations and vibrational spectra -->Vibrations and vibrational spectra | Fixture-covered | Vibration containers are exposed for covered frequency outputs; fixture inventory includes IR and Raman examples. | Feature counts prove fixture presence; only fields asserted in the structured parse contract are guaranteed structured. | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Charge and spin populations -->Charge and spin populations | Partial | Charge/spin population containers are exposed for covered ORCA outputs that include population sections. | Specific population schemes outside the tested fixtures may remain unstructured. | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Solvation -->Solvation | Partial | Solvent or solvation model information is exposed on file or frame models for covered CPCM/SMD solvation outputs. | Only fixture-backed implicit-solvation fields are advertised; detailed solvent energy decompositions may remain unstructured. | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Excited states and multireference results -->Excited states and multireference results | Partial | Electronic state containers are exposed for covered TDDFT, ADC2, EOM-CCSD, STEOM-CCSD, STEOM-DLPNO-CCSD, and ROCIS outputs; multireference result containers are part of the structured contract where present. | Coverage is fixture-driven and does not imply every ORCA excited-state or multireference print variant is normalized. | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Response properties -->Response properties | Partial | Polarizability is exposed as a structured field for covered outputs; fixture inventory also includes NMR and spin-spin coupling examples. | NMR and spin-spin coupling are currently fixture-anchored and may not yet be lifted into dedicated common fields. | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
