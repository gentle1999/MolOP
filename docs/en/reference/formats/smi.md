# SMILES

<!-- format-support:smi -->

| Item | Value |
| ---- | ----- |
| Format ID | `smi` |
| Extensions | `.smi`, `.txt` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates on read; graph required on write |

SMILES record IO with graph-derived charge/multiplicity and canonical SMILES writing.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Reader -->Reader | Supported | Non-empty SMILES records; first whitespace token parsed; 2D coordinates, graph, formal charge/radical, charge, and multiplicity derived by RDKit. | Trailing record names or columns are not part of the advertised parsed-data contract; input 3D coordinates cannot be preserved. | `tests/test_autoparser_smi_tmpfile.py::test_autoparser_smi_tmpfile`<br>`tests/test_structure_format_features.py::test_smi_reader_uses_first_token_and_writer_canonicalizes` |
| <!-- feature-area:Writer graph policy -->Writer graph policy | Supported | Canonical graph-derived SMILES writer with strict graph semantics. | Coords-only rendering is not the default writer contract. | `tests/test_structure_format_features.py::test_smi_reader_uses_first_token_and_writer_canonicalizes`<br>`tests/test_io_registry_and_batch_more.py::test_registry_builtin_default_graph_policies_match_format_semantics`<br>`tests/test_io_convert_registry_smoke.py::test_registry_convert_xyz_smoke` |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Gaussian output structures can be converted to SMILES through the registry. | Conversion depends on successful graph recovery or transformation. | `tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke` |
