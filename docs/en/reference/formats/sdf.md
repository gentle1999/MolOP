# SDF/MOL

<!-- format-support:sdf -->

| Item | Value |
| ---- | ----- |
| Format ID | `sdf` |
| Extensions | `.sdf`, `.sd`, `.mol` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates on read; graph required on write |

SDF/MOL structure IO backed by RDKit graph extraction and strict graph writer semantics.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Reader -->Reader | Supported | SDF/MOL blocks parsed into atom, coordinate, bond, formal charge, radical, total charge, and multiplicity fields. | Arbitrary SD data-field round-trip is not advertised by the tested contract. | `tests/test_autoparser_sdf_tmpfile.py::test_autoparser_sdf_tmpfile` |
| <!-- feature-area:Writer graph policy -->Writer graph policy | Supported | SDF file/frame writer requires graph-capable input, defaults to strict graph semantics, and supports RDKit/OpenBabel rendering engines. | Coords-only override is rejected when no coords-only SDF writer exists. | `tests/test_structure_format_features.py::test_sdf_frame_writer_supports_rdkit_openbabel_engines_and_rejects_unknown_engine`<br>`tests/test_io_registry_and_batch_more.py::test_registry_builtin_default_graph_policies_match_format_semantics`<br>`tests/test_io_registry_and_batch_more.py::test_registry_raises_when_coords_override_has_no_coords_writer` |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Parsed structures can be converted to SDF through the codec registry. | Conversion quality depends on recovered or transformed graph data. | `tests/test_io_convert_registry_smoke.py::test_registry_convert_xyz_smoke` |
