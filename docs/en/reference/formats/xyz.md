# XYZ

<!-- format-support:xyz -->

| Item | Value |
| ---- | ----- |
| Format ID | `xyz` |
| Extensions | `.xyz` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates |

Standard XYZ coordinate IO with charge/multiplicity comment support.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Reader -->Reader | Supported | Standard multi-frame XYZ files; charge and multiplicity in comments; file-level charge/multiplicity finalized from the first frame. | No graph recovery; malformed atom counts or incomplete frames fail during parse-phase mismatch handling. | `tests/test_parsers_smoke.py::test_autoparser_xyz_smoke`<br>`tests/test_io_registry_and_batch_more.py::test_base_file_parser_finalizes_file_charge_and_multiplicity_from_first_frame` |
| <!-- feature-area:Writer and conversion -->Writer and conversion | Supported | File/frame XYZ rendering, explicit comment override, and registry conversion from coordinate-bearing parsed files. | Writer uses coordinate semantics; graph-only metadata is not encoded by XYZ. | `tests/test_structure_format_features.py::test_xyz_frame_writer_uses_explicit_comment_override`<br>`tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke`<br>`tests/test_format_transform_output_dir.py::test_format_transform_batch_output_dir_writes_under_requested_dir` |
| <!-- feature-area:Format mismatch handling -->Format mismatch handling | Supported | Simple-format mismatch detection is deferred to normal parsing to avoid a second IO probe. | No separate file-level fingerprint is required for XYZ. | `tests/test_io_registry_and_batch_more.py::test_simple_coordinate_readers_defer_mismatch_to_parse_phase` |
