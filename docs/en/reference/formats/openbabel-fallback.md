# OpenBabel Fallback

<!-- format-support:openbabel-fallback -->

| Item | Value |
| ---- | ----- |
| Format ID | `openbabel-fallback` in the support matrix; runtime reader reports `openbabel` |
| Extensions | Unknown extension path, then OpenBabel-supported candidate formats |
| Read | Yes |
| Write | No |
| Registry role | Special fallback reader |
| Data level | Coordinates |

Fallback reader for unknown extensions when OpenBabel can parse the source file.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Unknown-extension fallback -->Unknown-extension fallback | Partial | Unknown extensions can be read through OpenBabel-compatible formats and converted to an XYZFile/XYZFileFrame model with detected format ID persistence and first-molecule coordinate preservation. | Only the first parsed molecule is converted; output is coordinate-level and OpenBabel availability controls format reach. | `tests/test_structure_format_features.py::test_openbabel_unknown_extension_fallback_keeps_first_molecule_coordinates`<br>`tests/test_io_convert_registry_smoke.py::test_unknown_extension_fallback_smoke`<br>`tests/test_filebatchmodeldisk_codec_filter.py::test_filter_by_codec_id_uses_detected_format_id` |
