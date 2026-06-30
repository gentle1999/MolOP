# Gaussian-Like Renderer

<!-- format-support:fakeg -->

| Item | Value |
| ---- | ----- |
| Format ID | `fakeg` |
| Extensions | `.fakeg` |
| Read | No |
| Write | Yes |
| Registry role | File writer |
| Data level | Parsed Gaussian output data |

Gaussian-like text rendering from parsed Gaussian output data.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Rendered Gaussian-like output sections -->Rendered Gaussian-like output sections | Partial | Renders Gaussian-like structure, charge/multiplicity, frequency, and thermochemistry sections from parsed Gaussian output data. | Compatibility renderer only; not a byte-for-byte Gaussian log reproduction. | `tests/test_g16log_parser_v3_render.py::test_g16log_v3_render_contains_expected_gaussian_sections`<br>`tests/test_g16log_parser_v3_render.py::test_g16log_file_model_can_render_fakeg_for_all_frames` |
| <!-- feature-area:Frequency and thermochemistry reparse -->Frequency and thermochemistry reparse | Supported | Rendered frequency and thermochemistry content can be parsed back into frame fields. | Frame-level fakeg rendering is intentionally unsupported. | `tests/test_g16log_parser_v3_render.py::test_g16log_file_format_transform_registers_fakeg_file_renderer_only`<br>`tests/test_g16log_parser_v3_render.py::test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry` |
