# CML Writer

<!-- format-support:cml -->

| Item | Value |
| ---- | ----- |
| Format ID | `cml` |
| Extensions | `.cml` |
| Read | No |
| Write | Yes |
| Registry role | File writer, frame writer |
| Data level | Graph |

Chemical Markup Language writer backed by RDKit or OpenBabel molecule rendering.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:File and frame writer -->File and frame writer | Supported | Selected-frame file rendering, negative frame index selection, separate selected-frame blocks, and single-frame rendering through the registry. | No CML reader is registered. | `tests/test_cml_codec.py::test_cml_writer_renders_selected_frames_in_one_file`<br>`tests/test_cml_codec.py::test_cml_writer_supports_negative_frame_id`<br>`tests/test_cml_codec.py::test_cml_writer_can_return_selected_frames_as_separate_blocks`<br>`tests/test_cml_codec.py::test_cml_frame_writer_renders_single_frame` |
| <!-- feature-area:Rendering engines and validation -->Rendering engines and validation | Supported | RDKit-backed default rendering, OpenBabel-backed rendering, frame-less file rejection, and unsupported-engine diagnostics. | Requires a recoverable RDKit or OpenBabel molecule; arbitrary XML/CML source preservation is not covered. | `tests/test_cml_codec.py::test_cml_writer_supports_openbabel_engine`<br>`tests/test_cml_codec.py::test_cml_writer_requires_frames`<br>`tests/test_cml_codec.py::test_cml_frame_writer_rejects_unknown_engine` |
