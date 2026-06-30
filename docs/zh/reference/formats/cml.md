# CML Writer

<!-- format-support:cml -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `cml` |
| 扩展名 | `.cml` |
| 读取 | 否 |
| 写入 | 是 |
| Registry 角色 | 文件 writer、帧 writer |
| 数据层级 | Graph |

Chemical Markup Language writer，基于 RDKit 或 OpenBabel 分子渲染。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:File and frame writer -->File and frame writer | 已支持 | 选定帧文件渲染、负数 frame index 选择、选定帧拆分为独立 block，以及通过 registry 渲染单帧。 | 没有注册 CML reader。 | `tests/test_cml_codec.py::test_cml_writer_renders_selected_frames_in_one_file`<br>`tests/test_cml_codec.py::test_cml_writer_supports_negative_frame_id`<br>`tests/test_cml_codec.py::test_cml_writer_can_return_selected_frames_as_separate_blocks`<br>`tests/test_cml_codec.py::test_cml_frame_writer_renders_single_frame` |
| <!-- feature-area:Rendering engines and validation -->Rendering engines and validation | 已支持 | 默认 RDKit 后端渲染、OpenBabel 后端渲染、无 frames 文件模型校验，以及不支持 engine 的诊断。 | 要求能够恢复 RDKit 或 OpenBabel 分子；不覆盖任意 XML/CML 源文本保留。 | `tests/test_cml_codec.py::test_cml_writer_supports_openbabel_engine`<br>`tests/test_cml_codec.py::test_cml_writer_requires_frames`<br>`tests/test_cml_codec.py::test_cml_frame_writer_rejects_unknown_engine` |
