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

下载[共享 ORCA 样例](../../../assets/examples/water_mp2.out)后运行：

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
cml_text = frame.format_transform("cml")
print(cml_text.splitlines()[0])
```

??? example "输出"

    ```text
    <?xml version="1.0" encoding="utf-8"?>
    ```

当前没有 CML reader。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:File and frame writer -->File and frame writer | 已支持 | 选定帧文件渲染、负数 frame index 选择、选定帧拆分为独立 block，以及通过 registry 渲染单帧。 | 没有注册 CML reader。 |
| <!-- feature-area:Rendering engines and validation -->Rendering engines and validation | 已支持 | 默认 RDKit 后端渲染、OpenBabel 后端渲染、无 frames 文件模型校验，以及不支持 engine 的诊断。 | 要求能够恢复 RDKit 或 OpenBabel 分子；不覆盖任意 XML/CML 源文本保留。 |
