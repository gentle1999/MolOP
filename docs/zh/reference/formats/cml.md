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
| <!-- feature-area:File and frame writer -->文件与帧写出 | 已支持 | 选定帧文件渲染、负数 frame index 选择、选定帧拆分为独立 block，以及通过 registry 渲染单帧。 | 没有注册 CML reader。 |
| <!-- feature-area:Rendering engines and validation -->渲染后端与格式差异 | 已支持 | 默认 RDKit 后端渲染、OpenBabel 后端渲染、无 frames 文件模型校验，以及不支持 engine 的诊断。 | 要求能够恢复 RDKit 或 OpenBabel 分子；不覆盖任意 XML/CML 源文本保留。 |

## 数据来源与写出约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[RDKit 文件接口](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html); [Open Babel 格式说明](https://openbabel.org/docs/FileFormats/Overview.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:File and frame writer -->文件与帧写出 | `选定帧的原子、键、电荷和坐标；所选后端的 XML 写出接口` | 仅写出：图和坐标来自已有或恢复后的模型，不是读取观测；帧选择及拼接由写出参数确定，未注册 CML 读取器。 |
| <!-- parsing-basis:Rendering engines and validation -->渲染后端与格式差异 | `RDKit MolToMrvBlock`; `Open Babel omol.write("cml")` | 明确实际后端：当前 RDKit 调用 MolToMrvBlock 生成 Marvin XML，Open Babel 调用 cml；不将两者声明为相同 XML 模式，；写出前需提供分子图及有效的 engine 参数。 |

<!-- parsing-basis:end -->
