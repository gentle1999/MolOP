# XYZ

<!-- format-support:xyz -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `xyz` |
| 扩展名 | `.xyz` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 坐标 |

标准 XYZ 坐标读写，支持从 comment 行读取电荷和自旋多重度。

```python
from molop import AutoParser

frame = AutoParser("molecule.xyz")[0][0]
print(len(frame.atoms), frame.coords.shape, frame.charge, frame.multiplicity)
```

??? example "输出"

    一个 3 原子中性单重态输出为 `3 (3, 3) 0 1`。

XYZ 不保存能量或完整分子图。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Reader -->读取 | 已支持 | 标准多帧 XYZ；comment 中的电荷和自旋多重度；文件级电荷与多重度采用第一帧的值。 | 不恢复分子图；非法原子数或不完整帧在解析阶段按格式不匹配处理。 |
| <!-- feature-area:Writer and conversion -->写出与转换 | 已支持 | 文件/帧 XYZ 渲染、显式 comment 覆盖，以及从带坐标的已解析文件通过 registry 转换为 XYZ。 | Writer 使用坐标语义；XYZ 不编码 graph-only 元数据。 |
| <!-- feature-area:Format mismatch handling -->帧边界与无效输入 | 已支持 | 读取时检查原子数、坐标列和帧是否完整。 | 不完整的帧会报告格式不匹配。 |

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[XYZ 格式约定](https://openbabel.org/docs/FileFormats/XYZ_cartesian_coordinates_format.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Reader -->读取 | 原子数行；注释行；`element X Y Z` 坐标行 | `atoms`、`coords`（Å）、`comment`；注释中的电荷/多重度是 MolOP 扩展，缺省为 0/1，不据此声称观测到键图。 |
| <!-- parsing-basis:Writer and conversion -->写出与转换 | 选定帧的 `atoms`、`coords` 与可选 `comment` | 由原子数、注释及 Å 坐标生成；显式/存储/自动注释的优先级由 MolOP 定义，不编码仅属于图结构的属性。 |
| <!-- parsing-basis:Format mismatch handling -->帧边界与无效输入 | `正整数原子数 N 与后续 N 行坐标` | 由原子数确定帧边界；不完整帧在解析时报告格式不匹配，不额外猜测文件类型。 |

<!-- parsing-basis:end -->
