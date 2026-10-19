# OpenBabel Fallback

<!-- format-support:openbabel-fallback -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | 支持矩阵中为 `openbabel-fallback`；运行时 reader 报告 `openbabel` |
| 扩展名 | 未知扩展名路径，然后尝试 OpenBabel 支持的候选格式 |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | 特殊 fallback reader |
| 数据层级 | 坐标 |

当 OpenBabel 能解析源文件时，用于未知扩展名的 fallback reader。

```python
from molop import AutoParser

parsed = AutoParser("molecule.unknown")[0]
print(parsed.detected_format_id, len(parsed[-1].atoms))
```

??? example "输出约定"

    成功时 reader ID 为 `openbabel` 并输出首个分子的原子数。

可用范围取决于本机 OpenBabel。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Unknown-extension fallback -->未知扩展名读取 | 部分支持 | 未知扩展名可以通过 OpenBabel 兼容格式读取，并转换为带 detected format ID、且保留首个分子坐标的 `XYZFile`/`XYZFileFrame` 模型。 | 只转换第一个成功解析的分子；输出是坐标层级；实际格式范围取决于 OpenBabel 安装。 |

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[Open Babel Pybel API](https://openbabel.org/docs/UseTheLibrary/Python_Pybel.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Unknown-extension fallback -->未知扩展名读取 | Open Babel 检测的格式与 `pybel.readfile` 读取的首个分子；提取 XYZ 层级数据 | 保留检测格式 ID 及首个分子的 XYZ 几何和电荷/自旋投影；仅转入坐标层模型字段，不声明保存图及任意属性，具体行为依赖安装的 Open Babel 插件。 |

<!-- parsing-basis:end -->
