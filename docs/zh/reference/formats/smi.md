# SMILES

<!-- format-support:smi -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `smi` |
| 扩展名 | `.smi`, `.txt` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 读取为坐标；写出要求 graph |

SMILES 记录读写，支持 graph 派生的电荷/自旋多重度和 canonical SMILES 写出。

```python
from molop import AutoParser

frame = AutoParser("molecules.smi")[0][0]
print(frame.to_canonical_SMILES(), frame.coords.shape)
```

??? example "输出"

    输入 `CCO ethanol` 时输出类似 `CCO (3, 3)`。

行尾名称不属于解析契约，三维坐标不会保留。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Reader -->读取 | 已支持 | 非空 SMILES 记录；每行第一个空白分隔 token；由 RDKit 派生 2D 坐标、graph、形式电荷/自由基、电荷和自旋多重度。 | 行尾名称或额外列不属于已声明解析数据契约；输入 3D 坐标无法保留。 |
| <!-- feature-area:Writer graph policy -->分子图与写出 | 已支持 | 使用 strict graph 语义写出 canonical graph-derived SMILES。 | Coords-only 渲染不是默认 writer 契约。 |
| <!-- feature-area:Registry conversion -->格式转换 | 已支持 | Gaussian 输出结构可以通过 registry 转换为 SMILES。 | 转换取决于成功的 graph 恢复或转换。 |

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[RDKit 文件接口](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Reader -->读取 | 每行第一个 SMILES 字段；RDKit `MolFromSmiles` 与 `Compute2DCoords` | 原子、键、形式电荷及自由基来自解码图；坐标是生成的二维绘图坐标，不是原文几何；尾随词不作为分子属性保存。 |
| <!-- parsing-basis:Writer graph policy -->分子图与写出 | `RDKit / Open Babel 分子图生成的规范 SMILES` | 规范化由后端算法执行，不保留输入词序；图恢复与立体化学取决于当前分子信息。 |
| <!-- parsing-basis:Registry conversion -->格式转换 | `已有分子图及目标格式写出参数` | 由图字段生成规范 SMILES，不编码 QM 能量、NMR、光谱及原始坐标。 |

<!-- parsing-basis:end -->
