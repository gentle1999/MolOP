# SDF/MOL

<!-- format-support:sdf -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `sdf` |
| 扩展名 | `.sdf`, `.sd`, `.mol` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、文件 writer、帧 writer |
| 数据层级 | 读取为坐标；写出要求 graph |

基于 RDKit graph 提取和 strict graph writer 语义的 SDF/MOL 结构读写。

```python
from molop import AutoParser

frame = AutoParser("molecule.sdf")[0][0]
print(frame.rdmol.GetNumAtoms(), frame.rdmol.GetNumBonds())
```

??? example "输出"

    输出原子数和键数，例如水分子为 `3 2`。

任意 SD data field 不保证无损 round-trip。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Reader -->读取 | 已支持 | SDF/MOL block 解析为原子、坐标、键、形式电荷、自由基、总电荷和自旋多重度字段。 | 任意 SD 属性不保证往返保留。 |
| <!-- feature-area:Writer graph policy -->分子图与写出 | 已支持 | SDF 文件/帧 writer 要求 graph-capable 输入，默认 strict graph 语义，并支持 RDKit/OpenBabel 渲染后端。 | 当不存在 coords-only SDF writer 时，coords-only override 会被拒绝。 |
| <!-- feature-area:Registry conversion -->格式转换 | 已支持 | 已解析结构可以通过 codec registry 转换为 SDF。 | 转换质量取决于已恢复或已转换的 graph 数据。 |

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[RDKit 文件接口](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Reader -->读取 | MOL 计数、原子与键区块；`M CHG`、`M RAD`；`$$$$` 分隔符 | 提取原子、Å 坐标、键、形式电荷/自由基及标题；总电荷求和、多重度=自由基数+1 属于 MolOP 投影，不是自旋观测值；不映射任意 SD 属性。 |
| <!-- parsing-basis:Writer graph policy -->分子图与写出 | 已有或恢复的分子图；RDKit `MolToMolBlock` / Open Babel SDF 写出接口 | 从带图结构的分子生成，原子/键编码由 RDKit/Open Babel 后端决定；仅坐标输入的拓扑恢复是另一项 MolOP 操作。 |
| <!-- parsing-basis:Registry conversion -->格式转换 | `选定帧的分子图、坐标与写出参数` | 转换依据现有图、坐标及写出参数，不产生新的 QM 属性观测依据。 |

<!-- parsing-basis:end -->
