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
| <!-- feature-area:Reader -->Reader | 已支持 | SDF/MOL block 解析为原子、坐标、键、形式电荷、自由基、总电荷和自旋多重度字段。 | 测试契约不声明任意 SD data field 的 round-trip。 |
| <!-- feature-area:Writer graph policy -->Writer graph policy | 已支持 | SDF 文件/帧 writer 要求 graph-capable 输入，默认 strict graph 语义，并支持 RDKit/OpenBabel 渲染后端。 | 当不存在 coords-only SDF writer 时，coords-only override 会被拒绝。 |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | 已解析结构可以通过 codec registry 转换为 SDF。 | 转换质量取决于已恢复或已转换的 graph 数据。 |
