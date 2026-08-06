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
| <!-- feature-area:Reader -->Reader | 已支持 | 非空 SMILES 记录；每行第一个空白分隔 token；由 RDKit 派生 2D 坐标、graph、形式电荷/自由基、电荷和自旋多重度。 | 行尾名称或额外列不属于已声明解析数据契约；输入 3D 坐标无法保留。 |
| <!-- feature-area:Writer graph policy -->Writer graph policy | 已支持 | 使用 strict graph 语义写出 canonical graph-derived SMILES。 | Coords-only 渲染不是默认 writer 契约。 |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | Gaussian 输出结构可以通过 registry 转换为 SMILES。 | 转换取决于成功的 graph 恢复或转换。 |
