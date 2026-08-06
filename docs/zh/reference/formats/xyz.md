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
| <!-- feature-area:Reader -->Reader | 已支持 | 标准多帧 XYZ；comment 中的电荷和自旋多重度；文件级电荷/自旋多重度从第一帧 finalize。 | 不恢复分子图；非法原子数或不完整帧在解析阶段按格式不匹配处理。 |
| <!-- feature-area:Writer and conversion -->Writer and conversion | 已支持 | 文件/帧 XYZ 渲染、显式 comment 覆盖，以及从带坐标的已解析文件通过 registry 转换为 XYZ。 | Writer 使用坐标语义；XYZ 不编码 graph-only 元数据。 |
| <!-- feature-area:Format mismatch handling -->Format mismatch handling | 已支持 | 简单格式的不匹配检测延迟到正常解析阶段，避免额外 IO probe。 | XYZ 不要求独立文件级指纹。 |
