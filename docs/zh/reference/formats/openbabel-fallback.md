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
| <!-- feature-area:Unknown-extension fallback -->Unknown-extension fallback | 部分支持 | 未知扩展名可以通过 OpenBabel 兼容格式读取，并转换为带 detected format ID、且保留首个分子坐标的 `XYZFile`/`XYZFileFrame` 模型。 | 只转换第一个成功解析的分子；输出是坐标层级；实际格式范围取决于 OpenBabel 安装。 |
