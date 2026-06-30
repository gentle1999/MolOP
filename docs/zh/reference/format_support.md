# 格式支持

<!-- format-support-overview -->

本页只保留 MolOP 内置格式支持概览。每个具体格式的支持范围见独立页面。可机器检查的
事实源是 `tests/format_feature_coverage/support_matrix.py`；
`tests/format_feature_coverage/test_support_matrix.py` 会校验该矩阵与已注册 codec
一致，并且每个特性都声明支持范围、限制和真实存在的测试节点。

| 分组 | 格式 | 概要 |
| ---- | ---- | ---- |
| 结构格式 | [`xyz`](formats/xyz.md), [`sdf`](formats/sdf.md), [`smi`](formats/smi.md), [`cml`](formats/cml.md) | 坐标/分子图结构读取与渲染。 |
| QM 输入格式 | [`gjf`](formats/gjf.md), [`orcainp`](formats/orcainp.md) | Gaussian 和 ORCA 输入解析；目前只有 Gaussian 输入支持渲染。 |
| QM 输出格式 | [`g16log`](formats/g16log.md), [`orcaout`](formats/orcaout.md), [`fakeg`](formats/fakeg.md) | Gaussian/ORCA 输出解析，以及从 Gaussian 输出数据渲染 Gaussian-like 文本。 |
| 特殊 reader | [OpenBabel fallback](formats/openbabel-fallback.md) | 通过 OpenBabel 兼容格式 fallback 读取未知扩展名文件。 |

| 格式 ID | 常用扩展名 | 读取 | 写入 | 支持详情 |
| ------- | ---------- | ---- | ---- | -------- |
| `xyz` | `.xyz` | 是 | 是 | [XYZ](formats/xyz.md) |
| `sdf` | `.sdf`, `.sd`, `.mol` | 是 | 是 | [SDF/MOL](formats/sdf.md) |
| `smi` | `.smi`, `.txt` | 是 | 是 | [SMILES](formats/smi.md) |
| `gjf` | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` | 是 | 是 | [Gaussian 输入](formats/gjf.md) |
| `g16log` | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` | 是 | 否 | [Gaussian 输出](formats/g16log.md) |
| `orcainp` | `.inp` | 是 | 否 | [ORCA 输入](formats/orcainp.md) |
| `orcaout` | `.out`, `.log`, `.orcaout` | 是 | 否 | [ORCA 输出](formats/orcaout.md) |
| `cml` | `.cml` | 否 | 是 | [CML writer](formats/cml.md) |
| `fakeg` | `.fakeg` | 否 | 是 | [Gaussian-like renderer](formats/fakeg.md) |
| OpenBabel fallback | 未知扩展名或 OpenBabel 支持的扩展名 | 是 | 否 | [OpenBabel fallback](formats/openbabel-fallback.md) |

!!! note
    扩展名并不总是唯一对应格式。例如 `.out` 和 `.log` 可能是 ORCA 输出，也可能是
    Gaussian 输出；`.gau` 也可能是 Gaussian 输入或输出。Reader 选择会先给出扩展名候选，
    再执行内容检查；不匹配的 reader 必须抛出 `FormatMismatchError`，以便继续尝试下一个候选。
