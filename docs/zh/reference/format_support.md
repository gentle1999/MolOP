# 格式支持

<!-- format-support-overview -->

本页概览 MolOP 内置格式支持。每个具体格式的独立页面会说明可读取或写出的信息、
支持范围和已知边界；这些覆盖声明会在测试中与已注册 codec 保持同步。

=== "按工作流选择"

    - 结构文件从 [`xyz`](formats/xyz.md)、[`sdf`](formats/sdf.md) 或 [`smi`](formats/smi.md) 开始。
    - 计算输入使用 [`gjf`](formats/gjf.md) 或 [`orcainp`](formats/orcainp.md)。
    - 计算结果使用 [`g16log`](formats/g16log.md)、[`g16fchk`](formats/g16fchk.md)、
      [`orcaout`](formats/orcaout.md) 或 [`xtbout`](formats/xtbout.md)。
    - 需要从已解析数据生成 Gaussian-like 文本时使用 [`fakeg`](formats/fakeg.md)。

=== "按扩展名选择"

    扩展名只用于提供候选 reader。`.out` 和 `.log` 可能属于不同 QM 输出格式，最终 reader
    仍由内容探测决定。完整的扩展名到 codec 映射见下表。

| 分组 | 格式 | 概要 |
| ---- | ---- | ---- |
| 结构格式 | [`xyz`](formats/xyz.md), [`sdf`](formats/sdf.md), [`smi`](formats/smi.md), [`cml`](formats/cml.md) | 坐标/分子图结构读取与渲染。 |
| QM 输入格式 | [`gjf`](formats/gjf.md), [`orcainp`](formats/orcainp.md) | Gaussian 和 ORCA 输入解析与渲染。 |
| QM 输出格式 | [`g16log`](formats/g16log.md), [`g16fchk`](formats/g16fchk.md), [`orcaout`](formats/orcaout.md), [`xtbout`](formats/xtbout.md), [`fakeg`](formats/fakeg.md) | Gaussian log/fchk、ORCA 和 xTB 输出解析，以及从 Gaussian 输出数据渲染 Gaussian-like 文本。 |
| 特殊 reader | [OpenBabel fallback](formats/openbabel-fallback.md) | 通过 OpenBabel 兼容格式 fallback 读取未知扩展名文件。 |

| 格式 ID | 常用扩展名 | 读取 | 写入 | 支持详情 |
| ------- | ---------- | ---- | ---- | -------- |
| `xyz` | `.xyz` | 是 | 是 | [XYZ](formats/xyz.md) |
| `sdf` | `.sdf`, `.sd`, `.mol` | 是 | 是 | [SDF/MOL](formats/sdf.md) |
| `smi` | `.smi`, `.txt` | 是 | 是 | [SMILES](formats/smi.md) |
| `gjf` | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` | 是 | 是 | [Gaussian 输入](formats/gjf.md) |
| `g16log` | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` | 是 | 否 | [Gaussian 输出](formats/g16log.md) |
| `g16fchk` | `.fchk`, `.fch`, `.fck` | 是 | 否 | [Gaussian 格式化检查点](formats/g16fchk.md) |
| `orcainp` | `.inp` | 是 | 是 | [ORCA 输入](formats/orcainp.md) |
| `orcaout` | `.out`, `.log`, `.orcaout` | 是 | 否 | [ORCA 输出](formats/orcaout.md) |
| `xtbout` | `.out`, `.log`, `.xtbout` | 是 | 否 | [xTB 输出](formats/xtbout.md) |
| `cml` | `.cml` | 否 | 是 | [CML writer](formats/cml.md) |
| `fakeg` | `.fakeg` | 否 | 是 | [Gaussian-like renderer](formats/fakeg.md) |
| OpenBabel fallback | 未知扩展名或 OpenBabel 支持的扩展名 | 是 | 否 | [OpenBabel fallback](formats/openbabel-fallback.md) |

!!! note
    扩展名并不总是唯一对应格式。例如 `.out` 和 `.log` 可能是 ORCA、xTB 或
    Gaussian 输出；`.gau` 也可能是 Gaussian 输入或输出。Reader 选择会先给出扩展名候选，
    再执行内容检查；不匹配的 reader 必须抛出 `FormatMismatchError`，以便继续尝试下一个候选。
