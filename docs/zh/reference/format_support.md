# 支持的格式

MolOP 支持多种化学文件格式的读取和写入。下表总结了内置的支持情况：

| 格式 ID   | 常用扩展名                             | 读取 | 写入 | 备注                                      |
| --------- | -------------------------------------- | ---- | ---- | ----------------------------------------- |
| `g16log`  | `.log`, `.out`, `.g16`, `.gal`, `.irc` | 是   | 否   | Gaussian 16 输出文件。                    |
| `xyz`     | `.xyz`                                 | 是   | 是   | 标准 XYZ 坐标文件。                       |
| `sdf`     | `.sdf`, `.sd`, `.mol`                  | 是   | 是   | MDL 结构数据文件。                        |
| `smi`     | `.smi`, `.txt`                         | 是   | 是   | SMILES 字符串。                           |
| `gjf`     | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` | 是   | 是   | Gaussian 输入文件。                       |
| `orcainp` | `.inp`                                 | 是   | 是   | ORCA 输入文件。                           |
| `cml`     | `.cml`                                 | 否 | 是   | 化学标记语言 (Chemical Markup Language)。 |
| `fakeg`   | `.fakeg`                               | 否   | 是   | 从已解析 Gaussian 数据渲染出的 Gaussian-like 输出。 |

!!! note
    writer 格式会延迟注册。安装 shell completion 后，可通过
    `molop parse PATTERN format-transform --format <TAB>` 查看当前环境可用的 writer 格式。
