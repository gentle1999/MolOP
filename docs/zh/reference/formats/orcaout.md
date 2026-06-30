# ORCA 输出

<!-- format-support:orcaout -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `orcaout` |
| 扩展名 | `.out`, `.log`, `.orcaout` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | 坐标和 QM 结果 |

ORCA 输出解析以用户可获得的结果属性为口径，例如能量、力、振动、布居、溶剂化、
激发态和响应性质。

| 特性 | 支持程度 | 支持范围 | 明确边界 | 测试证据 |
| ---- | -------- | -------- | -------- | -------- |
| <!-- feature-area:Software metadata, frames, and status -->Software metadata, frames, and status | fixture 覆盖 | 本地和 cclib ORCA fixtures 暴露软件名/版本、帧数和正常终止状态，并覆盖 legacy 与 modern ORCA 版本。 | 覆盖有意以 fixture 驱动；该集合外的异常或部分输出可能仍不支持。 | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_versions_cover_legacy_and_modern_orca`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Energies -->Energies | 部分支持 | 对声明了期望 final energy 的 fixtures，暴露 Hartree 单位的最终 single-point total energy。 | 只声明已测试的 ORCA 能量字段；方法专有 correction/decomposition 表可能仍是非结构化。 | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Forces and geometry optimization -->Forces and geometry optimization | 部分支持 | 对覆盖的 ORCA gradient 和 optimization 输出暴露 gradient/force 数组和几何优化状态。 | 覆盖范围跟随本地/cclib fixture 清单，不表示每个 ORCA optimizer diagnostic 都已结构化。 | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Vibrations and vibrational spectra -->Vibrations and vibrational spectra | fixture 覆盖 | 对覆盖的 frequency 输出暴露 vibration 容器；fixture 清单包含 IR 和 Raman 示例。 | Feature counts 证明 fixture 存在；只有 structured parse contract 明确断言的字段保证结构化。 | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Charge and spin populations -->Charge and spin populations | 部分支持 | 对包含 population 区段的覆盖 ORCA 输出暴露 charge/spin population 容器。 | 测试 fixture 之外的特定 population scheme 可能仍非结构化。 | `tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract` |
| <!-- feature-area:Solvation -->Solvation | 部分支持 | 对覆盖的 CPCM/SMD 溶剂化输出，在文件或帧模型上暴露 solvent 或 solvation model 信息。 | 只声明 fixture 支持的隐式溶剂化字段；详细 solvent energy decomposition 可能仍非结构化。 | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Excited states and multireference results -->Excited states and multireference results | 部分支持 | 对覆盖的 TDDFT、ADC2、EOM-CCSD、STEOM-CCSD、STEOM-DLPNO-CCSD 和 ROCIS 输出暴露 electronic state 容器；multireference result 容器属于存在时的结构化契约。 | 覆盖以 fixture 驱动，不表示每个 ORCA 激发态或多参考打印变体都已归一化。 | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
| <!-- feature-area:Response properties -->Response properties | 部分支持 | 对覆盖输出暴露结构化 polarizability；fixture 清单也包含 NMR 和 spin-spin coupling 示例。 | NMR 和 spin-spin coupling 目前是 fixture-anchored，可能尚未提升到专门公共字段。 | `tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad`<br>`tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract`<br>`tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts` |
