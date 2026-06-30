# ORCA 输入文件模型与解析指南

本文说明 MolOP 中 ORCA `.inp` 的结构化读取模型。ORCA writer 暂不注册；不要依赖 `orca_raw_preamble` / `orca_raw_postamble` 或 raw round-trip 兼容。

## 数据层次

ORCA 的 `$new_job` 在文件解析层切分为多个 `ORCAInpFileFrame`，与 Gaussian `--Link1--` 分帧思路一致：

```text
ORCAInpFile
  frames: list[ORCAInpFileFrame]

ORCAInpFileFrame
  comment_lines
  keyword_lines
  blocks
  geometry
  trailing_lines
  auxiliary_basis_set
  dispersion_correction
  has_mixed_basis
  output_print_settings
  atoms / coords / charge / multiplicity
```

`ORCAInpFileFrame` 必须直接承载坐标和主语义字段。`geometry` 是 frame 的结构化字段来源，而不是额外的 job 包装层。

## Frame 字段

`comment_lines`

: 保存顶层 `#` 注释，去掉前导 `#`。

`keyword_lines`

: 保存顶层 `!` 关键词行，去掉前导 `!`。frame validator 会把这些行合并到通用 `keywords` 字段，并尽量派生 `method`、`functional`、`basis_set`。ORCA 辅助基组关键词会写入 `auxiliary_basis_set`，例如 `def2-SVP/C`。ORCA 色散关键词也会写入 `dispersion_correction`，例如 `D3BJ` 或 `D4`；当识别到 DFT 泛函时，通用 `functional` 字段会附加色散后缀，与 Gaussian 行为一致，例如 `B3LYP-D3BJ`、`PBE0-D4`。

`blocks`

: 保存顶层 `%...` 块。每个 `ORCABlock` 包含：

- `name`：去掉 `%` 后的小写块名，例如 `pal`、`maxcore`、`scf`
- `raw_header`：原始 header 行
- `lines`：块体行，去掉行尾换行

`geometry`

: 保存 ORCA 几何输入结构。支持 `* xyz`、`* cart`、`* xyzfile`、`* gzmtfile` 和 `%coords` 中的 Cartesian 坐标。Cartesian 坐标会同步到 frame 的 `atoms` / `coords`。

`trailing_lines`

: 保存暂未归类的顶层非空文本，供调试和后续扩展使用。它不是 raw 兼容接口。

`dispersion_correction`

: 从 ORCA `!` 关键词行识别出的色散校正。该字段是 ORCA 输入 frame 的结构化字段，不是 `BaseQMInputFrame` 的通用字段。当同一 frame 识别到 DFT 泛函时，色散会同时作为 `FUNCTIONAL-DISPERSION` 后缀补充到通用 `functional` 字符串。

`has_mixed_basis`

: 当任一原子行带有 `newgto` 或 `newauxgto` 时为 `true`。

`output_print_settings`

: 从 `%output` 块中识别出的 `Print[ ... ] = ...` 设置。原始 `%output` 块仍保存在 `blocks` 和 `resources_raw`。

## Geometry 模型

`ORCAGeometry` 字段：

- `ctype`：`xyz`、`cart`、`cartesian`、`int`、`internal`、`gzmt`、`xyzfile` 或 `gzmtfile`
- `charge` / `multiplicity`
- `units`：例如 `bohr`
- `external_path`：外部几何文件路径
- `atoms`：Cartesian 原子列表
- `point_charges`：`Q q x y z` 点电荷
- `source`：`star` 或 `percent_coords`

`ORCAGeometryAtom` 字段：

- `symbol` / `atomic_number`
- `x` / `y` / `z`
- `is_dummy` / `is_ghost`
- `fragment_id`
- `frozen`
- `isotope`
- `nuclear_charge`
- `basis_set` / `auxiliary_basis_set`：来自原子行上的 `newgto` / `newauxgto`
- `basis_overrides`：原子级基组覆盖指令列表，保留指令类型、命名基组和 token

## 解析流程

文件解析器：

1. 读取 `.inp` 文件文本。
2. 在顶层 `$new_job` 处分帧。
3. 对每个 frame 调用 `ORCAInpFileFrameParser`。

frame 解析器：

1. 识别 `* ...` 几何段，或 `%coords ... end` 几何块。
2. 扫描顶层 `%` 块，生成 `ORCABlock`。
3. 过滤已被几何段和 `%` 块占用的行。
4. 从剩余行提取 `#` 注释和 `!` 关键词。
5. 组装 `ORCAInpFileFrame` 顶层字段。
6. frame validator 派生通用 QM input 字段。

## 单位与坐标

Cartesian 输入默认按 Å 处理。`units` 为 `bohr`、`a0`、`au` 等时，解析器把坐标和点电荷位置转换为 Å 后写入模型。

`xyzfile` / `gzmtfile` 只保留 `external_path`，不会读取外部文件，也不会伪造 `atoms` / `coords`。

## 开发约束

- 不新增 `jobs` 包装层；`$new_job` 对应多个 frame。
- 不恢复 `orca_raw_preamble` / `orca_raw_postamble`。
- 不注册 ORCA writer，直到有结构化 renderer。
- 程序特有语法优先写入 ORCA 结构化字段和公共 QM 结构化容器，再由 `project_common_qm_fields()` 投影到 `method`、`functional`、`basis_set` 等兼容字段。
- 新语法应先扩展结构化模型，再扩展 parser；不要把语义藏进 raw 字符串。
