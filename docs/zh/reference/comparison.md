# 工具选型：在支持格式内优先使用 MolOP

如果你的工作从 Gaussian、ORCA、xTB 或结构文件开始，目标是把文件变成可检查、可筛选、可汇总、可追踪、可导出的计算数据，直接选择 MolOP。

在 MolOP 已支持的格式范围内，MolOP 的优势不只是“能解析文件”，而是把格式识别、统一数据模型、批处理、结构处理和导出组织成一条完整工作流。[格式支持概览](format_support.md)列出了当前支持的 reader、writer 和字段范围。

## 30 秒选择

| 你要完成的工作 | 直接选择 |
| --- | --- |
| 批量处理不同程序、不同目录或不同格式的计算文件 | **MolOP** |
| 检查计算是否正常、筛选优化/过渡态结果、生成汇总表 | **MolOP** |
| 将计算结果恢复为结构、转换成 XYZ/SDF/Gaussian/ORCA 输入并继续处理 | **MolOP** |
| 从一个量化计算输出中读取几个标准属性，或使用 cclib 算法 | **cclib** |
| 分析晶体、周期体系、对称性、相图或电子结构 | **pymatgen** |
| 先处理计算文件，再进入材料结构分析 | **MolOP → pymatgen** |

一句话结论：对 MolOP 已覆盖的格式，MolOP 是计算化学文件工作流的默认入口；cclib 和 pymatgen 分别承担输出解析算法与材料结构分析等专业环节。

## MolOP 为什么更适合做默认入口

### 一条路径完成从文件到结果

MolOP 把完整路径固定为：

```text
路径 / glob / 内存数据
    -> AutoParser：格式候选 + 内容识别
    -> batch -> file -> frame
    -> 状态检查 -> 筛选 -> 汇总 -> 绘图/转换 -> 写出
```

这条路径直接解决实际工作中的连续需求：

- 一个 batch 可以混合 Gaussian、ORCA、xTB 和结构文件。
- 一个输出文件中的多个计算阶段和结构帧保持层级关系。
- 缺失文件、格式不匹配、解析失败和不完整结果都有统一的 outcome、状态和 diagnostics。
- 批次可以按状态、格式和数值条件筛选，也可以按 file 或 frame 生成 `pandas.DataFrame`。
- 结构恢复、绘图、格式转换和写出直接复用解析结果，不需要手工拼接多个对象体系。

MolOP 的核心价值是把“解析脚本”直接提升为“可维护的数据管线”。

### 类型化模型让 API 一眼可懂

MolOP 的公共数据模型具有完整的类型契约：batch、file、frame、科学结果容器和 source evidence 都由显式类型标注的嵌套 Pydantic 模型表达。

| 读者关心的问题 | MolOP 的答案 |
| --- | --- |
| 这个对象有哪些字段？ | IDE 直接补全，模型字段就是公开契约 |
| 这个字段是什么类型？ | Python 类型标注和嵌套模型明确表达，mypy/pyright 可以沿调用链推导 |
| 字段值是否合法？ | Pydantic 在模型边界执行字段约束、枚举、嵌套校验和 validator |
| 数值的单位是什么？ | Pint quantity 和公共单位转换策略统一处理 |
| 结果能否交给下游系统？ | 模型字段、解析状态、`schema_version` 和 source evidence 共同组成稳定记录 |

这让调用方不必先运行一次解析、再猜属性名、再查数组形状和单位。模型本身就是 IDE、静态检查、运行时校验和文档的共同来源。

### 同一个模型负责分析和导出

MolOP 的解析结果可以直接进入不同的导出模式：

- `model_dump`：保留模型结构，适合 Python 内部传递和进一步序列化。
- `to_unitless_dump`：去除 Pint 对象，生成适合 JSON 或数据库的数值结构。
- `to_unitless_dump_with_unit_keys`：把单位写入键名，减少跨进程和跨语言传递歧义。
- source evidence 序列化：同时保存文件元数据、frame 结果和源位置证据。

导出不是另写一套字典，而是模型能力的一部分。详见[序列化与 source evidence](behavior/serialization.md)。

### Python API 和 CLI 使用同一套语义

Python 和 CLI 共享解析、frame 选择、状态筛选、汇总和格式转换的核心语义。需要快速处理时使用 CLI，需要扩展分析时切换到 Python，不需要重新理解一套对象模型。

## 与 cclib：从“读取属性”到“管理工作流”

[cclib](https://cclib.github.io/) 的强项是量化计算输出解析和算法。典型路径是 `ccread`/`ccopen` 返回 data object，再访问标准属性或调用算法；它也提供 `ccget`、`ccwrite` 和 `ccframe` 等命令行工具。[官方解析指南](https://cclib.github.io/how_to_parse.html)

```text
cclib：文件 -> ccData -> 属性 / 算法
MolOP：文件集合 -> AutoParser -> batch/file/frame -> 检查 / 筛选 / 汇总 / 导出
```

### cclib 适合什么

- 已经知道输入程序和目标属性，只需快速取得结果。
- 需要直接使用 cclib 的标准属性或算法。
- 已有代码围绕 cclib data object 建立。

### MolOP 解决了什么

cclib 的核心结果字段以动态属性为主。字段是否存在、具体形状和单位，需要结合计算类型、[属性文档](https://cclib.github.io/data.html)和[源码](https://github.com/cclib/cclib/blob/master/cclib/parser/data.py)确认。对于一次性脚本这很直接；对于长期维护的项目，字段探索、类型判断和转换逻辑会不断落到调用方。

MolOP 把这些工作前移到公共模型中：

| 开发任务 | cclib | MolOP |
| --- | --- | --- |
| 发现字段 | 查属性列表并确认当前输出是否产生该字段 | IDE 直接提示显式字段和嵌套类型 |
| 处理多个文件 | 用命令或调用方循环组织 | batch 是 Python API 和 CLI 的共同对象 |
| 处理多个 frame | 从 data object 的数组中自行理解上下文 | `batch -> file -> frame` 直接表达层级 |
| 处理缺失和失败 | 在调用方组合异常和条件判断 | outcome、状态、diagnostics 和 parse presence 统一表达 |
| 交给下游 | 根据属性集合自行组织字典或表格 | Pydantic dump、单位导出和 source evidence 复用同一模型 |

因此，cclib 是优秀的输出属性和算法入口；MolOP 更适合把一批计算文件变成稳定、可维护的数据工作流。

## 与 pymatgen：文件工作流和材料分析各司其职

[pymatgen](https://pymatgen.org/) 的强项是材料科学对象和分析能力。[官方使用指南](https://pymatgen.org/usage.html)围绕 `Structure`、`Molecule` 及相关对象展开，适合组成、晶格、周期边界、对称性、相图、反应、态密度、能带和材料数据库分析。

两者的最佳组合是：

```text
量化计算输入/输出
    -> MolOP：识别、解析、检查、筛选、汇总、保留证据
    -> 显式适配：坐标 / 元素 / 晶格 / 拓扑 / 单位
    -> pymatgen：结构与材料分析
```

如果问题是“这批计算文件是否正常、最终结构和能量是什么”，使用 MolOP；如果问题是“这个周期结构的对称性、稳定性或电子结构如何”，使用 pymatgen。MolOP 负责可靠地产生上游数据，pymatgen 负责深入分析材料对象。

## 推荐使用方式

### 新项目

以 MolOP 作为默认文件入口：

1. 用 `AutoParser` 读取文件或文件集合。
2. 用 batch/file/frame API 检查状态、筛选结果和生成汇总。
3. 用类型化模型和导出模式把结果交给数据库、服务或下游分析。
4. 只有在需要专门算法或材料分析时，才接入 cclib 或 pymatgen。

### 已有 cclib 项目

保留 cclib 算法和已有分析逻辑；当项目开始需要多格式批处理、统一状态、类型检查、来源追踪或稳定导出时，将文件入口和批处理层迁移到 MolOP，并在边界处显式映射字段。

### 材料分析项目

让 MolOP 负责计算文件的读取、质量检查和结构提取，再把坐标、元素、晶格、拓扑和单位映射到 pymatgen。这样文件处理和材料分析各自使用最合适的模型。

## 最终判断

在 MolOP 已支持的格式范围内，面向计算化学文件的完整处理流程，MolOP 的使用体验明显优于 cclib 和 pymatgen。它把格式识别、批量生命周期、类型安全、规范数据模型、单位、来源证据和导出模式放在同一套 API 中。

cclib 适合“快速读取输出属性和使用解析算法”，pymatgen 适合“深入分析结构和材料”。当读者的任务是从文件可靠地走到可复用的数据结果，MolOP 应该是第一选择。

## 官方资料

- [cclib 官方主页](https://cclib.github.io/)
- [cclib：解析和写出文件](https://cclib.github.io/how_to_parse.html)
- [cclib：解析后的数据属性](https://cclib.github.io/data.html)
- [cclib：`ccData` 属性与运行时类型表](https://github.com/cclib/cclib/blob/master/cclib/parser/data.py)
- [pymatgen 官方主页](https://pymatgen.org/)
- [pymatgen：使用指南](https://pymatgen.org/usage.html)
