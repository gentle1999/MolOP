# Gaussian-like 渲染器

<!-- format-support:fakeg -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `fakeg` |
| 扩展名 | `.fakeg` |
| 读取 | 否 |
| 写入 | 是 |
| Registry 角色 | 文件 writer |
| 数据层级 | 已解析 Gaussian 输出数据 |

`fakeg` 将已解析的 Gaussian 输出数据渲染为 Gaussian-like 文本，适合人工检查、
兼容性测试和下游工具衔接。它不是 Gaussian log 的逐字复现器。

`format_transform("fakeg")` 遵循通用转换默认值 `frame=-1`；需要全文件
Gaussian-like 渲染时传入 `frame="all"`。

```python
from molop import AutoParser

gaussian_log = AutoParser(
    "calculation.log",
    parser_detection="g16log",
    n_jobs=1,
)[0]
rendered = gaussian_log.format_transform("fakeg", frame="all")
print(type(rendered).__name__)
```

??? example "输出类型"

    ```text
    str
    ```

`G16LogFile.render_fakeg()` 的默认 frame 是 `"all"`，与 `format_transform()` 的默认值不同。

| 能力 | 支持程度 | 可渲染内容 | 边界 |
| ---- | -------- | ---------- | ---- |
| <!-- feature-area:File-level Gaussian-like writer -->文件级 Gaussian-like writer | 部分支持 | 从已解析 Gaussian 输出数据写出 `.fakeg` 文件。 | 只支持文件级渲染；不提供帧级 fakeG writer；选中帧遵循 `frame`；不保证逐字复现 Gaussian log。 |
| <!-- feature-area:Structure and SCF energy rendering -->结构与 SCF 能量渲染 | 部分支持 | 从坐标和能量字段渲染归一化 orientation 与 SCF-cycle 文本。 | 渲染结果是语义化 Gaussian-like 文本，不是原始 log 文本。 |
| <!-- feature-area:Vibrational frequency rendering -->振动频率渲染 | 部分支持 | frequency、reduced mass、force constant、IR intensity 和逐模式位移区段。 | 当前只声明 frequency 与 IR 相关字段；不声明 Raman/VCD 区段。 |
| <!-- feature-area:Thermochemistry rendering -->热力学渲染 | 部分支持 | temperature、correction、energy、entropy、heat capacity、mass、inertia 和转动/振动 metadata。 | 这是归一化的热力学摘要，不是完整 Gaussian thermochemistry pretty-printer。 |
| <!-- feature-area:Reparseable frequency and thermochemistry output -->可再次解析的频率与热力学输出 | 已支持 | 渲染出的频率和热力学内容可以再次解析回帧字段。 | Round trip 证明的是支持的语义字段，不表示与 Gaussian 原始输出逐字等价。 |

## fakeG 能力覆盖矩阵

下表区分“从结构化字段重新生成”和“沿用 component 的 raw 文本”。只有前者才表示
fakeG 可以脱离原始 Gaussian log 合成该区段。

| fakeG 能力 | 数据来源 | 结构化重建 | 再解析 | 原文保真 | 当前边界 |
| ---------- | -------- | ---------- | ------ | -------- | -------- |
| Registry 文件 writer | 带坐标的 file 模型 | 支持 | 部分 | 不支持 | 注册为 file writer，扩展名 `.fakeg`；rich output 取决于输入 frame 是否已有 Gaussian 结果字段。 |
| Registry frame writer | 单个 frame | 不支持 | 不适用 | 不适用 | `frame.format_transform("fakeg")` 明确报错；registry 未注册 fakeG frame writer。 |
| `frame.render_fakeg()` | G16Log frame | 支持 | 部分 | 不支持 | 模型方法可直接渲染单 frame，与 registry frame writer 是不同接口。 |
| frame 选择 | `frame` selector | 支持 | 不适用 | 不适用 | `format_transform()` 默认 `-1`；`G16LogFile.render_fakeg()` 默认 `"all"`，支持 int/slice/list/all。 |
| 分离输出 | `embed_in_one_file=False` | 支持 | 部分 | 不支持 | `render_fakeg()` 可返回逐 frame 字符串列表；registry 写盘行为仍遵循通用 file writer 合同。 |
| 文件 header | version、options、keywords、title | 支持 | 部分 | 不支持 | 重建 Gaussian-like banner、Link0、route 和 title；缺失 route 时使用 `#p fakeg`。 |
| shared-memory CPU 行 | `%nprocshared` | 支持 | 部分 | 不支持 | 复用 Gaussian Link0 parser；非纯数字值只生成泛化说明。 |
| `Symbolic Z-matrix` header | 首个 frame 的 atoms/coords/charge/multiplicity | 支持 | 支持 | 不支持 | 名称仿照 Gaussian，但当前实际输出 Cartesian 原子坐标，不重建原始 Z-matrix。 |
| input/standard orientation | atoms、coords、standard_coords | 支持 | 支持 | 不支持 | 优先 standard orientation，否则 input orientation；center/type 和小数位按固定模板生成。 |
| SCF energy | `energies.reference_energy` | 支持 | 支持 | 不支持 | 统一生成 `SCF Done: E(SCF)` 和单 cycle 描述；不保留 functional label、cycle 数或 convergence 过程。 |
| post-HF energy | MP2-MP5、CCSD、CCSD(T) 字段 | 不支持 | 不支持 | 不支持 | 当前不会为这些字段生成对应 Gaussian post-HF energy 行。 |
| total spin | `total_spin` | 支持 | 部分 | 不支持 | 可生成 `S**2` 和 `S`；不重建 spin-contamination 诊断。 |
| MO 和 population | `molecular_orbitals`、`charge_spin_populations` | Raw/不支持 | 不保证 | 不支持 | 没有字段驱动的规范化 renderer；解析树存在 raw component 时可能沿用原片段。 |
| dipole / polarizability / multipoles | `polarizability` | Raw/不支持 | 不保证 | 不支持 | component tree 可携带 raw 文本，但从纯模型字段不能完整合成 response 区段。 |
| harmonic frequencies | `vibrations.frequencies` | 支持 | 支持 | 不支持 | 每批最多三 mode，统一使用 `A` symmetry label；不保留原始 symmetry assignment。 |
| reduced mass / force constant / IR | `vibrations` 对应数组 | 支持 | 支持 | 不支持 | 只在字段存在时输出；固定格式和精度。 |
| normal-mode displacement | `vibration_modes`、atoms | 支持 | 支持 | 不支持 | 完整 frequency block 可输出逐原子向量；单 child-node renderer 只给出有限预览。 |
| Raman / VCD / ROA | 无对应稳定字段 | 不支持 | 不支持 | 不支持 | frequency header 含 Gaussian 惯用 Raman 文案，但当前不输出 Raman activity 或 depolarization 数值。 |
| thermochemistry header | temperature、pressure | 支持 | 支持 | 不支持 | 只输出已有值，不推断缺失条件。 |
| mass、惯性矩和转动/振动温度 | `thermal_informations` | 支持 | 支持 | 不支持 | 采用 Gaussian-like 固定标签和精度，不复现原始表格布局。 |
| ZPVE / thermal corrections | ZPVE、TCE、TCH、TCG | 支持 | 支持 | 不支持 | 单位和值来自规范化容器。 |
| U0 / UT / H / G | thermal summary 字段 | 支持 | 支持 | 不支持 | 生成 Gaussian parser 可再次识别的 summary 行。 |
| entropy / heat capacity | `S`、`C_V` | 支持 | 支持 | 不支持 | 输出规范化 total value 和 thermochemistry 表头，不重建各分量。 |
| Cartesian forces | `forces` | Raw/不支持 | 不保证 | 不支持 | 无字段驱动的 force table renderer；有原 component raw 文本时可能保留片段。 |
| Cartesian Hessian | `hessian` | Raw/不支持 | 不保证 | 不支持 | 无字段驱动的完整 second-derivative renderer。 |
| optimization convergence | `geometry_optimization_status` | 部分 | 支持 | 不支持 | 文件级多帧渲染会生成 step number、四项 convergence；单 frame 主要依赖 component 内容。 |
| 多帧 optimization 选择 | optimization status + total energy + 无 vibrations | 支持 | 支持 | 不支持 | 只纳入满足条件的 optimization frames；其他选中 frame 可能不会出现在合并 body。 |
| frequency frame 选择 | 最后一个含 vibrations 的 frame | 支持 | 支持 | 不支持 | 多帧合并只选择最后一个 frequency frame，不保证逐 frame 完整输出。 |
| runtime | `running_time` | 部分 | 部分 | 不支持 | 单 frame 可生成 `Job cpu time`；多帧 file body 会移除逐 frame runtime 行。 |
| termination | renderer 合成 | 支持 | 支持 | 不支持 | 多帧 optimization/frequency 段后写入 `Normal termination`，不代表原文件实际正常结束。 |
| component tree 全树渲染 | parsed/synthetic components | 部分 | 部分 | 不支持 | 有专用 renderer 的节点从 payload 生成；其他节点可能返回 `raw_text`。 |
| 指定 component node | `component_tree.render_node(name)` | 支持 | 不适用 | 不支持 | 可渲染如 `l716.forceconstants`、`l716.vibration.mode[0]`、temperature 等节点。 |
| fakeG 再解析 | 生成的 Gaussian-like 文本 | 部分 | 支持 | 不适用 | 已验证 structure、optimization frame、frequency 和 thermochemistry；不保证所有源字段 round-trip。 |
| 原始 Gaussian log 复现 | 原始 log | 不支持 | 不适用 | 不支持 | 不保留原始 link 顺序、迭代细节、banner、诊断、时间戳或字节级排版。 |

fakeG 的设计目标是为测试和下游兼容生成“足够像 Gaussian 且可被 MolOP 再解析”的文本。
它不应作为 Gaussian 计算成功、原始终止状态或完整 provenance 的证据。
