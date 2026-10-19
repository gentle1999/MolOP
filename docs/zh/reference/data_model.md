# 数据容器与跨 QM 语义

MolOP 将不同程序的输入和输出映射到公共数据容器。读取 Gaussian、ORCA 或 xTB 后，
可以使用相同字段访问同类物理量，同时通过来源、计算条件和约定识别它们之间的差异。

本页介绍公共容器及跨程序比较时的解释规则。按性质查字段见[科学字段索引](model_fields.md)，
各格式实际可填充的字段及解析依据见[格式支持概览](format_support.md)。容器中定义了某个字段，
并不表示每个程序或计算任务都会提供它。

## 数据层级

`AutoParser(path)` 返回批量文件容器；`batch[0]` 是文件，`batch[0][-1]` 是该文件保留的最后一帧。
文件可以包含多个计算步骤（segment），每个步骤又可以有多个几何或结果帧。

| 层级 | 容器与入口 | 内容 |
| ---- | ---------- | ---- |
| 批量 | `FileBatchModelDisk` | 文件集合、选择、汇总和批量转换。 |
| 文件 | `BaseChemFile`；结构、输入及结果文件分别使用 `BaseCoordsFile`、`BaseQMInputFile`、`BaseCalcFile` | 文件格式、来源、步骤、文件级计算设置和帧集合。 |
| 帧 | `BaseChemFileFrame`；`BaseCoordsFrame`、`BaseQMInputFrame`、`BaseCalcFrame` | 某个几何或结果的结构、请求、科学性质及状态。 |
| 性质 | 例如 `frame.energies`、`frame.vibrations` | 一个领域的结构化数据，或包含逐项记录的集合。 |
| 通用表示 | `PropertyTable`、`Spectrum`、`PropertyBundle` 等 | 供表格、光谱及跨系统导出使用的数据视图。 |

比较多步骤文件时，先选择对应任务和帧；文件级汇总不能替代每帧的计算条件。
`file_frame_index` 保留源文件中的帧顺序，`segment_index` / `segment_frame_index` 标识步骤及步内帧。
`frame_id` 属于当前帧集合，筛选或重排后不宜作为来源身份。

## 公共数据容器

以下科学数据类型可从 `molop.io.base_models.DataClasses` 导入；通用表示类型位于
`molop.io.base_models.Bases`。表中列出用途和主要访问方式，完整签名见 [Python API](api/io.md)。

### 输入请求与计算设置

| 容器 | 入口或关系 | 用途 |
| ---- | ---------- | ---- |
| `QMModelChemistry` | `model_chemistry` | 方法类别、具体方法、参考方法、泛函、基组、色散、自旋处理及输入溶剂设置。 |
| `QMBasisSet` | `model_chemistry.basis_sets` 中的记录 | 按角色和作用范围区分轨道基组、辅助基组、ECP 及元素/原子级覆盖。 |
| `QMResourceRequest` | `resource_request` | CPU 数、内存和程序特有资源选项；保留原始指令。 |
| `QMTaskRequest` | `task_requests` 中的记录 | 单点、优化、频率、扫描等任务及导数阶数、目标态和性质请求。 |
| `ExcitedStateRequest` | `excited_state_requests` 中的记录 | 激发态方法、请求根数、目标根、自旋通道和 root following。 |
| `ActiveSpace` | 多参考请求或结果中的 `active_space` | 活性电子数、轨道数及活性、非活性、冻结轨道列表。 |
| `MultireferenceStateBlock` | `MultireferenceRequest.state_blocks` 中的记录 | 按多重度和对称性划分的根数、参考空间及活性空间。 |
| `MultireferenceRequest` | `multireference_requests` 中的记录 | 多参考方法、态区块、阈值和请求修正。 |
| `ExplicitSolventRequest` | ORCA 输入的 `explicit_solvent_requests` | 显式溶剂放置数量、液滴/团簇设置和溶剂文件等请求。 |

### 坐标与注释

| 容器或类型 | 入口或关系 | 用途 |
| ---------- | ---------- | ---- |
| `CoordinateContainer` | 有序坐标记录的通用集合 | 提供 `items`、索引、遍历，以及真实原子与 dummy/ghost 的区分。 |
| `AtomInInternalCoords` | 内坐标中的单原子记录 | 参考原子索引、键长、键角、二面角及 dummy/ghost 标记。 |
| `InternalCoords` | `AtomInInternalCoords` 集合 | 内坐标与笛卡尔坐标的转换；长度为 Å，角度为度。 |
| `CoordinateParameter` | 单个坐标参数或扫描变量 | 参数名、原始表达式、起止值和步数。 |
| `CoordinateParameters` | `CoordinateParameter` 集合 | 有序参数、扫描定义和可用数值映射。 |
| `Comment`、`CommentContainer` | 文件或帧的 `comments` | 注释文字、原文、来源行及有序注释集合；`CommentKind` 表示注释类别。 |

笛卡尔结构直接保存在帧的 `atoms`、`coords`、`charge`、`multiplicity` 中；
`coords` 的形状为 `(N, 3)`。原始坐标、内坐标转换结果、SMILES 生成的二维坐标应按来源区分。

### 能量、热化学、溶剂与状态

| 容器 | 常用入口 | 用途 |
| ---- | -------- | ---- |
| `Energies` | `frame.energies` | 参考、电子及方法分辨能量；`total_energy` 提供当前可用的总能量。 |
| `EnergyObservation` | `energies.observations` 中的记录 | 通过 `method`、`quantity_semantics`、`source_label` 区分总能量、相关修正和分量。 |
| `ThermalInformations` | `frame.thermal_informations` | 零点能、热修正、内能/焓/自由能、熵、热容，以及温压和低频处理条件。 |
| `ImplicitSolvation` | `solvent` | 输出中识别的溶剂、模型、原子半径方案及介电常数。 |
| `Status` | `status` | SCF 收敛和程序正常终止分别记录。 |
| `GeometryOptimizationStatus` | `geometry_optimization_status` | 优化完成状态、能量/力/位移判据、阈值，以及程序报告的收敛结论。 |

### 轨道、布居、键级与总自旋

| 容器 | 常用入口或关系 | 用途 |
| ---- | -------------- | ---- |
| `MolecularOrbitals`、`MoleculeOrbital` | `frame.molecular_orbitals` 及单轨道索引访问 | alpha/beta 能量、占据数、对称性、总占据数、轨道基底及前线轨道标记。 |
| `NaturalAtomicOrbitals`、`NaturalAtomicOrbital` | 自然原子轨道集合及其记录 | 轨道编号、原子、角动量、轨道类型、占据数和能量；集合可区分总量及 alpha/beta。 |
| `NaturalBondOrbitals`、`NaturalBondOrbital` | 自然键轨道集合及其记录 | 键轨道类型、成键原子、占据数、能量及离域信息；可区分自旋通道。 |
| `ChargeSpinPopulations`、`AtomicPopulationSeries` | `frame.charge_spin_populations` 及 `populations[name]` | 可扩展的逐原子电荷/自旋布居；每组记录保留方案、量类型、自旋通道与来源标签。 |
| `BondOrders`、`BondOrderMatrix` | `frame.bond_orders` 及 `bond_orders[name]` | 可扩展的原子对矩阵；每组记录保留方案、自旋、打印阈值和缺失原因。 |
| `TotalSpin` | `frame.total_spin` | `spin_square` 与 `spin_quantum_number`，分别表示自旋平方期望值和自旋量子数。 |

自然轨道等类型可以表达独立数据；是否由某个格式自动填充，以该格式页为准。

### 振动、响应与其他性质

| 容器 | 常用入口或关系 | 用途 |
| ---- | -------------- | ---- |
| `Vibrations`、`Vibration` | `frame.vibrations` 及单模式索引访问 | 频率、约化质量、力常数、IR 强度、振型、拉曼活性和退偏振比，以及模式身份与方法信息。 |
| `Polarizability` | `frame.polarizability` | 偶极矩、极化率、四极及更高多极矩和电子空间范围。 |
| `NMR`、`ShieldingTensor` | `frame.nmr` 及 `nmr.shielding_tensors` | 逐原子屏蔽张量、各向同性值、主值、各向异性，以及 K/J 耦合矩阵和分量。 |
| `SinglePointProperties` | `frame.single_point_properties` | 垂直 IP、EA、全局亲电指数、逐原子 Fukui 指数及 FOD 等专用性质。 |
| `Dispersions` | 独立色散性质记录 | `C6AA`、`C8AA`；区别于输入中请求的经验色散修正。 |

力和 Hessian 直接保存在 `frame.forces`、`frame.hessian` 中，其原子顺序、轴顺序和方向由
对应的 `*_atom_order`、`*_axis_order`、`*_orientation` 字段说明。

### 电子态与多参考结果

| 容器 | 常用入口或关系 | 用途 |
| ---- | -------------- | ---- |
| `ElectronicConfiguration` | `ElectronicState.configurations` 中的记录 | 电子组态、系数、权重、占据模式与轨道编号。 |
| `ElectronicState` | `ElectronicStates.states` 中的记录 | 总能量、激发能、根编号、自旋/对称性、振子强度、跃迁偶极及明确的基态身份。 |
| `ElectronicStates` | `frame.electronic_states` | 态集合及态表、跃迁和电子光谱视图。 |
| `MultireferenceResult` | `frame.multireference_result` | 多参考输出的方法、活性空间、电子态、修正及诊断。 |

### 通用表示与导出

| 容器 | 用途 | 主要字段 |
| ---- | ---- | -------- |
| `BaseDataClassWithUnit` | 带单位数据类型的公共基础；提供数值和摘要导出。 | `to_unitless_dump()`、`to_unitless_dump_with_unit_keys()` |
| `PropertyPoint`、`PropertySeries` | 单个观测与有序序列。 | `value`、`uncertainty`、`points`、轴标签及元数据 |
| `TensorProperty` | 带坐标系说明的通用张量。 | `tensor`、`frame`、`metadata` |
| `PropertyTable` | 保留每列单位的等长列表格。 | `columns`、`row_labels`、`column_names`、`row()` |
| `PropertyTransition` | 两态间的跃迁。 | 初末态编号、能量、波长、振子强度和跃迁偶极 |
| `SpectralBand`、`Spectrum` | 谱带与谱带/跃迁集合。 | 峰位、宽度、强度、归属及 `bands` / `transitions` |
| `PropertyBundle` | 按数据形状汇集性质。 | 标量、向量、张量、表、光谱、序列和跃迁映射 |

提供 `to_property_bundle()` 的性质容器可转换为通用表示；振动使用 `to_ir_spectrum()` 和
`to_raman_spectrum()` 生成光谱。导出保留可用元数据，不会补齐来源未提供的物理量。
序列化方式见[来源证据与序列化](behavior/serialization.md)。

### 来源与解析状态

| 记录或类型 | 入口或关系 | 用途 |
| ---------- | ---------- | ---- |
| `SourceSpan` | 帧或步骤的 `source_span` | 回溯原文的字节、字符和行区间。 |
| `SourceSegmentEvidence` | 文件的 `source_segments` | 计算步骤的来源、任务、状态、原始帧数及实际保留的帧。 |
| `ParserProvenance` | 文件的 `parser_provenance` | 解析器版本、相关依赖版本及生效配置。 |
| `ParseDiagnostic` | 文件或帧的诊断集合 | 诊断编号、级别、作用范围、关联字段及消息。 |
| `ParsePresence`、`ParseCompleteness` | `parse_presence`、`parse_completeness` | 区分内容是否出现及解析是否完整；它们与计算收敛状态含义不同。 |

这些类型位于 `molop.io.base_models.source`。需要额外来源证据时，使用
`capture_source_evidence=True`；原文回溯与数据导出的用法见[来源证据与序列化](behavior/serialization.md)。

## 跨 QM 格式的语义统一

### 计算设置与结果的关系

跨程序读取计算设置时，优先查看 `model_chemistry`、`resource_request` 和各类请求容器。
`method`、`functional`、`basis_set`、`keywords` 等平铺字段是便捷投影。
Gaussian 的 R/U/RO 参考类型保存在 `spin_treatment` 中；泛函和色散同时保留规范名称及修正标识，
例如 `functional="B3LYP-GD3BJ"` 与 `dispersion_correction="GD3BJ"`。

`Freq`、NMR、激发态或显式溶剂请求描述拟执行的任务。对应的振动、屏蔽、电子态和溶剂结构
需要输出数据支持。`basis_sets` 中的角色和作用范围也应参与比较，不能只比较主基组字符串。
ORCA 的 `%maxcore` 是每进程内存，其资源含义不能仅凭 MB 数值与其他程序直接等同。

### 单位、原子顺序与坐标方向

| 物理量 | 公共默认单位 | 统一规则 |
| ------ | ------------ | -------- |
| 坐标 | Å | Bohr 坐标转换为 Å；保持来源原子顺序。 |
| 电子能、轨道能 | Hartree（Eh） | 不合并总能量、参考能和相关修正。 |
| 热校正、热力学总量 | kcal/mol | 原文 Hartree 按每个粒子的能量解释后转为摩尔量。 |
| 熵、定容热容 | cal/(mol·K) | 温度 K，压力 atm；缺少条件时不推测。 |
| 力、笛卡尔 Hessian | Eh/Bohr、Eh/Bohr² | 梯度按 `F = -dE/dR` 转为力；已打印力不再次取负。 |
| 频率、约化质量、力常数 | cm⁻¹、amu、mdyne/Å | 模式身份和归一化另行保留。 |
| IR 强度、拉曼活性 | km/mol、Å⁴/amu | 保留各自物理量名称。 |
| 偶极矩、极化率 | Debye、Bohr³ | 张量的分量排列与坐标方向同时参与解释。 |
| 屏蔽、K/J 耦合 | ppm、Hz | 屏蔽与化学位移、K 与 J 分别保存。 |

带单位数值可用 `.to(...)` 或 `.m_as(...)` 转为需要的兼容单位。
单位一致不表示坐标系一致：`forces` 为 `(N, 3)`，Hessian 为 `(3N, 3N)`；比较前应核对
原子顺序、轴顺序以及 input/standard/source 方向。没有方向依据时保持未知，不自动旋转所有张量。
帧会核对布居数组长度、键级矩阵大小和位移中的原子数是否与自身原子一致。
各可选逐模式数组要么缺失，要么与全部已存频率对齐；不完整表格不能使数值错配到其他模式。

### 能量与热化学

`energies.total_energy` 按容器定义选择可用总能量：优先 `electronic_energy`，否则依次使用
CCSD(T)、CCSD、MP5 至 MP2、`reference_energy`。各方法字段继续保留；比较计算协议时，
还需核对方法和来源标签，不能只按该便捷字段判断方法一致。

热化学使用以下定义，其中 `Eel` 为电子能量：

- `U_0 = Eel + ZPVE`；
- `U_T = Eel + TCE`；
- `H_T = Eel + TCH`；
- `G_T = Eel + TCG`。

Gaussian log 区分校正与“电子能 + 热量”总和；FCHK 由总量相减得到的热校正保留
`derivations` 公式；ORCA 的增量焓修正不能直接当作 `TCH`，打印的熵能项需除以温度得到熵。
xTB 保留实际提供的热化学汇总及 RRHO/QRRHO 条件。

比较自由能时，同时核对温度、压力、`entropy_model`、频率缩放、低频模式排除阈值、
`rotor_reference_frequency` 和虚频处理。谐振子/转子插值参考频率与模式排除阈值是不同条件；
仅有 FCHK 总量时不补造熵模型。

### 轨道与自旋

alpha/beta 数组按轨道索引对齐，较短或缺席的通道补缺失值，保留较长通道的数据。
补齐本身不证明两个通道共享轨道谱；只有来源明确为共享受限轨道时才共享能量。

`orbital_basis` 区分共享空间轨道、自旋分辨轨道及未知情况；`total_occupancies` 保留空间轨道总占据数，
`occupancy_origin` 说明占据数是直接读取、按电子数构造还是依据受限/高自旋约定推导。
闭壳层总占据数可以等分；有高自旋依据的整数占据数可以拆分；其余无法确定的自旋占据数保持未知。
xTB 省略的轨道区间不能视为未占据，分数占据数本身也不能确定 HOMO/LUMO 或 SOMO。
打印的前线轨道标记保存在 `reported_homo_index`、`reported_lumo_index` 中。

`spin_square` 表示 `<S²>`。由它得到的
`spin_quantum_number = sqrt(<S²> + 1/4) - 1/2` 是推导的有效自旋，
对有自旋污染的波函数可能不同于理想自旋。输入多重度不替代 `<S²>` 观测。

### 布居与键级

布居通过 `population_items()` / `population_names` 遍历，键级通过
`bond_order_items()` / `bond_order_names` 遍历。按记录的 `scheme`、`quantity`（布居）和
`spin_channel` 比较，避免将不同方案或总量与 alpha/beta 分量混合。

`populations["mulliken_charges"].values` 是逐原子数组；
`bond_orders["mayer_bond_order"].values` 是原子对矩阵。Wiberg、Mayer、MO 及其他方案使用同一类记录，
方案名保持开放。阈值表未打印的原子对保留未知，`print_threshold` 和 `missing_reason` 说明原因；
实际打印的零及允许的负值保持原值。
向 RDKit 嵌入键级时，仅给现有化学键添加属性，分子图保持不变；
非成键原子对的数值保留在帧的矩阵中。

### 振动、IR 与拉曼

`mode_indices` 是内部零基编号，`source_mode_indices` 单独保留源模式号；未打印编号的模式使用 `None`。
频率、谱强度与位移按模式身份对齐。筛选模式时，两种编号一同保留。

振型数组使用 `(mode, atom, xyz)` 顺序。来源的 `normalization`、`mass_weighting`、`orientation`
和原子顺序随数据保留，不通过再次缩放强行统一不同程序的振型。位移使用 Å 表示是容器存储约定，
归一化向量的模长不能直接解释为实验绝对振幅。

`IR_intensities` 是 IR 强度，`raman_activities` 是拉曼活性。Gaussian 的 P/U 退偏振比分别位于
`raman_depolarization_ratios_polarized` / `raman_depolarization_ratios_unpolarized`；
未明确偏振约定的来源使用无后缀字段。退偏振比无量纲，拉曼活性不包含实验激光与温度因子。
xTB 的 `frequency_method` 与 `raman_activity_method` 可不同；未打印的位移和退偏振比保持缺失。

### NMR 与电响应

`ShieldingTensor.atom_index` 是源原子顺序的零基编号；`coupling_atom_indices` 明确耦合矩阵的轴。
K、J 及 FC/SD/PSO/DSO 分量独立保留。FCHK 缺少同位素依据时，不由 K 补造 J。

屏蔽张量与已打印的各向同性值、主值、各向异性分别保留；缺少标量时可由张量推导。
`anisotropy_convention` 标明各向异性约定：公共缺省推导采用 Haeberlen，特定格式可保留其程序约定。
屏蔽值需结合参考物才能用于化学位移，方向未知的张量也不能直接按分量与其他程序比较。
偶极与极化率按来源单位转换，紧凑张量的分量顺序以对应格式说明为准。

### 电子态、状态与缺失值

`ElectronicState.state_index` 是内部唯一的非负编号，`root` 保留程序根编号。
`is_ground_state=True` / `False` / `None` 分别表示明确基态、激发态和未知；
列表首项或 root 为零都不能单独证明基态。只保存激发态的集合仍可导出跃迁；
未保存基态记录时，跃迁内部初态编号可为空，来源标签仍随跃迁保留。

程序正常终止、SCF 收敛、几何优化完成分别记录。缺少状态证据时使用 `None`，不能当成 `False`，
也不能从单个能量或正常终止行推断所有任务成功。

未知轨道能量、键级与耦合矩阵项在数值数组中使用 `np.nan`，未知占据数和标签使用 `None`。
允许缺失的字段在 JSON 导出时将缺项表示为 `null`；`array_mode="ndarray"` 可保留数组 NaN。
打印零保持为零。容器整体缺失、数组为空和数组内某项未知是不同情况，比较与导出时应分别处理。

SMILES 的二维坐标、从坐标恢复的分子图及 fakeG 文本是生成数据。格式转换使用已有数据，
不会执行新的量化计算；公共容器统一访问与解释规则，也保留不同程序的可用信息和物理约定。
