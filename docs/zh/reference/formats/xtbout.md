# xTB 输出

<!-- format-support:xtbout -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `xtbout` |
| 扩展名 | `.out`, `.log`, `.xtbout` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | QM 结果；仅在输出实际打印坐标时提供坐标 |

## 快速读取

```python
from molop import AutoParser

frame = AutoParser("xtb.out", parser_detection="xtbout")[0][-1]
print(frame.qm_software_version, frame.is_normal)
print(frame.energies.total_energy.m_as("hartree"))
print(len(frame.atoms), frame.coords.shape)
```

??? example "输出约定"

    输出版本、终止状态、总能量和几何规模。如果 xTB 单点 stdout 只引用外部坐标而没有打印
    几何，最后一行可能是 `0 (0, 3)`。

MolOP 不会隐式读取相邻坐标文件。

MolOP 解析 xTB 5 和 `>=6` 的命令行标准输出，同时覆盖 xTB 5/6.1 的 legacy
打印族和 xTB `>=6` 的 modern 打印族。由于 xTB 通常直接从命令行接收坐标文件与选项，
因此不单独定义 xTB 输入文件格式。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Version scope, setup, tasks, and status -->版本、设置、任务与状态 | 部分支持 | 识别 xTB 5，并将 `>=6` 作为 modern 打印族；解析版本、程序调用、坐标文件名、方法、电荷、多重度、OMP 线程数、任务请求、终止状态、SCC 状态和总 wall time。 | 明确拒绝 5 之前的版本。版本早于 5 的输出不支持；xTB 5 仅部分支持，后续大版本的兼容性未确认。 |
| <!-- feature-area:Geometry availability -->可用几何 | 部分支持 | 解析 legacy Bohr `$coord`、legacy setup 坐标、现代 XYZ final structure 和内嵌 V2000 SDF final structure。 | xTB 单点 stdout 经常只引用外部坐标文件而不打印坐标。这类日志仍保留结果 frame，但 geometry 为空；parser 不会隐式读取相邻文件。 |
| <!-- feature-area:Energies and geometry optimization -->能量与几何优化 | 部分支持 | 暴露 Hartree 单位的最终总能量、带来源标签的能量证据、优化收敛状态、能量变化判据和 xTB gradient norm 字段。 | 不结构化能量分解表和完整优化轨迹；只表示 stdout 实际打印的最终结构。 |
| <!-- feature-area:Orbitals, populations, and dipole -->轨道、布居、偶极矩与键级 | 部分支持 | 解析 legacy/modern 轨道能量与占据数、GFN1 Mulliken/CM5 电荷、GFN2 Mulliken 电荷、分子偶极矩和已打印的转动常数。 | 支持阈值 Wiberg 矩阵，未打印的原子对保持未知；四极矩、色散性质表和省略轨道仍不完整。 |
| <!-- feature-area:Vibrations and thermochemistry -->振动与热化学 | 部分支持 | 解析投影后的物理频率、约化质量、IR 强度、零点能、总焓、总自由能、热容、熵、温度和已打印的分子质量。 | 支持明确 Å⁴/amu 单位的拉曼活性并保留频率/活性方法；不映射位移、Hessian sidecar、退偏振比及旧式未定义 (amu) 列。 |
| <!-- feature-area:xTB analysis properties -->xTB 分析性质 | 部分支持 | 在请求且打印时解析垂直 IP、垂直 EA、全局亲电指数和逐原子 Fukui 指数。 | FOD、metadynamics、MD、ONIOM、GFN-FF topology 以及 JSON/sidecar 输出不在当前 stdout 合同内。 |

## xTB 5 / xTB >=6 能力覆盖矩阵

下表描述 stdout reader 的字段级合同。状态含义如下：

- **支持**：可以读取对应版本族的所列输出形式。
- **部分**：只覆盖已知打印形式、最终值或有限任务语义。
- **未验证**：解析路径可能适用，但缺少该版本族的真实输出样本，不作为稳定承诺。
- **不支持**：当前不会生成对应结构化字段。
- **不适用**：不属于 stdout reader 的职责。

`xTB >=6 modern` 列的“支持”由现有 xTB 6 真实样本建立。xTB 7 及后续版本使用同一
modern 合同，但在获得真实输出前仍属于前向兼容，而不是独立的版本验证。

| xTB 输出能力 | xTB 5 legacy | xTB >=6 modern | 结构化结果 | 当前边界 |
| ------------ | ------------ | -------------- | ---------- | -------- |
| banner 与版本 | 部分 | 支持 | `qm_software`、`qm_software_version` | 接受主版本 `5` 和 `>=6`；xTB 5 只覆盖部分旧式布局，后续大版本兼容性未确认。 |
| 拼接的多次命令行运行 | 未验证 | 支持 | 每个 version banner 开始一个 segment 和结果 frame | 依赖每次运行打印版本 banner；不根据任意 shell 分隔文本拆分。 |
| program call 与坐标文件名 | 部分 | 支持 | `keywords`、`input_file_name`、模型化学 raw keywords | 只读取 stdout 中的 `program call` 和 `coordinate file` 行，不恢复 shell 环境或相对路径基准。 |
| GFN Hamiltonian | 部分 | 支持 | `method`、`model_chemistry.method`，方法族为 `SEMIEMPIRICAL` | 识别已打印的 GFN、GFN1、GFN2 和 GFN-FF 风格标签；不验证方法与版本兼容性。 |
| 电荷、多重度与 OMP 线程 | 部分 | 支持 | `charge`、`multiplicity`、`request_num_cpu`、`resource_request` | 多重度由 `--uhf` 或打印的未成对电子数加一得到；未打印时使用 xTB 常见默认值。 |
| SP、优化、频率与梯度请求 | 部分 | 支持 | `task_requests`，含 derivative order 和性质请求 | 从 program call 的 `--opt`、`--ohess`、`--hess`、`--grad`、`--vip`、`--vea`、`--vipea`、`--fukui`/`--vfukui` 推断；不覆盖所有 xTB CLI 组合。 |
| 隐式溶剂与温度 | 未验证 | 部分 | `solvent`、`temperature`、`electron_temperature` | 只解析已打印的 model、solvent、振动/溶剂温度和 electronic temperature；不重建完整 ALPB/GBSA 参数。 |
| legacy setup 坐标与 `$coord` final structure | 部分 | 部分 | `atoms`、Angstrom `coords`、坐标精度和 observed provenance | 输入为 Bohr 时转换到 Angstrom；不保留原始单位、空白和逐行文本。 |
| modern XYZ final structure | 未验证 | 支持 | `atoms`、Angstrom `coords`、坐标精度和 observed provenance | 只读取 stdout 中 `final structure:` 后的最终 XYZ，不保留优化中间结构。 |
| 内嵌 V2000 SDF final structure | 未验证 | 支持 | 原子与坐标 | 当前只取 atom block；键、formal charge、SDF property block 不在此 reader 中重建。 |
| 外部坐标、`xtbopt.xyz` 等 sidecar | 不支持 | 不支持 | 只保留已打印的坐标文件名 | 不隐式读取相邻文件，以维持 stdout source span 和哈希的单源合同。 |
| 最终总能量 | 部分 | 支持 | Hartree `energies.total_energy`；启用 source evidence 时附带 `EnergyObservation` | 选择最后一个已知 final-energy 打印形式；不把每个 SCC iteration 作为能量轨迹。 |
| SCC/能量分解表 | 不支持 | 不支持 | 无专用分解字段 | isotropic electrostatic、dispersion、repulsion 等分量当前不结构化。 |
| SCC 与终止状态 | 部分 | 支持 | `status.scf_converged`、`status.normal_terminated` | SCC 成功可由收敛行或最终能量推断；正常终止依赖 `finished run` 标记。 |
| CPU 与墙钟时间 | 部分 | 支持 | 秒单位独立的 `cpu_time` / `wall_time`；`running_time` 是墙钟时间别名 | 文件按时间类型分别累计选中任务段。优先打印的总计，总计块缺少的时间类型保持未知；没有总计时，各自采用最后一条匹配记录。不暴露详细模块计时和调度器时间。 |
| 优化收敛与阈值 | 部分 | 支持 | `geometry_optimization_status`、energy change/threshold、gradient norm/threshold | 只保留最终打印值；不结构化每一步的收敛历史。 |
| 完整优化轨迹 | 不支持 | 不支持 | 无逐步 frame 序列 | 一个命令行 run 当前对应一个结果 frame，而不是每个优化 step 一个 frame。 |
| 梯度向量 | 不支持 | 不支持 | 仅能记录 `gradient` 任务请求和最终 gradient norm | stdout 中的逐原子梯度和 `gradient` sidecar 尚未映射到结构化数组。 |
| 轨道能量与占据数 | 部分 | 支持 | 轨道能量、`total_occupancies` 及可确定的自旋占据数 | 读取 legacy `occ./eps` 与现代轨道表；省略区间不能恢复。总占据数保存在 `total_occupancies`，只有具备自旋依据时才拆分到 alpha/beta 通道。 |
| 原子电荷 | 未验证 | 部分 | GFN1 Mulliken/CM5、GFN2 Mulliken，并带 source label 和可扩展 population metadata | 不解析自旋密度、完整 population table 或所有 Hamiltonian 的电荷方案。 |
| 偶极矩与转动常数 | 未验证 | 部分 | Debye dipole 向量、GHz rotation constants | 不结构化偶极矩总模、四极矩、极化率张量和更高多极矩。 |
| 频率、约化质量与 IR 强度 | 未验证 | 支持 | `vibrations`，含 imaginary-mode 计数 | 只覆盖已打印的 projected physical frequencies；缺失列不会合成。 |
| 模式位移与 Hessian | 不支持 | 不支持 | 无对应结构化数组 | 不读取位移与 Hessian sidecar。 |
| 拉曼活性 | 未验证 | 部分支持 | `vibrations.raman_activities`，Å⁴/amu | 读取 xTB 6.7.1 PTB stdout 的活性表；不补造退偏振比。 |
| Wiberg 键级 | 未验证 | 部分支持 | `bond_orders["wiberg_bond_order"]` | 保留 AO 基、阈值、双向舍入策略与未知原子对。 |
| 热化学 | 未验证 | 支持 | ZPVE、enthalpy、free energy、`C_V`、entropy、temperature、molecular mass | 只保留汇总值；不结构化 rotor/interpolation 明细和各项热校正分解。 |
| VIP、VEA、GEI 与 Fukui | 未验证 | 支持 | `single_point_properties` 标量和逐原子 Fukui 数组 | 仅在对应分析被请求且 stdout 实际打印结果时存在。 |
| FOD、MD、metadynamics、ONIOM、GFN-FF topology | 不支持 | 不支持 | 无专用字段 | 当前 stdout 合同不覆盖这些分析或工作流。 |
| JSON、molden、Hessian 等辅助输出 | 不支持 | 不支持 | 无 sidecar 合并 | parser 只消费传入的 stdout 文本，不扫描工作目录。 |
| xTB 输入与输出 writer | 不适用 | 不适用 | Reader only | MolOP 不定义 xTB 输入格式，也不尝试从结果模型重建命令行或 stdout。 |

## 坐标边界

xTB 单点任务通常只打印坐标文件路径。MolOP 不会静默读取这个相邻文件，因为这会把
第二个源文件混入只描述 stdout 的 source span 和哈希合同。优化输出实际打印 final
structure 时，正常提供坐标字段。

## 版本边界

parser 接受 major version `5`，并把所有 `>=6` 版本交给 modern 解析合同。对于识别为
xTB、但 major version 小于 `5` 的输出，会明确报告版本不受支持，而不是继续交给
Gaussian 或 ORCA reader。

PTB Raman 结果用 `frequency_method`、`raman_activity_method` 区分频率与活性方法，
单模式和频谱导出也保留这些字段。stdout 解析器保持主输出精度，不静默用 `vibspectrum`
或仿 Gaussian `g98.out` 替换数值。详见 [公共数据模型](../data_model.md)。

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[xTB 命令行手册](https://github.com/grimme-lab/xtb/blob/v6.7.1/man/xtb.1.adoc); [xTB 分子性质](https://xtb-docs.readthedocs.io/en/latest/properties.html); [xTB PTB](https://xtb-docs.readthedocs.io/en/latest/ptb.html)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Version scope, setup, tasks, and status -->版本、设置、任务与状态 | `Version banners`; `program call`; `Hamiltonian`; `charge/spin`; `OMP`; `solvent`; `T(el)`; `normal/converged/failure labels`; `total/component wall-time` | 映射版本、方法、任务、资源、溶剂请求及条件；多重度由未配对电子设置加一推导，不是 <S²> 观测。版本 banner 确定一次运行；支持 xTB 5 的部分旧式布局和 xTB 6 的现代布局；后续大版本的兼容性未确认。 `cpu_time`、`wall_time` 分别读取打印的 CPU/墙钟秒数，`running_time` 是墙钟时间别名；存在总计块时，其中缺少的物理量保持未知，否则各自采用最后一条记录；电子温度与振动温度区分。 |
| <!-- parsing-basis:Geometry availability -->可用几何 | `final structure: XYZ / $coord / V2000 atom block`; `legacy setup coordinates` | XYZ/SDF 坐标为 Å，旧式 Bohr 坐标转为 Å；仅打印的最终几何算作观测，文件名或布居表不是坐标，不自动读取相邻文件。 |
| <!-- parsing-basis:Energies and geometry optimization -->能量与几何优化 | `TOTAL ENERGY / total energy / total E`; `energy/gradient changes and limits`; `optimization convergence` | 使用最后匹配的 Eh 总能量；优化值、阈值与收敛状态独立于程序正常终止；标量梯度范数不当作逐原子力数组。 |
| <!-- parsing-basis:Orbitals, populations, and dipole -->轨道、布居、偶极矩与键级 | `occ./eps and modern orbital tables`; `HOMO/LUMO labels/ellipsis`; `GFN1 Mulliken/CM5 or GFN2 q`; `dipole full`; `rotational constants`; `Wiberg/Mayer heading with Wiberg subheading` | 保留总占据、省略轨道及打印前线标记，自旋拆分需依据；电荷方案保留标签，full 偶极分量按原子单位转换；转动波数转 GHz。Wiberg AO 矩阵保留阈值、未知对及双向舍入首次值策略，完整 WBO 行号可在无坐标时确定矩阵轴。 |
| <!-- parsing-basis:Vibrations and thermochemistry -->振动与热化学 | `Projected eigval`; `indexed masses/IR/Raman`; `Hessian/intensity method labels`; `ZPVE/TOTAL ENTHALPY/FREE ENERGY`; `TOT CV/S`; `thermodynamic settings` | 采用 cm⁻¹、amu、km/mol 及明确 Å⁴/amu 活性；不映射旧式未定义 (amu) 数组和 g98 零占位。按 MolOP 策略去除开头 5/6 个近零外部模式并对齐索引。保留 GFN2-xTB/PTB 方法差别、K 温度、RRHO/QRRHO 条件及 HO/转子参考，不补造位移、退偏振比或压力。 热能原文为 Eh/particle，默认转换为 kcal/mol，C_V/S 为 cal/mol/K。 |
| <!-- parsing-basis:xTB analysis properties -->xTB 分析性质 | `Vertical IP / EA / global electrophilicity index`; `f(+)/f(-)/f(0) atom rows` | `single_point_properties` 的 VIP/VEA/GEI 转换前为 eV/particle；Fukui 无量纲数组按打印原子序；仅把实际打印的表视为观测，不据此映射 FOD、MD 或 metadynamics。 |

<!-- parsing-basis:end -->
