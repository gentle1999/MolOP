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

## 快速读取

下载[共享 ORCA 样例](../../../assets/examples/water_mp2.out)后运行：

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", parser_detection="orcaout", n_jobs=1)[0][-1]
print(frame.qm_software, frame.method, frame.is_normal)
print(frame.energies.total_energy.m_as("hartree"))
```

??? example "输出"

    ```text
    ORCA MP2 True
    -74.999374598107
    ```

ORCA 输出解析以用户可获得的结果属性为口径，例如能量、力、振动、布居、溶剂化、
激发态和响应性质。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Molecular orbitals -->分子轨道 | 部分支持 | 结构化提取轨道能量、分数占据数、alpha/beta 自旋通道及打印的对称性标签。 | 按源轨道编号对齐两个自旋通道；未打印的能量在带单位数组中补 `np.nan`，JSON 导出为 `null`；占据数及对称性标签补 `None`。尚不提供轨道系数。 |
| <!-- feature-area:Thermochemistry -->热化学 | 部分支持 | 提取最后一个温度区块的 ZPVE、热修正、内能、焓、Gibbs 自由能、熵、分子质量及转动常数；温度和压力写入帧字段。 | 将熵的能量项除以温度得到熵；尚不提供热容及其他配分函数诊断。 |
| <!-- feature-area:Software metadata, frames, and status -->程序信息、帧与状态 | 部分支持 | 读取软件名、版本、任务与几何帧、正常终止状态及耗时。 | 只支持已识别的程序标题、坐标和终止标记；截断文件可能缺少状态或属性。 |
| <!-- feature-area:Energies -->能量 | 部分支持 | 分别保留 Hartree 单位的最终总能量、参考能量及已支持的相关能量。 | 方法专用的修正项与能量分解表可能未映射；相关修正不当作总能量。 |
| <!-- feature-area:Forces and geometry optimization -->力与几何优化 | 部分支持 | 由完整笛卡尔梯度生成力数组，读取几何优化判据、阈值及完成状态。 | 不提供全部优化诊断；缺失或截断梯度不生成完整力数组。 |
| <!-- feature-area:Vibrations and vibrational spectra -->振动与振动光谱 | 部分支持 | 结构化提取频率、源模式编号、Cartesian 正常模式、旧版及新版 IR 强度、拉曼活性及退偏振比。 | 各模式数组需要完整行，并按已保留频率的模式编号对齐；拉曼活性不包含实验激光及温度因子。 |
| <!-- feature-area:Charge and spin populations -->电荷、自旋布居与键级 | 部分支持 | 解析 Mulliken、Loewdin、Hirshfeld 原子电荷及自旋表、已打印的 Mayer 键级及总自旋平方期望值。 | Mayer 矩阵将未打印项保留为 NaN（JSON null），并保存打印阈值；尚不提供 reduced orbital population、Mayer 价态诊断及其他未支持的布居方案。 |
| <!-- feature-area:Solvation -->溶剂化 | 部分支持 | 读取已识别的 CPCM/SMD 溶剂模型、名称及打印参数。 | 不提供完整溶剂能量分解；溶剂名称和输入请求不代表已读取自由能校正。 |
| <!-- feature-area:Excited states and multireference results -->激发态与多参考结果 | 部分支持 | 读取已支持的 TDDFT、ADC2、EOM-CCSD、STEOM-CCSD、STEOM-DLPNO-CCSD 和 ROCIS 态能量及电偶极跃迁性质。 | SOC、磁偶极与速度形式不并入同一电偶极结果；多参考输入请求不代表输出结果已被读取。 |
| <!-- feature-area:Response properties -->电响应、NMR 与总自旋 | 部分支持 | 结构化提取 polarizability、NMR 屏蔽张量及各向同性自旋耦合 J 矩阵，包括屏蔽主值、各向异性及 FC/SD/PSO/DSO 的 J 值。 | 尚不提供完整 Cartesian 耦合张量或 SD/FC 交叉项；耦合矩阵保留源原子编号，未打印项使用 NaN（JSON null）。 |

拉曼数据位于 `frame.vibrations.raman_activities`（`angstrom^4/amu`）和
`raman_depolarization_ratios`（无量纲）。`to_raman_spectrum()` 返回以拉曼活性为谱带强度的频谱对象。
单位及拉曼活性含义遵循 [ORCA 手册](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/properties.html#raman-spectra)。

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

表中参考资料采用 ORCA 6.0/6.1 手册与教程。它们说明输入语法、物理量及打印单位；
实际读取范围以本页为准，旧版输出的表头和列序可能不同。

| 属性与参考资料 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Molecular orbitals --><!-- source:ORCA-OUTPUT --><!-- source:ORCA-METHOD -->分子轨道（[输出控制][output]、[波函数类型][hftype]） | `ORBITAL ENERGIES`; `SPIN UP/DOWN ORBITALS`; `occupation/E(Eh)/E(eV)/Irrep`; `omitted-orbital notices` | 保留打印占据、能量单位、参考/自旋身份及标签；共享受限轨道需有依据，开壳层分数总占据不等于均分自旋；缺失能量用 NaN、占据用 None。 |
| <!-- parsing-basis:Thermochemistry --><!-- source:ORCA-THERMO -->热化学（[热化学][thermo]、[热力学教程][thermo-tutorial]） | 最后一个 `THERMOCHEMISTRY AT` 区块；电子能、零点能、热能、焓、Gibbs 能及熵能项 | 总量与校正分别保存，原文 Eh/particle 转为 kcal/mol；U₀/校正相对 E(el) 推导，熵由打印的熵能项除以 T；保留温压、质量及转动信息，不将增量校正混作总校正。 |
| <!-- parsing-basis:Software metadata, frames, and status --><!-- source:ORCA-INPUT --><!-- source:ORCA-COORD --><!-- source:ORCA-RESOURCE --><!-- source:ORCA-UNITS -->程序信息、帧与状态（[输入结构][input]、[坐标][coordinates]、[并行与内存][parallel]、[单位][units]） | `Program Version`; `echoed INPUT FILE`; `coordinate/SCF/optimization boundaries`; `ORCA TERMINATED NORMALLY`; `TOTAL RUN TIME` | 分别读取程序版本、任务/资源与状态；优先打印 Å 坐标，AU 坐标按 Bohr 转换但 MASS 列仍为 amu；每个任务保留各自的帧，正常终止不自动证明全部收敛字段。 `wall_time` 及其别名 `running_time` 将打印的日/时/分/秒/毫秒转为秒；文件按选中任务段累计墙钟时间。缺少 `cpu_time` 时按墙钟时间 × 进程数回退估算，记录 `cpu_time_source="estimated"` 与 `cpu_time_num_processes`。优先运行时 MPI 进程数，其次 `%pal nprocs` / `PALn`；完整串行输入默认 1。墙钟或进程数未知时，CPU 时间仍为未知。这是参考 [cclib ORCA 实现][cpu-estimate] 的估算，不是实测 CPU 计时。 |
| <!-- parsing-basis:Energies --><!-- source:ORCA-ENERGY -->能量（[能量与梯度][energy]、[MP2][mp2]） | `FINAL SINGLE POINT ENERGY`; `E(SCF)/MPn/CC reference/correlation/total labels` | 分别映射 Eh 的最终、参考、相关总能量，保留来源观测并区分相关校正和总量；由方法身份限制字段映射。 |
| <!-- parsing-basis:Forces and geometry optimization --><!-- source:ORCA-OPT -->力与几何优化（[几何优化][optimization]、[数值梯度][gradient]） | `CARTESIAN GRADIENT`; `geometry-convergence values/limits`; `optimization-converged markers` | 力为梯度负值（Eh/Bohr），保留原子顺序；判据、阈值与完成状态分别记录；截断/缺失梯度不补零。 |
| <!-- parsing-basis:Vibrations and vibrational spectra --><!-- source:ORCA-VIB -->振动与振动光谱（[IR 与 Raman][vibrations]） | `VIBRATIONAL FREQUENCIES`; `NORMAL MODES`; `IR SPECTRUM Int/T**2`; `RAMAN SPECTRUM activity/depolarization` | 按打印模式号对齐 cm⁻¹ 频率、位移和 km/mol IR；新版读取 `Int`（km/mol），不将其 `T**2` 原子单位列当作 IR 强度；旧版按旧表头读取强度。活性为 Å⁴/amu，当前未分类偏振约定保留无后缀比值字段；保留归一化/轴信息，不完整可选表不补造。 |
| <!-- parsing-basis:Charge and spin populations --><!-- source:ORCA-POP -->电荷、自旋布居与键级（[布居与 Mayer 分析][population]） | `MULLIKEN/LOEWDIN/HIRSHFELD atomic charge/spin tables`; `Mayer bond orders larger than threshold` | 电荷/自旋布居保留方案与通道，Mayer 记录保留打印阈值；原子对按来源原子号排列，未打印对为 NaN 而非零。 |
| <!-- parsing-basis:Solvation --><!-- source:ORCA-SOLVENT -->溶剂化（[溶剂化][solvation]） | `CPCM/SMD/solvation section labels and echoed input solvent options` | `solvent` 来自已识别的模型/名称及打印元数据；输入溶剂请求与输出能量校正区分，溶剂名不等于自由能。 |
| <!-- parsing-basis:Excited states and multireference results --><!-- source:ORCA-EXCITED -->激发态与多参考结果（[激发态][excited]） | `Covered TD/ADC/EOM/STEOM/ROCIS root tables`; `electric-dipole absorption spectrum` | 以 eV 映射已覆盖根/态能量，电偶极跃迁矢量由原子单位转 Debye，并明确态身份；SOC、磁偶极、速度表及未支持 MRCI/CASSCF 结果不混入同一个电偶极态集合，多参考输入请求不证明结果已解析。 |
| <!-- parsing-basis:Response properties --><!-- source:ORCA-ELECTRIC --><!-- source:ORCA-NMR --><!-- source:ORCA-SPIN -->电响应、NMR 与总自旋（[电学性质][electric]、[NMR][nmr]、[SCF][scf]、[UHF 自旋诊断][spin]） | `DIPOLE MOMENT/raw polarizability`; `Total shielding tensor (ppm)`; `J[...] / isotropic coupling summary`; `<S**2>` | 偶极矢量由 ea₀ 转为 Debye并保留极化率张量；NMR 屏蔽为 ppm，J 为 Hz，保留原子索引及独立贡献；<S²> 为观测，S 由 sqrt(<S²>+1/4)-1/2 推导，不单凭多重度生成。 |

<!-- parsing-basis:end -->

[input]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html
[coordinates]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html
[parallel]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html
[output]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/output.html
[hftype]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hftype.html
[energy]: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/energygradients.html
[mp2]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/MP2.html
[population]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html
[scf]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html
[spin]: https://www.faccts.de/docs/orca/6.1/tutorials/prop/single_point.html
[vibrations]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html
[thermo]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html
[thermo-tutorial]: https://www.faccts.de/docs/orca/6.1/tutorials/prop/thermo.html
[nmr]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nmr.html
[electric]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/electric.html
[optimization]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html
[gradient]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalgradients.html
[solvation]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html
[excited]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/excitedstates.html
[units]: https://www.faccts.de/docs/orca/6.1/manual/contents/orcaarchitecture/conversionfactors.html

[cpu-estimate]: https://github.com/cclib/cclib/blob/master/cclib/parser/orcaparser.py
