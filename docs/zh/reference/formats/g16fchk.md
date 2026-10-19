# Gaussian 格式化检查点

<!-- format-support:g16fchk -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `g16fchk` |
| 扩展名 | `.fchk`, `.fch`, `.fck` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | 单 frame Gaussian QM 结果和坐标 |

## 快速读取

```python
from molop import AutoParser

frame = AutoParser("molecule.fchk", parser_detection="g16fchk")[0][0]
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy if frame.energies else None)
print(frame.charge_spin_populations.population_names if frame.charge_spin_populations else [])
```

??? example "输出约定"

    输出依次给出原子/坐标规模、可用总能量和实际存在的布居方案。

fchk 只产生一个当前几何 frame，不包含 log 中的完整优化轨迹。

MolOP 读取 Gaussian `formchk` 生成的文本格式化检查点。解析器按 fchk 固定宽度
record 和 `N=` 数组长度解码，不依赖字段相邻顺序；未知 record 会被跳过。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:File recognition and fixture inventory -->文本记录与文件识别 | 部分支持 | 读取固定宽度的标量和数组记录，生成当前几何的一个结果帧。 | 不接受二进制 `.chk`；未映射的记录不生成结构化属性。 |
| <!-- feature-area:Metadata and route semantics -->元数据与 Route 设置 | 部分支持 | 解析标题、Gaussian 版本、Route、电荷、多重度、模型化学、任务、色散、溶剂和开壳层处理；第二行摘要用于补齐 Route 中未识别的 method/basis。 | Route 仍受共享 Gaussian 语法覆盖范围限制；fchk 不含完整 Link 0 输入和运行时间。 |
| <!-- feature-area:Geometry and Cartesian derivatives -->几何与笛卡尔导数 | 部分支持 | 解析原子序数、Bohr 坐标、由 Cartesian gradient 取负得到的 forces，以及由 packed force constants 展开的对称 Hessian。 | fchk 只表示当前几何，不提供完整优化轨迹；坐标和导数保持 source order，不反演原输入方向变换。 |
| <!-- feature-area:Energies, orbitals, and spin -->能量、轨道与自旋 | 部分支持 | 解析 SCF/reference、MP2、MP3、MP4、CCSD、CCSD(T) 和最终总能量，以及 alpha/beta 轨道能、占据数和 S-squared。 | 不暴露 MO coefficient、密度矩阵、自然轨道和未映射的方法专用能量标签。 |
| <!-- feature-area:Populations and electric response -->布居与电响应 | 部分支持 | 解析与原子数一致的 Mulliken、NPA、ESP、APT、Hirshfeld/CM5 及可用 spin population records，并解析偶极矩、packed polarizability 和四极矩。 | 只暴露 fchk 实际包含的 record；键级、超精细张量和更高阶电响应数组尚未映射。 |
| <!-- feature-area:Vibrations, thermochemistry, and status -->振动、热化学与状态 | 部分支持 | 解析频率、约化质量、力常数、IR 强度、简正模式、thermal energy/enthalpy/free energy、Job Status 和优化完成状态。 | 完整极化率导数记录可支持拉曼活性及 P/U 退偏振比；TCE/TCH/TCG 保留推导式，温度/压力仅在打印时保留。ROA/VCD 和收敛历史暂不映射。 |
| <!-- feature-area:NMR shielding and spin-spin coupling -->NMR 屏蔽与自旋耦合 | 部分支持 | 解析逐原子 shielding tensor 及 isotropic、anisotropy、principal values，FC、SD、PSO、DSO 约化耦合矩阵，以及 Hz 单位的 total K 矩阵。 | 需要对应 NMR record；fchk 缺少重建 J 所需的同位素信息，也不映射参考化学位移和 EPR 数据。 |

## 能力覆盖矩阵

| fchk 能力 | 读取 | 结构化结果 | 当前边界 |
| --------- | ---- | ---------- | -------- |
| 固定宽度 scalar record | 支持 | `I`、`R`、`C`、`L`、`H` 值 | 未知 label 会跳过。 |
| `N=` 数组 record | 支持 | 按声明长度读取数值或 12 字符块 | 截断数组会报告格式不匹配。 |
| 标题、Route 与 Gaussian 版本 | 支持 | `title_card`、`keywords`、`qm_software_version` | 不保留 record 的原始排版。 |
| 电荷与多重度 | 支持 | `charge`、`multiplicity` | 使用 fchk 标量值。 |
| method、functional、basis 与任务 | 部分 | `model_chemistry`、`task_requests` | Route 优先，第二行摘要补齐缺失 method/basis；未知 Route token 仍保留在 diagnostics。 |
| 文件与 frame 数量 | 支持 | 一个 segment、一个 terminal/SP frame | fchk 不表示 Gaussian log 的多 Link/优化步时间序列。 |
| 原子和当前 Cartesian 坐标 | 支持 | Angstrom `coords`、source-order `atoms` | fchk Bohr 坐标会转换到 Angstrom。 |
| Cartesian gradient | 支持 | Hartree/Bohr `forces` | forces 是 gradient 的逐元素负值。 |
| Cartesian force constants | 支持 | `(3N, 3N)` 对称 `hessian` | 只接受合法 packed lower-triangle 长度。 |
| SCF、MP2-4、CCSD、CCSD(T)、总能量 | 支持 | `Energies` 和可选来源证据 | 未映射的能量 label 不自动猜测语义。 |
| alpha/beta 轨道能与占据 | 支持 | `MolecularOrbitals` | 闭壳层缺少 beta 数组时复制 alpha 能量，并按 beta electron 数生成占据。 |
| MO coefficients 与密度矩阵 | 不支持 | 无 | 大型原始数组会安全跳过。 |
| 原子布居 records | 部分 | `ChargeSpinPopulations`，含可扩展 ESP/NPA spin series | 仅当数组长度与原子数一致时映射 Mulliken、NPA、ESP、APT、Hirshfeld/CM5 和已知 spin labels。 |
| S-squared | 支持 | `TotalSpin` | spin quantum number 由 S-squared 推导。 |
| dipole、polarizability、quadrupole | 部分 | `Polarizability` | 保留 fchk packed 顺序；不展开完整张量命名。 |
| 频率、质量、力常数与 IR | 支持 | `Vibrations` | 从 `Vib-E2` 的已知前四个 mode 分组读取。 |
| 简正模式位移 | 支持 | 每个 mode 的 `(N, 3)` 数组 | 标记为 source-program normalization，mass weighting 未知。 |
| thermal totals | 部分 | `U_T`、`H_T`、`G_T` | 不提供完整 ZPVE、熵和热容分解。 |
| NMR magnetic shielding | 部分 | `nmr.gauge`、逐原子 `shielding_tensors`、isotropic、anisotropy、principal values | 读取 `NMR shielding`；张量单位转换为 ppm，Cartesian orientation 当前标记为 `unknown`。逐原子 shielding 可通过 `qm_embedded_rdmol(embed_nmr=True)` 写入 RDKit atom properties。 |
| NMR reduced spin-spin coupling | 部分 | `spin_spin_coupling_k`、`spin_spin_coupling_k_components`、`coupling_atom_indices` | 读取四组 packed FC/SD/PSO/DSO K 贡献并求和得到 Hz 单位 total K；fchk 缺少同位素信息，因此 `spin_spin_coupling_j` 和 J 分量保持为空；coupling 矩阵仍保留为 frame-level 原子对数据。 |
| Job Status 与优化完成 | 部分 | `Status`、`GeometryOptimizationStatus` | 优化完成由 opt 请求和正常 Job Status 共同推断，不含逐步收敛值。 |
| source span 与哈希 | 支持 | 单一全文 segment/frame provenance | 不合并 `.chk`、log 或其他 sidecar。 |
| basis、ECP、MO/density 原始数组 | 可跳过 | 无专用公共字段 | 能扫描不代表已结构化支持。 |
| fchk/chk writer | 不支持 | Reader only | 不从 MolOP 模型重建 fchk，也不处理二进制 checkpoint。 |

## 格式边界

`.fchk` 是 Gaussian formatted checkpoint 文本，不是二进制 `.chk`。解析器不会调用
Gaussian `formchk`，也不会从同目录读取 log 或 checkpoint；所有 provenance 只指向传入的
fchk 文本。

完整极化率导数记录可确认 Raman 数据存在，此时提取拉曼活性及独立 P/U 退偏振比。
热校正 `TCE`、`TCH`、`TCG` 的推导式保存在 `derivations` 中。
容器语义与数据约定见 [公共数据模型](../data_model.md)。

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

参考资料：[Gaussian 格式化检查点](https://gaussian.com/formchk/)。

| 属性 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:File recognition and fixture inventory -->文本记录与文件识别 | 固定宽度标量/数组记录；`I/R/C/L/H` 类型；`N=` 数量；文件头 | 按声明类型和数量读取记录，支持 D 指数；仅已映射标签获得语义，不读取二进制 .chk。 |
| <!-- parsing-basis:Metadata and route semantics -->元数据与 Route 设置 | 标题及第二行的任务/方法/基组摘要；`Route`、`Gaussian Version`、`Charge`、`Multiplicity` | 提取元数据与输入请求，第二行摘要补充 route 缺失信息；不从 fchk 推断运行时间或完整 Link0。 |
| <!-- parsing-basis:Geometry and Cartesian derivatives -->几何与笛卡尔导数 | `Atomic numbers`; `Real atomic weights`; `Current cartesian coordinates`; `Cartesian Gradient/Force Constants` | 坐标 Bohr→Å，质量 amu；力为梯度负值（Eh/Bohr）；下三角 Hessian 对称展开为 Eh/Bohr²，轴顺序保持来源原子序。 |
| <!-- parsing-basis:Energies, orbitals, and spin -->能量、轨道与自旋 | `SCF/MPn/Cluster/Total Energy`; `Alpha/Beta Orbital Energies`; `electron counts`; `S**2` | 能量按标签分别保存为 Eh；占据由自旋电子数推导，仅有受限参考依据时共享缺失 beta 谱；S 由 sqrt(<S²>+1/4)-1/2 推导。 |
| <!-- parsing-basis:Populations and electric response -->布居与电响应 | `Mulliken/NPA/ESP/APT/Hirshfeld/CM5 charge/spin arrays`; `Dipole Moment`; `Polarizability`; `Quadrupole Moment` | 布居数组按原子数对齐并保留方案及自旋标签；偶极 ea₀、四极矩 ea₀² 分别转为 Debye、Debye·Å，极化率采用 a₀³ 并保持紧凑排列；不将其解释为键级矩阵。 |
| <!-- parsing-basis:Vibrations, thermochemistry, and status -->振动、热化学与状态 | `Number of Normal Modes`; `Vib-E2/Vib-Modes`; `Polarizability Derivatives`; `Thermal Energy/Enthalpy/Free Energy`; `Temperature/Pressure`; `Job Status` | Vib-E2 第 0–3 块为频率（cm⁻¹）、质量（amu）、力常数（mdyne/Å）、IR（km/mol）；4–6 块为活性（Å⁴/amu）、P/U，要求完整极化率导数。振型按（mode,atom,xyz）重排。热量原文为 Eh/particle，默认转换为 kcal/mol，TCE/TCH/TCG 为总量减 Total Energy 并记录推导式；状态依据标签，不重建优化历史。 |
| <!-- parsing-basis:NMR shielding and spin-spin coupling -->NMR 屏蔽与自旋耦合 | `NMR shielding`; `Isotropic Spin-Spin Couplings`; `NMR route gauge` | 屏蔽张量转置后用 α²×10⁶ 转为 ppm；各向同性、各向异性及主值由张量推导。FC/SD/PSO/DSO 四个紧凑块映射 K（Hz）并求和；缺少同位素依据时不补造 J，原子索引从 0 开始。 |

<!-- parsing-basis:end -->
