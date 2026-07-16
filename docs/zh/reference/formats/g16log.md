# Gaussian 输出

<!-- format-support:g16log -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `g16log` |
| 扩展名 | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` |
| 读取 | 是 |
| 写入 | 否 |
| Registry 角色 | Reader |
| 数据层级 | 坐标和 QM 结果 |

## 快速读取

```python
from molop import AutoParser

frame = AutoParser("calculation.log", parser_detection="g16log")[0][-1]
print(frame.is_normal)
if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))
```

输出为终止状态和最终总能量，例如 `True` 与 `-317.592366596`；具体数值取决于文件。
其他结果入口见[按科学性质查找字段](../model_fields.md)。

MolOP 读取 Gaussian 输出文件，并提取后处理常用的结构、能量、热力学、振动、
轨道、布居、梯度、NMR、响应性质、archive 和终止状态信息。覆盖范围按计算化学用户熟悉的
Gaussian 输出内容组织；“部分支持”表示已覆盖常见或测试样例中的打印形式，但不等价于完整
Gaussian 输出语法。

| 能力 | 支持程度 | 可提取信息 | 边界 |
| ---- | -------- | ---------- | ---- |
| <!-- feature-area:Gaussian job metadata -->Gaussian 任务元数据 | 样例覆盖 | 软件版本、route 文本、title、电荷、自旋多重度，以及存在时累计的运行时间。 | 字段会被归一化；不保留原始 log 的行顺序或字节级排版。 |
| <!-- feature-area:Title card -->Title card | 样例覆盖 | Gaussian title card 会作为归一化的任务标题暴露。 | 只声明标题文本；不保留周边原始 Gaussian 格式。 |
| <!-- feature-area:Route keywords and calculation setup -->Route keywords 与计算设置 | 部分支持 | 原始 route 文本，以及可识别的模型化学、任务类型、色散、溶剂化、布居请求和 HF 设置。 | 测试外的方法专有关键字可能仍保留为 raw 文本。 |
| <!-- feature-area:Input and standard orientations -->Input/standard orientation | 样例覆盖 | Gaussian 打印的 input orientation 与 standard orientation 中的原子和坐标。 | Distance matrix 和 stoichiometry 打印块不作为独立结构化字段声明。 |
| <!-- feature-area:Rotational constants -->Rotational constants | 样例覆盖 | Gaussian 输出中的转动常数，以帧级频率量暴露。 | 不单独保留 principal-axis 诊断文本。 |
| <!-- feature-area:SCF and electronic energies -->SCF 与电子能量 | 部分支持 | SCF、reference/electronic 和方法相关能量；主 log 中的能量优先于 archive-tail 重复值。 | 只声明已测试能量字段；不保证所有 post-HF 或 correction energy 表都结构化。 |
| <!-- feature-area:Molecular orbitals and population analysis -->分子轨道与布居分析 | 部分支持 | 轨道能级/占据数，Mulliken charge/spin、APT、Lowdin、Hirshfeld/CM5、NPA、ESP charge series，electronic spatial extent、多极矩，以及存在时的早期 polarizability。 | 不解析轨道系数矩阵、NBO 轨道/键级细节，以及没有可识别逐原子表格的方案。 |
| <!-- feature-area:Vibrational frequencies and IR intensities -->振动频率与 IR 强度 | 部分支持 | 频率、约化质量、力常数、IR 强度、逐模式位移向量和虚频标记。 | Raman、VCD 等其他谱学变体只有在后续明确覆盖时才声明支持。 |
| <!-- feature-area:Thermochemistry -->热力学 | 部分支持 | 温度、压力、分子质量、惯性矩、转动对称数、转动/振动温度、转动常数、ZPVE、热能、焓、Gibbs 自由能、熵和热容。 | 覆盖主要面向 frequency 类型 Gaussian 输出和已有样例。 |
| <!-- feature-area:Dipole and polarizability -->偶极矩与 polarizability | 部分支持 | Gaussian 已覆盖响应性质区段中的 dipole 和 polarizability 值。 | 只声明样例覆盖的 response 字段；不保证所有响应性质打印变体都结构化。 |
| <!-- feature-area:NMR shielding and spin-spin coupling -->NMR 屏蔽与自旋-自旋耦合 | 部分支持 | 解析 Gaussian GIAO 的逐原子 shielding tensor、源 isotropic/anisotropy、principal values，以及 Hz 单位的 total K/J 和 FC、SD、PSO、DSO 贡献矩阵。 | 当前覆盖标准 SCF GIAO 打印族；参考化学位移、EPR 和其他 gauge 尚未结构化。 |
| <!-- feature-area:Cartesian gradients and forces -->Cartesian gradients/forces | 部分支持 | Gaussian forces 区段中的 Cartesian force 数组。 | 只声明归一化 force 数组；辅助 force diagnostics 不单独记录。 |
| <!-- feature-area:Cartesian Hessian -->Cartesian Hessian | 部分支持 | Gaussian second-derivative 区段中的 Cartesian Hessian 数据。 | 契约是归一化 Hessian 字段，不是每个 second-derivative 诊断项。 |
| <!-- feature-area:Geometry optimization convergence -->几何优化收敛 | 部分支持 | Berny 优化摘要，包括收敛阈值、force/displacement、energy change 和 optimized-state 标记。 | 测试样例外的 optimizer diagnostics 可能仍是非结构化。 |
| <!-- feature-area:Gaussian archive section -->Gaussian archive section | 部分支持 | Archive-tail 中的 metadata、坐标、能量、热力学、polarizability 和 Hessian fallback/补充字段。 | 已从主 log 解析出的字段优先；archive-tail 不会凭空生成 live status、temperature 或 pressure。 |
| <!-- feature-area:Termination status -->终止状态 | 部分支持 | Gaussian 终止证据属于 segment，并参与文件级聚合，不复制到 frame。 | Frame status 只保留 frame-local SCF 证据；结构化解析诊断之外的失败分类不做推断。 |
| <!-- feature-area:CPU and elapsed time -->CPU 与 elapsed time | 样例覆盖 | Job CPU / elapsed-time 风格记录会累计到运行时间字段。 | Per-link timing rows 不作为独立 timing record 暴露。 |
| <!-- feature-area:Link1 multi-step jobs -->Link1 多步任务 | 部分支持 | Link1 section 元数据会传播到后续帧，使多步 Gaussian 任务保留逐帧上下文。 | 低层 link 边界行不作为用户级记录暴露。 |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | 已解析 Gaussian 输出可以通过 registry 转换为坐标格式和 graph 格式。 | 转换质量取决于结构和 graph 的成功恢复。 |

## Gaussian 输出能力覆盖矩阵

`g16log` 是结构化 reader，不提供 Gaussian log writer。下表中的“fakeG”表示对应字段能否
参与 [Gaussian-like renderer](fakeg.md)，不表示能够重建原始 Gaussian 输出。

| Gaussian 输出内容 | 解析 | 结构化载体 | fakeG | 原文保真 | 当前边界 |
| ----------------- | ---- | ---------- | ----- | -------- | -------- |
| Gaussian 软件和版本 | 样例覆盖 | `qm_software`、`qm_software_version` | 支持 | 不支持 | 覆盖已有 Gaussian 版本打印形式；无法从缺失 header 的截断文件推断准确版本。 |
| Link0/options | 部分 | file/frame `options` 文本 | 支持 | 不支持 | 保存已打印的 Link0 选项；fakeG 可从 `%nprocshared` 生成 CPU 说明，但不承诺完整资源请求结构化。 |
| route 和 title | 样例覆盖 | `keywords`、`title_card`、`semantic_route` | 支持 | 不支持 | route 语义与 GJF 共用；标题和 route 均为归一化文本，不保留装饰线。 |
| method、functional、basis、任务 | 部分 | `model_chemistry`、`task_requests` 及兼容字段 | 部分 | 不适用 | 只结构化共用 route parser 已识别的关键词；结果性质不会仅因 route 请求而被伪造。 |
| charge 和 multiplicity | 支持 | frame/file 顶层字段 | 支持 | 不支持 | 从 Gaussian charge/multiplicity 行或 archive metadata 提取，主日志值优先。 |
| Link1 / 多 segment 上下文 | 部分 | segment/frame 索引和传播后的任务元数据 | 部分 | 不支持 | 后续 frame 继承适用的 section metadata；低层 Link enter/leave 不作为公共业务记录。 |
| input orientation | 样例覆盖 | `atoms`、`coords` | 支持 | 不支持 | 解析最后适用的 input orientation；不保存表格空白和 center number。 |
| standard orientation | 样例覆盖 | `standard_coords`、变换矩阵 | 支持 | 不支持 | 可计算 input/standard 刚体变换；缺少 input orientation 时可能用 standard orientation 回填坐标。 |
| 仅结构快速解析 | 支持 | atoms、coords 和基础 frame 数据 | 不适用 | 不适用 | `only_extract_structure` 在 orientation 后停止，不提取能量、频率、热化学等昂贵字段。 |
| rotational constants | 样例覆盖 | `rotation_constants` | Raw/有限 | 不支持 | 以 GHz 数组暴露；fakeG 没有完整的规范化 rotational-constant renderer。 |
| SCF/reference energy | 支持 | `energies.reference_energy`、energy observations | 支持 | 不支持 | 从 `SCF Done` 提取；fakeG 统一写成 `E(SCF)`，不保留原方法标签和 cycle 详情。 |
| MP2、MP3、MP4、MP5 | 部分 | `energies.mp*_energy` | 不支持 | 不支持 | 覆盖当前正则识别的 `EUMP*`/MP5 打印；并非所有 MP correction 和 spin-component 表。 |
| CCSD、CCSD(T) | 部分 | `energies.ccsd_energy`、`ccsd_t_energy` | 不支持 | 不支持 | 覆盖已测试总能量行；其他 coupled-cluster 变体和 correction table 不保证。 |
| archive energy fallback | 部分 | 合并后的 `energies` | 间接 | 不支持 | 主 log observation 优先，archive 只补充缺失/重复能量，不覆盖更强的 live 证据。 |
| 总自旋 | 部分 | `total_spin.spin_square`、`spin_quantum_number` | 支持 | 不支持 | 只覆盖标准 `S**2`/`S` 打印形式。 |
| MO 能级、占据和对称性 | 部分 | `molecular_orbitals` | Raw/不支持 | 不支持 | 支持 alpha/beta、开壳层和连写数值回归；不解析完整 MO coefficient matrix。 |
| Mulliken / spin population | 部分 | `charge_spin_populations` | Raw/不支持 | 不支持 | 支持独立表和开壳层合并 charge/spin 表；表头、总和及打印精度不保留。 |
| APT / Lowdin population | 部分 | `charge_spin_populations` | Raw/不支持 | 不支持 | 只覆盖已有样例的标准表格。 |
| Hirshfeld / CM5 | 部分 | `populations["hirshfeld_charges"]`、`populations["hirshfeld_spins"]`、`populations["cm5_charges"]` | Raw/不支持 | 不支持 | 提取 Hirshfeld charge/spin 和 CM5 charge series。 |
| NPA / ESP 原子电荷 | 部分 | `populations["npa_charges"]`、`populations["esp_charges"]` | Raw/不支持 | 不支持 | 保留 frame 中最后一份完整 NPA summary 和 ESP 原子电荷表；不解析详细 NBO 轨道与键级。 |
| electronic spatial extent | 部分 | `polarizability.electronic_spatial_extent` | Raw/不支持 | 不支持 | 只覆盖 population 区段中的标准 scalar 行。 |
| dipole 和高阶多极矩 | 部分 | `polarizability` 的 dipole、quadrupole、traceless quadrupole、octapole、hexadecapole | Raw/不支持 | 不支持 | 覆盖已测试 field-independent/response 打印；不保证全部单位和频率依赖变体。 |
| polarizability | 部分 | isotropic、anisotropic、tensor 字段 | Raw/不支持 | 不支持 | 合并 population、response 和 archive 来源；后出现的明确 response 数据可覆盖早期近似值。 |
| harmonic frequencies | 部分 | `vibrations.frequencies`、虚频标记 | 支持 | 不支持 | 解析连写数值；不支持 anharmonic/VPT2 结果表。 |
| reduced mass、force constant、IR | 部分 | `vibrations` 对应数组 | 支持 | 不支持 | 与频率按 mode 对齐；Raman activity、VCD/ROA 不在当前结构化合同。 |
| normal-mode displacement | 部分 | `vibrations.vibration_modes` 和 axis metadata | 支持 | 不支持 | 规范为 mode/atom/Cartesian；normalization 和 mass weighting 当前标为 `unknown`。 |
| temperature 和 pressure | 部分 | frame `temperature`、`pressure` | 支持 | 不支持 | 来自 frequency/thermochemistry 标准行；archive 不凭空补造缺失的温压。 |
| molecular mass、惯性矩和转动数据 | 部分 | `thermal_informations` | 支持 | 不支持 | 覆盖 mass、moments、symmetry number、rotational temperatures/constants。 |
| vibrational temperatures | 部分 | `thermal_informations.vibrational_temperatures` 及 mode indices | 支持 | 不支持 | 优先映射正频 mode；数量无法一致时只保留可证明的数组关系。 |
| ZPVE 和 thermal corrections | 部分 | `thermal_informations` 的 ZPVE/TCE/TCH/TCG | 支持 | 不支持 | 使用 Hartree/particle 归一化；不包含所有分项 partition-function 诊断。 |
| U0、UT、H、G、entropy、Cv | 部分 | `thermal_informations` | 支持 | 不支持 | 覆盖 Gaussian 标准 thermochemistry summary；不是完整热化学原文模型。 |
| Cartesian forces | 部分 | `forces` 和 axis/order/orientation metadata | Raw/不支持 | 不支持 | 数组形状固定为 `(N, 3)`、source atom order；orientation 当前通常为 `unknown`。 |
| Cartesian Hessian | 部分 | `hessian` 和 axis/order/orientation metadata | Raw/不支持 | 不支持 | 主 second-derivative block 优先，archive 可回填；矩阵要求 `(3N, 3N)`。 |
| Berny convergence | 部分 | `geometry_optimization_status` | 支持 | 不支持 | 提取 force/displacement、threshold、energy change 和 optimized 标记；其他 optimizer 形式可能未结构化。 |
| SCF 与 termination status | 部分 | frame-local `status`、segment/file 聚合 status | 合成 | 不支持 | SCF 证据属于 frame；normal/error termination 属于 segment/file。fakeG 终止行不是原始状态复现。 |
| CPU / elapsed time | 样例覆盖 | 累计 `running_time` | 部分 | 不支持 | file 汇总 segment 时间；不暴露逐 Link timing records，多帧 fakeG 会移除 frame runtime 行。 |
| archive tail metadata/coords | 部分 | 缺失字段 fallback | 间接 | 不支持 | archive 可补充 metadata、坐标、能量、热化学、polarizability 和 Hessian；live 数据优先。 |
| TD excited-state 结果 | 不支持 | 无稳定公共结果合同 | 不支持 | 不支持 | route 可识别 TD 请求，但 excitation energies、oscillator strengths 等结果尚未结构化。 |
| NMR magnetic shielding | 部分 | `nmr.gauge`、逐原子 `shielding_tensors`、isotropic、anisotropy、principal values | 不支持 | 不支持 | 覆盖标准 SCF GIAO `(ppm)` 打印；张量按 source atom index 对齐，Cartesian orientation 当前标记为 `unknown`。逐原子 shielding 可通过 `qm_embedded_rdmol(embed_nmr=True)` 写入 RDKit atom properties。 |
| NMR spin-spin coupling | 部分 | `spin_spin_coupling_k`、`spin_spin_coupling_j`、`spin_spin_coupling_k_components`、`spin_spin_coupling_j_components`、`coupling_atom_indices` | 不支持 | 不支持 | 将分块下三角 total 和 FC/SD/PSO/DSO K/J 表展开为 Hz 单位对称方阵；coupling 矩阵仍保留为 frame-level 原子对数据。 |
| EPR、参考化学位移、NBO bond order 等专用结果 | 不支持 | 无稳定公共结果合同 | 不支持 | 不支持 | route 请求不代表结果字段已覆盖；需要独立样例和公共语义合同。 |
| 坐标/graph 格式转换 | 支持 | frame structure / recovered graph | 不适用 | 不适用 | 坐标转换依赖成功提取结构；graph 转换还依赖键级和电荷状态恢复质量。 |
| 原始 log 无损重写 | 不支持 | component raw snippets 仅供检查/fallback | 不支持 | 不支持 | parser 是字段提取器，不是完整 Gaussian 输出 CST；不提供 byte-for-byte writer。 |

`capture_source_evidence=True` 会保留能量、SCF 和优化收敛等来源证据，适合审计解析结论；
默认模式更偏向结果提取。无论是否捕获证据，`g16log` 都不承诺保存完整原始输出结构。
