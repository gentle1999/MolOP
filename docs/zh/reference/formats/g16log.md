<!--
 * @Author: TMJ
 * @Date: 2026-07-09 21:59:43
 * @LastEditors: TMJ
 * @LastEditTime: 2026-07-09 22:29:14
 * @Description: 请填写简介
-->
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

MolOP 读取 Gaussian 输出文件，并提取后处理常用的结构、能量、热力学、振动、
轨道、布居、梯度、响应性质、archive 和终止状态信息。覆盖范围按计算化学用户熟悉的
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
| <!-- feature-area:Molecular orbitals and population analysis -->分子轨道与布居分析 | 部分支持 | 轨道能级/占据数、电荷/自旋布居、electronic spatial extent、多极矩，以及存在时的早期 polarizability。 | 不声明轨道系数矩阵或每一种 population analysis 变体都已覆盖。 |
| <!-- feature-area:Vibrational frequencies and IR intensities -->振动频率与 IR 强度 | 部分支持 | 频率、约化质量、力常数、IR 强度、逐模式位移向量和虚频标记。 | Raman、VCD 等其他谱学变体只有在后续明确覆盖时才声明支持。 |
| <!-- feature-area:Thermochemistry -->热力学 | 部分支持 | 温度、压力、分子质量、惯性矩、转动对称数、转动/振动温度、转动常数、ZPVE、热能、焓、Gibbs 自由能、熵和热容。 | 覆盖主要面向 frequency 类型 Gaussian 输出和已有样例。 |
| <!-- feature-area:Dipole and polarizability -->偶极矩与 polarizability | 部分支持 | Gaussian 已覆盖响应性质区段中的 dipole 和 polarizability 值。 | 只声明样例覆盖的 response 字段；不保证所有响应性质打印变体都结构化。 |
| <!-- feature-area:Cartesian gradients and forces -->Cartesian gradients/forces | 部分支持 | Gaussian forces 区段中的 Cartesian force 数组。 | 只声明归一化 force 数组；辅助 force diagnostics 不单独记录。 |
| <!-- feature-area:Cartesian Hessian -->Cartesian Hessian | 部分支持 | Gaussian second-derivative 区段中的 Cartesian Hessian 数据。 | 契约是归一化 Hessian 字段，不是每个 second-derivative 诊断项。 |
| <!-- feature-area:Geometry optimization convergence -->几何优化收敛 | 部分支持 | Berny 优化摘要，包括收敛阈值、force/displacement、energy change 和 optimized-state 标记。 | 测试样例外的 optimizer diagnostics 可能仍是非结构化。 |
| <!-- feature-area:Gaussian archive section -->Gaussian archive section | 部分支持 | Archive-tail 中的 metadata、坐标、能量、热力学、polarizability 和 Hessian fallback/补充字段。 | 已从主 log 解析出的字段优先；archive-tail 不会凭空生成 live status、temperature 或 pressure。 |
| <!-- feature-area:Termination status -->终止状态 | 部分支持 | Gaussian 终止证据属于 segment，并参与文件级聚合，不复制到 frame。 | Frame status 只保留 frame-local SCF 证据；结构化解析诊断之外的失败分类不做推断。 |
| <!-- feature-area:CPU and elapsed time -->CPU 与 elapsed time | 样例覆盖 | Job CPU / elapsed-time 风格记录会累计到运行时间字段。 | Per-link timing rows 不作为独立 timing record 暴露。 |
| <!-- feature-area:Link1 multi-step jobs -->Link1 多步任务 | 部分支持 | Link1 section 元数据会传播到后续帧，使多步 Gaussian 任务保留逐帧上下文。 | 低层 link 边界行不作为用户级记录暴露。 |
| <!-- feature-area:Registry conversion -->Registry conversion | 已支持 | 已解析 Gaussian 输出可以通过 registry 转换为坐标格式和 graph 格式。 | 转换质量取决于结构和 graph 的成功恢复。 |
