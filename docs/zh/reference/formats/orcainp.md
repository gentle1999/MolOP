# ORCA 输入

<!-- format-support:orcainp -->

| 项目 | 值 |
| ---- | -- |
| 格式 ID | `orcainp` |
| 扩展名 | `.inp` |
| 读取 | 是 |
| 写入 | 是 |
| Registry 角色 | Reader、file writer、frame writer |
| 数据层级 | 坐标和 QM 输入语义 |

读取 ORCA 坐标、多步骤任务、计算方法及结构化输入请求。

带坐标的 frame 可以通过显式 ORCA 计算设置生成输入文件：

下载[共享 ORCA 样例](../../../assets/examples/water_mp2.out)后运行：

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
rendered = frame.format_transform(
    "orcainp",
    keywords="wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP",
    nprocs=16,
    maxcore=4000,
    blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
)
```

`rendered` 是包含 ORCA simple input line、`%pal/%maxcore/%scf` 和 `* xyz` 坐标块的字符串。

??? example "输出"
    ```text
    ! wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP

    %pal
      nprocs 16
    end

    %maxcore 4000

    %scf
      MaxIter 300
      STABPerform true
    end

    * xyz 0 1
    O      1.7849140000     1.2624220000     0.5119850000
    H      2.6482370000     1.0729290000     0.1316310000
    H      1.1831680000     1.2568160000    -0.2388350000
    *
    ```

跨格式转换必须提供 `keywords`；MolOP 不会复用其他 QM 软件的关键字语法。
`blocks` 也接受原始 `%block` 文本。

| 特性 | 支持程度 | 支持范围 | 明确边界 |
| ---- | -------- | -------- | -------- |
| <!-- feature-area:Geometry and job splitting -->几何与多任务 | 部分支持 | 直接坐标、点电荷、外部 `xyzfile`/`pdbfile` 引用，以及 `$new_job` 分割。 | 外部几何作为引用解析；解析器不要求读取外部文件。 |
| <!-- feature-area:Model chemistry and options -->计算方法与选项 | 部分支持 | Method、functional、basis、辅助基组、色散作为 functional 后缀、混合基组、print 设置、PARAS 变量和扫描坐标。 | 未支持的 ORCA block 选项可能保留在 raw/resource 字段，而不是专门语义容器中。 |
| <!-- feature-area:ORCA 6 manual fixture families -->ORCA 6 输入语法 | 部分支持 | ORCA 6 手册中的单点、SCF stability、优化、频率、激发态、MRCI 和 Solvator 输入示例。 | 完整坐标需要内联几何声明；外部几何保持路径引用。 |
| <!-- feature-area:Excited-state and multireference requests -->激发态与多参考请求 | 部分支持 | 结构化激发态和多参考任务语义，包括 MRCI 多步骤任务。 | 覆盖示例之外的 ORCA 激发态和多参考关键字族可能仍为 raw。 |
| <!-- feature-area:Optimization, coordinates, frequency, and Solvator structures -->优化、频率与 Solvator 设置 | 部分支持 | 优化约束、内坐标、fragments、NEB、频率 restart 和 Solvator 示例。 | 只解析下表列出的语法与选项，未识别的内容保留原文。 |
| <!-- feature-area:Writer availability -->输入文件写出 | 部分支持 | 从带坐标的模型规范化渲染 file/frame，支持显式 keywords、PAL、maxcore、任意 `%block`、笛卡尔坐标和 `.inp` 输出路径。 | 尚不规范化非笛卡尔内联几何和任意原文语句顺序；跨格式转换必须显式提供 ORCA keywords。 |

## ORCA 6.1 能力覆盖矩阵

下表区分四种能力，避免把“能够保留 raw 文本”误认为“已经结构化支持”：

- **支持**：可以读取或规范化写出所列语法。
- **部分**：只覆盖所列形式，不能承诺覆盖 ORCA 6.1 的全部变体。
- **Raw**：可以通过原始 `%block` 文本传递，但 MolOP 不解释或校验其中语义。
- **不支持**：当前会拒绝、丢失信息，或者不能从模型重新生成。

| ORCA 输入能力 | 解析 | 结构化语义 | 规范化渲染 | 原文保真 | 当前边界 |
| ------------- | ---- | ---------- | ---------- | -------- | -------- |
| `!` 简单关键词行 | 支持 | 部分 | 支持 | 不支持 | 可保留关键词文本；method、functional 和 basis 仅识别内置词表，输出时重新排版。 |
| 顶层 `#` 注释 | 部分 | 不适用 | 支持 | 不支持 | 保存独立注释行；第二个 `#` 或行尾结束注释，闭合注释后的代码继续解析，带引号文件名中的 `#` 保持字面值；不保留注释的全局位置。 |
| `%maxcore`、`%pal` | 支持 | 支持 | 支持 | 部分 | 支持 `maxcore`、`nprocs` 结构化覆盖；输出采用规范格式。 |
| 普通非嵌套 `%block ... end` | 支持 | 部分 | Raw | 部分 | 未建模的 block 可保存在 `ORCABlock.raw_text`，但不保证其与其他语句的原始相对位置。 |
| 单行 `%` 指令，如 `%moinp`、`%base` | 部分 | 不支持 | Raw | 部分 | 通用 block 扫描可以捕获；除资源字段外通常没有专用模型或参数校验。 |
| 嵌套 block，如 `%scf/SOSCF`、`%basis/NewGTO` | 部分 | 部分 | 部分 | 不支持 | 嵌套识别目前依赖有限白名单；未识别的内层 `end` 可能截断外层 block。 |
| 重复 block 和输入优先级 | 部分 | 不支持 | 部分 | 不支持 | 解析保留 block 列表顺序；按名称覆盖 block 时会移除旧实例并把新实例追加到末尾。 |
| 未分类顶层语句 | 部分 | 不支持 | 部分 | 不支持 | 存入 `trailing_lines`，渲染时统一移动到几何之后。 |
| `* xyz`、`* cart`、`* cartesian` | 支持 | 支持 | 支持 | 不支持 | 坐标规范化为 `* xyz`；空白、有效数字和原始坐标类型不保留。 |
| `%coords` Cartesian 几何 | 支持 | 支持 | 支持 | 不支持 | 解析后规范化为 `* xyz`，不重新生成 `%coords`。 |
| 点电荷 `Q q x y z` | 支持 | 支持 | 支持 | 不支持 | 输出电荷使用固定六位小数，坐标使用 writer 的精度设置。 |
| ghost、dummy、fragment、frozen、原子级基组 | 部分 | 部分 | 部分 | 不支持 | 支持下表列出的形式；原始 token、引号和排列不保证保留。 |
| isotope `M=`、nuclear charge `Z=` | 支持 | 支持 | 不支持 | 不支持 | parser 写入原子模型，但 Cartesian writer 尚未输出这两个属性。 |
| `* int`、`* internal`、`* gzmt` | 部分 | 部分 | 不支持 | 不支持 | 可解析已覆盖的数值内坐标并投影 Cartesian 坐标；writer 会抛出 `NotImplementedError`。 |
| `* xyzfile`、`* gzmtfile`、`* pdbfile` | 支持 | 支持 | 支持 | 部分 | 只保存并重新输出路径引用，不读取或验证外部文件。 |
| `%paras` 和坐标参数表达式 | 部分 | 部分 | 部分 | 不支持 | 支持简单参数、范围和 `name +/- offset`；几何渲染时表达式会被解析后的数值替代。 |
| `$new_job` 多任务 | 支持 | 支持 | 支持 | 不支持 | 每个 job 映射为一个 frame；输出统一使用 `$new_job` 分隔符。 |
| `%Compound` 工作流 | 部分 | 不支持 | 不支持 | 不支持 | 无顶层几何的合法 Compound-only 输入会被格式检测拒绝，writer 也要求关键词和几何。 |
| 模型化学与常见任务 | 部分 | 部分 | 部分 | 不适用 | 覆盖常见 SP、OPT、TS/NEB、FREQ、IRC、gradient、DFT、部分波函数方法和色散关键词。 |
| 激发态、多参考、Solvator | 部分 | 部分 | Raw | 部分 | parser 提供结构化请求；writer 依赖已解析的 raw block 或调用方显式传入 block，不从请求模型完整反向生成。 |
| 完整 `%basis`、ECP、AuxJ/AuxJK/AuxC/CABS | 部分 | 不支持 | Raw | 部分 | 简单关键词和原子级覆盖可识别；完整 block 内容不做字段级建模与校验。 |
| QMMM、MD、GOAT、DOCKER、RESP、Raman 等 ORCA 6.1 功能 | 部分 | 不支持 | Raw | 部分 | 非嵌套 raw block 可以传递；不承诺嵌套语法、任务兼容性或字段级构建。 |
| ORCA 版本和兼容性校验 | 不支持 | 不支持 | 不支持 | 不适用 | 当前版本记为 `Any`，不检查关键词、方法、基组和 ORCA 版本之间的兼容性。 |

当前 writer 是 **canonical writer**，适合从坐标和显式计算设置构建新输入；它不是
**source-preserving writer**。需要读取现有复杂输入、只修改一个字段并保持其他文本不变时，
不应依赖当前 round-trip 行为。

参考：[ORCA 6.1 输入结构](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html)、
[坐标](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html)、
[基组](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html)和
[Compound](https://www.faccts.de/docs/orca/6.1/manual/contents/workflowsautomatization/compound.html)。

## 属性来源与读取约定

<!-- parsing-basis:start -->

下表说明文件中的信息如何对应 MolOP 数据。公共单位、编号和缺失值说明见
[数据容器与跨 QM 语义](../data_model.md)。

表中参考资料采用 ORCA 6.0/6.1 手册与教程。它们说明输入语法、物理量及打印单位；
实际读取范围以本页为准，旧版输出的表头和列序可能不同。

| 属性与参考资料 | 文件标记或数据来源 | 读取或写出约定 |
| ---- | ------------------ | -------------- |
| <!-- parsing-basis:Geometry and job splitting --><!-- source:ORCA-COORD --><!-- source:ORCA-INPUT --><!-- source:ORCA-UNITS -->几何与多任务（[坐标][coordinates]、[输入结构][input]、[单位][units]） | `* xyz/internal/xyzfile`; `%coords`; `$new_job`; `charge/multiplicity`; `unit declarations` | 可解的内联几何映射为原子/Å 坐标；内坐标与参数保留类型化输入语义；外部路径保持引用，任务分隔符确定来源区间。 |
| <!-- parsing-basis:Model chemistry and options --><!-- source:ORCA-METHOD --><!-- source:ORCA-RESOURCE --><!-- source:ORCA-OUTPUT --><!-- source:ORCA-SOLVENT -->计算方法与选项（[波函数类型][hftype]、[并行与内存][parallel]、[输出控制][output]、[溶剂化][solvation]） | `! keyword lines`; `%scf/%basis/%pal`; `%maxcore`; `%output`; `PARAS variables` | 按有限语义映射方法、自旋、基组/辅助基、色散、溶剂、资源和打印请求；maxcore 为每进程内存，未知块选项保留原文，不补造结果字段。 |
| <!-- parsing-basis:ORCA 6 manual fixture families -->ORCA 6 输入语法（[输入结构][input]） | `!` 关键字、`%...end` 输入块与坐标声明 | 读取已支持的关键字、输入块和几何语法；输入文件仅描述计算请求，能量、光谱和收敛状态需从输出文件读取。 |
| <!-- parsing-basis:Excited-state and multireference requests --><!-- source:ORCA-EXCITED -->激发态与多参考请求（[激发态][excited]、[多参考输入][mrci]） | `%tddft/%cis/%rocis/%casscf/%mrci and covered method keywords` | 在已映射范围内保存激发态/多参考请求的根数、态、参考及活性空间；根号属于请求标签，不是已观测本征态。 |
| <!-- parsing-basis:Optimization, coordinates, frequency, and Solvator structures --><!-- source:ORCA-OPT --><!-- source:ORCA-VIB -->优化、频率与 Solvator 设置（[几何优化][optimization]、[数值梯度][gradient]、[IR 与 Raman][vibrations]） | `%geom Constraints/Scan/Fragments`; `NEB`; `%freq restart`; `%solvator` | 保留支持的任务、约束、扫描和环境输入及嵌套块边界；原始坐标与请求位移区分，未识别的选项保留原文。 |
| <!-- parsing-basis:Writer availability -->输入文件写出（[输入结构][input]、[坐标][coordinates]、[并行与内存][parallel]） | `Common geometry/model fields`; `explicit ! keywords`; `PAL/maxcore and %blocks` | 生成规范笛卡尔 ORCA 输入；跨格式关键词需显式给定，不将 Gaussian 关键词静默解释为 ORCA 语义，不恢复任意原文顺序。 |

<!-- parsing-basis:end -->

[input]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html
[coordinates]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html
[parallel]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html
[output]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/output.html
[hftype]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hftype.html
[mrci]: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/mrci.html
[vibrations]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html
[optimization]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html
[gradient]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalgradients.html
[solvation]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html
[excited]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/excitedstates.html
[units]: https://www.faccts.de/docs/orca/6.1/manual/contents/orcaarchitecture/conversionfactors.html
