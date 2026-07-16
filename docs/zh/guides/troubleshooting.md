# 常见问题

从安装、格式检测、字段缺失和导出四个方向定位常见问题。

## `AutoParser` 返回空 batch

先确认路径实际匹配文件：

```python
from pathlib import Path

paths = list(Path("results").glob("*.out"))
print(len(paths), paths[:3])
```

若文件存在，单进程并显式指定格式重试：

```python
from molop import AutoParser

batch = AutoParser(
    "results/job.out",
    parser_detection="orcaout",
    n_jobs=1,
)
print(len(batch))
```

预期成功时输出 `1`。仍为 `0` 时查看终端和 `molop.log` 中的格式不支持、内容探测或解析错误。

## 文件扩展名与内容不一致

`.out` 可能来自 ORCA 或 xTB，`.log` 也可能不是 Gaussian。使用对应格式 ID：

```python
AutoParser("orca.out", parser_detection="orcaout")
AutoParser("xtb.out", parser_detection="xtbout")
AutoParser("gaussian.log", parser_detection="g16log")
```

如果显式格式仍失败，先确认文件不是被截断的输入、调度器日志或空文件。

## 结果字段是 `None`

`None` 通常表示源文件没有该结果或当前打印形式尚未结构化，不表示数值为零。

```python
frame = batch[0][-1]
if frame.vibrations is None:
    print("该 frame 没有结构化频率数据")
```

检查事项：

- 是否取到了真正包含结果的 frame。
- 计算是否请求并打印了该性质。
- 对应格式页是否声明支持该字段。
- 是否误用了 `only_extract_structure=True`。

## `is_normal` 或 `is_optimized` 是 `None`

MolOP 保留“未知”状态。输入文件、截断输出或没有终止证据的片段可能无法给出真假判断。批量
决策时单独记录未知值，不要用 `bool(None)` 把它静默变成失败。

## 输出目录错误

Python API 的 `output_dir` 必须存在：

```python
from pathlib import Path

Path("output").mkdir(parents=True, exist_ok=True)
batch.format_transform(
    "xyz", output_dir="output", write_to_disk=True
)
```

CLI 的 `--output-dir` 会创建目录。Python API 只有 `write_to_disk=True` 才写盘。

## SDF/SMILES 转换失败

图级格式需要可用拓扑：

```python
frame = batch[0][-1]
print(frame.rdmol)
print(frame.topology_reconstruction_status)
```

若输出为 `None` 和 `failed`，先核对电荷、多重度和几何，再参考[结构恢复](structure-recovery.md)。

## 单位如何转换

不要直接读取内部 magnitude 后猜单位：

```python
energy = frame.energies.total_energy
print(energy.m_as("hartree"))
print(energy.m_as("eV"))
```

`hartree` 到 `eV` 是单粒子能量转换。若需要 `kcal/mol`，必须显式引入每摩尔语义，不能直接
对普通电子能 quantity 调用 `.m_as("kcal/mol")`。

## 怎样提交可复现问题

提供以下信息：

- `molop --version` 与 `python --version`。
- 最小输入文件；敏感文件可裁剪，但需保留触发问题的完整区段。
- 使用的代码或完整 CLI 命令。
- `parser_detection`、`n_jobs` 和错误文本。
- 预期字段及对应原始输出行。

在 [GitHub Issues](https://github.com/gentle1999/MolOP/issues) 提交问题。
