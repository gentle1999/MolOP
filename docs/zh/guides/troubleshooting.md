# 常见问题

先选择一个有代表性的文件并设置 `n_jobs=1`。这样可以保持解析诊断的顺序，也能区分格式问题
与并行执行问题。

## `AutoParser` 返回空 batch

先确认路径或 glob 确实对应普通文件，再指定预期格式 ID 重试一次：

```python
from pathlib import Path

from molop import AutoParser

path = Path("results/job.out")
print(path.is_file())

batch = AutoParser(path, parser_detection="orcaout", n_jobs=1)
print(len(batch), batch[0].detected_format_id if batch else None)
```

??? example "有效 ORCA 输出的预期结果"
    ```text
    True
    1 orcaout
    ```

glob 没有匹配、文件不存在、格式 ID 未注册或内容探测失败时，不会产生已解析文件，并会输出诊断
信息。显式指定格式后仍然失败，应检查文件是否为空、被截断，或实际是输入文件、调度器日志而非
计算输出。

## 扩展名与内容不一致

`.out` 和 `.log` 等扩展名并不唯一。排查时可使用规范格式 ID 选择 reader：

```python
orca = AutoParser("orca.out", parser_detection="orcaout", n_jobs=1)
xtb = AutoParser("xtb.out", parser_detection="xtbout", n_jobs=1)
gaussian = AutoParser("gaussian.log", parser_detection="g16log", n_jobs=1)
```

[格式支持概览](../reference/format_support.md)列出了已注册的 ID 和扩展名。显式检测只会缩小
候选 reader 范围，不会让不完整的内容变成有效输入。

## 结果字段是 `None`

`None` 表示所选 frame 没有该性质的结构化值，不表示数值为零。应核对 frame、计算请求、源文本
以及对应格式页：

```python
frame = batch[0][-1]
if frame.vibrations is None:
    print("当前 frame 没有结构化频率数据")
```

??? example "缺少频率数据时的输出"
    ```text
    当前 frame 没有结构化频率数据
    ```

`only_extract_structure=True` 会按约定跳过许多非结构字段。输入文件、截断输出和缺少终止证据的
任务也可能令 `is_normal` 或 `is_optimized` 为 `None`；应保留这一未知状态，不要用
`bool(...)` 将其强制转换。

## 写入输出目录失败

批量 Python API 要求目录已经存在，并且只有 `write_to_disk=True` 才会写盘：

```python
from pathlib import Path

Path("structures").mkdir(parents=True, exist_ok=True)
batch.format_transform(
    "xyz",
    output_dir="structures",
    write_to_disk=True,
)
```

??? example "生成文件"
    ```text
    structures/job.xyz
    ```

CLI 会按需创建 `--output-dir`。`write_to_disk=False` 时，Python 返回渲染文本并忽略
`output_dir`。

## SDF、SMILES 或 CML 转换失败

图级格式需要可用的分子图。重试 writer 前先检查拓扑恢复状态：

```python
frame = batch[0][-1]
print(frame.rdmol is not None)
print(frame.topology_reconstruction_status)
```

??? example "拓扑恢复失败时的输出"
    ```text
    False
    failed
    ```

修改恢复策略前，先核对电荷、多重度和几何。状态契约与配置选项见
[结构恢复](structure-recovery.md)。

## 导入 RDKit 和 Open Babel 时发生原生崩溃

如果先导入 Open Babel，再导入 RDKit，可能触发原生 segmentation fault。使用 MolOP 时先导入
MolOP；如果代码需要直接导入两个原生库，请先加载 RDKit：

```python
import molop
from rdkit import Chem
from openbabel import pybel
```

MolOP 会在加载依赖 Open Babel 的模块前先初始化 RDKit。修改导入顺序后应重启 Python 进程；
第一次发生冲突后，无法在同一进程内修复原生库状态。

## 提交可复现问题

提供 MolOP 与 Python 版本、最小复现文件、完整命令或 Python 代码、`parser_detection`、
`n_jobs`、完整错误文本，以及预期字段对应的源文件行。在
[GitHub Issues](https://github.com/gentle1999/MolOP/issues) 提交。
