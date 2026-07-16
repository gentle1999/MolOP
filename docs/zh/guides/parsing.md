# 解析文件

用 `AutoParser` 读取单个文件、glob 或混合路径列表。

## 最短示例

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
for parsed_file in batch:
    print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

返回的 `FileBatchModelDisk` 只包含成功解析且至少有一个 frame 的文件。
对共享样例，输出为：

```text
water_mp2.out 1 orcaout
```

## 输入方式

```python
one = AutoParser("calculation.log")
many = AutoParser("results/*.out")
mixed = AutoParser([
    "gaussian/*.log",
    "orca/job.out",
    "structures/*.xyz",
])
```

glob 由 MolOP 展开，不需要在 shell 中预先展开。输入会转成绝对路径、去重并排序。

## 选择 frame

```python
parsed_file = batch[0]
first = parsed_file[0]
final = parsed_file[-1]

for frame in parsed_file:
    print(frame.frame_id, frame.energies.total_energy if frame.energies else None)
```

只关心最终结果且希望降低内存占用时：

```python
batch = AutoParser("results/*.log", only_last_frame=True)
```

`only_last_frame=True` 改变解析结果中保留的 frame，不等同于解析全部 frame 后再取 `[-1]`；
需要完整优化轨迹时不要开启。

## 自动与显式格式检测

默认 `parser_detection="auto"` 根据扩展名选择候选 reader，并用内容探测排除不匹配项。
扩展名含义不明确时显式指定格式 ID：

```python
orca = AutoParser("job.out", parser_detection="orcaout")
xtb = AutoParser("xtb.out", parser_detection="xtbout")
fchk = AutoParser("molecule.fchk", parser_detection="g16fchk")
```

可用格式 ID 见[格式支持概览](../reference/format_support.md)。

## 并行与结构优先

```python
batch = AutoParser("results/*", n_jobs=4)

structures = AutoParser(
    "results/*.log",
    n_jobs=4,
    only_extract_structure=True,
)
```

`n_jobs=-1` 使用 MolOP 配置允许的最大并行度；小批量或排查问题时使用 `n_jobs=1`。
`only_extract_structure=True` 跳过非结构结果，不适合能量或热化学提取。

## 发现失败文件

```python
from pathlib import Path
from molop import AutoParser

inputs = sorted(Path("results").glob("*.out"))
batch = AutoParser(inputs, n_jobs=1)
parsed = {Path(path).resolve() for path in batch.file_paths}
failed = [path for path in inputs if path.resolve() not in parsed]

for path in failed:
    print("未加入 batch:", path)
```

全部成功时没有输出；失败时每个遗漏路径一行，例如：

```text
未加入 batch: results/truncated.out
```

不存在、格式不支持或没有可用 frame 的文件不会进入 batch，并会写入 MolOP 日志。排查时先用
`n_jobs=1`，再按[常见问题](troubleshooting.md)检查扩展名、编码和格式 ID。

## 深入了解

- [读取计算结果](results.md)
- [批量汇总](batch.md)
- [格式支持概览](../reference/format_support.md)
- [Source evidence 与序列化](../reference/api_contracts.md)
