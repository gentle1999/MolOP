# 解析文件

用 `AutoParser` 读取单个文件、glob 或混合路径列表。

如果只处理单个文件，不需要 batch 准备或 worker 调度，可使用 `AutoFileParser`：

```python
from molop import AutoFileParser

parsed_file = AutoFileParser("water_mp2.out")
print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

该接口会根据扩展名和源内容自动检测格式，返回值本身就是 file 级对象。

## 最短示例

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
for parsed_file in batch:
    print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

返回的 `FileBatchModelDisk` 只包含成功解析且至少有一个 frame 的文件。
对共享样例，输出为：

??? example "输出"

    ```text
    water_mp2.out 1 orcaout
    ```

如果调用方已经持有源内容，可以直接使用内存入口，不需要创建临时文件：

```python
from molop import AutoBytesParser, AutoTextParser

parsed_text = AutoTextParser("1\nwater\nH 0.0 0.0 0.0\n", parser_detection="xyz")
parsed_bytes = AutoBytesParser(raw_bytes, parser_detection="xyz")
```

`AutoMemoryParser` 同时接受文本和 bytes。传入其中的字符串始终按文本处理，不会当作路径；路径请使用
`AutoFileParser`。

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

print([(frame.frame_id, len(frame.atoms)) for frame in parsed_file])
```

??? example "输出"

    ```text
    [(0, 3)]
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
batch = AutoParser("results/*")

structures = AutoParser(
    "results/*.log",
    n_jobs=-1,
    only_extract_structure=True,
)
```

`AutoParser` 和所有 batch 操作默认使用 `n_jobs=-1`。MolOP 会将它解析为当前进程可用 CPU 与
全局 `molopconfig.max_jobs` 上限中较小的值。默认 `max_jobs=None` 时，会跟随系统能够暴露的
调度器/容器配额和 CPU affinity；需要进程级上限时设置正整数。传入正整数 `n_jobs` 可指定
worker 数量；小批量或排查问题时使用 `n_jobs=1`。CLI 的 parse 和操作命令也都独立默认使用
`-1`。

文件解析本身不使用 MolGR 安全上限：`AutoParser` 仍使用完整的 `effective_max_jobs`。只有进入
MolGR 的分子图相关操作才使用更严格的可用 CPU 三分之二上限。

`only_extract_structure=True` 会跳过非结构结果，不适合提取能量或热化学。

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

??? example "存在遗漏文件时的输出"

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
