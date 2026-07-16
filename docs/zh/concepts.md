# 理解 Batch、File 与 Frame

使用 MolOP 只需要掌握三个层级。

```text
FileBatchModelDisk
  ├─ 文件 0
  │   ├─ frame 0
  │   └─ frame 1
  └─ 文件 1
      └─ frame 0
```

## Batch

`AutoParser(...)` 的返回值。它保存多个已解析文件，并提供汇总、筛选、分组和批量转换。

```python
batch = AutoParser("results/*.log")
print(len(batch), batch.file_names)
```

输出是文件数量和文件名列表，例如：

```text
2 ['reactant.log', 'product.log']
```

## File

`batch[0]` 表示一个源文件。一个文件可能只有一个单点 frame，也可能包含优化轨迹、频率步骤
或多个计算 segment。

```python
parsed_file = batch[0]
print(parsed_file.filename, len(parsed_file), parsed_file.detected_format_id)
```

示例输出：

```text
reactant.log 12 g16log
```

## Frame

frame 是一次结构/结果快照。大多数用户任务读取最后一帧：

```python
final = parsed_file[-1]
print(final.frame_id, final.charge, final.multiplicity)
```

需要轨迹时遍历全部 frame，不要假设每一帧都包含相同结果容器。

## 最后一帧不总是“最好结构”

错误终止、Link1、多 segment 或频率任务可能使最后一帧与优化收敛 frame 不同。筛选优化结果
优先使用 `filter_state("opt")` 和 `parsed_file.closest_optimized_frame`，再结合原始输出复核。

## 深入了解

Registry、codec、source lifecycle 和内部模型边界属于开发者概念，见
[架构概览](developer/overview.md)和 [API 契约](reference/api_contracts.md)。
