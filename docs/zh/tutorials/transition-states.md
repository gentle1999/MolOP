# 过渡态分析

本流程适用于包含频率结果的过渡态候选输出。MolOP 可以筛选候选 frame、检查虚频，并按几何生成
前后体候选；它不会优化这些结构，也不会证明反应路径成立。

## 筛选候选 frame

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
normal = batch.filter_state("normal", n_jobs=-1)
ts_batch = normal.filter_state("ts", n_jobs=-1)

print("解析文件：", len(batch))
print("正常结束：", len(normal))
print("过渡态文件：", len(ts_batch))
```

??? example "筛选结果（由输入决定）"

    ```text
    解析文件：<数量>
    正常结束：<数量>
    过渡态文件：<数量>
    ```

计数取决于输入目录。一个文件可能包含多个 frame，因此不要假设 `parsed_file[-1]` 一定是
`filter_state("ts")` 选中的 frame，应显式定位：

```python
for parsed_file in ts_batch:
    for frame in parsed_file:
        if frame.is_TS:
            imaginary = frame.vibrations.num_imaginary if frame.vibrations else None
            print(parsed_file.filename, frame.frame_id, imaginary)
```

??? example "候选 frame 输出形状（由输入决定）"

    ```text
    <源文件名> <frame_id> <虚频数或 None>
    ```

`is_TS` 是程序化候选判断。`None` 表示源文件证据不足，与 `False` 不同。

输出中的文件数、frame ID 和虚频数取决于输入数据。

## 检查虚频

```python
ts_frame = next(
    (frame for parsed_file in ts_batch for frame in parsed_file if frame.is_TS),
    None,
)
if ts_frame is None:
    raise ValueError("输入中没有可用的过渡态 frame")
vibrations = ts_frame.vibrations

if vibrations is None:
    raise ValueError("过渡态候选没有结构化频率结果")

print("虚频数：", vibrations.num_imaginary)
print("频率（cm^-1）：", vibrations.frequencies.m_as("cm^-1"))
```

??? example "虚频检查输出（由输入决定）"

    ```text
    虚频数：<数量>
    频率（cm^-1）：<频率数组>
    ```

需要结合频率正负、位移方向和化学含义复核。虚频数量本身不能证明反应路径正确。

## 沿虚频生成结构

`ts_vibration` 沿第一个振动模式位移坐标，并返回候选 `Molecule`。它要求 `frame.is_TS` 为真，且
可能因基本拥挤检查而跳过部分几何：

```python
candidates = ts_frame.ts_vibration(ratio=1.75, steps=7)
print("候选几何数：", len(candidates))
```

??? example "虚频位移候选（图像输出）"

    ![沿虚频模式生成并重建的候选结构](../../assets/examples/ts_imaginary_mode.svg)

    图像由一个具有单个虚频的真实解析 frame 生成。每个面板对应 `ts_vibration(...)`
    返回候选中成功恢复出分子图的一项；二维布局用于比较连接关系，不代表优化后的反应路径。

## 推断前后体候选与键变化

```python
try:
    reactant, product = ts_frame.possible_pre_post_ts(show_3D=True)
    print(reactant.GetNumAtoms(), product.GetNumAtoms())
except ValueError as exc:
    print("前后体推断失败：", exc)

difference = ts_frame.to_diff_rdmol()
if difference is not None:
    print("差异图键数：", difference.GetNumBonds())
```

??? example "前后体与虚键差异图（图像输出）"

    ![前体候选、虚键差异图和后体候选](../../assets/examples/ts_endpoints_difference.svg)

    红色标出前后体候选之间发生变化的键。中央面板来自 `to_diff_rdmol()`；其中的零阶键
    是 MolOP 的可机读差异标记，生成脚本只在用于绘图的副本中把它们转换为可见线段。

`possible_pre_post_ts` 返回基于几何的前后体候选，不是优化后的真实反应物和产物。
`to_diff_rdmol` 在无法推断出支持的断键差异时可能返回 `None`。自动反应流程使用这些结果前，
仍需检查原子映射、连接关系和原始计算结果。

## 相关操作

- 优化轨迹使用 `filter_state("opt")` 和 `parsed_file.closest_optimized_frame`。
- 需要指定振动模式时使用 `frame.vibrate(...)`，不要固定使用第一模式。
- 选中结构后可通过[格式转换与导出](../guides/conversion.md)写出。
