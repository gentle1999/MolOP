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

`ts_vibration` 沿唯一虚频模式位移坐标，并返回候选 `Molecule`。它要求 `frame.is_TS` 为真，且可能
因基本拥挤检查而跳过部分几何：

```python
candidates = ts_frame.ts_vibration(ratio=1.75, steps=7)
print("候选几何数：", len(candidates))
```

??? example "虚频位移候选（图像输出）"

    ![沿虚频模式生成并重建的候选结构](../../assets/examples/ts_imaginary_mode.svg)

    图像由一个具有单个虚频的真实解析 frame 生成。每个面板对应 `ts_vibration(...)`
    返回候选中成功恢复出分子图的一项；二维布局用于比较连接关系，不代表优化后的反应路径。

## 绘制动图

对任一解析文件调用 `draw_animation(...)`，即可将可绘制 frame 输出为 GIF 或动画 SVG，适用于优化、IRC
和扫描轨迹。无法重建 RDKit 分子图的 frame 会被跳过；只有所有 frame 都不可绘制时才报错。默认 legend
保留原始 frame ID；过渡态 frame 会标记 `TS`，有总能量时还会加入能量。

```python
trajectory = ts_batch[0]
trajectory.draw_animation(file_path="trajectory.gif", duration=120, size=(640, 480))
trajectory.draw_animation(
    image_format="svg",
    file_path="trajectory.svg",
    duration=120,
    size=(640, 480),
)
```

指定振动模式时使用 `draw_vibration_animation(...)`；`draw_ts_vibration_animation(...)` 会自动选择唯一虚频。
两者的默认 legend 都包含模式序号、频率和候选位置。传入 `legends=[...]` 可自行指定标签。

```python
ts_frame.draw_vibration_animation(vibration_id=0, file_path="mode-0.gif", steps=9)
ts_frame.draw_ts_vibration_animation(file_path="ts-imaginary-mode.gif", steps=9)
```

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

默认情况下，前后体推断会在虚频振动的正、负两侧分别采样 7 个振幅，振幅在
`min_ratio=0.75` 到 `max_ratio=1.75` 之间等距分布。两侧各自按重建拓扑的出现频次投票，
并为胜出的拓扑保留其最大采样振幅下的三维构象。随后，碎片数更多的候选被判为前体；
碎片数相同时保持负空间、正空间的顺序。`steps` 用于设置每侧的振幅采样数量。

使用 `save_pre_post_ts(...)` 可直接将一对候选写成 XYZ（默认）或 SDF 文件；SDF 会保留重建的分子图和
推断得到的三维构象：

```python
pre_path, post_path = ts_frame.save_pre_post_ts("ts-endpoints", prefix="candidate")
sdf_pre_path, sdf_post_path = ts_frame.save_pre_post_ts(
    "ts-endpoints", prefix="candidate", format="sdf"
)
print(pre_path, post_path)
```

??? example "端点文件路径"

    ```text
    ts-endpoints/candidate_pre.xyz ts-endpoints/candidate_post.xyz
    ```

文件保留推断出的三维候选几何，不是已经优化的真实反应物和产物。

省略 `prefix` 时，带磁盘来源信息的 frame 会将源文件 stem 和 frame ID 加入生成的文件名；仅存在于内存中的
frame 则回退为 `ts_frame_<id>`。

在解析后的计算文件上调用同名方法，会导出其中全部 TS frame；在批处理对象上调用时，每个支持该操作的
计算文件会写入以完整文件 stem 命名的独立目录，不支持 TS 前后体导出的文件会记录警告并跳过。若不同
目录中的源文件具有相同 stem，MolOP 会追加稳定的源路径摘要，保证输出目录彼此独立：

```python
file_endpoints = trajectory.save_pre_post_ts("ts-endpoints", format="sdf")
batch_endpoints = ts_batch.save_pre_post_ts("ts-endpoints-batch", format="sdf", n_jobs=4)
```

两种返回值均以原始 frame ID 为键；批处理返回值额外以源文件路径作为最外层键。

## 相关操作

- 优化轨迹使用 `filter_state("opt")` 和 `parsed_file.closest_optimized_frame`。
- 需要指定振动模式时使用 `frame.vibrate(...)`，不要固定使用第一模式。
- 选中结构后可通过[格式转换与导出](../guides/conversion.md)写出。
