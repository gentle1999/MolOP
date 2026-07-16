# 过滤与选择

从 batch 中选择正常结束、成功优化、过渡态或满足特定电荷与格式条件的文件。

## 状态筛选

```python
from molop import AutoParser

batch = AutoParser("results/*")
normal = batch.filter_state("normal")
optimized = normal.filter_state("opt")
transition_states = normal.filter_state("ts")

print(len(batch), len(normal), len(optimized), len(transition_states))
```

输出四个整数，依次为总文件数、正常结束数、优化成功数和过渡态数。筛选返回新 batch，不会
修改原对象。

状态值：

| 状态 | 保留条件 |
| --- | --- |
| `normal` | 文件中的 frame 均未显示错误，且格式提供正常状态证据 |
| `error` | 至少一个 frame 显示错误 |
| `opt` | 存在优化状态，且最接近优化完成的 frame 已收敛 |
| `ts` | 至少一个 frame 的 `is_TS` 为真 |
| `thermal` | 至少一个 frame 有热化学结果 |
| `no-img` | 有频率结果，且所有相关 frame 都没有虚频 |

反向选择：

```python
not_normal = batch.filter_state("normal", negate=True)
```

## 用共享样例核对

```python
batch = AutoParser("water_mp2.out", n_jobs=1)
print(len(batch), len(batch.filter_state("normal")))
```

输出：

```text
1 1
```

## 按电荷、多重度或扩展名

```python
neutral = batch.filter_value("charge", 0)
doublets = batch.filter_value("multiplicity", 2)
logs = batch.filter_value("format", ".log")
```

`format` 比较源文件扩展名。按实际检测 reader 选择时使用 codec ID：

```python
orca_outputs = batch.filter_by_codec_id("orcaout")
gaussian_logs = batch.filter_by_codec_id("g16log")
```

## 自定义科学条件

```python
def has_one_imaginary_frequency(parsed_file):
    frame = parsed_file[-1]
    return bool(frame.vibrations and frame.vibrations.num_imaginary == 1)

one_imaginary = batch.filter_custom(has_one_imaginary_frequency)
```

## CLI 等价操作

```bash
molop -q parse "results/*.log" \
  filter-state --state normal \
  filter-state --state ts \
  to-summary-df --full --out transition_states.csv
```

```bash
molop -q parse "results/*" \
  filter-by-codec --codec-id orcaout \
  filter-value --target charge --value 0 \
  to-summary-df --out neutral_orca.csv
```

## 注意

状态值可能为 `None`，表示源文件没有足够证据。严格科研筛选时应把“未知”与明确 `False`
分开，并抽查原始输出。

## 下一步

- [筛选优化结果与过渡态](../tutorials/select-results.md)
- [格式转换与导出](conversion.md)
- [CLI 常用任务](cli-recipes.md)
