# 导出结构和下一步计算输入

从正常结束的最终 frame 同时生成 XYZ、Gaussian input 和 ORCA input。

## Python

```python
from pathlib import Path
from molop import AutoParser

batch = AutoParser("results/*")
selected = batch.filter_state("normal")

for directory in ("xyz", "gjf", "orca"):
    Path(directory).mkdir(exist_ok=True)

selected.format_transform(
    "xyz", output_dir="xyz", write_to_disk=True
)
selected.format_transform(
    "gjf",
    output_dir="gjf",
    write_to_disk=True,
    route_section="#p B3LYP/6-31G(d) opt",
)
selected.format_transform(
    "orcainp",
    output_dir="orca",
    write_to_disk=True,
    keywords="B3LYP def2-SVP Opt",
)
```

## 输出

```text
xyz/
  job.xyz
gjf/
  job.gjf
orca/
  job.inp
```

每个源文件按主名生成一个最终 frame 文件。提交计算前检查 route/keywords、资源、charge、
multiplicity 和溶剂设置；转换不会替你决定科研参数。

## CLI 导出 XYZ

```bash
molop -q parse "results/*" \
  filter-state --state normal \
  format-transform --format xyz --output-dir xyz
```

Gaussian/ORCA writer 的动态选项较多，批量科研工作流通常用 Python 明确传参更易审计。
