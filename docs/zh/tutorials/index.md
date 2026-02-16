# 教程

本节提供 MolOP 的使用路径导航与要点说明。

推荐阅读顺序：

- [01-Gaussian 解析与检查（Notebook）](../examples/01-gaussian-parse-and-inspect.ipynb)
- [02-批量汇总、过滤与选择（Notebook）](../examples/02-batch-summary-filter-select.ipynb)
- [03-转换与导出（Notebook）](../examples/03-transform-and-export.ipynb)

## 1. 文件解析 (IO Parsing)

### 前提条件

- 已安装 MolOP。
- 拥有计算化学输出文件（如 `.log`）或坐标文件（如 `.xyz`, `.sdf`）。

### 代码示例

```python
from molop.io import AutoParser
# 解析单个文件
files = AutoParser("example.log")
file = files[0]
# 获取最后一帧
frame = file[-1]
print(f"能量: {frame.energies.total_energy}")
```

### 预期输出

- `files` 是一个 `FileBatchModelDisk` 对象。
- `frame` 包含 `energies`, `coords`, `atoms` 等字段。

## 2. 批处理 (Batch Processing)

### 代码示例

```python
from molop.io import AutoParser
# 使用通配符并行解析
batch = AutoParser("*.log", n_jobs=-1)
# 过滤成功优化的结构
opt_batch = batch.filter_state("opt")
print(f"成功解析 {len(opt_batch)} 个优化结构")
```

### 预期输出

- 返回过滤后的 `FileBatchModelDisk` 对象。

## 3. 结构恢复 (Structure Reconstruction)

### 前提条件

- 只有原子坐标（如来自 XYZ 文件），需要推导成键信息。

### 代码示例

```python
from molop.io import AutoParser
batch = AutoParser("molecule.xyz")
# 自动触发结构恢复算法
rdmol = batch[0][0].to_rdmol
```

### 预期输出

- 返回一个带有成键信息和键级的 RDKit `Mol` 对象。

## 4. 格式转换 (Transforms / Format Conversion)

### 前提条件

- 已解析的批处理对象。

### 代码示例

```python
from molop.io import AutoParser
batch = AutoParser("*.log")
# 批量转换为 SDF 并保存（具体参数以 Notebook 03 为准）
batch.format_transform(format="sdf", output_dir="./output", frameID=-1)
```

### 预期输出

- `./output` 目录下生成对应的 `.sdf` 文件。

## 5. CLI 使用示例

### 前提条件

- 终端已配置 `molop` 命令。

### 命令示例

```bash
# 生成摘要 CSV
uv run molop -q summary "tests/test_files/g16log/2-TS1-Opt.log" --out tutorial_ts_summary.csv --format csv --mode frame --frame -1 --n-jobs 1

# 将分子文件转换为另一种格式
uv run molop -q transform "tests/test_files/orca/h2_grad_orca.inp" --to sdf --output-dir ./tutorial_transform_out --frame -1 --embed --parser-detection orcainp --n-jobs 1
```

## 6. 插件与编解码器扩展 (Plugin/Codec Extension)

### 前提条件

- 需要支持自定义的私有格式。

### 代码示例

在 `src/molop/io/logic/coords_parsers/MyParser.py` 中定义：

```python
def register(registry):
    @registry.reader_factory(format_id="myfmt", extensions={".myfmt"}, priority=0)
    def my_reader_factory():
        return MyReader()

class MyReader:
    format_id = "myfmt"
    extensions = frozenset({".myfmt"})
    priority = 0
    def read(self, path, **kwargs):
        # 实现解析逻辑
        ...
```

### 预期输出

- `AutoParser` 现在可以识别 `.myfmt` 后缀。

### 注意

- 注册函数必须命名为 `register`，详见 [核心概念](../concepts.md)。
