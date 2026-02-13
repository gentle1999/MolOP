# 快速上手

## 目标

通过一个可复现的小例子，快速了解 MolOP 的核心工作流：解析 -> 摘要。

本页的具体行为与输出以 Notebook `examples/01-gaussian-parse-and-inspect.ipynb` 为绝对基准；如果本页与 Notebook 有冲突，以 Notebook 为准。

## 前置条件

- 已安装 MolOP（参见[安装指南](installation.md)）。
- 可以访问仓库的测试文件（如果从源码运行）。

## 步骤

### 1. 解析 Gaussian 输出文件并生成摘要

使用 `AutoParser` 读取 Gaussian `.log` 文件，并通过 `to_summary_df()` 得到紧凑、可核对的摘要信息。

```python
import os
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")

from molop.io import AutoParser

# 通过搜索 pyproject.toml 定位仓库根目录
def get_repo_root():
    current = Path(os.getcwd())
    for parent in [current] + list(current.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    return current

repo_root = get_repo_root()
file_path = repo_root / "tests/test_files/g16log/2-TS1-Opt.log"

# 解析文件
# AutoParser 返回一个 FileBatchModelDisk 对象
batch = AutoParser(file_path.as_posix(), n_jobs=1)
parsed_file = batch[0]

# 紧凑摘要（默认 brief=True）
df_brief = parsed_file.to_summary_df()
df_brief
```

### 2. 获取更丰富的摘要

Notebook 中展示了 `brief=False` 的更完整摘要：

```python
df_full = parsed_file.to_summary_df(brief=False)
df_full
```

## 相关链接

- [安装指南](installation.md)
- [01-Gaussian 解析与检查（Notebook）](../examples/01-gaussian-parse-and-inspect.ipynb)
- [教程概览](../tutorials/index.md)
- [API 参考](../reference/api.md)
