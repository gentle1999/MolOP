# Quickstart

## Goal

Get a quick, reproducible tour of MolOP's core workflow: parse -> summarize.

This page is an entry-level guide; the canonical baseline (behavior/output/API) is the notebook `examples/01-gaussian-parse-and-inspect.ipynb`. If there is any conflict, the notebook wins.

## Prerequisites

- MolOP installed (see [Installation](installation.md)).
- Access to the repository's test files (if running from the source).

## Steps

### 1. Parse a Gaussian Output File and generate a summary

Use `AutoParser` to read a Gaussian `.log` file and call `to_summary_df()` to obtain a compact, verifiable summary.

```python
import os
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")

from molop.io import AutoParser

# Locate the repository root by searching for pyproject.toml
def get_repo_root():
    current = Path(os.getcwd())
    for parent in [current] + list(current.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    return current

repo_root = get_repo_root()
file_path = repo_root / "tests/test_files/g16log/2-TS1-Opt.log"

# Parse the file
# AutoParser returns a FileBatchModelDisk object
batch = AutoParser(file_path.as_posix(), n_jobs=1)
parsed_file = batch[0]

# Compact summary (default brief=True)
df_brief = parsed_file.to_summary_df()
df_brief
```

### 2. Get a richer summary

The notebook demonstrates a more complete summary via `brief=False`:

```python
df_full = parsed_file.to_summary_df(brief=False)
df_full
```

## Related Links

- [Installation](installation.md)
- [01-Gaussian Parse and Inspect (Notebook)](../examples/01-gaussian-parse-and-inspect.ipynb)
- [Tutorials Overview](../tutorials/index.md)
- [API Reference](../reference/api.md)
