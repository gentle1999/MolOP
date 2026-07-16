# 安装

安装 MolOP 并确认 Python API 与命令行入口可用。

## 准备

- Python 3.10 或更高版本。
- 能从 GitHub 安装 Python 包。
- 建议使用独立虚拟环境，避免与已有 RDKit 环境冲突。

MolOP 当前未发布到 PyPI 或 Conda，终端用户需要从 GitHub 安装。

## 安装

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install git+https://github.com/gentle1999/MolOP.git
```

Windows PowerShell 使用 `.venv\Scripts\Activate.ps1` 激活环境。

## 验证

```bash
python -c "import molop; print(molop.__version__)"
molop --version
molop --help
```

`molop --help` 应显示 `parse` 命令。版本来自 Git tag；直接从未标记的源码运行时可能显示开发
版本标识。

## 常见问题

### Python 版本不兼容

先检查 `python --version`。同一终端中的 `python` 和 `python -m pip` 必须指向刚激活的
Python 3.10+ 环境。

### RDKit 安装失败

RDKit 是 MolOP 的运行时依赖。优先使用支持当前 Python 和操作系统的环境；若 pip 无可用
wheel，可先在 Conda 环境中安装 RDKit，再从 GitHub 安装 MolOP。

### OpenBabel 是否必需

专用 Gaussian、ORCA、xTB、XYZ、SDF 和 SMILES reader 不依赖 OpenBabel fallback。只有读取
未知扩展名或明确选择 OpenBabel 渲染后端时才需要可用的 OpenBabel Python 绑定。

### 如何安装开发环境

源码检出、`uv sync` 和测试命令属于贡献者流程，见[开发环境与质量门禁](../developer/quality.md)。

## 下一步

使用共享样例完成[5 分钟上手](quickstart.md)。
