<!--
 * @Author: TMJ
 * @Date: 2026-02-11 16:48:41
 * @LastEditors: TMJ
 * @LastEditTime: 2026-02-11 20:23:31
 * @Description: 请填写简介
-->

# 安装指南

## 目标

安装 MolOP 及其依赖项，开始您的计算化学工作流程。

## 前置条件

- **Python 3.10+**: 请确保您安装了兼容的 Python 版本。
- **RDKit & OpenBabel**: 这些是运行时实现完整功能（分子图重构、格式转换等）所必需的。

## 步骤

### 选项 1：面向终端用户（通过 pip）

直接从 GitHub 仓库安装最新版本：

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

### 选项 2：面向开发人员（通过 uv）

如果您想为 MolOP 贡献代码或运行测试套件，我们建议使用 `uv`：

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

## 预期输出

安装完成后，您应该能够运行以下验证命令而不会报错。

### 验证版本

```bash
python -c "import molop; print(molop.__version__)"
```

_预期：类似于 `0.1.0` 的版本字符串。_

### 验证核心组件

```bash
python -c "from molop.io import AutoParser; print(AutoParser)"
```

_预期：`<function AutoParser at ...>`。_

## 相关链接

- [快速上手](quickstart.md)
- [GitHub 仓库](https://github.com/gentle1999/MolOP)
