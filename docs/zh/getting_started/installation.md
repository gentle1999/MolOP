# 安装

```bash
pip install molop
```

## 验证

```bash
python -c "import molop; print(molop.__version__)"
molop --version
molop --help
```

`molop --help` 应显示 `parse` 命令。直接从未标记的源码运行时，版本可能带有开发版本标识。

验证输出的稳定形状为：

??? example "验证输出"

    ```text
    <version>
    molop, version <version>
    Usage: molop [OPTIONS] COMMAND [ARGS]...
    ...
      parse       Parse files into a FileBatchModelDisk state, then run...
    ```

两个 `<version>` 应一致；具体版本号和完整 help 文本随发布版本变化。

## 常见问题

???+ note "从源码检出建立开发环境"
    需要可编辑源码环境时，请使用[开发环境与质量门禁](../developer/quality.md)。面向终端用户的
    示例假设 `molop` 命令已经位于 `PATH`。

### 导入 RDKit 和 Open Babel 时发生原生崩溃

如果先导入 Open Babel，再导入 RDKit，可能触发原生 segmentation fault。使用 MolOP 时先导入
MolOP；如果直接导入这两个原生库，请先加载 RDKit：

```python
import molop
from rdkit import Chem
from openbabel import pybel
```

MolOP 会在初始化 Open Babel 前先加载 RDKit，以避免原生库冲突。

## 开发环境

源码检出、`uv sync` 和测试命令属于贡献者流程，见[开发环境与质量门禁](../developer/quality.md)。

## 下一步

使用共享样例完成[5 分钟上手](quickstart.md)。
