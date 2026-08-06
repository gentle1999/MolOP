# 文档贡献

用户文档以任务为入口，中英文核心页面必须成对维护。

## 页面顺序

1. 说明任务和输入。
2. 给出 5-20 行最短可运行示例。
3. 紧跟真实输出、返回形状或生成文件树。
4. 说明常用变体和实际限制。
5. 链接到格式矩阵或精确契约。

终端用户示例使用 `molop`，只有开发环境命令使用 `uv run molop`。快速上手不得依赖
`tests/test_files/`、源码根目录或测试环境变量。

## 事实来源

- API 名称和默认值：当前 Python 源码与类型提示。
- CLI：`molop --help` 和具体子命令 help。
- 格式能力：`tests/format_feature_coverage/support_matrix.py`。
- 核心示例：`docs/assets/examples/` 与 `tests/test_documentation_examples.py`。

不要手写与能力矩阵重复的“支持/不支持”表。

## 输出渲染

Markdown 示例中的输出必须作为独立的可折叠内容呈现。文档站点使用
Material details：

````markdown
??? example "输出"
    ```text
    ...
    ```
````

命令、输入代码和解释文字放在折叠块外。打印值、返回文本、生成的文件树以及用于核对的其他
产物都遵守这条规则。README 使用原生 HTML `<details>`，因为 GitHub 不渲染 Material admonition。

Notebook 的输出由另一条流程负责：只编辑代码单元和 Markdown 单元，然后由 CI 在构建 MkDocs
前使用 `nbconvert --execute` 执行 Notebook。`mkdocs-jupyter` 以 `execute: false` 读取保存的
Notebook 并渲染 HTML；不要手动修改 Notebook 的 `outputs`、执行计数或生成的 HTML。

需要在 Markdown 页面直接展示某个 Notebook 单元格的 HTML 输出时，使用单元格 `id` 声明：

```markdown
<!-- notebook-output: examples/02-batch-summary-filter-select.ipynb#batch-summary -->
```

`scripts/mkdocs_hooks.py` 会读取已执行 Notebook 的 `text/html` 输出并嵌入当前语言页面。目标
单元格必须有稳定 `id`，且 CI 执行后必须产生 HTML；不要把 DataFrame HTML 复制到 Markdown。

## 新页面

同一改动中完成：

- `docs/zh/...` 中文页面。
- `docs/en/...` 英文语义本地化页面。
- `mkdocs.yml` 导航。
- 可执行代码对应测试。

## 校验

```bash
uv run rumdl check README.md README.zh.md $(rg --files docs -g '*.md')
uv run pytest -q tests/test_documentation_examples.py --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```
