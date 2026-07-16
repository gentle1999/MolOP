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

## 新页面

同一改动中完成：

- `docs/zh/...` 中文页面。
- `docs/en/...` 英文语义本地化页面。
- `mkdocs.yml` 导航。
- 可执行代码对应测试。

## 校验

```bash
uv run rumdl check docs README.md README.zh.md
uv run pytest -q tests/test_documentation_examples.py --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```
