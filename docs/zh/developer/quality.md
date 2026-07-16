# 开发环境与质量门禁

本页面向源码贡献者，不是终端用户安装流程。

## 建立环境

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
uv sync
```

## 常用检查

```bash
uv run ruff format . --check
uv run ruff check .
uv run mypy src
uv run pyright
uv run pytest --no-cov -q
```

仓库聚合门禁：

```bash
make check
```

## 文档检查

```bash
uv run rumdl check docs README.md README.zh.md
uv run pytest -q tests/test_documentation_examples.py --no-cov
uv run pytest -q tests/format_feature_coverage --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```

格式能力页面由 hooks/测试维护。修改 reader/writer 后需要运行格式覆盖测试，而不只是构建站点。

## 提交前

```bash
git diff --check
git status --short
```

确认没有把本地日志、构建目录或无关用户改动加入提交。
