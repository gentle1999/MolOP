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
uv run rumdl check README.md README.zh.md $(rg --files docs -g '*.md')
uv run pytest -q tests/test_documentation_examples.py --no-cov
uv run pytest -q tests/format_feature_coverage --no-cov
NO_MKDOCS_2_WARNING=1 uv run mkdocs build --strict
```

格式能力页面由 hooks/测试维护。修改 reader/writer 后需要运行格式覆盖测试，而不只是构建站点。

## 发布流程

发布需要在 GitHub 上手动创建。将发布改动合并到 `main` 后，打开[新建 release 页面](https://github.com/gentle1999/MolOP/releases/new)，
创建一个从 `main` 指向的 `v*` tag（例如 `v0.2.13`），填写 release 内容并发布。

生成的 `v*` tag 事件会触发 `.github/workflows/ci.yaml` 和 `.github/workflows/docs-deploy.yml`。质量、构建、文档和测试门禁全部通过后，
CI 会把构建产物上传到已存在的 GitHub Release，并通过 Trusted Publishing 发布到 PyPI。文档 workflow 会发布版本化文档并更新 `latest` 别名。
发布前确认两个 workflow 都已成功。

## 提交前

```bash
git diff --check
git status --short
```

确认没有把本地日志、构建目录或无关用户改动加入提交。
