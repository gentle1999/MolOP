# 文档与 MolOP 版本

本网站从 MolOP 当前源码构建，并保留稳定发布 tag 的文档快照。使用页面顶部版本选择器可以在最新稳定版、此前发布版和 `main` 开发版之间切换。页面顶部横幅给出当前快照对应的源码版本、Git ref、commit 和最近发布版本；报告问题或核对 API 时，应同时提供这些信息。

## 两种文档状态

| 状态 | 含义 |
| --- | --- |
| 开发版文档 | 源码位于 `main` 的未发布提交或本地修改，可能包含 PyPI 最新版本尚未提供的 API。 |
| 发布版文档 | 构建源码精确对应一个稳定发布 tag，且工作区没有额外修改。 |

稳定版本快照使用 `/<版本>/` 路径，`latest` 指向最新稳定发布版，`main` 和 `dev` 指向当前主线。根地址默认进入 `latest`；因此旧的无版本根路径仍会落到最新稳定文档，版本化链接则保持不变。

版本快照从稳定版本 tag 触发构建，并更新 `latest`。预发布 tag 不进入此版本系统。历史稳定 tag 在版本系统启用前没有自动生成快照；需要通过文档部署 workflow 的 `workflow_dispatch` 选择对应 tag 回填。

回填时将 `source_ref` 设置为目标稳定 tag，例如 `v0.2.2`；只有目标确实是最新稳定发布版时才把 `update_latest` 设为 `true`，否则保持 `false`，避免历史版本覆盖 `latest`。

## 核对本地版本

```bash
python -c "import molop; print(molop.__version__)"
```

如果输出与横幅中的源码版本不同，先查看横幅列出的最近发布版本：

- 使用 PyPI 发布版时，以对应发布 tag 的 API 为准。
- 需要验证主线新增 API 时，安装横幅中显示的 commit，而不是假定最新 PyPI 包已经包含该功能。

```bash
uv add "molop @ git+https://github.com/gentle1999/MolOP.git@COMMIT"
```

将 `COMMIT` 替换为横幅中的 commit。生产环境和可复现实验应固定发布版本或完整 commit，不要依赖浮动的 `main`。

## 插件兼容性

第三方 reader/writer 插件应在 `pyproject.toml` 中声明经过测试的 MolOP 版本范围。若插件依赖开发版文档中的新接口，应先固定主线 commit 进行开发测试；发布插件前，等待该接口进入正式 MolOP 版本并更新依赖下限。

插件接入方式见 [Reader/Writer 插件开发](../extensions/plugins.md)。
