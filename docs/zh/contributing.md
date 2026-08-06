# 贡献 MolOP

感谢参与 MolOP。先根据任务选择开发者文档：

- [开发者概览](developer/overview.md)
- [Reader/Writer 插件开发](developer/plugins.md)
- [完整 Parser 契约](developer/parser-contract.md)
- [文档贡献](developer/documentation.md)
- [开发环境与质量门禁](developer/quality.md)

提交 issue 时提供最小复现文件、完整命令、版本和预期结果。提交代码时保持改动范围清晰，并为
新增或变化的公共行为补充测试与中英文文档。

## 提交 pull request 前

1. 运行变更行为对应的针对性测试。
2. 运行相关的格式覆盖或文档检查。
3. 当源码变更跨越多个模块时运行 `make check`。
4. 检查 `git diff --check`，确认没有加入生成文件、本地日志或无关改动。

新增 reader 或 writer 时，先阅读[完整 Parser 契约](developer/parser-contract.md)和
[插件开发指南](developer/plugins.md)。修改文档时同步维护中英文页面，并遵循
[文档贡献规则](developer/documentation.md)。
