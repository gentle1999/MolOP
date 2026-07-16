# Reader/Writer 插件开发

新增格式前先确认现有专用 codec 或 OpenBabel fallback 不能满足需求。

## Reader 最小职责

1. 声明稳定 `format_id`、扩展名、优先级和内容 probe。
2. 通过共享解析生命周期读取源文件。
3. 生成至少一个有效 File/Frame 模型。
4. 对可选科学字段记录 presence、diagnostic 和必要 source evidence。
5. 为真实或最小 fixture 增加解析、错误和格式探测测试。

## Writer 最小职责

1. 声明目标格式 ID、扩展名和能力层级。
2. 实现 frame 选择、合并/拆分和 render/write 行为。
3. 明确坐标级、图级和格式专用信息边界。
4. 在 registry 中声明动态选项及 CLI help。
5. 增加字符串渲染、写盘、命名和失败行为测试。

## 注册位置

内置 codec 由 `src/molop/io/codecs/catalog.py` 和对应格式包注册。不要在文档中维护第二份格式
清单；能力状态由 `tests/format_feature_coverage/support_matrix.py` 及格式页生成区段维护。

## 必读契约

完整类层级、注册顺序、probe、source span、解析诊断、transform 和测试约束保存在
[Parser 契约](parser-contract.md)。

实现完成后运行[质量门禁](quality.md)，并同时更新中英文格式页及用户示例。
