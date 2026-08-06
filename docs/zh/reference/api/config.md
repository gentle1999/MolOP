# molop.config

本页面提供了 `molop.config` 模块的 API 参考。

模块暴露进程级 `molopconfig` 对象，用于控制进度显示、并行任务上限、日志和结构恢复选项。
单次调用的行为优先使用显式函数参数；只有需要影响后续调用的全局策略时才修改全局配置。

::: molop.config
