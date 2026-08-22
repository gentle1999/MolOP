# molop.io

本页面提供了 `molop.io` 模块的 API 参考。

单文件解析可使用轻量入口 `AutoFileParser`：它自动检测格式并直接返回一个 file 级对象，
不会构造 batch parser，也不会调度 worker 进程。批量解析仍使用 `AutoParser`，它接受单个路径、
glob 或路径可迭代对象并返回 `FileBatchModelDisk`。需要显式控制批量解析或路径模式拆分时，再使用
`FileBatchParserDisk` 和 `split_path_pattern` 等底层辅助 API。

::: molop.io
    options:
      members:
        - AutoFileParser
        - FileBatchParserDisk
        - split_path_pattern
