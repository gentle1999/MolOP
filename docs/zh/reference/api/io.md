# molop.io

本页面提供了 `molop.io` 模块的 API 参考。

公共 IO 入口是 `AutoParser`：它接受单个路径、glob 或路径可迭代对象，并返回
`FileBatchModelDisk`。需要显式控制批量解析或路径模式拆分时，再使用 `FileBatchParserDisk` 和
`split_path_pattern` 等底层辅助 API。

::: molop.io
    options:
      members:
        - FileBatchParserDisk
        - split_path_pattern
