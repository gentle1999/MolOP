# molop.io

This page provides the API reference for the `molop.io` module.

The public IO entry point is `AutoParser`, which accepts a path, glob, or
iterable of paths and returns a `FileBatchModelDisk`. `FileBatchParserDisk` and
`split_path_pattern` are lower-level helpers for callers that need explicit
batch parsing or path-pattern handling.

::: molop.io
    options:
      members:
        - FileBatchParserDisk
        - split_path_pattern
