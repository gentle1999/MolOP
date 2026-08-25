# molop.io

This page provides the API reference for the `molop.io` module.

`AutoFileParser` is the lightweight single-file entry point. It detects the
format and returns one file model directly, without constructing a batch
parser or scheduling worker processes. `AutoMemoryParser`, `AutoTextParser`,
and `AutoBytesParser` provide the same file-level result for already-loaded
text, bytes, or binary/text streams. A string passed to `AutoMemoryParser` is
source text, not a path. `AutoParser` remains the batch entry point for a
path, glob, or iterable of paths and returns a
`FileBatchModelDisk`. `FileBatchParserDisk` and `split_path_pattern` are
lower-level helpers for callers that need explicit batch parsing or
path-pattern handling.

::: molop.io
    options:
      members:
        - AutoFileParser
        - AutoMemoryParser
        - AutoTextParser
        - AutoBytesParser
        - FileBatchParserDisk
        - split_path_pattern
