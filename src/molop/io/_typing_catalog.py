"""Runtime-safe typing catalog.

This module intentionally avoids importing concrete file/frame model classes at
runtime (to keep imports light and avoid import cycles).

Static type checkers should consume the sibling stub file:
`src/molop/io/_typing_catalog.pyi`.
"""

from __future__ import annotations

from typing import Any, TypeAlias


# At runtime we don't need precise unions here.
FileDiskObj: TypeAlias = Any
FrameDiskObj: TypeAlias = Any
