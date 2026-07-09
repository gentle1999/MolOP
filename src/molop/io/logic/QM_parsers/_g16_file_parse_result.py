from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any


@dataclass(slots=True)
class G16FileParseContext:
    """Mutable local scan context for Gaussian file metadata extraction."""

    content: str

    def split(self, pattern: Any) -> str:
        focus_content, self.content = pattern.split_content(self.content)
        return focus_content


@dataclass(slots=True)
class G16FileParseResult:
    """Canonical Gaussian file metadata before file-model validation."""

    fields: dict[str, Any] = field(default_factory=lambda: {"qm_software": "Gaussian"})

    def set(self, key: str, value: Any) -> None:
        self.fields[key] = value

    def update(self, values: Mapping[str, Any]) -> None:
        self.fields.update(values)

    def model_data(self) -> dict[str, Any]:
        return dict(self.fields)
