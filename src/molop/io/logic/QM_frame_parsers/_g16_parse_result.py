from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any


def _payload_mapping(value: Any) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return {key: item for key, item in value.items() if item is not None}
    if hasattr(value, "model_dump"):
        return value.model_dump(
            exclude_unset=True,
            exclude_none=True,
            exclude_computed_fields=True,
        )
    return {}


def _has_payload_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    try:
        return len(value) > 0  # type: ignore[arg-type]
    except Exception:
        return True


def _merge_payload(
    existing: Any,
    incoming: Mapping[str, Any],
    *,
    overwrite: bool = False,
) -> dict[str, Any]:
    merged = _payload_mapping(existing)
    for key, value in incoming.items():
        if value is None:
            continue
        if overwrite or not _has_payload_value(merged.get(key)):
            merged[key] = value
    return merged


@dataclass(slots=True)
class G16FrameParseResult:
    """Canonical state-machine output before model validation."""

    fields: dict[str, Any] = field(default_factory=lambda: {"qm_software": "Gaussian"})

    def set(self, key: str, value: Any) -> None:
        self.fields[key] = value

    def update(self, values: Mapping[str, Any]) -> None:
        for key, value in values.items():
            self.fields[key] = value

    def set_missing_from(self, values: Mapping[str, Any]) -> None:
        for key, value in values.items():
            if not _has_payload_value(self.fields.get(key)):
                self.fields[key] = value

    def merge_payload_field(
        self,
        key: str,
        incoming: Mapping[str, Any],
        *,
        overwrite: bool = False,
    ) -> None:
        if key in self.fields:
            self.fields[key] = _merge_payload(self.fields[key], incoming, overwrite=overwrite)
        else:
            self.fields[key] = _payload_mapping(incoming)

    def has_value(self, key: str) -> bool:
        return _has_payload_value(self.fields.get(key))

    def model_data(self) -> dict[str, Any]:
        return dict(self.fields)
