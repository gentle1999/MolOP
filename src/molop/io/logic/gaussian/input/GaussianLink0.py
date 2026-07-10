from __future__ import annotations

from pint.facets.plain import PlainQuantity
from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.unit import atom_ureg


class GaussianLink0Command(BaseDataClassWithUnit):
    key: str = Field(description="Gaussian Link0 keyword key")
    value: str | None = Field(description="Gaussian Link0 keyword value")

    def _render(self, **kwargs) -> str:
        return f"%{self.key}={self.value}\n" if self.value else f"%{self.key}\n"

    def memory_request(self) -> PlainQuantity | None:
        if self.key.lower() == "mem":
            if self.value and self.value.endswith("B"):
                return atom_ureg.Quantity(self.value)
            if self.value and self.value.endswith("W"):
                return atom_ureg.Quantity(self.value.replace("W", "B")) * 8
        return None

    def cpu_request(self) -> int | None:
        key = self.key.lower()
        if key == "cpu" and self.value:
            return len(self.value.split(","))
        if key in {"nproc", "nprocshared"} and self.value:
            return int(self.value)
        return None

    def shared_memory_cpu_value(self) -> str | None:
        if self.key.lower() != "nprocshared" or self.value is None:
            return None
        return self.value.strip()


class GaussianLink0Commands(BaseDataClassWithUnit):
    link0_keywords: list[GaussianLink0Command] = Field(
        default_factory=list, description="Gaussian Link0 keywords"
    )

    @classmethod
    def from_dict(cls, data: dict[str, str | None]) -> GaussianLink0Commands:
        link0_keywords: list[GaussianLink0Command] = []
        for key, value in data.items():
            link0_keywords.append(GaussianLink0Command(key=key, value=value))
        return cls(link0_keywords=link0_keywords)

    @classmethod
    def from_str(cls, data: str) -> GaussianLink0Commands:
        from molop.io.logic.gaussian.input.GaussianInputParsing import (
            parse_gjf_link0_commands,
        )

        return parse_gjf_link0_commands(data)

    def _render(self, **kwargs) -> str:
        return "".join([link0._render() for link0 in self.link0_keywords])

    def memory_request(self) -> PlainQuantity:
        for link0 in self.link0_keywords:
            if (memory_request := link0.memory_request()) is not None:
                return memory_request
        return atom_ureg.Quantity("800MB")

    def explicit_cpu_request(self) -> int | None:
        for link0 in self.link0_keywords:
            if (cpu_request := link0.cpu_request()) is not None:
                return cpu_request
        return None

    def cpu_request(self) -> int:
        return self.explicit_cpu_request() or 1

    def shared_memory_cpu_value(self) -> str | None:
        for link0 in self.link0_keywords:
            if (value := link0.shared_memory_cpu_value()) is not None:
                return value
        return None

    def shared_memory_cpu_render_line(self) -> str | None:
        nproc = self.shared_memory_cpu_value()
        if nproc is None:
            return None
        return (
            f" Will use up to {int(nproc):4d} processors via shared memory."
            if nproc.isdigit()
            else " Will use up to shared-memory processors."
        )

    def add_link0_keyword(self, key: str, value: str | None = None) -> None:
        self.link0_keywords.append(GaussianLink0Command(key=key, value=value))


def render_gaussian_link0_shared_memory_line(raw_link0: str) -> str | None:
    return GaussianLink0Commands.from_str(raw_link0).shared_memory_cpu_render_line()


GJFLink0 = GaussianLink0Command
GJFLink0Commands = GaussianLink0Commands
