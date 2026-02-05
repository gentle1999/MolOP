import subprocess
import sys


def _run_python(code: str) -> str:
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout.strip()


def test_lazy_import_boundary_codec_registry_import_only() -> None:
    out = _run_python(
        "import sys; import molop.io.codec_registry as cr; "
        "assert 'molop.io.codecs.builtin_readers' not in sys.modules; "
        "assert 'molop.io.codecs.builtin_writers' not in sys.modules; "
        "assert 'molop.io.codecs.openbabel_reader' not in sys.modules; "
        "print('ok')"
    )
    assert out == "ok"


def test_lazy_import_boundary_codecs_package_import_only() -> None:
    out = _run_python(
        "import sys; import molop.io.codecs; "
        "assert 'molop.io.codecs.builtin_readers' not in sys.modules; "
        "assert 'molop.io.codecs.builtin_writers' not in sys.modules; "
        "assert 'molop.io.codecs.openbabel_reader' not in sys.modules; "
        "print('ok')"
    )
    assert out == "ok"


def test_registry_deterministic_order_same_priority() -> None:
    from dataclasses import dataclass
    from pathlib import Path
    from typing import Any

    from molop.io.codec_registry import Registry
    from molop.io.codec_types import ParseResult, StructureLevel

    @dataclass
    class DummyReader:
        format_id: str
        extensions: frozenset[str]
        priority: int

        def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
            return ParseResult(
                value=object(),
                level=StructureLevel.COORDS,
                detected_format=self.format_id,
            )

    reg = Registry(autoload_defaults=False)

    @reg.reader_factory(format_id="a", extensions={".x"}, priority=0)
    def make_a():
        return DummyReader("a", frozenset({".x"}), 0)

    @reg.reader_factory(format_id="b", extensions={".x"}, priority=0)
    def make_b():
        return DummyReader("b", frozenset({".x"}), 0)

    codecs = reg.select_reader("file.x")
    assert [c.format_id for c in codecs] == ["a", "b"]


def test_select_reader_dedup_hint_and_extension() -> None:
    from dataclasses import dataclass
    from pathlib import Path
    from typing import Any

    from molop.io.codec_registry import Registry
    from molop.io.codec_types import ParseResult, StructureLevel

    @dataclass
    class DummyReader:
        format_id: str
        extensions: frozenset[str]
        priority: int

        def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
            return ParseResult(
                value=object(),
                level=StructureLevel.COORDS,
                detected_format=self.format_id,
            )

    reg = Registry(autoload_defaults=False)

    @reg.reader_factory(format_id="foo", extensions={".x"}, priority=0)
    def make_foo():
        return DummyReader("foo", frozenset({".x"}), 0)

    codecs = reg.select_reader("file.x", hint_format="foo")
    assert len(codecs) == 1
    assert codecs[0].format_id == "foo"


def test_default_activation_idempotent() -> None:
    from molop.io.codec_registry import default_registry

    default_registry.ensure_default_codecs_registered()
    before = sum(len(v) for v in default_registry._readers_by_format.values())
    default_registry.ensure_default_codecs_registered()
    after = sum(len(v) for v in default_registry._readers_by_format.values())
    assert before == after


def test_lazy_import_boundary_entry_points_not_called_on_import() -> None:
    out = _run_python(
        "import importlib.metadata as md; "
        "md.entry_points = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('called')); "
        "import molop.io.codecs; print('ok')"
    )
    assert out == "ok"


def test_default_activation_calls_entry_points() -> None:
    out = _run_python(
        """
import importlib.metadata as md


class _Eps:
    def select(self, **kw):
        print("select", kw.get("group"))
        return []


def _entry_points(*a, **k):
    print("entry_points")
    return _Eps()


md.entry_points = _entry_points

import molop.io.codec_registry as cr

cr.default_registry.ensure_default_codecs_registered()
print("ok")
"""
    )
    # This is intentionally a substring check because activation may emit
    # additional output in the future.
    assert "entry_points" in out
    assert "select molop.codecs" in out
    assert out.endswith("ok")


def test_catalog_builtin_scan_orders_and_dedups(monkeypatch) -> None:
    from molop.io.codec_registry import Registry
    from molop.io.codecs import catalog

    monkeypatch.setattr(catalog, "_loaded_builtin_registries", set())
    calls: list[str] = []

    monkeypatch.setattr(
        catalog,
        "_discover_builtin_codec_module_names",
        lambda: [
            "molop.io.codecs.readers.b",
            "molop.io.codecs.readers.a",
            "molop.io.codecs.readers.b",
        ],
    )
    monkeypatch.setattr(catalog, "_register_module_codecs", lambda name, reg: calls.append(name))
    monkeypatch.setattr(
        catalog, "_register_openbabel_fallback", lambda reg: calls.append("openbabel")
    )

    reg = Registry(autoload_defaults=False)
    catalog.load_builtin_codecs(reg)
    assert calls == ["molop.io.codecs.readers.a", "molop.io.codecs.readers.b", "openbabel"]


def test_catalog_builtin_scan_missing_register_fails(monkeypatch) -> None:
    import pytest

    from molop.io.codec_registry import Registry
    from molop.io.codecs import catalog

    monkeypatch.setattr(catalog, "_loaded_builtin_registries", set())
    monkeypatch.setattr(catalog, "_discover_builtin_codec_module_names", lambda: ["x.y.z"])
    monkeypatch.setattr(
        catalog,
        "_register_module_codecs",
        lambda name, reg: (_ for _ in ()).throw(AttributeError("no register")),
    )

    reg = Registry(autoload_defaults=False)
    with pytest.raises(AttributeError):
        catalog.load_builtin_codecs(reg)


def test_catalog_plugins_sorted_and_warn_on_failure(monkeypatch) -> None:
    import pytest

    from molop.io.codec_registry import Registry
    from molop.io.codecs import catalog

    monkeypatch.setattr(catalog, "_loaded_plugin_registries", set())
    events: list[str] = []

    def _make_register(label: str):
        def register(reg):
            events.append(f"register:{label}")

        return register

    class _EP:
        def __init__(self, name: str, value: str, register):
            self.name = name
            self.value = value
            self._register = register

        def load(self):
            events.append(f"load:{self.name}")
            return self._register

    class _BadEP:
        name = "bad"
        value = "oops"

        def load(self):
            raise RuntimeError("boom")

    monkeypatch.setattr(
        catalog,
        "_iter_entry_points",
        lambda group: [
            _EP("b", "2", _make_register("b")),
            _BadEP(),
            _EP("a", "1", _make_register("a")),
        ],
    )

    reg = Registry(autoload_defaults=False)
    with pytest.warns(catalog.CodecPluginWarning):
        catalog.load_plugin_codecs(reg)

    assert events == ["load:a", "register:a", "load:b", "register:b"]
