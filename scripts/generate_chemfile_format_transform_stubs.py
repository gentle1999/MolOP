from __future__ import annotations

import argparse
import ast
import difflib
import subprocess
from dataclasses import dataclass
from pathlib import Path


_CATALOG_PATH = Path("src/molop/io/codecs/catalog.py")


@dataclass(frozen=True)
class _Param:
    name: str
    annotation: str
    default: str | None
    has_default: bool


@dataclass(frozen=True)
class _RenderSpec:
    format_id: str
    params: tuple[_Param, ...]
    doc: str | None
    source: str


@dataclass(frozen=True)
class _WriterSpec:
    format_id: str
    writer_module: str
    frame_cls_name: str | None


def _iter_py_files(dir_path: Path) -> list[Path]:
    return sorted(
        [
            p
            for p in dir_path.glob("*.py")
            if p.is_file() and p.name != "__init__.py" and not p.stem.startswith("_")
        ],
        key=lambda p: p.as_posix(),
    )


def _module_name(src_root: Path, py_file: Path) -> str:
    rel = py_file.relative_to(src_root)
    return ".".join(rel.with_suffix("").parts)


def _extract_import_map(tree: ast.Module) -> dict[str, str]:
    out: dict[str, str] = {}
    for node in tree.body:
        if not isinstance(node, ast.ImportFrom) or not node.module:
            continue
        for alias in node.names:
            out[alias.asname or alias.name] = node.module
    return out


def _module_to_dir(src_root: Path, module_name: str) -> Path:
    return src_root / Path(*module_name.split("."))


def _parse_catalog_tree(repo_root: Path) -> ast.Module:
    catalog_path = repo_root / _CATALOG_PATH
    return ast.parse(catalog_path.read_text(encoding="utf-8"), filename=str(catalog_path))


def _collect_function_defs(tree: ast.Module) -> dict[str, ast.FunctionDef | ast.AsyncFunctionDef]:
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _builtin_scan_package_names(repo_root: Path) -> tuple[str, ...]:
    tree = _parse_catalog_tree(repo_root)
    for node in tree.body:
        value: ast.expr | None = None
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.target.id == "_BUILTIN_SCAN_PACKAGE_NAMES":
                value = node.value
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id == "_BUILTIN_SCAN_PACKAGE_NAMES":
                value = node.value

        if not isinstance(value, ast.Tuple):
            continue

        values: list[str] = []
        for elt in value.elts:
            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                values.append(elt.value)
        return tuple(values)
    return ()


def _explicit_builtin_module_names(repo_root: Path) -> set[str]:
    tree = _parse_catalog_tree(repo_root)
    fn_defs = _collect_function_defs(tree)
    load_builtin = fn_defs.get("load_builtin_codecs")
    if load_builtin is None:
        return set()

    helper_fns: set[str] = set()
    for node in ast.walk(load_builtin):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Name):
            continue
        name = node.func.id
        if name.startswith("_register_") and name != "_register_module_codecs":
            helper_fns.add(name)

    module_names: set[str] = set()
    for helper_name in helper_fns:
        helper_fn = fn_defs.get(helper_name)
        if helper_fn is None:
            continue
        for node in ast.walk(helper_fn):
            if isinstance(node, ast.ImportFrom) and node.module:
                module_names.add(node.module)
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    module_names.add(alias.name)
    return {m for m in module_names if m.startswith("molop.io")}


def _discover_writer_modules(repo_root: Path, src_root: Path) -> list[str]:
    module_names: set[str] = set()

    for pkg_name in _builtin_scan_package_names(repo_root):
        pkg_dir = _module_to_dir(src_root, pkg_name)
        if not pkg_dir.is_dir():
            continue
        for py in _iter_py_files(pkg_dir):
            module_names.add(_module_name(src_root, py))

    for explicit_mod in _explicit_builtin_module_names(repo_root):
        explicit_path = src_root / Path(*explicit_mod.split(".")).with_suffix(".py")
        if (
            explicit_path.is_file()
            and explicit_path.name != "__init__.py"
            and not explicit_path.stem.startswith("_")
        ):
            module_names.add(explicit_mod)

    return sorted(module_names)


def _extract_writer_specs(src_root: Path, writer_modules: list[str]) -> list[_WriterSpec]:
    specs: list[_WriterSpec] = []
    for mod in writer_modules:
        py = (src_root / Path(*mod.split(".")).with_suffix(".py")).resolve()
        tree = ast.parse(py.read_text(encoding="utf-8"), filename=str(py))

        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue

            format_id: str | None = None
            for dec in node.decorator_list:
                if not isinstance(dec, ast.Call):
                    continue
                if not isinstance(dec.func, ast.Attribute) or dec.func.attr != "writer_factory":
                    continue
                for kw in dec.keywords:
                    if (
                        kw.arg == "format_id"
                        and isinstance(kw.value, ast.Constant)
                        and isinstance(kw.value.value, str)
                    ):
                        format_id = kw.value.value

            if not format_id:
                continue

            frame_cls_name: str | None = None
            for inner in ast.walk(node):
                if not isinstance(inner, ast.Call):
                    continue
                if not isinstance(inner.func, ast.Name) or inner.func.id != "FileRendererWriter":
                    continue
                for kw in inner.keywords:
                    if kw.arg == "frame_cls" and isinstance(kw.value, ast.Name):
                        frame_cls_name = kw.value.id

            specs.append(
                _WriterSpec(
                    format_id=format_id,
                    writer_module=mod,
                    frame_cls_name=frame_cls_name,
                )
            )
    return specs


def _resolve_module_to_path(repo_root: Path, module: str) -> Path:
    return repo_root / "src" / Path(*module.split(".")).with_suffix(".py")


def _class_map(tree: ast.Module) -> dict[str, ast.ClassDef]:
    out: dict[str, ast.ClassDef] = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            out[node.name] = node
    return out


def _base_name(expr: ast.expr) -> str | None:
    if isinstance(expr, ast.Name):
        return expr.id
    return None


def _find_render_def(
    classes: dict[str, ast.ClassDef], class_name: str, seen: set[str]
) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    if class_name in seen:
        return None
    seen.add(class_name)
    c = classes.get(class_name)
    if c is None:
        return None
    for stmt in c.body:
        if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)) and stmt.name == "_render":
            return stmt
    for base in c.bases:
        b = _base_name(base)
        if b and b in classes:
            hit = _find_render_def(classes, b, seen)
            if hit is not None:
                return hit
        # Base classes outside this module AST are intentionally ignored.
        # Safe fallback: return None so codegen emits no format-specific params.
    return None


def _extract_params(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> tuple[_Param, ...]:
    params: list[_Param] = []

    args = [*fn.args.posonlyargs, *fn.args.args]
    if args and args[0].arg == "self":
        args = args[1:]

    defaults = fn.args.defaults
    default_offset = len(args) - len(defaults)
    for idx, arg in enumerate(args):
        ann = ast.unparse(arg.annotation) if arg.annotation is not None else "Any"
        default_node = defaults[idx - default_offset] if idx >= default_offset else None
        has_default = default_node is not None
        default = ast.unparse(default_node) if default_node is not None else None
        params.append(
            _Param(name=arg.arg, annotation=ann, default=default, has_default=has_default)
        )

    for arg, default_node in zip(fn.args.kwonlyargs, fn.args.kw_defaults, strict=True):
        ann = ast.unparse(arg.annotation) if arg.annotation is not None else "Any"
        has_default = default_node is not None
        default = ast.unparse(default_node) if default_node is not None else None
        params.append(
            _Param(name=arg.arg, annotation=ann, default=default, has_default=has_default)
        )

    return tuple(params)


def _normalize_doc(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n").strip("\n")


def _format_with_ruff(text: str) -> str:
    res = subprocess.run(
        ["uv", "run", "ruff", "format", "-", "--stdin-filename", "stub.pyi"],
        input=text,
        text=True,
        capture_output=True,
    )
    return res.stdout if res.returncode == 0 else text


def _escape_doc(text: str) -> str:
    return _normalize_doc(text).replace('"""', r"\"\"\"")


def _render_param_decl(p: _Param) -> str:
    if p.has_default:
        return f"{p.name}: {p.annotation} = ...,"
    return f"{p.name}: {p.annotation},"


def _render_doc_block(doc_text: str, indent: str = "        ") -> list[str]:
    escaped = _escape_doc(doc_text)
    body = escaped.splitlines()
    out = [f'{indent}"""{body[0]}']
    out.extend(f"{indent}{line}" for line in body[1:])
    out.append(f'{indent}"""')
    return out


def _file_base_args() -> tuple[tuple[str, str], ...]:
    return (
        ("format", "Target output format literal for this overload."),
        (
            "frameID",
            'Frame selector forwarded to rendering; accepts index, sequence, slice, or `"all"`.',
        ),
        ("file_path", "Optional output path. If omitted, a format-derived path is used."),
        (
            "embed_in_one_file",
            "Whether selected frames are embedded into one file when supported by the format.",
        ),
        (
            "graph_policy",
            'Molecular graph policy used before rendering (for example, `"prefer"`).',
        ),
        ("**kwargs", "Additional writer options accepted by the implementation."),
    )


def _batch_base_args() -> tuple[tuple[str, str], ...]:
    return (
        ("format", "Target output format literal for this overload."),
        (
            "output_dir",
            "Output directory for generated files; defaults to batch behavior when omitted.",
        ),
        ("frameID", "Frame selector forwarded to per-file rendering."),
        (
            "embed_in_one_file",
            "Whether selected frames are embedded into one file when supported by the format.",
        ),
        ("n_jobs", "Parallelism used by batch transform."),
        (
            "graph_policy",
            'Molecular graph policy used before rendering (for example, `"prefer"`).',
        ),
        ("**kwargs", "Additional writer options accepted by the implementation."),
    )


def _format_specific_lines(spec: _RenderSpec) -> list[str]:
    if not spec.params:
        return ["    None."]

    out: list[str] = []
    for p in spec.params:
        default_hint = " (optional)" if p.has_default else ""
        out.append(f"    {p.name}: `{p.annotation}`{default_hint}.")
    return out


def _render_literal_overload_doc(
    *,
    format_id: str,
    spec: _RenderSpec,
    base_args: tuple[tuple[str, str], ...],
    return_type: str,
    return_desc: str,
) -> str:
    lines = [f"Transform using the `{format_id}` writer overload.", "", "Args:"]
    lines.extend(f"    {name}: {desc}" for name, desc in base_args)
    lines.extend(["", "Format-specific Args:"])
    lines.extend(_format_specific_lines(spec))
    lines.extend(
        [
            "",
            "Returns:",
            f"    ({return_type}): {return_desc}",
            "",
            "Source:",
            f"    {spec.source}",
        ]
    )
    return "\n".join(lines)


def _render_combined_overload_doc(
    *,
    base_args: tuple[tuple[str, str], ...],
    specs: list[_RenderSpec],
    return_type: str,
    return_desc: str,
) -> str:
    lines = [
        "Transform using a format string and format-specific writer options.",
        "",
        "Args:",
    ]
    lines.extend(f"    {name}: {desc}" for name, desc in base_args)
    lines.extend(["", "Format-specific Args by format:"])
    for spec in specs:
        lines.append(f"    {spec.format_id}:")
        if spec.params:
            for p in spec.params:
                default_hint = " (optional)" if p.has_default else ""
                lines.append(f"        {p.name}: `{p.annotation}`{default_hint}.")
        else:
            lines.append("        None.")

    lines.extend(["", "Returns:", f"    ({return_type}): {return_desc}", "", "Source:"])
    lines.extend(f"    {spec.format_id}: {spec.source}" for spec in specs)
    return "\n".join(lines)


def _render_spec_from_writer(repo_root: Path, writer: _WriterSpec) -> _RenderSpec:
    writer_path = _resolve_module_to_path(repo_root, writer.writer_module)
    writer_tree = ast.parse(writer_path.read_text(encoding="utf-8"), filename=str(writer_path))
    import_map = _extract_import_map(writer_tree)

    rel_writer_path = writer_path.relative_to(repo_root).as_posix()

    frame_cls_name = writer.frame_cls_name
    if not frame_cls_name:
        return _RenderSpec(
            format_id=writer.format_id,
            params=(),
            doc=None,
            source=f"{rel_writer_path}",
        )

    frame_mod = import_map.get(frame_cls_name, writer.writer_module)
    frame_path = _resolve_module_to_path(repo_root, frame_mod)
    frame_tree = ast.parse(frame_path.read_text(encoding="utf-8"), filename=str(frame_path))
    classes = _class_map(frame_tree)
    render_def = _find_render_def(classes, frame_cls_name, seen=set())

    rel_frame_path = frame_path.relative_to(repo_root).as_posix()

    if render_def is None:
        return _RenderSpec(
            format_id=writer.format_id,
            params=(),
            doc=None,
            source=f"{rel_frame_path}",
        )

    doc = ast.get_docstring(render_def)
    return _RenderSpec(
        format_id=writer.format_id,
        params=_extract_params(render_def),
        doc=_normalize_doc(doc) if doc else None,
        source=f"{rel_frame_path}:{render_def.lineno}",
    )


def _stub_prelude(lines: list[str], title: str, imports: list[str]) -> None:
    lines.extend(
        [
            f'"""{title} (generated).',
            "",
            "DO NOT EDIT MANUALLY.",
            "Generated by: scripts/generate_chemfile_format_transform_stubs.py",
            '"""',
            "",
            "from __future__ import annotations",
            "",
        ]
    )
    lines.extend(imports)
    lines.extend(["", ""])


def _render_file_stub(specs: list[_RenderSpec]) -> str:
    specs = sorted(specs, key=lambda s: s.format_id)
    base_args = _file_base_args()
    return_type = "str | list[str]"
    return_desc = "Output path string or list of output path strings."
    combined_doc = _render_combined_overload_doc(
        base_args=base_args,
        specs=specs,
        return_type=return_type,
        return_desc=return_desc,
    )
    lines: list[str] = []
    _stub_prelude(
        lines,
        "Static file format_transform overload stubs",
        [
            "import os",
            "from collections.abc import Sequence",
            "from typing import Any, Literal, overload",
            "",
            "from molop.io.codec_types import GraphPolicy",
        ],
    )
    lines.append("class FormatTransformMixin:")
    for s in specs:
        lines.append("    @overload")
        lines.append("    def format_transform(")
        lines.append("        self,")
        lines.append(f'        format: Literal["{s.format_id}"],')
        lines.append('        frameID: Sequence[int] | int | Literal["all"] | slice = -1,')
        lines.append("        file_path: os.PathLike | str | None = None,")
        lines.append("        embed_in_one_file: bool = True,")
        lines.append("        *,")
        lines.append('        graph_policy: GraphPolicy = "prefer",')
        for p in s.params:
            lines.append(f"        {_render_param_decl(p)}")
        lines.append("        **kwargs: Any,")
        lines.append("    ) -> str | list[str]:")
        lines.extend(
            _render_doc_block(
                _render_literal_overload_doc(
                    format_id=s.format_id,
                    spec=s,
                    base_args=base_args,
                    return_type=return_type,
                    return_desc=return_desc,
                )
            )
        )
        lines.append("        ...")
    lines.append("    @overload")
    lines.append("    def format_transform(")
    lines.append("        self,")
    lines.append("        format: str,")
    lines.append('        frameID: Sequence[int] | int | Literal["all"] | slice = -1,')
    lines.append("        file_path: os.PathLike | str | None = None,")
    lines.append("        embed_in_one_file: bool = True,")
    lines.append("        *,")
    lines.append('        graph_policy: GraphPolicy = "prefer",')
    lines.append("        **kwargs: Any,")
    lines.append("    ) -> str | list[str]:")
    lines.extend(_render_doc_block(combined_doc))
    lines.append("        ...")
    lines.append("")
    return "\n".join(lines)


def _render_batch_stub(specs: list[_RenderSpec]) -> str:
    specs = sorted(specs, key=lambda s: s.format_id)
    base_args = _batch_base_args()
    return_type = "dict[str, str | list[str]]"
    return_desc = "Mapping from source path to generated output path(s)."
    combined_doc = _render_combined_overload_doc(
        base_args=base_args,
        specs=specs,
        return_type=return_type,
        return_desc=return_desc,
    )
    lines: list[str] = []
    _stub_prelude(
        lines,
        "Static batch format_transform overload stubs",
        [
            "from collections.abc import Sequence",
            "from typing import Any, Literal, overload",
            "",
            "from molop.io.codec_types import GraphPolicy",
        ],
    )
    lines.append("class BatchFormatTransformMixin:")
    for s in specs:
        lines.append("    @overload")
        lines.append("    def format_transform(")
        lines.append("        self,")
        lines.append(f'        format: Literal["{s.format_id}"],')
        lines.append("        output_dir: str | None = None,")
        lines.append('        frameID: int | Literal["all"] | Sequence[int] = -1,')
        lines.append("        embed_in_one_file: bool = True,")
        lines.append("        n_jobs: int = 1,")
        lines.append("        *,")
        lines.append('        graph_policy: GraphPolicy = "prefer",')
        for p in s.params:
            lines.append(f"        {_render_param_decl(p)}")
        lines.append("        **kwargs: Any,")
        lines.append("    ) -> dict[str, str | list[str]]:")
        lines.extend(
            _render_doc_block(
                _render_literal_overload_doc(
                    format_id=s.format_id,
                    spec=s,
                    base_args=base_args,
                    return_type=return_type,
                    return_desc=return_desc,
                )
            )
        )
        lines.append("        ...")
    lines.append("    @overload")
    lines.append("    def format_transform(")
    lines.append("        self,")
    lines.append("        format: str,")
    lines.append("        output_dir: str | None = None,")
    lines.append('        frameID: int | Literal["all"] | Sequence[int] = -1,')
    lines.append("        embed_in_one_file: bool = True,")
    lines.append("        n_jobs: int = 1,")
    lines.append("        *,")
    lines.append('        graph_policy: GraphPolicy = "prefer",')
    lines.append("        **kwargs: Any,")
    lines.append("    ) -> dict[str, str | list[str]]:")
    lines.extend(_render_doc_block(combined_doc))
    lines.append("        ...")
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    src_root = repo_root / "src"
    io_root = src_root / "molop" / "io"

    writer_modules = _discover_writer_modules(repo_root, src_root)
    writer_specs = sorted(
        _extract_writer_specs(src_root, writer_modules),
        key=lambda x: (x.writer_module, x.format_id),
    )
    by_format: dict[str, _WriterSpec] = {}
    for w in writer_specs:
        by_format.setdefault(w.format_id, w)

    render_specs = [
        _render_spec_from_writer(repo_root, w)
        for w in sorted(by_format.values(), key=lambda x: x.format_id)
    ]

    outputs: dict[Path, str] = {
        (io_root / "base_models" / "_format_transform.pyi").resolve(): _format_with_ruff(
            _render_file_stub(render_specs)
        ),
        (io_root / "_batch_format_transform.pyi").resolve(): _format_with_ruff(
            _render_batch_stub(render_specs)
        ),
    }

    dirty: list[Path] = []
    for path, text in outputs.items():
        existing = path.read_text(encoding="utf-8") if path.exists() else ""
        if existing != text:
            dirty.append(path)
            if args.check:
                diff = "\n".join(
                    difflib.unified_diff(
                        existing.splitlines(),
                        text.splitlines(),
                        fromfile=str(path),
                        tofile="(generated)",
                        lineterm="",
                    )
                )
                print(diff)

    if args.check:
        if not dirty:
            for p in sorted(outputs):
                print(f"OK: {p} is up to date")
            return 0
        for p in dirty:
            print(
                "ERROR: "
                f"{p} is out of date. Run: uv run python scripts/generate_chemfile_format_transform_stubs.py"
            )
        return 1

    for path, text in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
