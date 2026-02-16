from __future__ import annotations

import argparse
import ast
import difflib
import subprocess
from dataclasses import dataclass
from pathlib import Path


_CATALOG_PATH = Path("src/molop/io/codecs/catalog.py")
_APP_PATH = Path("src/molop/cli/app.py")
_APP_STUB_PATH = Path("src/molop/cli/app.pyi")

_TRANSFORM_BASE_ARG_DOCS: tuple[tuple[str, str], ...] = (
    ("pattern", "File pattern to match."),
    ("to", "Target writer format."),
    ("output_dir", "Directory to save transformed files."),
    ("frame", "Frame selection (for example, `all`, `-1`, `1,2`)."),
    ("embed", "Whether to embed multiple frames in one file."),
    ("parser_detection", "Parser detection mode."),
    ("n_jobs", "Number of parallel jobs."),
)


@dataclass(frozen=True)
class _Param:
    name: str
    annotation: str
    default: str | None
    has_default: bool


@dataclass(frozen=True)
class _WriterSpec:
    format_id: str
    writer_module: str
    frame_cls_name: str | None


@dataclass(frozen=True)
class _RenderSpec:
    format_id: str
    params: tuple[_Param, ...]
    source: str


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


def _module_to_dir(src_root: Path, module_name: str) -> Path:
    return src_root / Path(*module_name.split("."))


def _collect_function_defs(tree: ast.Module) -> dict[str, ast.FunctionDef | ast.AsyncFunctionDef]:
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _parse_catalog_tree(repo_root: Path) -> ast.Module:
    catalog_path = repo_root / _CATALOG_PATH
    return ast.parse(catalog_path.read_text(encoding="utf-8"), filename=str(catalog_path))


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


def _extract_import_map(tree: ast.Module) -> dict[str, str]:
    out: dict[str, str] = {}
    for node in tree.body:
        if not isinstance(node, ast.ImportFrom) or not node.module:
            continue
        for alias in node.names:
            out[alias.asname or alias.name] = node.module
    return out


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
    class_def = classes.get(class_name)
    if class_def is None:
        return None
    for stmt in class_def.body:
        if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)) and stmt.name == "_render":
            return stmt
    for base in class_def.bases:
        base_name = _base_name(base)
        if base_name and base_name in classes:
            hit = _find_render_def(classes, base_name, seen)
            if hit is not None:
                return hit
    return None


def _simplify_annotation(node: ast.expr | None) -> str:
    if node is None:
        return "Any"
    if isinstance(node, ast.Subscript):
        value = node.value
        if isinstance(value, ast.Name) and value.id == "Annotated":
            sub = node.slice
            if isinstance(sub, ast.Tuple) and sub.elts:
                return ast.unparse(sub.elts[0])
            return ast.unparse(sub)
    return ast.unparse(node)


def _extract_params(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> tuple[_Param, ...]:
    params: list[_Param] = []

    args = [*fn.args.posonlyargs, *fn.args.args]
    if args and args[0].arg == "self":
        args = args[1:]

    defaults = fn.args.defaults
    default_offset = len(args) - len(defaults)
    for idx, arg in enumerate(args):
        ann = _simplify_annotation(arg.annotation)
        default_node = defaults[idx - default_offset] if idx >= default_offset else None
        has_default = default_node is not None
        default = ast.unparse(default_node) if default_node is not None else None
        params.append(
            _Param(name=arg.arg, annotation=ann, default=default, has_default=has_default)
        )

    for arg, default_node in zip(fn.args.kwonlyargs, fn.args.kw_defaults, strict=True):
        ann = _simplify_annotation(arg.annotation)
        has_default = default_node is not None
        default = ast.unparse(default_node) if default_node is not None else None
        params.append(
            _Param(name=arg.arg, annotation=ann, default=default, has_default=has_default)
        )

    return tuple(params)


def _render_spec_from_writer(repo_root: Path, writer: _WriterSpec) -> _RenderSpec:
    writer_path = _resolve_module_to_path(repo_root, writer.writer_module)
    writer_tree = ast.parse(writer_path.read_text(encoding="utf-8"), filename=str(writer_path))
    import_map = _extract_import_map(writer_tree)

    frame_cls_name = writer.frame_cls_name
    if not frame_cls_name:
        source = writer_path.relative_to(repo_root).as_posix()
        return _RenderSpec(format_id=writer.format_id, params=(), source=source)

    frame_mod = import_map.get(frame_cls_name, writer.writer_module)
    frame_path = _resolve_module_to_path(repo_root, frame_mod)
    frame_tree = ast.parse(frame_path.read_text(encoding="utf-8"), filename=str(frame_path))
    classes = _class_map(frame_tree)
    render_def = _find_render_def(classes, frame_cls_name, seen=set())
    if render_def is None:
        source = frame_path.relative_to(repo_root).as_posix()
        return _RenderSpec(format_id=writer.format_id, params=(), source=source)

    source = f"{frame_path.relative_to(repo_root).as_posix()}:{render_def.lineno}"
    return _RenderSpec(
        format_id=writer.format_id, params=_extract_params(render_def), source=source
    )


def _stub_param_decl(param: _Param) -> str:
    if param.has_default:
        return f"{param.name}: {param.annotation} = ...,"
    return f"{param.name}: {param.annotation},"


def _render_doc_block(doc_text: str, indent: str = "    ") -> list[str]:
    escaped = doc_text.replace('"""', r"\"\"\"")
    body = escaped.splitlines()
    out = [f'{indent}"""{body[0]}']
    out.extend(f"{indent}{line}" for line in body[1:])
    out.append(f'{indent}"""')
    return out


def _render_transform_overload_doc(spec: _RenderSpec) -> str:
    lines = [f"Transform files using the `{spec.format_id}` writer.", "", "Args:"]
    lines.extend(f"    {name}: {desc}" for name, desc in _TRANSFORM_BASE_ARG_DOCS)
    lines.extend(["", "Format-specific Args:"])
    if spec.params:
        for param in spec.params:
            optional = " (optional)" if param.has_default else ""
            lines.append(f"    {param.name}: `{param.annotation}`{optional}.")
    else:
        lines.append("    None.")
    lines.extend(["", "Source:", f"    {spec.source}"])
    return "\n".join(lines)


def _render_transform_fallback_doc(specs: list[_RenderSpec]) -> str:
    lines = ["Transform files using a dynamic writer format.", "", "Args:"]
    lines.extend(f"    {name}: {desc}" for name, desc in _TRANSFORM_BASE_ARG_DOCS)
    lines.extend(["", "Format-specific Args:"])
    for spec in specs:
        lines.append(f"    {spec.format_id}:")
        if spec.params:
            for param in spec.params:
                optional = " (optional)" if param.has_default else ""
                lines.append(f"        {param.name}: `{param.annotation}`{optional}.")
        else:
            lines.append("        None.")
    lines.extend(["", "Source:"])
    lines.extend(f"    {spec.format_id}: {spec.source}" for spec in specs)
    return "\n".join(lines)


def _is_name(node: ast.expr, name: str) -> bool:
    return isinstance(node, ast.Name) and node.id == name


def _is_app_attr(node: ast.expr, attr: str) -> bool:
    return isinstance(node, ast.Attribute) and _is_name(node.value, "app") and node.attr == attr


def _has_decorator(func: ast.FunctionDef | ast.AsyncFunctionDef, attr: str) -> bool:
    for dec in func.decorator_list:
        if _is_app_attr(dec, attr):
            return True
        if isinstance(dec, ast.Call) and _is_app_attr(dec.func, attr):
            return True
    return False


def _is_public_name(name: str) -> bool:
    return not name.startswith("_")


def _function_signature_lines(
    fn: ast.FunctionDef | ast.AsyncFunctionDef,
    *,
    return_override: str | None = None,
) -> list[str]:
    lines = [f"def {fn.name}("]

    posonly = list(fn.args.posonlyargs)
    normal = list(fn.args.args)
    defaults = list(fn.args.defaults)
    all_pos = [*posonly, *normal]
    default_offset = len(all_pos) - len(defaults)

    for idx, arg in enumerate(all_pos):
        ann = _simplify_annotation(arg.annotation)
        default_node = defaults[idx - default_offset] if idx >= default_offset else None
        if default_node is None:
            lines.append(f"    {arg.arg}: {ann},")
        else:
            lines.append(f"    {arg.arg}: {ann} = {ast.unparse(default_node)},")
        if idx + 1 == len(posonly) and posonly:
            lines.append("    /,")

    if fn.args.vararg is not None:
        vararg_ann = _simplify_annotation(fn.args.vararg.annotation)
        lines.append(f"    *{fn.args.vararg.arg}: {vararg_ann},")
    elif fn.args.kwonlyargs:
        lines.append("    *,")

    for arg, default_node in zip(fn.args.kwonlyargs, fn.args.kw_defaults, strict=True):
        ann = _simplify_annotation(arg.annotation)
        if default_node is None:
            lines.append(f"    {arg.arg}: {ann},")
        else:
            lines.append(f"    {arg.arg}: {ann} = {ast.unparse(default_node)},")

    if fn.args.kwarg is not None:
        kwarg_ann = _simplify_annotation(fn.args.kwarg.annotation)
        lines.append(f"    **{fn.args.kwarg.arg}: {kwarg_ann},")

    ret = return_override if return_override is not None else _simplify_annotation(fn.returns)
    lines.append(f") -> {ret}: ...")
    return lines


def _render_transform_overloads(
    specs: list[_RenderSpec],
) -> list[str]:
    lines: list[str] = []
    for spec in specs:
        lines.append("@overload")
        lines.append("def transform(")
        lines.append("    pattern: str,")
        lines.append(f'    to: Literal["{spec.format_id}"],')
        lines.append("    output_dir: Path,")
        lines.append('    frame: str = "-1",')
        lines.append("    embed: bool = True,")
        lines.append('    parser_detection: str = "auto",')
        lines.append("    n_jobs: int = -1,")
        if spec.params:
            lines.append("    *,")
            for param in spec.params:
                lines.append(f"    {_stub_param_decl(param)}")
        lines.append("    **kwargs: Any,")
        lines.append(") -> None:")
        lines.extend(_render_doc_block(_render_transform_overload_doc(spec)))
        lines.append("    ...")
        lines.append("")

    lines.append("@overload")
    lines.append("def transform(")
    lines.append("    pattern: str,")
    lines.append("    to: str,")
    lines.append("    output_dir: Path,")
    lines.append('    frame: str = "-1",')
    lines.append("    embed: bool = True,")
    lines.append('    parser_detection: str = "auto",')
    lines.append("    n_jobs: int = -1,")
    lines.append("    **kwargs: Any,")
    lines.append(") -> None:")
    lines.extend(_render_doc_block(_render_transform_fallback_doc(specs)))
    lines.append("    ...")
    return lines


def _collect_public_names(tree: ast.Module) -> set[str]:
    public: set[str] = set()
    for node in tree.body:
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ) and _is_public_name(node.name):
            public.add(node.name)
            continue
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and _is_public_name(target.id):
                    public.add(target.id)
        if (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and _is_public_name(node.target.id)
        ):
            public.add(node.target.id)
    return public


def _render_module_stub(app_tree: ast.Module, specs: list[_RenderSpec]) -> str:
    lines = [
        '"""Static Typer CLI stubs for transform overloads (generated).',
        "",
        "DO NOT EDIT MANUALLY.",
        "Generated by: scripts/generate_cli_transform_stubs.py",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "from pathlib import Path",
        "from typing import Any, Literal, overload",
        "",
        "import typer",
        "",
        "app: typer.Typer",
        "",
    ]

    public_names = _collect_public_names(app_tree)

    for node in app_tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _is_public_name(node.name):
            if node.name == "transform" and _has_decorator(node, "command"):
                lines.extend(_render_transform_overloads(specs))
            else:
                lines.extend(_function_signature_lines(node, return_override="None"))
            lines.append("")
            continue

        if isinstance(node, ast.Assign):
            for target in node.targets:
                if (
                    isinstance(target, ast.Name)
                    and _is_public_name(target.id)
                    and target.id != "app"
                ):
                    lines.append(f"{target.id}: Any")
                    lines.append("")
            continue

        if (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and _is_public_name(node.target.id)
            and node.target.id != "app"
        ):
            ann = _simplify_annotation(node.annotation)
            lines.append(f"{node.target.id}: {ann}")
            lines.append("")

    declared = {
        node.name
        for node in app_tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _is_public_name(node.name)
    }
    missing = sorted(public_names - declared - {"app"})
    for name in missing:
        lines.append(f"{name}: Any")
        lines.append("")

    return "\n".join(lines)


def _format_with_ruff(text: str) -> str:
    result = subprocess.run(
        ["uv", "run", "ruff", "format", "-", "--stdin-filename", "src/molop/cli/app.pyi"],
        input=text,
        text=True,
        capture_output=True,
    )
    return result.stdout if result.returncode == 0 else text


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    src_root = repo_root / "src"

    writer_modules = _discover_writer_modules(repo_root, src_root)
    writer_specs = sorted(
        _extract_writer_specs(src_root, writer_modules),
        key=lambda spec: (spec.writer_module, spec.format_id),
    )

    by_format: dict[str, _WriterSpec] = {}
    for spec in writer_specs:
        by_format.setdefault(spec.format_id, spec)

    render_specs = [
        _render_spec_from_writer(repo_root, spec)
        for spec in sorted(by_format.values(), key=lambda s: s.format_id)
    ]

    app_tree = ast.parse(
        (repo_root / _APP_PATH).read_text(encoding="utf-8"), filename=str(_APP_PATH)
    )
    output_text = _format_with_ruff(_render_module_stub(app_tree, render_specs))

    out_path = (repo_root / _APP_STUB_PATH).resolve()
    existing = out_path.read_text(encoding="utf-8") if out_path.exists() else ""
    if existing == output_text:
        if args.check:
            print(f"OK: {out_path} is up to date")
        return 0

    if args.check:
        diff = "\n".join(
            difflib.unified_diff(
                existing.splitlines(),
                output_text.splitlines(),
                fromfile=str(out_path),
                tofile="(generated)",
                lineterm="",
            )
        )
        if diff:
            print(diff)
        print(
            "ERROR: "
            f"{out_path} is out of date. Run: uv run python scripts/generate_cli_transform_stubs.py"
        )
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(output_text, encoding="utf-8")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
