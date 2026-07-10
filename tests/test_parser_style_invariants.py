import ast
import importlib
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from molop.io.base_models.SearchPattern import MolOPPattern


ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
TESTS_ROOT = ROOT / "tests"
LOGIC_ROOT = ROOT / "src" / "molop" / "io" / "logic"
BASE_MODELS_ROOT = ROOT / "src" / "molop" / "io" / "base_models"
GAUSSIAN_LOGIC_ROOT = LOGIC_ROOT / "gaussian"
ORCA_LOGIC_ROOT = LOGIC_ROOT / "orca"
COORDS_LOGIC_ROOT = LOGIC_ROOT / "coords"
GAUSSIAN_INPUT_ROOT = GAUSSIAN_LOGIC_ROOT / "input"
GAUSSIAN_LOG_ROOT = GAUSSIAN_LOGIC_ROOT / "log"
ORCA_INPUT_ROOT = ORCA_LOGIC_ROOT / "input"
ORCA_LOG_ROOT = ORCA_LOGIC_ROOT / "log"
FORMAT_SPECIFIC_BASE_MODEL_FILES = (
    "GaussianInput.py",
    "GaussianInputParsing.py",
    "GaussianInputPatterns.py",
    "GaussianLink0.py",
    "GaussianRoute.py",
    "GaussianRouteParsing.py",
    "ORCA.py",
)
FORBIDDEN_LAYER_DIRECTORIES = (
    "QM_frame_models",
    "QM_frame_parsers",
    "QM_models",
    "QM_parsers",
    "coords_frame_models",
    "coords_frame_parsers",
    "coords_models",
    "coords_parsers",
    "qminput_frame_models",
    "qminput_frame_parsers",
    "qminput_models",
    "qminput_parsers",
)
FORBIDDEN_ROOT_LOGIC_HELPERS = (
    "gaussian_common.py",
    "gaussian_patterns.py",
    "gaussian_route_models.py",
)
FORBIDDEN_INTERMEDIATE_GAUSSIAN_HELPERS = (
    "src/molop/io/logic/qminput_frame_parsers/_gaussian_common.py",
    "src/molop/io/logic/qminput_frame_parsers/_gaussian_route_models.py",
    "src/molop/io/logic/qminput_frame_parsers/_gaussian_route_parser.py",
    "src/molop/io/logic/qminput_frame_parsers/_gjf_patterns.py",
    "src/molop/io/logic/qminput_frame_parsers/_gjf_sections.py",
)
FORBIDDEN_GAUSSIAN_HELPER_IMPORTS = (
    "molop.io.logic.gaussian_common",
    "molop.io.logic.gaussian_patterns",
    "molop.io.logic.gaussian_route_models",
    "molop.io.logic.qminput_frame_parsers._gaussian_common",
    "molop.io.logic.qminput_frame_parsers._gaussian_route_models",
    "molop.io.logic.qminput_frame_parsers._gaussian_route_parser",
    "molop.io.logic.qminput_frame_parsers._gjf_patterns",
    "molop.io.logic.qminput_frame_parsers._gjf_sections",
)


def _project_python_files() -> list[Path]:
    return sorted(
        path
        for root in (SRC_ROOT, TESTS_ROOT)
        for path in root.rglob("*.py")
        if "__pycache__" not in path.parts
    )


def _logic_python_files() -> list[Path]:
    return sorted(path for path in LOGIC_ROOT.rglob("*.py") if "__pycache__" not in path.parts)


def _base_model_python_files() -> list[Path]:
    return sorted(
        path for path in BASE_MODELS_ROOT.rglob("*.py") if "__pycache__" not in path.parts
    )


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def _module_name_from_path(path: Path) -> str:
    return path.relative_to(ROOT / "src").with_suffix("").as_posix().replace("/", ".")


def _molop_pattern_module_names() -> list[str]:
    return [
        _module_name_from_path(path)
        for path in (*_logic_python_files(), *_base_model_python_files())
        if path.name != "SearchPattern.py"
        if "MolOPPattern" in _read(path)
    ]


def _collect_molop_patterns(
    owner: Any, label: str, seen: set[int]
) -> list[tuple[str, MolOPPattern]]:
    owner_id = id(owner)
    if owner_id in seen:
        return []
    seen.add(owner_id)

    if isinstance(owner, MolOPPattern):
        return [(label, owner)]
    if isinstance(owner, Mapping):
        patterns: list[tuple[str, MolOPPattern]] = []
        for key, value in owner.items():
            patterns.extend(_collect_molop_patterns(value, f"{label}[{key!r}]", seen))
        return patterns
    if isinstance(owner, (list, tuple, set, frozenset)):
        patterns = []
        for idx, value in enumerate(owner):
            patterns.extend(_collect_molop_patterns(value, f"{label}[{idx}]", seen))
        return patterns
    if (
        isinstance(owner, type)
        or owner.__class__.__module__.startswith("molop.io.logic")
        or owner.__class__.__module__.startswith("molop.io.base_models")
    ):
        patterns = []
        for attr_name, value in vars(owner).items():
            if attr_name.startswith("_"):
                continue
            patterns.extend(_collect_molop_patterns(value, f"{label}.{attr_name}", seen))
        return patterns
    return []


def _iter_named_patterns() -> list[tuple[str, str, MolOPPattern]]:
    patterns: list[tuple[str, str, MolOPPattern]] = []
    for module_name in _molop_pattern_module_names():
        module = importlib.import_module(module_name)
        for owner_name, owner in vars(module).items():
            if owner_name.startswith("_"):
                continue
            for pattern_name, pattern in _collect_molop_patterns(owner, owner_name, set()):
                patterns.append((module_name, pattern_name, pattern))
    return patterns


def _logic_frame_parser_files() -> list[Path]:
    return sorted(
        path
        for path in LOGIC_ROOT.rglob("*FrameParser.py")
        if path.parent.name == "frame_parsers" and "__pycache__" not in path.parts
    )


def _state_machine_frame_parser_files() -> list[Path]:
    state_machine_roots = (
        GAUSSIAN_INPUT_ROOT / "frame_parsers",
        GAUSSIAN_LOG_ROOT / "frame_parsers",
        ORCA_INPUT_ROOT / "frame_parsers",
        ORCA_LOG_ROOT / "frame_parsers",
    )
    return sorted(
        path
        for root in state_machine_roots
        for path in root.rglob("*FileFrameParser.py")
        if "__pycache__" not in path.parts
    )


def _logic_file_parser_files() -> list[Path]:
    return sorted(
        path
        for path in LOGIC_ROOT.rglob("*FileParser.py")
        if path.parent.name == "parsers" and "__pycache__" not in path.parts
    )


def _state_machine_file_parser_files() -> list[Path]:
    return sorted(
        path
        for root in (GAUSSIAN_LOG_ROOT / "parsers", ORCA_LOG_ROOT / "parsers")
        for path in root.rglob("*FileParser.py")
        if "__pycache__" not in path.parts
    )


def _state_machine_file_split_parser_files() -> list[Path]:
    return sorted(
        path
        for root in (GAUSSIAN_INPUT_ROOT / "parsers", ORCA_INPUT_ROOT / "parsers")
        for path in root.rglob("*FileParser.py")
        if "__pycache__" not in path.parts
    )


def _logic_model_files() -> list[Path]:
    return sorted(
        path
        for path in LOGIC_ROOT.rglob("*.py")
        if path.parent.name in {"frame_models", "models"} and "__pycache__" not in path.parts
    )


def _base_model_files() -> list[Path]:
    return sorted(
        path
        for path in (SRC_ROOT / "molop" / "io" / "base_models").rglob("*.py")
        if "__pycache__" not in path.parts
    )


def _decorator_name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Call):
        return _decorator_name(node.func)
    return ""


def _is_parser_helper_module(module: str | None) -> bool:
    return module is not None and (".parsers" in module or "_parsers" in module)


def _iter_parser_helper_imports(
    node: ast.AST, function_stack: tuple[str, ...] = ()
) -> list[tuple[ast.ImportFrom, tuple[str, ...]]]:
    scoped_imports: list[tuple[ast.ImportFrom, tuple[str, ...]]] = []
    if isinstance(node, ast.FunctionDef):
        function_stack = (*function_stack, node.name)
    if isinstance(node, ast.ImportFrom) and _is_parser_helper_module(node.module):
        scoped_imports.append((node, function_stack))
    for child in ast.iter_child_nodes(node):
        scoped_imports.extend(_iter_parser_helper_imports(child, function_stack))
    return scoped_imports


def test_logic_layer_does_not_compile_regex_directly() -> None:
    banned_patterns = [
        re.compile(r"^\s*import\s+(?:re|regex)\b", re.MULTILINE),
        re.compile(r"^\s*from\s+(?:re|regex)\s+import\b", re.MULTILINE),
        re.compile(r"\b(?:re|regex)\.compile\("),
    ]
    violations: list[str] = []

    for path in _logic_python_files():
        text = _read(path)
        for pattern in banned_patterns:
            if pattern.search(text):
                violations.append(_rel(path))
                break

    assert violations == []


def test_logic_layer_uses_named_regex_groups() -> None:
    banned_patterns = [
        re.compile(r"\.group\(\s*\d+"),
        re.compile(r"\.groups\(\s*\)"),
        re.compile(r"\b(?:match|matched)\s*\[\s*\d+\s*\]"),
    ]
    violations: list[str] = []

    for path in _logic_python_files():
        text = _read(path)
        for pattern in banned_patterns:
            if pattern.search(text):
                violations.append(_rel(path))
                break

    assert violations == []


def test_logic_layer_uses_single_molop_pattern_class() -> None:
    violations: list[str] = []

    for path in _logic_python_files():
        text = _read(path)
        if "MolOPPatternV2" in text:
            violations.append(_rel(path))

    assert violations == []


def test_molop_pattern_v2_is_only_a_compatibility_alias() -> None:
    from molop.io.base_models.SearchPattern import MolOPPatternV2

    assert MolOPPatternV2 is MolOPPattern
    assert "class MolOPPatternV2" not in _read(BASE_MODELS_ROOT / "SearchPattern.py")


def test_base_models_do_not_depend_on_logic_layer() -> None:
    violations: list[str] = []

    for path in _base_model_python_files():
        tree = ast.parse(_read(path), filename=_rel(path))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.ImportFrom)
                and node.module is not None
                and node.module.startswith("molop.io.logic")
            ):
                violations.append(f"{_rel(path)} imports {node.module}")
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("molop.io.logic"):
                        violations.append(f"{_rel(path)} imports {alias.name}")

    assert violations == []


def test_format_specific_containers_do_not_live_in_base_models() -> None:
    violations = [
        _rel(BASE_MODELS_ROOT / filename)
        for filename in FORMAT_SPECIFIC_BASE_MODEL_FILES
        if (BASE_MODELS_ROOT / filename).exists()
    ]

    for required_path in (
        GAUSSIAN_INPUT_ROOT / "GaussianInput.py",
        GAUSSIAN_INPUT_ROOT / "GaussianRoute.py",
        ORCA_LOGIC_ROOT / "common.py",
    ):
        if not required_path.exists():
            violations.append(f"{_rel(required_path)} missing")

    assert violations == []


def test_molop_patterns_do_not_define_positional_capture_groups() -> None:
    violations: list[str] = []

    for module_name, pattern_name, pattern in _iter_named_patterns():
        compiled = pattern.content_pattern_compiled
        if compiled is None:
            continue
        if compiled.groups != len(compiled.groupindex):
            violations.append(
                f"{module_name}:{pattern_name} has {compiled.groups} captures "
                f"but {len(compiled.groupindex)} named groups"
            )

    assert violations == []


def test_retired_gaussian_helpers_stay_removed() -> None:
    violations: list[str] = []

    for filename in FORBIDDEN_ROOT_LOGIC_HELPERS:
        if (LOGIC_ROOT / filename).exists():
            violations.append(f"retired root helper exists: {filename}")
    for filename in FORBIDDEN_INTERMEDIATE_GAUSSIAN_HELPERS:
        if (ROOT / filename).exists():
            violations.append(f"retired intermediate helper exists: {filename}")

    banned_imports = set(FORBIDDEN_GAUSSIAN_HELPER_IMPORTS)
    for path in _project_python_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module in banned_imports:
                violations.append(f"{_rel(path)} imports retired Gaussian helper")
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name in banned_imports:
                        violations.append(f"{_rel(path)} imports retired Gaussian helper")

    assert violations == []


def test_logic_root_keeps_format_specific_helpers_out() -> None:
    violations = [
        _rel(path)
        for path in LOGIC_ROOT.glob("*.py")
        if path.name != "__init__.py" and "__pycache__" not in path.parts
    ]
    for dirname in FORBIDDEN_LAYER_DIRECTORIES:
        path = LOGIC_ROOT / dirname
        if path.exists():
            violations.append(f"{_rel(path)} should stay retired")

    assert violations == []


def test_frame_parsers_use_model_parse_result_boundary() -> None:
    violations: list[str] = []

    for path in _logic_frame_parser_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or not node.name.endswith("ParserMixin"):
                continue
            methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
            parse_block = methods.get("_parse_block_to_result")
            parse_frame = methods.get("_parse_frame")
            if parse_block is None:
                violations.append(f"{_rel(path)}:{node.name} missing _parse_block_to_result")
                continue
            returns = ast.unparse(parse_block.returns) if parse_block.returns is not None else ""
            if returns != "ModelParseResult":
                violations.append(
                    f"{_rel(path)}:{node.name} _parse_block_to_result returns {returns!r}"
                )
            if parse_frame is None:
                violations.append(f"{_rel(path)}:{node.name} missing _parse_frame")
                continue
            source = ast.get_source_segment(text, parse_frame) or ""
            if "._parse_block_to_result(" not in source or ".model_data()" not in source:
                violations.append(
                    f"{_rel(path)}:{node.name} _parse_frame must project ModelParseResult"
                )

    assert violations == []


def test_complex_frame_parsers_use_explicit_phase_machine() -> None:
    violations: list[str] = []

    for path in _state_machine_frame_parser_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        phase_classes = [
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef)
            and node.name.endswith("ParsePhase")
            and any(ast.unparse(base) == "Enum" for base in node.bases)
        ]
        if not phase_classes:
            violations.append(f"{_rel(path)} missing *ParsePhase Enum")
            continue
        if not any(
            isinstance(item, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "DONE" for target in item.targets)
            for phase_class in phase_classes
            for item in phase_class.body
        ):
            violations.append(f"{_rel(path)} *ParsePhase missing DONE member")

        function_phase_runners = [
            node.name
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name.startswith("_run_")
            and node.name.endswith("_phase")
        ]
        result_functions = [
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name.startswith("parse_")
            and node.name.endswith("_result")
        ]
        if result_functions and function_phase_runners:
            if len(function_phase_runners) < 2:
                violations.append(f"{_rel(path)} has too few phase runners")
            for result_function in result_functions:
                source = ast.get_source_segment(text, result_function) or ""
                if "while phase is not" not in source or ".DONE" not in source:
                    violations.append(
                        f"{_rel(path)}:{result_function.name} does not loop until DONE"
                    )
                if not any(f"{name}(" in source for name in function_phase_runners):
                    violations.append(
                        f"{_rel(path)}:{result_function.name} does not dispatch phase runners"
                    )
                if "raise AssertionError" not in source:
                    violations.append(
                        f"{_rel(path)}:{result_function.name} does not reject unexpected phases"
                    )
            continue

        found_mixin = False
        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or not node.name.endswith("ParserMixin"):
                continue
            found_mixin = True
            methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
            parse_block = methods.get("_parse_block_to_result")
            if parse_block is None:
                continue
            source = ast.get_source_segment(text, parse_block) or ""
            run_phase_methods = [
                name for name in methods if name.startswith("_run_") and name.endswith("_phase")
            ]
            if len(run_phase_methods) < 2:
                violations.append(f"{_rel(path)}:{node.name} has too few phase runners")
            if "while phase is not" not in source or ".DONE" not in source:
                violations.append(f"{_rel(path)}:{node.name} does not loop until DONE")
            if not any(f"self.{name}(" in source for name in run_phase_methods):
                violations.append(f"{_rel(path)}:{node.name} does not dispatch phase runners")
            if "raise AssertionError" not in source:
                violations.append(f"{_rel(path)}:{node.name} does not reject unexpected phases")
        if not found_mixin:
            violations.append(f"{_rel(path)} missing ParserMixin or parse_*_result function")

    assert violations == []


def test_file_parsers_use_metadata_result_boundary() -> None:
    violations: list[str] = []

    for path in _logic_file_parser_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or not node.name.endswith("ParserMixin"):
                continue
            methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
            parse_metadata_result = methods.get("_parse_metadata_result")
            parse_metadata = methods.get("_parse_metadata")
            if parse_metadata_result is None:
                violations.append(f"{_rel(path)}:{node.name} missing _parse_metadata_result")
                continue
            returns = (
                ast.unparse(parse_metadata_result.returns)
                if parse_metadata_result.returns is not None
                else ""
            )
            if returns != "ModelParseResult":
                violations.append(
                    f"{_rel(path)}:{node.name} _parse_metadata_result returns {returns!r}"
                )
            if parse_metadata is None:
                violations.append(f"{_rel(path)}:{node.name} missing _parse_metadata")
                continue
            source = ast.get_source_segment(text, parse_metadata) or ""
            if "._parse_metadata_result(" not in source or ".model_data()" not in source:
                violations.append(
                    f"{_rel(path)}:{node.name} _parse_metadata must project ModelParseResult"
                )

    assert violations == []


def test_complex_file_parsers_use_explicit_metadata_phase_machine() -> None:
    violations: list[str] = []

    for path in _state_machine_file_parser_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        phase_classes = [
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef)
            and node.name.endswith("MetadataParsePhase")
            and any(ast.unparse(base) == "Enum" for base in node.bases)
        ]
        if not phase_classes:
            violations.append(f"{_rel(path)} missing *MetadataParsePhase Enum")
            continue
        if not any(
            isinstance(item, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "DONE" for target in item.targets)
            for phase_class in phase_classes
            for item in phase_class.body
        ):
            violations.append(f"{_rel(path)} *MetadataParsePhase missing DONE member")

        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or not node.name.endswith("ParserMixin"):
                continue
            methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
            parse_metadata_result = methods.get("_parse_metadata_result")
            if parse_metadata_result is None:
                continue
            source = ast.get_source_segment(text, parse_metadata_result) or ""
            run_phase_methods = [
                name
                for name in methods
                if name.startswith("_run_") and name.endswith("_metadata_phase")
            ]
            if len(run_phase_methods) < 2:
                violations.append(f"{_rel(path)}:{node.name} has too few metadata phase runners")
            if "while phase is not" not in source or ".DONE" not in source:
                violations.append(f"{_rel(path)}:{node.name} does not loop metadata until DONE")
            if not any(f"self.{name}(" in source for name in run_phase_methods):
                violations.append(f"{_rel(path)}:{node.name} does not dispatch metadata runners")
            if "raise AssertionError" not in source:
                violations.append(
                    f"{_rel(path)}:{node.name} does not reject unexpected metadata phases"
                )

    assert violations == []


def test_complex_file_splitters_use_explicit_phase_machine() -> None:
    violations: list[str] = []

    for path in _state_machine_file_split_parser_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        phase_classes = [
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef)
            and node.name.endswith("FileSplitPhase")
            and any(ast.unparse(base) == "Enum" for base in node.bases)
        ]
        if not phase_classes:
            violations.append(f"{_rel(path)} missing *FileSplitPhase Enum")
            continue
        if not any(
            isinstance(item, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "DONE" for target in item.targets)
            for phase_class in phase_classes
            for item in phase_class.body
        ):
            violations.append(f"{_rel(path)} *FileSplitPhase missing DONE member")

        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or not node.name.endswith("ParserMixin"):
                continue
            methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
            split_file = methods.get("_split_file")
            if split_file is None:
                continue
            source = ast.get_source_segment(text, split_file) or ""
            run_phase_methods = [
                name
                for name in methods
                if name.startswith("_run_") and name.endswith("_split_phase")
            ]
            if len(run_phase_methods) < 2:
                violations.append(f"{_rel(path)}:{node.name} has too few split phase runners")
            if "while phase is not" not in source or ".DONE" not in source:
                violations.append(f"{_rel(path)}:{node.name} does not loop split until DONE")
            if not any(f"self.{name}(" in source for name in run_phase_methods):
                violations.append(f"{_rel(path)}:{node.name} does not dispatch split runners")
            if "raise AssertionError" not in source:
                violations.append(
                    f"{_rel(path)}:{node.name} does not reject unexpected split phases"
                )

    assert violations == []


def test_model_validators_do_not_parse_raw_sections() -> None:
    violations: list[str] = []

    for path in _logic_model_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            decorator_names = {_decorator_name(decorator) for decorator in node.decorator_list}
            if not decorator_names.intersection({"model_validator", "field_validator"}):
                continue

            for child in ast.walk(node):
                if isinstance(child, ast.ImportFrom) and _is_parser_helper_module(child.module):
                    violations.append(f"{_rel(path)}:{node.name} imports parser helper")
                if (
                    isinstance(child, ast.Call)
                    and isinstance(child.func, ast.Attribute)
                    and child.func.attr == "from_str"
                ):
                    violations.append(f"{_rel(path)}:{node.name} calls from_str")

    assert violations == []


def test_model_parser_helper_imports_stay_explicit_factories() -> None:
    violations: list[str] = []

    for path in _logic_model_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node, function_stack in _iter_parser_helper_imports(tree):
            if not function_stack or function_stack[-1] != "from_str":
                location = f"{_rel(path)}:{node.lineno}"
                scope = ".".join(function_stack) if function_stack else "<module>"
                violations.append(f"{location} parser helper import in {scope}")

    assert violations == []


def test_base_model_parser_helper_imports_stay_explicit_factories() -> None:
    violations: list[str] = []

    for path in _base_model_files():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        for node, function_stack in _iter_parser_helper_imports(tree):
            if not function_stack or function_stack[-1] != "from_str":
                location = f"{_rel(path)}:{node.lineno}"
                scope = ".".join(function_stack) if function_stack else "<module>"
                violations.append(f"{location} parser helper import in {scope}")

    assert violations == []


def test_gaussian_models_use_shared_semantic_projection_helper() -> None:
    violations: list[str] = []
    gaussian_route_model = GAUSSIAN_INPUT_ROOT / "GaussianRoute.py"
    link0_model = GAUSSIAN_INPUT_ROOT / "GaussianLink0.py"
    g16_models = [
        GAUSSIAN_LOG_ROOT / "frame_models" / "G16LogFileFrame.py",
        GAUSSIAN_LOG_ROOT / "models" / "G16LogFile.py",
    ]

    for path in _logic_model_files():
        text = _read(path)
        if "populate_common_gaussian_qm_containers" in text:
            violations.append(f"{_rel(path)} calls low-level Gaussian container projection")

    base_text = _read(gaussian_route_model)
    if "class GaussianRouteSemanticFieldsMixin" not in base_text:
        violations.append(
            f"{_rel(gaussian_route_model)} missing Gaussian route semantic fields mixin"
        )
    if (
        "populate_gaussian_legacy_qm_fields_from_semantic(self, self.semantic_route)"
        not in base_text
    ):
        violations.append(
            f"{_rel(gaussian_route_model)} common Gaussian mixin does not normalize fields"
        )
    link0_text = _read(link0_model) if link0_model.exists() else ""
    if "def render_gaussian_link0_shared_memory_line" not in link0_text:
        violations.append(f"{_rel(link0_model)} missing shared Link0 render projection helper")

    for path in g16_models:
        text = _read(path)
        if "GaussianRouteSemanticFieldsMixin" not in text:
            violations.append(f"{_rel(path)} does not inherit shared Gaussian route mixin")
        if "populate_gaussian_legacy_qm_fields_from_semantic" in text:
            violations.append(f"{_rel(path)} calls Gaussian route projection directly")
        if path.name == "G16LogFile.py":
            if "render_gaussian_link0_shared_memory_line" not in text:
                violations.append(f"{_rel(path)} does not reuse shared Link0 projection helper")
            if "%nprocshared=" in text or 'split("=", 1)' in text:
                violations.append(f"{_rel(path)} parses Link0 options directly")
        for duplicated_fragment in (
            'qm_software: str = Field(default="Gaussian")',
            'options: str = Field(default="", description="options comment")',
            'title_card: str = Field(default="", description="title card")',
            'job_type: str = Field(default="", description="Job type")',
            "semantic_route: GaussianRouteSemantic = Field",
            "def dieze_tag",
        ):
            if duplicated_fragment in text:
                violations.append(f"{_rel(path)} duplicates Gaussian route semantic fields")

    gjf_frame_model = GAUSSIAN_INPUT_ROOT / "frame_models" / "GJFFileFrame.py"
    if "populate_gaussian_legacy_qm_fields_from_semantic" not in _read(gjf_frame_model):
        violations.append(f"{_rel(gjf_frame_model)} does not project GJF route_section semantics")

    assert violations == []


def test_gjf_frame_model_uses_gaussian_geometry_payload_helper() -> None:
    violations: list[str] = []
    frame_model = GAUSSIAN_INPUT_ROOT / "frame_models" / "GJFFileFrame.py"
    gaussian_input_model = GAUSSIAN_INPUT_ROOT / "GaussianInput.py"
    frame_text = _read(frame_model)
    base_text = _read(gaussian_input_model)

    if "normalize_gjf_frame_geometry_payload" not in frame_text:
        violations.append(f"{_rel(frame_model)} does not call Gaussian geometry payload helper")
    if "def build_gjf_molecule_specifications_from_frame_payload" not in base_text:
        violations.append(f"{_rel(gaussian_input_model)} missing GJF molecule payload builder")

    banned_fragments = {
        "_build_molecule_specifications_from_payload": "still owns molecule payload builder",
        "import numpy as np": "still owns NumPy frame geometry normalization",
        "from rdkit import Chem": "still owns RDKit element lookup",
        "molop.unit import atom_ureg": "still owns unit-bearing geometry normalization",
        "pt = Chem.GetPeriodicTable()": "still caches RDKit periodic table",
    }
    for fragment, message in banned_fragments.items():
        if fragment in frame_text:
            violations.append(f"{_rel(frame_model)} {message}")

    assert violations == []


def test_gjf_frame_render_overrides_use_gaussian_helper() -> None:
    violations: list[str] = []
    frame_model = GAUSSIAN_INPUT_ROOT / "frame_models" / "GJFFileFrame.py"
    gaussian_input_model = GAUSSIAN_INPUT_ROOT / "GaussianInput.py"
    frame_text = _read(frame_model)
    base_text = _read(gaussian_input_model)
    tree = ast.parse(frame_text, filename=_rel(frame_model))

    if "adapt_gjf_writer_payload" not in frame_text:
        violations.append(f"{_rel(frame_model)} does not call GJF writer payload helper")
    if "def adapt_gjf_writer_payload" not in base_text:
        violations.append(f"{_rel(gaussian_input_model)} missing GJF writer payload helper")
    if "resolve_gjf_render_parts" not in frame_text:
        violations.append(f"{_rel(frame_model)} does not call GJF render override helper")
    if "def resolve_gjf_render_parts" not in base_text:
        violations.append(f"{_rel(gaussian_input_model)} missing GJF render override helper")

    for node in ast.walk(tree):
        if not isinstance(node, ast.ClassDef) or node.name != "GJFFileFrameMixin":
            continue
        methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
        adapt_method = methods.get("adapt_writer_payload")
        if adapt_method is None:
            violations.append(f"{_rel(frame_model)}:GJFFileFrameMixin missing adapt_writer_payload")
        else:
            adapt_source = ast.get_source_segment(frame_text, adapt_method) or ""
            if ".from_str(" in adapt_source:
                violations.append(f"{_rel(frame_model)}:adapt_writer_payload parses raw sections")
        render_method = methods.get("_render")
        if render_method is None:
            violations.append(f"{_rel(frame_model)}:GJFFileFrameMixin missing _render")
            break
        render_source = ast.get_source_segment(frame_text, render_method) or ""
        if ".from_str(" in render_source:
            violations.append(f"{_rel(frame_model)}:_render parses raw sections directly")
        if "GJFLink0Commands.from_dict(" in render_source:
            violations.append(f"{_rel(frame_model)}:_render normalizes Link0 dict directly")
        break

    assert violations == []


def test_g16_component_tree_stays_model_derived_view() -> None:
    violations: list[str] = []
    frame_parser = GAUSSIAN_LOG_ROOT / "frame_parsers" / "G16LogFileFrameParser.py"
    frame_model = GAUSSIAN_LOG_ROOT / "frame_models" / "G16LogFileFrame.py"
    parser_text = _read(frame_parser)
    model_text = _read(frame_model)
    tree = ast.parse(model_text, filename=_rel(frame_model))

    if "G16ComponentTreeBuilder" in parser_text or "G16Components" in parser_text:
        violations.append(f"{_rel(frame_parser)} depends on G16 component tree construction")

    if "G16ComponentTreeBuilder.from_frame_data(self)" not in model_text:
        violations.append(f"{_rel(frame_model)} does not build component tree from frame data")

    for node in ast.walk(tree):
        if not isinstance(node, ast.ClassDef) or node.name != "G16LogFileFrameMixin":
            continue
        methods = {item.name: item for item in node.body if isinstance(item, ast.FunctionDef)}
        component_tree = methods.get("component_tree")
        render_fakeg = methods.get("render_fakeg")
        if component_tree is None:
            violations.append(f"{_rel(frame_model)} missing component_tree property")
        if render_fakeg is None:
            violations.append(f"{_rel(frame_model)} missing render_fakeg")
            break
        render_source = ast.get_source_segment(model_text, render_fakeg) or ""
        if "G16ComponentTreeBuilder.from_frame_data" in render_source:
            violations.append(f"{_rel(frame_model)}:render_fakeg bypasses component_tree view")
        if "self.component_tree.render_fakeg" not in render_source:
            violations.append(f"{_rel(frame_model)}:render_fakeg does not use component_tree")
        break

    assert violations == []


def test_orca_models_use_shared_common_projection_helper() -> None:
    violations: list[str] = []
    orca_model = ORCA_LOGIC_ROOT / "common.py"
    output_models = [
        ORCA_LOG_ROOT / "frame_models" / "ORCALogFileFrame.py",
        ORCA_LOG_ROOT / "models" / "ORCALogFile.py",
    ]
    input_models = [
        ORCA_INPUT_ROOT / "frame_models" / "ORCAInpFileFrame.py",
    ]

    base_text = _read(orca_model)
    if "class ORCACommonQMFieldsMixin" not in base_text:
        violations.append(f"{_rel(orca_model)} missing ORCA common QM fields mixin")
    if "class ORCAOutputQMFieldsMixin(ORCACommonQMFieldsMixin)" not in base_text:
        violations.append(f"{_rel(orca_model)} missing ORCA output QM fields mixin")
    if "populate_common_orca_qm_fields(self)" not in base_text:
        violations.append(f"{_rel(orca_model)} common ORCA mixin does not normalize fields")

    for path in output_models:
        text = _read(path)
        if "ORCAOutputQMFieldsMixin" not in text:
            violations.append(f"{_rel(path)} does not inherit shared ORCA output mixin")
        if "backfill_common_qm_containers_from_legacy(" in text:
            violations.append(f"{_rel(path)} calls common QM backfill directly")
        if "project_common_qm_fields(" in text:
            violations.append(f"{_rel(path)} calls common QM projection directly")
        if "populate_common_orca_qm_fields" in text:
            violations.append(f"{_rel(path)} calls ORCA common projection directly")
        for duplicated_fragment in (
            "input_file_name: str = Field",
            "auxiliary_basis_set: str = Field",
            "dispersion_correction: str = Field",
            "def _normalize_common_orca_fields",
        ):
            if duplicated_fragment in text:
                violations.append(f"{_rel(path)} duplicates ORCA common fields")

    for path in input_models:
        text = _read(path)
        if "ORCACommonQMFieldsMixin" not in text:
            violations.append(f"{_rel(path)} does not inherit shared ORCA common mixin")
        for duplicated_fragment in (
            "auxiliary_basis_set: str = Field",
            "dispersion_correction: str = Field",
            "def _normalize_common_orca_fields",
        ):
            if duplicated_fragment in text:
                violations.append(f"{_rel(path)} duplicates ORCA common fields")

    assert violations == []


def test_orca_input_frame_model_uses_orca_geometry_projection_helper() -> None:
    violations: list[str] = []
    frame_model = ORCA_INPUT_ROOT / "frame_models" / "ORCAInpFileFrame.py"
    orca_model = ORCA_LOGIC_ROOT / "common.py"
    frame_text = _read(frame_model)
    base_text = _read(orca_model)

    if "project_orca_geometry_to_qm_frame" not in frame_text:
        violations.append(f"{_rel(frame_model)} does not call ORCA geometry projection helper")
    if "def project_orca_geometry_to_qm_frame" not in base_text:
        violations.append(f"{_rel(orca_model)} missing ORCA geometry projection helper")

    banned_fragments = {
        "import numpy as np": "still owns NumPy coordinate projection",
        "from rdkit import Chem": "still owns RDKit element lookup",
        "molop.unit import atom_ureg": "still owns unit-bearing geometry projection",
        "GetPeriodicTable": "still resolves atomic numbers directly",
        "internal_coords.to_cartesian_coords": "still converts internal geometry directly",
        "geometry.real_atoms": "still scans ORCA geometry atoms directly",
    }
    for fragment, message in banned_fragments.items():
        if fragment in frame_text:
            violations.append(f"{_rel(frame_model)} {message}")

    assert violations == []


def test_orca_input_parser_imports_only_frame_classes_from_frame_model() -> None:
    violations: list[str] = []
    allowed_frame_model_imports = {"ORCAInpFileFrameDisk", "ORCAInpFileFrameMemory"}

    for path in (ORCA_INPUT_ROOT / "frame_parsers").rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(_read(path), filename=_rel(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.module != "molop.io.logic.orca.input.frame_models.ORCAInpFileFrame":
                continue
            imported_names = {alias.name for alias in node.names}
            disallowed = imported_names - allowed_frame_model_imports
            if disallowed:
                violations.append(
                    f"{_rel(path)} imports ORCA format containers from frame model: "
                    f"{sorted(disallowed)}"
                )

    assert violations == []


def test_orca_input_parser_uses_single_semantic_payload_helper() -> None:
    path = ORCA_INPUT_ROOT / "frame_parsers" / "ORCAInpFileFrameParser.py"
    text = _read(path)
    tree = ast.parse(text, filename=_rel(path))
    violations: list[str] = []

    if "build_orca_input_semantic_payload(" not in text:
        violations.append(f"{_rel(path)} does not call ORCA semantic payload helper")

    disallowed_imports = {
        "_build_excited_state_semantic",
        "_build_multi_reference_semantic",
        "_build_orca_excited_state_requests",
        "_build_orca_explicit_solvent_requests",
        "_build_orca_model_chemistry",
        "_build_orca_multireference_requests",
        "_build_orca_task_requests",
        "_derive_keyword_semantics",
        "_derive_semantic_payload",
        "_parse_orca_explicit_solvent_semantics",
        "_resolve_method_for_excited_state",
    }
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.module != "molop.io.logic.orca.input.frame_parsers._orca_inp_semantics":
            continue
        imported_names = {alias.name for alias in node.names}
        disallowed = imported_names.intersection(disallowed_imports)
        if disallowed:
            violations.append(
                f"{_rel(path)} imports ORCA semantic helper internals: {sorted(disallowed)}"
            )

    assert violations == []


def test_orca_input_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    parser = ORCA_INPUT_ROOT / "frame_parsers" / "ORCAInpFileFrameParser.py"
    result_helper = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_result.py"
    retired_extractor = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_extractors.py"
    block_helpers = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_blocks.py"
    geometry_helpers = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_geometry.py"
    resource_helpers = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_resources.py"
    required_calls = {
        "extract_star_geometry",
        "extract_percent_coords_geometry",
        "extract_block_spans",
        "geometry_from_section",
        "render_resource_blocks",
        "parse_request_num_cpu",
        "parse_request_memory",
        "parse_output_print_settings",
    }

    for helper_path in (block_helpers, geometry_helpers, resource_helpers):
        if not helper_path.exists():
            violations.append(f"{_rel(helper_path)} missing")
    if result_helper.exists():
        violations.append(f"{_rel(result_helper)} should stay retired")
    if retired_extractor.exists():
        violations.append(f"{_rel(retired_extractor)} should stay retired")

    text = _read(parser)
    if "parse_orca_input_frame_result" not in text:
        violations.append(f"{_rel(parser)} missing ORCA input result helper")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_result" in text:
        violations.append(f"{_rel(parser)} imports ORCA input result compatibility facade")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_extractors" in text:
        violations.append(f"{_rel(parser)} imports ORCA input compatibility extractor")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_patterns" in text:
        violations.append(f"{_rel(parser)} imports ORCA input patterns directly")
    if "molop.io.base_models.DataClasses" in text:
        violations.append(f"{_rel(parser)} imports ORCA input low-level data classes directly")

    helper_defs = {
        "_parse_parameter_blocks",
        "_parse_cartesian_coordinate_lines",
        "_parse_internal_coordinate_lines",
        "_extract_star_geometry",
        "_extract_percent_coords_geometry",
        "_extract_block_spans",
        "_render_resource_blocks",
        "_parse_request_num_cpu",
        "_parse_request_memory",
        "_parse_output_print_settings",
    }
    for helper_name in helper_defs:
        if f"def {helper_name}" in text:
            violations.append(f"{_rel(parser)} still defines {helper_name}")

    missing_calls = sorted(name for name in required_calls if name not in text)
    if missing_calls:
        violations.append(f"{_rel(parser)} missing ORCA input extractor calls: {missing_calls}")

    assert violations == []


def test_orca_input_file_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    parser = ORCA_INPUT_ROOT / "parsers" / "ORCAInpFileParser.py"
    retired_extractor = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_extractors.py"
    retired_file_split_helpers = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_file_split.py"
    file_split_helpers = ORCA_INPUT_ROOT / "parsers" / "_orca_inp_file_extractors.py"
    required_calls = {
        "build_orca_input_file_lines",
        "ensure_orca_input_content",
        "split_orca_input_frames",
    }

    if not file_split_helpers.exists():
        violations.append(f"{_rel(file_split_helpers)} missing")
    if retired_extractor.exists():
        violations.append(f"{_rel(retired_extractor)} should stay retired")
    if retired_file_split_helpers.exists():
        violations.append(f"{_rel(retired_file_split_helpers)} should stay retired")

    text = _read(parser)
    helper_text = _read(file_split_helpers) if file_split_helpers.exists() else ""
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_extractors" in text:
        violations.append(f"{_rel(parser)} imports ORCA input compatibility extractor")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_file_split" in text:
        violations.append(
            f"{_rel(parser)} imports ORCA input file splitter from frame parser layer"
        )
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_patterns" in text:
        violations.append(f"{_rel(parser)} imports ORCA input patterns directly")
    if "molop.io.logic.orca.input.frame_parsers" in helper_text:
        violations.append(f"{_rel(file_split_helpers)} imports frame parser helpers")

    helper_defs = {
        "_ensure_orca_input_content",
        "_is_coords_open_line",
        "_is_new_job_delimiter",
    }
    for helper_name in helper_defs:
        if f"def {helper_name}" in text:
            violations.append(f"{_rel(parser)} still defines {helper_name}")

    missing_calls = sorted(name for name in required_calls if name not in text)
    if missing_calls:
        violations.append(f"{_rel(parser)} missing ORCA input extractor calls: {missing_calls}")

    assert violations == []


def test_gjf_input_parser_imports_only_frame_classes_from_frame_model() -> None:
    violations: list[str] = []
    allowed_frame_model_imports = {"GJFFileFrameDisk", "GJFFileFrameMemory"}

    for path in (GAUSSIAN_INPUT_ROOT / "frame_parsers").rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(_read(path), filename=_rel(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.module != "molop.io.logic.gaussian.input.frame_models.GJFFileFrame":
                continue
            imported_names = {alias.name for alias in node.names}
            disallowed = imported_names - allowed_frame_model_imports
            if disallowed:
                violations.append(
                    f"{_rel(path)} imports GJF format containers from frame model: "
                    f"{sorted(disallowed)}"
                )

    assert violations == []


def test_gjf_input_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    parser = GAUSSIAN_INPUT_ROOT / "frame_parsers" / "GJFFileFrameParser.py"
    extractor = GAUSSIAN_INPUT_ROOT / "frame_parsers" / "_gjf_extractors.py"
    required_calls = {
        "build_gjf_context_lines",
        "find_first_gjf_charge_multiplicity_index",
        "find_first_gjf_route_index",
        "find_gjf_molecule_block_end",
        "find_next_gjf_blank_index",
        "gjf_route_flags",
        "join_gjf_lines",
        "skip_gjf_blank_indices",
    }

    if not extractor.exists():
        violations.append(f"{_rel(extractor)} missing")

    text = _read(parser)
    if "molop.io.logic.gaussian.input.frame_parsers._gjf_patterns" in text:
        violations.append(f"{_rel(parser)} imports GJF input patterns directly")

    helper_defs = {
        "_normalize_free_separators",
        "_strip_comment_content",
        "_is_blank",
        "_join_lines",
        "_find_first_route_index",
        "_find_next_blank_index",
        "_skip_blank_indices",
        "_is_charge_multiplicity_line",
        "_find_first_charge_multiplicity_index",
        "_is_zmat_variable_label",
        "_is_zmat_variable_assignment",
        "_route_flags",
        "_find_molecule_block_end",
    }
    for helper_name in helper_defs:
        if f"def {helper_name}" in text:
            violations.append(f"{_rel(parser)} still defines {helper_name}")

    missing_calls = sorted(name for name in required_calls if name not in text)
    if missing_calls:
        violations.append(f"{_rel(parser)} missing GJF extractor calls: {missing_calls}")

    assert violations == []


def test_gjf_input_file_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    parser = GAUSSIAN_INPUT_ROOT / "parsers" / "GJFFileParser.py"
    frame_extractor = GAUSSIAN_INPUT_ROOT / "frame_parsers" / "_gjf_extractors.py"
    file_extractor = GAUSSIAN_INPUT_ROOT / "parsers" / "_gjf_file_extractors.py"
    required_calls = {
        "ensure_gjf_content",
        "expand_gjf_at_includes",
        "find_gjf_link1_matches",
        "validate_gjf_link1_boundaries",
        "split_gjf_link1_frames",
    }

    if not frame_extractor.exists():
        violations.append(f"{_rel(frame_extractor)} missing")
    if not file_extractor.exists():
        violations.append(f"{_rel(file_extractor)} missing")

    text = _read(parser)
    file_extractor_text = _read(file_extractor) if file_extractor.exists() else ""
    frame_extractor_text = _read(frame_extractor) if frame_extractor.exists() else ""
    if "molop.io.logic.gaussian.input.frame_parsers._gjf_extractors" in text:
        violations.append(f"{_rel(parser)} imports GJF file extractors from frame parser layer")
    if "molop.io.logic.gaussian.input.frame_parsers._gjf_patterns" in text:
        violations.append(f"{_rel(parser)} imports GJF input patterns directly")

    helper_defs = {
        "_ensure_gjf_content",
        "_expand_at_includes",
    }
    for helper_name in helper_defs:
        if f"def {helper_name}" in text:
            violations.append(f"{_rel(parser)} still defines {helper_name}")

    missing_calls = sorted(name for name in required_calls if name not in text)
    if missing_calls:
        violations.append(f"{_rel(parser)} missing GJF extractor calls: {missing_calls}")
    missing_file_helpers = sorted(
        name for name in required_calls if f"def {name}" not in file_extractor_text
    )
    if missing_file_helpers:
        violations.append(
            f"{_rel(file_extractor)} missing GJF file extractor helpers: {missing_file_helpers}"
        )
    misplaced_file_helpers = sorted(
        name for name in required_calls if f"def {name}" in frame_extractor_text
    )
    if misplaced_file_helpers:
        violations.append(
            f"{_rel(frame_extractor)} keeps file-level helpers: {misplaced_file_helpers}"
        )

    assert violations == []


def test_orca_input_token_helpers_stay_out_of_semantics_module() -> None:
    violations: list[str] = []
    token_helper = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_tokens.py"
    helper_files = [
        ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_blocks.py",
        ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_geometry.py",
        ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_resources.py",
    ]
    forbidden_private_imports = {
        "_parse_float",
        "_parse_int",
        "_strip_inline_comment",
        "_split_key_value",
    }

    if not token_helper.exists():
        violations.append(f"{_rel(token_helper)} missing")

    token_text = _read(token_helper) if token_helper.exists() else ""
    for helper_name in (
        "parse_orca_float",
        "parse_orca_int",
        "strip_orca_inline_comment",
        "split_orca_key_value",
    ):
        if f"def {helper_name}" not in token_text:
            violations.append(f"{_rel(token_helper)} missing {helper_name}")

    for path in helper_files:
        tree = ast.parse(_read(path), filename=_rel(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.module != "molop.io.logic.orca.input.frame_parsers._orca_inp_semantics":
                continue
            imported_names = {alias.name for alias in node.names}
            disallowed = imported_names & forbidden_private_imports
            if disallowed:
                violations.append(
                    f"{_rel(path)} imports token helpers from semantics: {sorted(disallowed)}"
                )

    assert violations == []


def test_orca_output_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    frame_parser = ORCA_LOG_ROOT / "frame_parsers" / "ORCALogFileFrameParser.py"
    file_parser = ORCA_LOG_ROOT / "parsers" / "ORCALogFileParser.py"
    frame_extractor = ORCA_LOG_ROOT / "frame_parsers" / "_orca_extractors.py"
    file_extractor = ORCA_LOG_ROOT / "parsers" / "_orca_log_file_extractors.py"
    shared_extractor = ORCA_LOG_ROOT / "parsers" / "_orca_log_shared.py"
    required_frame_extractors = {
        "extract_orca_coords",
        "extract_orca_electronic_states",
        "extract_orca_energies",
        "extract_orca_forces",
        "extract_orca_geometry_optimization_status",
        "extract_orca_polarizability",
        "extract_orca_populations",
        "extract_orca_running_time",
        "extract_orca_solvent",
        "extract_orca_status",
        "extract_orca_vibrations",
    }

    if not frame_extractor.exists():
        violations.append(f"{_rel(frame_extractor)} missing")
    if not file_extractor.exists():
        violations.append(f"{_rel(file_extractor)} missing")
    if not shared_extractor.exists():
        violations.append(f"{_rel(shared_extractor)} missing")

    frame_text = _read(frame_parser)
    if "molop.io.logic.orca.log.parsers._orca_log_patterns" in frame_text:
        violations.append(f"{_rel(frame_parser)} imports ORCA patterns directly")
    if "molop.io.base_models.DataClasses" in frame_text:
        violations.append(f"{_rel(frame_parser)} imports ORCA output data classes directly")
    missing_frame_extractors = sorted(
        name for name in required_frame_extractors if name not in frame_text
    )
    if missing_frame_extractors:
        violations.append(
            f"{_rel(frame_parser)} missing ORCA extractor calls: {missing_frame_extractors}"
        )

    file_text = _read(file_parser)
    shared_text = _read(shared_extractor) if shared_extractor.exists() else ""
    for helper_name in ("_parse_status", "_parse_running_time"):
        if f"def {helper_name}" in file_text:
            violations.append(f"{_rel(file_parser)} duplicates {helper_name}")
    for extractor_name in ("extract_orca_status", "extract_orca_running_time"):
        if extractor_name not in file_text:
            violations.append(f"{_rel(file_parser)} does not reuse {extractor_name}")
        if f"def {extractor_name}" not in shared_text:
            violations.append(f"{_rel(shared_extractor)} missing {extractor_name}")
    if "molop.io.logic.orca.log.parsers._orca_log_patterns" in file_text:
        violations.append(f"{_rel(file_parser)} imports ORCA patterns directly")
    if "molop.io.logic.orca.log.frame_parsers._orca_extractors" in file_text:
        violations.append(f"{_rel(file_parser)} imports ORCA metadata from frame parser layer")
    file_helper_defs = {
        "_ensure_orca_output",
        "_parse_version",
        "_extract_printed_input",
        "_parse_printed_input_metadata",
        "_split_orca_frames",
        "_first_frame_value",
        "_last_frame_value",
    }
    for helper_name in file_helper_defs:
        if f"def {helper_name}" in file_text:
            violations.append(f"{_rel(file_parser)} still defines {helper_name}")
    for extractor_name in (
        "ensure_orca_output_content",
        "extract_orca_output_version",
        "extract_orca_printed_input",
        "parse_orca_printed_input_metadata",
        "split_orca_output_frames",
    ):
        if extractor_name not in file_text:
            violations.append(f"{_rel(file_parser)} does not call {extractor_name}")

    assert violations == []


def test_g16_output_file_parser_uses_extractor_boundary() -> None:
    violations: list[str] = []
    file_parser = GAUSSIAN_LOG_ROOT / "parsers" / "G16LogFileParser.py"
    file_extractor = GAUSSIAN_LOG_ROOT / "parsers" / "_g16_log_file_extractors.py"
    required_extractors = {
        "ensure_g16_output_content",
        "extract_g16_charge_multiplicity",
        "extract_g16_keywords",
        "extract_g16_options",
        "extract_g16_running_time",
        "extract_g16_solvent",
        "extract_g16_standard_orientation_transformation_matrix",
        "extract_g16_temperature_and_pressure",
        "extract_g16_termination_status",
        "extract_g16_title",
        "extract_g16_version",
        "split_g16_section_frames",
        "split_g16_sections",
    }

    if not file_extractor.exists():
        violations.append(f"{_rel(file_extractor)} missing")

    file_text = _read(file_parser)
    extractor_text = _read(file_extractor) if file_extractor.exists() else ""
    if "molop.io.logic.gaussian.log.parsers._g16_log_patterns" in file_text:
        violations.append(f"{_rel(file_parser)} imports G16 patterns directly")
    if "molop.io.base_models.DataClasses" in file_text:
        violations.append(f"{_rel(file_parser)} imports G16 metadata data classes directly")
    if "import numpy" in file_text:
        violations.append(f"{_rel(file_parser)} imports numpy for raw metadata extraction")
    if "find_rigid_transform" in file_text:
        violations.append(f"{_rel(file_parser)} owns G16 coordinate transform extraction")

    file_helper_defs = {
        "_ensure_gaussian_output",
        "_first_frame_value",
        "_last_frame_value",
        "_parse_charge_multiplicity",
        "_parse_coordinates",
        "_parse_keywords",
        "_parse_options",
        "_parse_running_time",
        "_parse_solvent",
        "_parse_standard_orientation_transformation_matrix",
        "_parse_temperature_and_pressure",
        "_parse_termination_status",
        "_parse_title",
        "_parse_version",
    }
    for helper_name in file_helper_defs:
        if f"def {helper_name}" in file_text:
            violations.append(f"{_rel(file_parser)} still defines {helper_name}")

    missing_calls = sorted(name for name in required_extractors if name not in file_text)
    if missing_calls:
        violations.append(f"{_rel(file_parser)} missing G16 extractor calls: {missing_calls}")
    missing_helpers = sorted(
        name for name in required_extractors if f"def {name}" not in extractor_text
    )
    if missing_helpers:
        violations.append(f"{_rel(file_extractor)} missing G16 extractors: {missing_helpers}")

    assert violations == []


def test_orca_printed_input_metadata_uses_orca_projection_helper() -> None:
    violations: list[str] = []
    extractor = ORCA_LOG_ROOT / "parsers" / "_orca_log_file_extractors.py"
    metadata_helper = ORCA_INPUT_ROOT / "parsers" / "_orca_inp_metadata.py"
    frame_parser = ORCA_INPUT_ROOT / "frame_parsers" / "ORCAInpFileFrameParser.py"
    result_helper = ORCA_INPUT_ROOT / "frame_parsers" / "_orca_inp_result.py"
    orca_model = ORCA_LOGIC_ROOT / "common.py"
    extractor_text = _read(extractor)
    metadata_helper_text = _read(metadata_helper) if metadata_helper.exists() else ""
    frame_parser_text = _read(frame_parser)
    base_text = _read(orca_model)

    if not metadata_helper.exists():
        violations.append(f"{_rel(metadata_helper)} missing")
    if result_helper.exists():
        violations.append(f"{_rel(result_helper)} should stay retired")
    if "def project_orca_printed_input_metadata" not in base_text:
        violations.append(f"{_rel(orca_model)} missing ORCA printed-input metadata helper")
    if "project_orca_printed_input_metadata" not in metadata_helper_text:
        violations.append(f"{_rel(metadata_helper)} does not use ORCA printed-input projection")
    if "def parse_orca_input_frame_result" not in frame_parser_text:
        violations.append(f"{_rel(frame_parser)} missing ORCA input result helper")
    if "parse_orca_input_frame_result" not in metadata_helper_text:
        violations.append(f"{_rel(metadata_helper)} does not use ORCA input result helper")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_result" in metadata_helper_text:
        violations.append(f"{_rel(metadata_helper)} imports retired ORCA input result facade")
    if "parse_orca_input_metadata" not in extractor_text:
        violations.append(f"{_rel(extractor)} does not call ORCA input metadata facade")
    if "molop.io.logic.orca.input.frame_parsers._orca_inp_result" in extractor_text:
        violations.append(f"{_rel(extractor)} imports ORCA input frame result helper directly")
    if "ORCAInpFileFrameParser" in extractor_text:
        violations.append(f"{_rel(extractor)} imports concrete ORCA input frame parser module")

    banned_fragments = {
        "ImplicitSolvation": "constructs ORCA solvent metadata directly",
        "ORCAInpFileFrameParserMemory": "instantiates ORCA input frame parser directly",
        "metadata[field]": "owns ORCA printed-input field projection",
        "excited_state_requests": "owns ORCA printed-input field list",
        "multireference_requests": "owns ORCA printed-input field list",
    }
    for fragment, message in banned_fragments.items():
        for path, text in ((extractor, extractor_text), (metadata_helper, metadata_helper_text)):
            if fragment in text:
                violations.append(f"{_rel(path)} {message}")

    assert violations == []


def test_coords_file_parsers_use_extractor_boundary() -> None:
    violations: list[str] = []
    parser_extractors = {
        COORDS_LOGIC_ROOT / "parsers" / "XYZFileParser.py": "split_xyz_frames",
        COORDS_LOGIC_ROOT / "parsers" / "SDFFileParser.py": "split_sdf_frames",
        COORDS_LOGIC_ROOT / "parsers" / "SMIFileParser.py": "split_smi_frames",
    }
    frame_extractor = COORDS_LOGIC_ROOT / "frame_parsers" / "_coords_extractors.py"
    file_extractor = COORDS_LOGIC_ROOT / "parsers" / "_coords_file_extractors.py"

    if not frame_extractor.exists():
        violations.append(f"{_rel(frame_extractor)} missing")
    if not file_extractor.exists():
        violations.append(f"{_rel(file_extractor)} missing")

    frame_extractor_text = _read(frame_extractor) if frame_extractor.exists() else ""
    file_extractor_text = _read(file_extractor) if file_extractor.exists() else ""

    for parser, extractor_name in parser_extractors.items():
        text = _read(parser)
        if extractor_name not in text:
            violations.append(f"{_rel(parser)} does not call {extractor_name}")
        if "molop.io.logic.coords.frame_parsers._coords_extractors" in text:
            violations.append(
                f"{_rel(parser)} imports coords file splitters from frame parser layer"
            )
        if "molop.io.codec_exceptions import FormatMismatchError" in text:
            violations.append(f"{_rel(parser)} still owns split mismatch errors")
        if "from rdkit import Chem" in text:
            violations.append(f"{_rel(parser)} still owns RDKit split implementation")
        if f"def {extractor_name}" not in file_extractor_text:
            violations.append(f"{_rel(file_extractor)} missing {extractor_name}")
        if f"def {extractor_name}" in frame_extractor_text:
            violations.append(f"{_rel(frame_extractor)} keeps file-level {extractor_name}")

    assert violations == []


def test_coords_frame_parsers_use_extractor_boundary() -> None:
    violations: list[str] = []
    parser_extractors = {
        COORDS_LOGIC_ROOT / "frame_parsers" / "XYZFileFrameParser.py": "extract_xyz_frame_payload",
        COORDS_LOGIC_ROOT / "frame_parsers" / "SDFFileFrameParser.py": "extract_sdf_frame_payload",
        COORDS_LOGIC_ROOT / "frame_parsers" / "SMIFileFrameParser.py": "extract_smi_frame_payload",
    }
    extractor = COORDS_LOGIC_ROOT / "frame_parsers" / "_coords_extractors.py"

    if not extractor.exists():
        violations.append(f"{_rel(extractor)} missing")

    banned_fragments = {
        "from rdkit import Chem": "still owns RDKit frame implementation",
        "import numpy as np": "still owns NumPy coordinate extraction",
        "molop.io.codec_exceptions import FormatMismatchError": "still owns frame mismatch errors",
        "molop.structure.StructureTransformation": "still owns structure projection helpers",
        "molop.unit import atom_ureg": "still owns unit-bearing payload construction",
        "molop.io.logic.coords.frame_parsers._xyz_patterns": "imports XYZ patterns directly",
    }
    for parser, extractor_name in parser_extractors.items():
        text = _read(parser)
        if extractor_name not in text:
            violations.append(f"{_rel(parser)} does not call {extractor_name}")
        for fragment, message in banned_fragments.items():
            if fragment in text:
                violations.append(f"{_rel(parser)} {message}")

    assert violations == []


def test_coords_frame_models_use_format_render_helpers() -> None:
    violations: list[str] = []
    frame_helpers = {
        COORDS_LOGIC_ROOT / "frame_models" / "XYZFileFrame.py": "render_xyz_frame",
        COORDS_LOGIC_ROOT / "frame_models" / "SDFFileFrame.py": "render_sdf_frame",
        COORDS_LOGIC_ROOT / "frame_models" / "SMIFileFrame.py": "render_smi_frame",
    }
    helper_module = COORDS_LOGIC_ROOT / "frame_models" / "_coords_renderers.py"
    helper_text = _read(helper_module) if helper_module.exists() else ""
    base_model = BASE_MODELS_ROOT / "Molecule.py"
    base_text = _read(base_model)

    if not helper_module.exists():
        violations.append(f"{_rel(helper_module)} missing")

    for helper_name in frame_helpers.values():
        if f"def {helper_name}" not in helper_text:
            violations.append(f"{_rel(helper_module)} missing {helper_name}")
        if f"def {helper_name}" in base_text:
            violations.append(f"{_rel(base_model)} defines format-specific {helper_name}")

    banned_fragments = {
        "from rdkit import Chem": "still owns RDKit rendering",
        "Chem.MolToMolBlock": "still renders SDF directly",
        '.write("sdf")': "still renders OpenBabel SDF directly",
        "to_canonical_SMILES()": "still renders SMI directly",
    }
    for frame_model, helper_name in frame_helpers.items():
        text = _read(frame_model)
        if helper_name not in text:
            violations.append(f"{_rel(frame_model)} does not call {helper_name}")
        if "molop.io.base_models.Molecule import render_" in text:
            violations.append(f"{_rel(frame_model)} imports format renderer from base_models")
        for fragment, message in banned_fragments.items():
            if fragment in text:
                violations.append(f"{_rel(frame_model)} {message}")

    assert violations == []


def test_simple_file_models_use_base_frame_rendering() -> None:
    violations: list[str] = []
    simple_file_models = {
        COORDS_LOGIC_ROOT / "models" / "XYZFile.py": 'file_frame_separator = "\\n"',
        COORDS_LOGIC_ROOT / "models" / "SDFFile.py": 'file_frame_separator = "$$$$\\n"',
        COORDS_LOGIC_ROOT / "models" / "SMIFile.py": 'file_frame_separator = "\\n"',
        GAUSSIAN_INPUT_ROOT / "models" / "GJFFile.py": 'file_frame_separator = "\\n--Link1--\\n"',
    }
    base_mixin = SRC_ROOT / "molop" / "io" / "base_models" / "Mixins.py"
    base_text = _read(base_mixin)

    if "file_frame_separator" not in base_text:
        violations.append(f"{_rel(base_mixin)} missing file frame separator")
    if "def _render_frames_in_one_file" not in base_text:
        violations.append(f"{_rel(base_mixin)} missing shared single-file renderer")
    if "def _render_frames" not in base_text:
        violations.append(f"{_rel(base_mixin)} missing shared frame renderer")

    for path, separator_source in simple_file_models.items():
        text = _read(path)
        tree = ast.parse(text, filename=_rel(path))
        if "file_frame_separator" not in text:
            violations.append(f"{_rel(path)} does not declare file_frame_separator")
        if separator_source not in text:
            violations.append(f"{_rel(path)} missing expected frame separator")
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name in {
                "_render_frames_in_one_file",
                "_render_frames",
            }:
                violations.append(f"{_rel(path)} reimplements {node.name}")
        banned_fragments = {
            "frame._render": "still renders frames directly",
            "frame.frame_id in frame_ids": "still filters frames directly",
            "_HasRenderableFrames": "still depends on renderable frame protocol directly",
        }
        for fragment, message in banned_fragments.items():
            if fragment in text:
                violations.append(f"{_rel(path)} {message}")

    assert violations == []
