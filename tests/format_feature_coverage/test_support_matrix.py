from __future__ import annotations

import ast
from pathlib import Path

from molop.io.codec_registry import Registry

from .support_matrix import (
    BUILTIN_READER_FORMATS,
    BUILTIN_WRITER_FORMATS,
    FORMAT_FEATURE_COVERAGE,
    FORMAT_GROUPS,
    FORMAT_SUPPORT,
    SPECIAL_CODEC_IDS,
    SupportLevel,
)


VALID_SUPPORT_LEVELS: set[SupportLevel] = {
    "supported",
    "partial",
    "fixture-covered",
    "intentionally-unsupported",
}
INTERNAL_IMPLEMENTATION_TERMS = ("component tree", "component-tree", "payload", "segmentation")


def _test_functions_by_file() -> dict[str, set[str]]:
    tests_root = Path(__file__).resolve().parents[1]
    result: dict[str, set[str]] = {}
    for path in tests_root.glob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        result[path.name] = {
            node.name
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name.startswith("test_")
        }
    return result


def test_format_feature_coverage_matches_builtin_codec_registry() -> None:
    registry = Registry()
    registry.ensure_default_codecs_registered()

    assert set(registry._readers_by_format) == set(BUILTIN_READER_FORMATS)
    assert set(registry._writers_by_format) == set(BUILTIN_WRITER_FORMATS)
    assert set(FORMAT_FEATURE_COVERAGE) == set(
        BUILTIN_READER_FORMATS | BUILTIN_WRITER_FORMATS | SPECIAL_CODEC_IDS
    )


def test_format_support_is_grouped_by_format_type_without_duplicates() -> None:
    grouped_format_ids = [format_id for formats in FORMAT_GROUPS.values() for format_id in formats]

    assert set(grouped_format_ids) == set(FORMAT_SUPPORT)
    assert len(grouped_format_ids) == len(set(grouped_format_ids))
    for group, format_ids in FORMAT_GROUPS.items():
        assert format_ids, f"{group} must include at least one format."
        for format_id in format_ids:
            assert FORMAT_SUPPORT[format_id].group == group


def test_every_format_feature_declares_scope_limits_level_and_test_support() -> None:
    functions_by_file = _test_functions_by_file()

    for format_id, report in FORMAT_SUPPORT.items():
        assert report.features, f"{format_id} must list at least one covered feature."
        feature_areas = [feature.area for feature in report.features]
        assert len(feature_areas) == len(set(feature_areas)), (
            f"{format_id} has duplicate feature areas."
        )
        for feature in report.features:
            assert feature.support in VALID_SUPPORT_LEVELS, (
                f"{format_id}: invalid support level {feature.support!r}."
            )
            assert feature.scope.strip(), f"{format_id}: {feature.area!r} has no scope."
            assert feature.limitations.strip(), f"{format_id}: {feature.area!r} has no limits."
            assert feature.tests, f"{format_id}: {feature.area!r} has no test support."
            for node_id in feature.tests:
                path, sep, test_name = node_id.partition("::")
                assert sep, f"{format_id}: invalid test node id {node_id!r}."
                assert path.startswith("tests/"), f"{format_id}: test node must live under tests/."
                test_file = Path(path).name
                assert test_file in functions_by_file, f"{format_id}: missing test file {path!r}."
                assert test_name in functions_by_file[test_file], (
                    f"{format_id}: missing test function {node_id!r}."
                )


def test_qm_output_feature_rows_describe_user_visible_properties() -> None:
    for format_id in ("g16log", "orcaout", "fakeg"):
        for feature in FORMAT_SUPPORT[format_id].features:
            text = f"{feature.area} {feature.scope} {feature.limitations}".lower()
            for term in INTERNAL_IMPLEMENTATION_TERMS:
                assert term not in text, (
                    f"{format_id}: {feature.area!r} describes parser internals instead of "
                    f"user-visible log properties: {term!r}."
                )


def test_format_support_docs_have_overview_and_per_format_pages() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    for locale in ("en", "zh"):
        overview = repo_root / "docs" / locale / "reference" / "format_support.md"
        assert overview.exists(), f"Missing {locale} format support overview."
        overview_text = overview.read_text(encoding="utf-8")
        assert "<!-- format-support-overview -->" in overview_text

        for format_id, report in FORMAT_SUPPORT.items():
            page = repo_root / "docs" / locale / "reference" / "formats" / f"{format_id}.md"
            assert page.exists(), f"Missing {locale} format page for {format_id}."
            text = page.read_text(encoding="utf-8")
            assert f"<!-- format-support:{format_id} -->" in text
            for feature in report.features:
                assert f"<!-- feature-area:{feature.area} -->" in text, (
                    f"{locale} {format_id} page does not document {feature.area!r}."
                )
