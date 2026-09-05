from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "resolve_docs_version.py"
WORKFLOW_DIR = Path(__file__).parents[1] / ".github" / "workflows"


def _resolve(source_ref: str, update_latest: bool = False) -> dict[str, str]:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), source_ref, str(update_latest).lower()],
        capture_output=True,
        check=False,
        text=True,
    )
    if result.returncode != 0:
        raise ValueError(result.stderr)
    return dict(line.split("=", maxsplit=1) for line in result.stdout.splitlines())


@pytest.mark.parametrize(
    ("source_ref", "expected_version"),
    [
        ("v1.2.3", "1.2.3"),
        ("v1.2.3rc1", "1.2.3rc1"),
        ("v1.2.3.post1", "1.2.3.post1"),
        ("v1.2.3.dev4", "1.2.3.dev4"),
        ("v0.2.1dev", "0.2.1dev"),
    ],
)
def test_resolve_docs_version_accepts_release_tags(source_ref: str, expected_version: str) -> None:
    resolved = _resolve(source_ref, update_latest=True)

    assert resolved == {
        "snapshot_kind": "tag",
        "version": expected_version,
        "source_version": expected_version,
        "channel": "release",
        "update_latest": "true",
    }


def test_resolve_docs_version_handles_main() -> None:
    resolved = _resolve("main", update_latest=True)

    assert resolved == {
        "snapshot_kind": "development",
        "version": "main",
        "source_version": "",
        "channel": "development",
        "update_latest": "false",
    }


@pytest.mark.parametrize("source_ref", ["v1.2", "1.2.3", "v1.2.3.post", "feature/test"])
def test_resolve_docs_version_rejects_unsupported_refs(source_ref: str) -> None:
    with pytest.raises(ValueError, match="source_ref must be main"):
        _resolve(source_ref)


def test_github_actions_use_version_tags() -> None:
    uses_pattern = re.compile(r"^\s+uses:\s+[^@]+@(?:v\d+(?:\.\d+)*|release/v\d+)$")
    for workflow in WORKFLOW_DIR.iterdir():
        if workflow.suffix not in {".yml", ".yaml"}:
            continue
        for line in workflow.read_text(encoding="utf-8").splitlines():
            if " uses: " in line:
                assert uses_pattern.match(line), (
                    f"Action must use a version tag in {workflow}: {line}"
                )


def test_release_deployment_requires_a_manually_created_github_release() -> None:
    workflow_text = (WORKFLOW_DIR / "ci.yaml").read_text(encoding="utf-8")

    assert "gh release view" in workflow_text
    assert "gh release upload" in workflow_text
    assert "Create and publish it manually" in workflow_text
    assert "gh release create" not in workflow_text
