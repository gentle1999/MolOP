from __future__ import annotations

from pathlib import Path

import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]


ORCA_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "orca"


@pytest.mark.parametrize("fixture_path", sorted(ORCA_FIXTURE_DIR.glob("*.inp")))
def test_orcainp_fixture_parse_and_raw_roundtrip(fixture_path: Path, tmp_path: Path) -> None:
    original = fixture_path.read_text(encoding="utf-8")

    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]
    assert len(file_model) > 0

    rendered = file_model.format_transform("orcainp", frameID="all")
    assert isinstance(rendered, str)
    assert rendered.strip()

    original_first = next((line for line in original.splitlines() if line.strip()), "")
    rendered_first = next((line for line in rendered.splitlines() if line.strip()), "")
    assert rendered_first == original_first
    for marker in ("%pal", "%maxcore", "%geom", "%scf", "%output"):
        if marker in original:
            assert marker in rendered

    assert "* xyz" in rendered
    assert rendered.strip().endswith("*")

    # Idempotence: parse rendered content and ensure it re-renders identically.
    roundtrip_path = tmp_path / "roundtrip.inp"
    roundtrip_path.write_text(rendered, encoding="utf-8")
    rt_file = AutoParser(str(roundtrip_path))[0]
    rendered2 = rt_file.format_transform("orcainp", frameID="all")
    assert isinstance(rendered2, str)
    assert rendered2.strip() == rendered.strip()
