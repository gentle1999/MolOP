from pathlib import Path

import pytest

from molop.io import AutoParser


_G16GJF_FIXTURE_DIR = Path(__file__).resolve().parent / "test_files" / "g16gjf"
_G16GJF_FIXTURES = sorted(_G16GJF_FIXTURE_DIR.glob("*.gjf"))
_VALID_G16GJF_FIXTURES = [fixture for fixture in _G16GJF_FIXTURES if "_invalid" not in fixture.stem]


def test_g16gjf_fixture_inventory_is_not_empty() -> None:
    assert _G16GJF_FIXTURES


def test_g16gjf_valid_fixture_inventory_is_not_empty() -> None:
    assert _VALID_G16GJF_FIXTURES


@pytest.mark.parametrize("fixture_path", _VALID_G16GJF_FIXTURES, ids=lambda p: p.name)
def test_g16gjf_all_fixtures_are_parseable_and_renderable(fixture_path: Path) -> None:
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0, f"No file model parsed for fixture: {fixture_path.name}"

    file_model = batch[0]
    assert len(file_model) > 0, f"No frames parsed for fixture: {fixture_path.name}"

    rendered = file_model.format_transform("gjf")
    assert isinstance(rendered, str), (
        f"Expected single GJF block string for fixture: {fixture_path.name}"
    )
    assert rendered.strip(), f"Empty rendered GJF block for fixture: {fixture_path.name}"
