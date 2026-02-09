from pathlib import Path
from typing import Any, cast

from molop.io import AutoParser  # type: ignore[reportMissingImports]


def _fixture_xyz_path() -> Path:
    return Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"


def test_detected_format_id_persisted_on_parsed_file() -> None:
    batch = AutoParser(str(_fixture_xyz_path()))
    assert len(batch) > 0
    file_model = batch[0]

    # New behavior: parsing should persist reader format_id onto disk objects.
    assert getattr(file_model, "detected_format_id", None) == "xyz"


def test_filter_by_codec_id_uses_detected_format_id(tmp_path: Path) -> None:
    xyz_content = _fixture_xyz_path().read_text()
    unknown_suffix_path = tmp_path / "fixture.foo"
    unknown_suffix_path.write_text(xyz_content)

    batch = AutoParser(str(unknown_suffix_path))
    assert len(batch) > 0
    file_model = batch[0]

    # For unknown suffix, parsing should still detect XYZ and persist it.
    assert getattr(file_model, "detected_format_id", None) == "xyz"

    # New API: batch filter based on detected codec id.
    assert hasattr(batch, "filter_by_codec_id")
    batch_any = cast(Any, batch)
    kept = batch_any.filter_by_codec_id("xyz")
    assert len(kept) == 1
