import tempfile
from pathlib import Path

from molop.io import AutoParser  # type: ignore[reportMissingImports]


def _ensure_str(rendered: str | list[str]) -> str:
    assert isinstance(rendered, str)
    return rendered


def _assert_non_empty_str(block: str) -> None:
    assert isinstance(block, str)
    assert block.strip()


def _assert_xyz_block(xyz_block: str) -> None:
    _assert_non_empty_str(xyz_block)
    first_line = xyz_block.strip().splitlines()[0]
    assert int(first_line) > 0


def _assert_smi_block(smi_block: str) -> None:
    _assert_non_empty_str(smi_block)
    assert "\n" not in smi_block.strip()


def test_registry_convert_xyz_smoke() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    sdf_block = _ensure_str(file_model.format_transform("sdf", graph_policy="prefer"))
    _assert_non_empty_str(sdf_block)
    assert "M  END" in sdf_block

    smi_block = _ensure_str(file_model.format_transform("smi", graph_policy="prefer"))
    _assert_smi_block(smi_block)

    gjf_block = _ensure_str(
        file_model.format_transform(
            "gjf",
            graph_policy="prefer",
            options="%mem=1GB",
            route="#p hf/3-21g",
            title_card="MolOP",
            suffix="end",
        )
    )
    _assert_non_empty_str(gjf_block)

    cml_block = _ensure_str(file_model.format_transform("cml", graph_policy="prefer"))
    _assert_non_empty_str(cml_block)


def test_registry_convert_g16log_smoke() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    xyz_block = _ensure_str(file_model.format_transform("xyz", graph_policy="prefer"))
    _assert_xyz_block(xyz_block)

    smi_block = _ensure_str(file_model.format_transform("smi", graph_policy="prefer"))
    _assert_smi_block(smi_block)


def test_unknown_extension_fallback_smoke() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    xyz_content = fixture_path.read_text()

    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_path = Path(tmp_dir) / "fixture.foo"
        temp_path.write_text(xyz_content)

        batch = AutoParser(str(temp_path))
        assert len(batch) > 0
        file_model = batch[0]

        xyz_block = _ensure_str(file_model.format_transform("xyz", graph_policy="prefer"))
        _assert_xyz_block(xyz_block)
