from pathlib import Path

import pytest

from molop.io import AutoParser
from molop.io.codec_exceptions import UnsupportedFormatError


def test_format_transform_unsupported_format_raises_error(tmp_path):
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    file_model = batch[0]

    with pytest.raises(UnsupportedFormatError):
        file_model.format_transform("nonexistent_format")


def test_format_transform_file_path_writes_under_requested_dir(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    cwd_dir = tmp_path / "cwd"
    out_dir.mkdir()
    cwd_dir.mkdir()

    monkeypatch.chdir(cwd_dir)

    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    file_model = batch[0]

    file_model.format_transform(
        "xyz",
        file_path=str(out_dir / "out.xyz"),
        graph_policy="prefer",
    )

    assert (out_dir / "out.xyz").exists()
    assert not (cwd_dir / "out.xyz").exists()


def test_format_transform_batch_output_dir_writes_under_requested_dir(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    cwd_dir = tmp_path / "cwd"
    out_dir.mkdir()
    cwd_dir.mkdir()

    monkeypatch.chdir(cwd_dir)

    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))

    batch.format_transform("xyz", output_dir=str(out_dir), n_jobs=1)

    expected = out_dir / "dsgdb9nsd_004015-7.xyz"
    leak = cwd_dir / "dsgdb9nsd_004015-7.xyz"
    assert expected.exists()
    assert not leak.exists()


def test_format_transform_multi_output_stays_in_output_dir(tmp_path, monkeypatch):
    out_dir = tmp_path / "out"
    cwd_dir = tmp_path / "cwd"
    out_dir.mkdir()
    cwd_dir.mkdir()

    monkeypatch.chdir(cwd_dir)

    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))

    batch.format_transform(
        "xyz",
        output_dir=str(out_dir),
        embed_in_one_file=False,
        frameID="all",
        n_jobs=1,
    )

    assert len(list(out_dir.glob("*.xyz"))) >= 2
    assert len(list(cwd_dir.glob("*.xyz"))) == 0


def test_format_transform_gjf_chk_propagation_single_file(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))
    file_model = batch[0]

    out_file = out_dir / "custom.gjf"
    file_model.format_transform(
        "gjf",
        file_path=str(out_file),
        chk=True,
    )

    assert out_file.exists()
    content = out_file.read_text()
    assert "%chk=custom.chk" in content


def test_format_transform_gjf_chk_propagation_multi_file(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))

    batch.format_transform(
        "gjf",
        output_dir=str(out_dir),
        embed_in_one_file=False,
        frameID=[0, 1],
        chk=True,
        n_jobs=1,
    )

    gjf0 = out_dir / "1000.gjf"
    gjf1 = out_dir / "1001.gjf"
    assert gjf0.exists()
    assert gjf1.exists()

    assert "%chk=1000.chk" in gjf0.read_text()
    assert "%chk=1001.chk" in gjf1.read_text()
