from pathlib import Path

from molop.io import AutoParser


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
