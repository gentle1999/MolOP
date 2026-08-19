import importlib
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


def test_batch_graph_transform_prewarms_once_before_file_workers(monkeypatch):
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    original_iterator = molecule_module.iter_xyz_to_rdmol_batch
    calls = 0

    def counting_iterator(requests, **kwargs):
        nonlocal calls
        calls += 1
        yield from original_iterator(requests, **kwargs)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", counting_iterator)

    rendered = batch.format_transform("smi", n_jobs=1)

    assert calls == 1
    assert isinstance(rendered[str(fixture_path)], str)


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
        write_to_disk=True,
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

    batch.format_transform("xyz", output_dir=str(out_dir), write_to_disk=True, n_jobs=1)

    expected = out_dir / "dsgdb9nsd_004015-7.xyz"
    leak = cwd_dir / "dsgdb9nsd_004015-7.xyz"
    assert expected.exists()
    assert not leak.exists()


def test_pure_filename_and_batch_output_preserve_hash_before_last_suffix(tmp_path):
    fixture_path = tmp_path / "source.abc123.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    batch = AutoParser(fixture_path)

    assert batch[0].pure_filename == "source.abc123"

    batch.format_transform(
        "gjf",
        output_dir=str(out_dir),
        write_to_disk=True,
        graph_policy="prefer",
        n_jobs=1,
    )

    assert (out_dir / "source.abc123.gjf").is_file()
    assert not (out_dir / "source.gjf").exists()


def test_format_transform_explicit_multi_suffix_output_preserves_hash(tmp_path):
    fixture_path = tmp_path / "source.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    file_model = AutoParser(fixture_path)[0]

    file_model.format_transform(
        "gjf",
        file_path=tmp_path / "target.abc123.placeholder",
        write_to_disk=True,
        graph_policy="prefer",
    )

    assert (tmp_path / "target.abc123.gjf").is_file()
    assert not (tmp_path / "target.gjf").exists()


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
        write_to_disk=True,
        frame="all",
        n_jobs=1,
    )

    assert len(list(out_dir.glob("*.xyz"))) >= 2
    assert len(list(cwd_dir.glob("*.xyz"))) == 0


def test_format_transform_multi_output_preserves_hash_before_last_suffix(tmp_path):
    source_fixture = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    fixture_path = tmp_path / "source.abc123.log"
    fixture_path.write_bytes(source_fixture.read_bytes())
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    batch = AutoParser(fixture_path)

    batch.format_transform(
        "xyz",
        output_dir=str(out_dir),
        embed_in_one_file=False,
        write_to_disk=True,
        frame=[0, 1],
        n_jobs=1,
    )

    assert (out_dir / "source.abc123000.xyz").is_file()
    assert (out_dir / "source.abc123001.xyz").is_file()
    assert not (out_dir / "source000.xyz").exists()


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
        write_to_disk=True,
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
        write_to_disk=True,
        frame=[0, 1],
        chk=True,
        n_jobs=1,
    )

    gjf0 = out_dir / "1000.gjf"
    gjf1 = out_dir / "1001.gjf"
    assert gjf0.exists()
    assert gjf1.exists()

    assert "%chk=1000.chk" in gjf0.read_text()
    assert "%chk=1001.chk" in gjf1.read_text()


def test_format_transform_write_to_disk_defaults_to_source_directory(tmp_path):
    fixture_path = tmp_path / "source.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    batch = AutoParser(str(fixture_path))
    file_model = batch[0]

    rendered = file_model.format_transform("gjf", write_to_disk=True, graph_policy="prefer")

    expected = tmp_path / "source.gjf"
    assert expected.exists()
    assert expected.read_text(encoding="utf-8") == rendered


def test_format_transform_file_path_is_ignored_when_not_writing(tmp_path):
    out_file = tmp_path / "ignored_name.gjf"
    fixture_path = Path(__file__).resolve().parent / "test_files" / "g16log" / "1.log"
    batch = AutoParser(str(fixture_path))
    file_model = batch[0]

    rendered = file_model.format_transform("gjf", file_path=str(out_file), chk=True)

    assert not out_file.exists()
    assert "%chk=ignored_name.chk" not in rendered


def test_frame_format_transform_write_to_disk_defaults_to_source_directory(tmp_path):
    fixture_path = tmp_path / "frame_source.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    batch = AutoParser(str(fixture_path))
    frame = batch[0][0]

    rendered = frame.format_transform("gjf", write_to_disk=True, graph_policy="prefer")

    expected = tmp_path / "frame_source.gjf"
    assert expected.exists()
    assert expected.read_text(encoding="utf-8") == rendered


def test_batch_format_transform_write_to_disk_defaults_to_source_directory(tmp_path):
    fixture_path = tmp_path / "batch_source.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    batch = AutoParser(str(fixture_path))

    result = batch.format_transform("gjf", write_to_disk=True, graph_policy="prefer", n_jobs=1)

    expected = tmp_path / "batch_source.gjf"
    assert expected.exists()
    assert result[str(fixture_path)] == expected.read_text(encoding="utf-8")


def test_batch_gjf_connectivity_uses_parent_prewarmed_topology(tmp_path):
    fixture_path = tmp_path / "batch_connectivity.xyz"
    fixture_path.write_text(
        "2\ncomment\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    batch = AutoParser(str(fixture_path))

    result = batch.format_transform(
        "gjf",
        add_gjf_connectivity=True,
        n_jobs=2,
    )

    assert "1 2 1.0" in result[str(fixture_path)]


@pytest.mark.parametrize("format_id", ["smi", "sdf", "cml"])
def test_batch_graph_transform_handles_disconnected_single_atom_in_workers(
    tmp_path,
    format_id,
):
    fixture_path = tmp_path / "single_atom.xyz"
    fixture_path.write_text(
        "1\ncomment\nHe 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    batch = AutoParser(str(fixture_path))

    result = batch.format_transform(format_id, n_jobs=2)

    assert result[str(fixture_path)]


def test_batch_format_transform_output_dir_is_ignored_when_not_writing(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))

    result = batch.format_transform("xyz", output_dir=str(out_dir), n_jobs=1)

    assert result
    assert not list(out_dir.iterdir())
