from molop.io import AutoParser  # type: ignore[reportMissingImports]


def test_autoparser_smi_tmpfile(tmp_path):
    path = tmp_path / "example.smi"
    path.write_text("CCO\n", encoding="utf-8")
    batch = AutoParser(str(path))
    assert len(batch) > 0
    assert len(batch[0]) > 0
    last_frame = batch[0][-1]
    atoms = getattr(last_frame, "atoms", None)
    coords = getattr(last_frame, "coords", None)
    assert atoms is not None and len(atoms) > 0
    assert coords is not None
