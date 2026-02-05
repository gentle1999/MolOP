from pathlib import Path

from molop.io import AutoParser  # type: ignore[reportMissingImports]


def _assert_last_frame_has_data(file_model):
    assert len(file_model) > 0
    last_frame = file_model[-1]
    atoms = getattr(last_frame, "atoms", None)
    coords = getattr(last_frame, "coords", None)
    energies = getattr(last_frame, "energies", None)
    has_atoms = atoms is not None and len(atoms) > 0
    has_coords = coords is not None
    has_energies = energies is not None
    assert has_atoms or has_coords or has_energies


def test_autoparser_xyz_smoke():
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_130336-5.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    _assert_last_frame_has_data(batch[0])


def test_autoparser_g16log_smoke():
    fixture_path = (
        Path(__file__).resolve().parent / "test_files" / "g16log" / "dsgdb9nsd_131941-4+.log"
    )
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    _assert_last_frame_has_data(batch[0])
