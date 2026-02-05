from pathlib import Path

from rdkit import Chem

from molop.io import AutoParser  # type: ignore[reportMissingImports]


def test_autoparser_sdf_tmpfile(tmp_path: Path) -> None:
    path = tmp_path / "example.sdf"
    mol = Chem.MolFromSmiles("CCO")
    mol_block = Chem.MolToMolBlock(mol)
    path.write_text(f"{mol_block}\n$$$$\n", encoding="utf-8")
    batch = AutoParser(str(path))
    assert len(batch) > 0
    assert len(batch[0]) > 0
    last_frame = batch[0][-1]
    atoms = getattr(last_frame, "atoms", None)
    coords = getattr(last_frame, "coords", None)
    assert atoms is not None and len(atoms) > 0
    assert coords is not None
