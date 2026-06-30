from __future__ import annotations

from pathlib import Path
from typing import cast

import pytest
from rdkit import Chem

from molop.io import AutoParser


def test_xyz_frame_writer_uses_explicit_comment_override(tmp_path: Path) -> None:
    path = tmp_path / "h2.xyz"
    path.write_text(
        "2\ncharge 1 multiplicity 2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n",
        encoding="utf-8",
    )
    frame = AutoParser(str(path), n_jobs=1)[0][0]

    rendered = cast(str, frame.format_transform("xyz", comment="custom comment"))

    assert rendered.splitlines()[0] == "2"
    assert rendered.splitlines()[1] == "comment custom comment"


def test_sdf_frame_writer_supports_rdkit_openbabel_engines_and_rejects_unknown_engine(
    tmp_path: Path,
) -> None:
    path = tmp_path / "ethanol.sdf"
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    path.write_text(f"{Chem.MolToMolBlock(mol)}\n$$$$\n", encoding="utf-8")
    frame = AutoParser(str(path), n_jobs=1)[0][-1]

    rdkit_block = cast(str, frame.format_transform("sdf", engine="rdkit"))
    openbabel_block = cast(str, frame.format_transform("sdf", engine="openbabel"))

    assert "M  END" in rdkit_block
    assert "M  END" in openbabel_block
    assert "OpenBabel" in openbabel_block
    with pytest.raises(ValueError, match="Unsupported engine"):
        frame.format_transform("sdf", engine="not-a-real-engine")


def test_smi_reader_uses_first_token_and_writer_canonicalizes(tmp_path: Path) -> None:
    path = tmp_path / "ethanol.smi"
    path.write_text("OCC ethanol name ignored\n", encoding="utf-8")
    file_model = AutoParser(str(path), n_jobs=1)[0]

    frame = file_model[-1]
    rendered = cast(str, file_model.format_transform("smi"))

    assert frame.to_canonical_SMILES() == "CCO"
    assert rendered == "CCO"


def test_openbabel_unknown_extension_fallback_keeps_first_molecule_coordinates(
    tmp_path: Path,
) -> None:
    path = tmp_path / "multi_frame_unknown.foo"
    path.write_text(
        "\n".join(
            [
                "2",
                "first molecule",
                "H 0.0 0.0 0.0",
                "H 0.0 0.0 0.7",
                "3",
                "second molecule",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 0.7",
                "H 0.7 0.0 0.0",
            ]
        ),
        encoding="utf-8",
    )

    file_model = AutoParser(str(path), n_jobs=1)[0]
    frame = file_model[0]

    assert file_model.detected_format_id == "xyz"
    assert len(file_model) == 1
    assert frame.atoms == [1, 1]
    assert tuple(frame.coords.shape) == (2, 3)
    assert frame.coords.to("angstrom").magnitude[1, 2] == pytest.approx(0.7)
