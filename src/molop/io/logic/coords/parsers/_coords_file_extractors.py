from __future__ import annotations

from rdkit import Chem

from molop.io.codec_exceptions import FormatMismatchError


def split_xyz_frames(file_content: str) -> list[str]:
    lines = file_content.splitlines()
    anchor = 0
    frames: list[str] = []
    while anchor < len(lines):
        try:
            num_atoms = int(lines[anchor].strip())
        except ValueError:
            anchor += 1
            continue
        if num_atoms <= 0:
            raise FormatMismatchError("Not an XYZ file: atom count must be positive.")
        if anchor + num_atoms + 2 > len(lines):
            raise FormatMismatchError("Not an XYZ file: incomplete frame.")
        frames.append("\n".join(lines[anchor : anchor + num_atoms + 2]))
        anchor += num_atoms + 2
    if not frames:
        raise FormatMismatchError("Not an XYZ file: no atom-count header found.")
    return frames


def split_sdf_frames(file_content: str) -> list[str]:
    suppl = Chem.SDMolSupplier()
    suppl.SetData(file_content, removeHs=False, sanitize=False)
    frames = [Chem.MolToMolBlock(mol) for mol in suppl if mol is not None]
    if not frames:
        raise FormatMismatchError("Not an SDF/MOL file: no valid mol block found.")
    return frames


def split_smi_frames(file_content: str) -> list[str]:
    lines = [line for line in file_content.splitlines() if line.strip()]
    if not lines:
        raise FormatMismatchError("Not a SMILES file: empty file.")
    return lines
