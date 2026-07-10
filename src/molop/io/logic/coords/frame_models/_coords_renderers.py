from typing import Any, Literal

from rdkit import Chem


def render_xyz_frame(target: Any, *, comment: str | None = None, stored_comment: str = "") -> str:
    comment_to_use = (
        comment or stored_comment or f"charge {target.charge} multiplicity {target.multiplicity}\n"
    )
    return (
        f"{len(target.atoms)}\n"
        + f"comment {comment_to_use.strip()}\n"
        + "\n".join(
            [
                f"{atom:10s}{x:18.10f}{y:18.10f}{z:18.10f}"
                for atom, (x, y, z) in zip(target.atom_symbols, target.coords.m, strict=True)
            ]
        )
    )


def render_sdf_frame(target: Any, *, engine: Literal["rdkit", "openbabel"] = "rdkit") -> str:
    if engine == "rdkit":
        rdmol = getattr(target, "qm_embedded_rdmol", target.rdmol)
        return Chem.MolToMolBlock(rdmol) if rdmol else ""
    if engine == "openbabel":
        return str(target.omol.write("sdf")) if target.omol else ""
    raise ValueError(f"Unsupported engine: {engine}")


def render_smi_frame(target: Any) -> str:
    return target.to_canonical_SMILES()
