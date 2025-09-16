"""
Author: TMJ
Date: 2025-09-14 21:54:49
LastEditors: TMJ
LastEditTime: 2025-09-14 22:30:35
Description: 请填写简介
"""

import io

import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def draw_molecule_with_dof_effect(
    mol: Chem.Mol,
    size: tuple[int, int] = (800, 800),
    legend: str = "",
    use_svg: bool = True,
    return_image: bool = True,
    show_atom_idx: bool = False,
    alpha: float = 0.4,
    keep_key_atom_colors: bool = False,
):
    if not mol or mol.GetNumConformers() == 0:
        raise ValueError("Invalid input molecule or no conformer available.")
    conf = mol.GetConformer()
    z_coords = np.array([conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())])
    if z_coords.size > 1 and z_coords.max() != z_coords.min():
        norm_z = (z_coords - z_coords.min()) / (z_coords.max() - z_coords.min())
        proximity = 1.0 - norm_z
    else:
        proximity = np.full_like(z_coords, 0.5, dtype=float)
    light_color = np.array((0.95, 0.95, 0.95, alpha), dtype=np.float32)
    dark_color = np.array((0.0, 0.0, 0.0, 1), dtype=np.float32)
    if use_svg:
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    dopts = drawer.drawOptions()
    dopts.continuousHighlight = False
    dopts.circleAtoms = False
    dopts.addAtomIndices = show_atom_idx
    atom_colors = {
        i: light_color + proximity[i] * (dark_color - light_color)
        for i in range(mol.GetNumAtoms())
        if not keep_key_atom_colors
        or mol.GetAtomWithIdx(i).GetAtomicNum() not in (7, 8, 9, 15, 16, 17, 35, 53)
    }
    bond_colors = {
        i: (
            atom_colors[mol.GetBondWithIdx(i).GetBeginAtomIdx()]
            + atom_colors[mol.GetBondWithIdx(i).GetEndAtomIdx()]
        )
        / 2
        for i in range(mol.GetNumBonds())
        if not keep_key_atom_colors
        or mol.GetBondWithIdx(i).GetBeginAtom().GetAtomicNum()
        not in (7, 8, 9, 15, 16, 17, 35, 53)
        and mol.GetBondWithIdx(i).GetEndAtom().GetAtomicNum()
        not in (7, 8, 9, 15, 16, 17, 35, 53)
    }
    drawer.DrawMolecule(
        mol,
        legend=legend,
        highlightAtoms=tuple(atom_colors.keys()),
        highlightAtomColors={k: tuple(v.tolist()) for k, v in atom_colors.items()},
        highlightBonds=tuple(bond_colors.keys()),
        highlightBondColors={k: tuple(v.tolist()) for k, v in bond_colors.items()},
    )
    drawer.FinishDrawing()
    if isinstance(drawer, rdMolDraw2D.MolDraw2DSVG):
        if return_image:
            try:
                from IPython.display import SVG
            except ImportError as e:
                raise ImportError(
                    "To return SVG image, please install IPython package."
                ) from e
            return SVG(drawer.GetDrawingText())
        return drawer.GetDrawingText()
    elif isinstance(drawer, rdMolDraw2D.MolDraw2DCairo):
        if return_image:
            return Image.open(io.BytesIO(drawer.GetDrawingText()))
        return drawer.GetDrawingText()
