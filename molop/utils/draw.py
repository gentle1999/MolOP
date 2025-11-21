"""
Author: TMJ
Date: 2025-09-14 21:54:49
LastEditors: TMJ
LastEditTime: 2025-11-16 21:56:04
Description: 请填写简介
"""

import io
from functools import lru_cache
from typing import overload

import numpy as np
from IPython.core.display import SVG
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdDepictor import Compute2DCoords

try:
    from IPython.display import SVG, display
except ImportError:
    display = print

DEFAULT_COLORS = {
    1: (0.55, 0.55, 0.55),  # H
    6: (0.2, 0.2, 0.2),  # C
    7: (0, 0, 1),  # N
    8: (1, 0, 0),  # O
    9: (0.2, 0.8, 0.8),  # F
    15: (1, 0.5, 0),  # P
    16: (0.8, 0.8, 0),  # S
    17: (0, 0.8, 0),  # Cl
    35: (0.5, 0.3, 0.1),  # Br
    53: (0.63, 0.12, 0.94),  # I
    0: (0.5, 0.5, 0.5),
}


@lru_cache(maxsize=1024)
def _get_atom_dof_color(
    base_color: tuple[float, float, float],
    proximity: float,
    min_alpha: float = 0.4,
    fog_color: np.ndarray = np.array([0.95, 0.95, 0.95]),
) -> tuple[float, float, float, float]:
    """
    Calculate the color of an atom based on its proximity to the camera.

    Parameters:
        base_color (tuple[float, float, float]): The base color of the atom (RGB, 0-1 range).
        proximity (float): The normalized proximity of the atom to the camera (0.0 = farthest, 1.0 = closest).
        min_alpha (float, optional): The minimum alpha value for the farthest atom (default is 0.4).
        fog_color (np.ndarray, optional): The color used to simulate fog (default is a light gray).

    Returns:
        tuple[float, float, float, float]: The final RGBA color of the atom.
    """
    base_rgb = np.array(base_color)
    dark_color_rgba = np.array([*base_rgb, 1.0])
    light_rgb = base_rgb * 0.2 + fog_color * 0.8
    light_color_rgba = np.array([*light_rgb, min_alpha])
    final_color = light_color_rgba + proximity * (dark_color_rgba - light_color_rgba)

    return tuple(final_color.tolist())


@overload
def draw_molecule_with_dof_effect(
    mol: Chem.Mol,
    size: tuple[int, int] = (800, 800),
    legend: str = "",
    use_svg: bool = True,
    return_image: bool = False,
    min_alpha: float = 0.4,
    keep_key_atom_colors: bool = True,
    addAtomIndices: bool = False,
    addBondIndices: bool = False,
) -> str: ...
@overload
def draw_molecule_with_dof_effect(
    mol: Chem.Mol,
    size: tuple[int, int] = (800, 800),
    legend: str = "",
    use_svg: bool = False,
    return_image: bool = True,
    min_alpha: float = 0.4,
    keep_key_atom_colors: bool = True,
    addAtomIndices: bool = False,
    addBondIndices: bool = False,
) -> Image.Image: ...
@overload
def draw_molecule_with_dof_effect(
    mol: Chem.Mol,
    size: tuple[int, int] = (800, 800),
    legend: str = "",
    use_svg: bool = True,
    return_image: bool = True,
    min_alpha: float = 0.4,
    keep_key_atom_colors: bool = True,
    addAtomIndices: bool = False,
    addBondIndices: bool = False,
) -> "SVG": ...
def draw_molecule_with_dof_effect(
    mol: Chem.Mol,
    size: tuple[int, int] = (800, 800),
    legend: str = "",
    use_svg: bool = True,
    return_image: bool = True,
    min_alpha: float = 0.4,
    keep_key_atom_colors: bool = True,
    addAtomIndices: bool = False,
    addBondIndices: bool = False,
):
    """
    Draw a molecule with depth-of-field effect based on Z-coordinates of atoms.

    Parameters:
        mol (Chem.Mol): The RDKit molecule object to draw.
        size (tuple[int, int], optional): The size of the output image (default is (800, 800)).
        legend (str, optional): The legend to display below the molecule (default is "").
        use_svg (bool, optional): Whether to use SVG format for the output (default is True).
        return_image (bool, optional): Whether to return the image as a PIL Image object (default is True).
        min_alpha (float, optional): The minimum alpha value for the farthest atom (default is 0.4).
        keep_key_atom_colors (bool, optional): Whether to keep the original colors of key atoms (default is True).
        addAtomIndices (bool, optional): Whether to add atom indices (default is False).
        addBondIndices (bool, optional): Whether to add bond indices (default is False).

    Returns:
        PIL.Image.Image or None: The drawn molecule image, or None if return_image is False.
    """
    if not mol:
        raise ValueError("Invalid RDKit molecule.")
    if mol.GetNumConformers() == 0:
        Compute2DCoords(mol)

    conf = mol.GetConformer()
    z_coords = np.array([conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())])

    if z_coords.size > 1 and z_coords.max() != z_coords.min():
        norm_z = (z_coords - z_coords.min()) / (z_coords.max() - z_coords.min())
        proximity = 1.0 - norm_z
    else:
        proximity = np.full_like(z_coords, 0.5, dtype=float)

    if use_svg:
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])

    dopts = drawer.drawOptions()
    dopts.continuousHighlight = False
    dopts.circleAtoms = False
    dopts.addAtomIndices = addAtomIndices
    dopts.addBondIndices = addBondIndices
    default_atom_colors = DEFAULT_COLORS.copy()
    carbon_color = default_atom_colors.get(6, (0.2, 0.2, 0.2))

    highlight_atom_colors = {}
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        atomic_num = atom.GetAtomicNum()
        base_color = default_atom_colors.get(atomic_num, carbon_color)
        highlight_atom_colors[i] = (
            _get_atom_dof_color(
                base_color=base_color, proximity=proximity[i], min_alpha=min_alpha
            )
            if keep_key_atom_colors
            else _get_atom_dof_color(
                base_color=carbon_color, proximity=proximity[i], min_alpha=min_alpha
            )
        )

    highlight_bond_colors = {}
    for i in range(mol.GetNumBonds()):
        bond = mol.GetBondWithIdx(i)
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        color1 = np.array(
            _get_atom_dof_color(
                base_color=carbon_color,
                proximity=proximity[atom1_idx],
                min_alpha=min_alpha,
            )
        )
        color2 = np.array(
            _get_atom_dof_color(
                base_color=carbon_color,
                proximity=proximity[atom2_idx],
                min_alpha=min_alpha,
            )
        )

        bond_color = (color1 + color2) / 2
        highlight_bond_colors[i] = tuple(bond_color.tolist())

    drawer.DrawMolecule(
        mol,
        legend=legend,
        highlightAtoms=tuple(range(mol.GetNumAtoms())),
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=tuple(range(mol.GetNumBonds())),
        highlightBondColors=highlight_bond_colors,
    )
    drawer.FinishDrawing()

    if isinstance(drawer, rdMolDraw2D.MolDraw2DSVG):
        if return_image:
            try:
                from IPython.display import SVG
            except ImportError as e:
                raise ImportError(
                    "To return SVG image, please install IPython package, or set return_image=False."
                ) from e
            return SVG(drawer.GetDrawingText())
        return drawer.GetDrawingText()
    elif isinstance(drawer, rdMolDraw2D.MolDraw2DCairo):
        if return_image:
            return Image.open(io.BytesIO(drawer.GetDrawingText()))
        return drawer.GetDrawingText()
