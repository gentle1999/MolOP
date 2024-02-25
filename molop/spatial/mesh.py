"""
Author: TMJ
Date: 2024-02-23 10:18:51
LastEditors: TMJ
LastEditTime: 2024-02-25 18:15:28
Description: 请填写简介
"""
"""
Author: TMJ
Date: 2024-02-23 10:18:51
LastEditors: TMJ
LastEditTime: 2024-02-25 14:55:18
Description: 请填写简介
"""
from typing import Literal, Sequence, Union, Tuple

import matplotlib as mpl
import numpy as np
import trimesh
from rdkit import Chem
from trimesh.visual import ColorVisuals

from molop.io.types import BLOCKTYPES, QMBLOCKTYPES

pt = Chem.GetPeriodicTable()
k = 8.9875517923e9
qe = 1.602176634 * 10**-19
radius_func = {
    "vdw": pt.GetRvdw,
    "covalent": pt.GetRcovalent,
    "atomic": pt.GetRb0,
}

epsilon_sigma_dict = {
    "H": (4.47789, 0.552357),
    "He": (0.0009421, 0.49890),
    "Li": (1.04969, 2.2807),
    "Be": (0.572942, 1.71053),
    "B": (2.96703, 1.51453),
    "C": (6.36953, 1.35417),
    "N": (9.75379, 1.26508),
    "O": (5.12647, 1.17599),
    "F": (1.60592, 1.01562),
    "Ne": (0.0036471, 1.03344),
    "Na": (0.736745, 2.95778),
    "Mg": (0.0785788, 2.51233),
    "Al": (2.70067, 2.15597),
    "Si": (3.17431, 1.9778),
    "P": (5.0305, 1.90652),
    "S": (4.36927, 1.87089),
    "Cl": (4.48328, 1.81743),
    "Ar": (0.0123529, 1.88871),
    "K": (0.551799, 3.61705),
    "Ca": (0.132679, 3.13596),
    "Sc": (1.6508, 3.02906),
    "Ti": (1.18027, 2.85088),
    "V": (2.75249, 2.72615),
    "Cr": (1.53679, 2.4767),
    "Mn": (0.599888, 2.4767),
    "Fe": (1.18442, 2.35197),
    "Co": (1.27769, 2.24506),
    "Ni": (2.07572, 2.20943),
    "Cu": (2.04463, 2.35197),
    "Zn": (0.191546, 2.17379),
    "Ga": (1.0642, 2.17379),
    "Ge": (2.70171, 2.13816),
    "As": (3.9599, 2.12034),
    "Se": (3.38677, 2.13816),
    "Br": (1.97063, 2.13816),
    "Kr": (0.0173276, 2.06689),
    "Rb": (0.468265, 3.91995),
    "Sr": (0.133923, 3.47451),
    "Y": (2.75975, 3.38542),
    "Zr": (3.05201, 3.11815),
    "Nb": (5.2782, 2.92215),
    "Mo": (4.47499, 2.74397),
    "Tc": (3.38159, 2.61924),
    "Ru": (1.96172, 2.60142),
    "Rh": (2.40582, 2.53015),
    "Pd": (1.37097, 2.4767),
    "Ag": (1.64976, 2.58361),
    "Cd": (0.0377447, 2.56579),
    "In": (0.811314, 2.53015),
    "Sn": (1.90057, 2.4767),
    "Sb": (3.08828, 2.4767),
    "Te": (2.63123, 2.45888),
    "I": (1.53938, 2.4767),
    "Xe": (0.023888, 2.49452),
    "Cs": (0.416642, 4.34759),
    "Ba": (1.9, 3.83086),
    "La": (2.49961, 3.68832),
    "Ce": (2.57008, 3.63487),
    "Pr": (1.29946, 3.61705),
    "Nd": (0.819605, 3.58141),
    "Pm": (3.24134, 3.54578),
    "Sm": (0.521122, 3.52796),
    "Eu": (0.429918, 3.52796),
    "Gd": (2.09956, 3.49232),
    "Tb": (1.39999, 3.45669),
    "Dy": (0.690055, 3.42105),
    "Ho": (0.690055, 3.42105),
    "Er": (0.738766, 3.3676),
    "Tm": (0.521122, 3.38542),
    "Yb": (0.130399, 3.33196),
    "Lu": (1.43315, 3.33196),
    "Hf": (3.36086, 3.11815),
    "Ta": (4.00343, 3.02906),
    "W": (6.86389, 2.88651),
    "Re": (4.43871, 2.69051),
    "Os": (4.26253, 2.56579),
    "Ir": (3.70287, 2.51233),
    "Pt": (3.1401, 2.42324),
    "Au": (2.3058, 2.42324),
    "Hg": (0.045414, 2.35197),
    "Tl": (0.577087, 2.58361),
    "Pb": (0.858988, 2.60142),
    "Bi": (2.07987, 2.63706),
    "Po": (1.89953, 2.49452),
    "At": (1.385442, 2.6727),
    "Rn": (0.0214992, 2.6727),
    "Fr": (0.3749778, 4.63267),
    "Ra": (1.71, 3.93777),
    "Ac": (2.249649, 3.83086),
    "Th": (2.313072, 3.6705),
    "Pa": (1.169514, 3.56359),
    "U": (0.7376445, 3.49232),
    "Np": (2.917206, 3.38542),
    "Pu": (0.4690098, 3.33196),
    "Am": (0.3869262, 3.20724),
    "Cm": (1.889604, 3.01124),
    "Bk": (1.259991, 2.99342),
    "Cf": (0.6210495, 2.99342),
    "Es": (0.6210495, 2.93997),
    "Fm": (0.6648894, 2.9756),
    "Md": (0.4690098, 3.08251),
    "No": (0.1173591, 3.13596),
    "Lr": (1.289835, 2.86869),
    "Rf": (3.024774, 2.79742),
    "Db": (3.603087, 2.65488),
    "Sg": (6.177501, 2.54797),
    "Bh": (3.994839, 2.51233),
    "Hs": (3.836277, 2.38761),
    "Mt": (3.332583, 2.29852),
    "Ds": (2.82609, 2.2807),
    "Rg": (2.07522, 2.15597),
    "Cn": (0.0408726, 2.17379),
    "Nh": (0.5193783, 2.42324),
    "Fl": (0.7730892, 2.54797),
    "Mc": (1.871883, 2.88651),
    "Lv": (1.709577, 3.11815),
    "Ts": (1.2468978, 2.93997),
    "Og": (0.0193493, 2.79742),
}


def clac_epsilon_sigma(atom1: str, atom2: str) -> Tuple[float, float]:
    if atom1 == atom2:
        epsilon,sigma = epsilon_sigma_dict[atom1]
    else:
        epsilon = np.sqrt(epsilon_sigma_dict[atom1][0] * epsilon_sigma_dict[atom2][0])
        sigma = (epsilon_sigma_dict[atom1][1] + epsilon_sigma_dict[atom2][1]) / 2
    return epsilon, sigma 


def build_molecule_mesh(
    block: BLOCKTYPES,
    radius: Union[Literal["vdw", "covalent", "atomic"], Sequence[float]] = "vdw",
    subdivision_iterations: int = 2,
):
    """
    Build a mesh from a molecule block
    Parameters:
        block (BLOCKTYPES): A molecule block
        radius (Union[Literal["vdw", "covalent", "atomic"], Sequence[float]]): The radius of the spheres
        subdivision_iterations (int): The number of subdivision iterations
    Returns:
        trimesh.Trimesh: The mesh of the molecule
    """
    if isinstance(radius, str):
        radii = [radius_func[radius](atom) for atom in block.atoms]
    elif isinstance(radius, Sequence):
        assert len(radius) == len(
            block.atoms
        ), f"The number of radii {len(radius)} must be equal to the number of atoms {len(block.atoms)}"
        radii = radius
    centers = np.array(block.dimensionless_coords)
    spheres = [
        trimesh.primitives.Sphere(center=center, radius=radi, subdivisions=1)
        for center, radi in zip(centers, radii)
    ]
    merged_sphere: trimesh.base.Trimesh = trimesh.boolean.union(spheres)
    merged_sphere.remove_duplicate_faces()
    merged_sphere.remove_unreferenced_vertices()
    trimesh.repair.fill_holes(merged_sphere)
    vertices, faces = trimesh.remesh.subdivide_loop(
        vertices=merged_sphere.vertices,
        faces=merged_sphere.faces,
        iterations=subdivision_iterations,
    )
    molecule_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    trimesh.repair.fill_holes(molecule_mesh)
    return molecule_mesh


class MoleculeMesh:
    def __init__(self, block: BLOCKTYPES) -> None:
        self._block: BLOCKTYPES = block
        self._mesh: trimesh.Trimesh = None

    def build_mesh(
        self,
        radius: Union[Literal["vdw", "covalent", "atomic"], Sequence[float]] = "vdw",
        subdivision_iterations: int = 2,
    ) -> None:
        self._mesh = build_molecule_mesh(
            self._block,
            radius=radius,
            subdivision_iterations=subdivision_iterations,
        )
        return self._mesh

    @property
    def mesh(self) -> trimesh.Trimesh:
        if self._mesh is None:
            return self.build_mesh()
        else:
            return self._mesh

    def __repr__(self) -> str:
        return f"MoleculeMesh({self._block})"

    def colorize(
        self,
        values: Union[Literal["vdw", "electric_potential"], np.ndarray] = "vdw",
        cmap: str = "coolwarm",
    ):
        if isinstance(values, str):
            if values == "vdw":
                _values = self.vdw_potential()
            elif values == "electric_potential":
                _values = self.electric_potential()
        else:
            _values = values
        
        colormap = mpl.colormaps[cmap]
        vertex_colors = colormap(_values)[np.newaxis, :, :3][0]
        face_colors = trimesh.visual.color.vertex_to_face_color(
            vertex_colors, self.mesh.faces
        )
        self.mesh.visual = ColorVisuals(
            vertex_colors=vertex_colors, face_colors=face_colors
        )
        return self.mesh

    def vdw_potential(self, probe="Ar"):
        coords = np.array(self._block.dimensionless_coords)
        epsilon_sigma = np.array(
            [clac_epsilon_sigma(probe, atom) for atom in self._block.atoms]
        )
        epsilon = epsilon_sigma[:, 0]
        sigma = epsilon_sigma[:, 1]

        def get_vdw_potential(pos):
            diff = coords - pos
            dis = np.linalg.norm(diff, axis=1)
            return np.sum(
                4 * epsilon * ((sigma / dis) ** 12 - (sigma / dis) ** 6)
            )

        vertices = self.mesh.vertices
        return np.apply_along_axis(get_vdw_potential, 1, vertices)

    def electric_potential(self):
        coords = np.array(self._block.dimensionless_coords)
        charges = np.array(self._block.partial_charges)

        def get_esp(pos):
            diff = (coords - pos) * 1e-10
            dis = np.linalg.norm(diff, axis=1)
            charge = charges * qe
            return np.sum(k * charge / dis)

        vertices = self.mesh.vertices
        return np.apply_along_axis(get_esp, 1, vertices)
