"""
Author: TMJ
Date: 2024-02-23 10:18:51
LastEditors: TMJ
LastEditTime: 2024-02-25 18:15:28
Description: 请填写简介
"""
from typing import Dict, List, Literal, Sequence, Tuple, Union

import matplotlib as mpl
import numpy as np
import trimesh
from rdkit import Chem
from sklearn.preprocessing import MinMaxScaler, RobustScaler, StandardScaler
from trimesh.visual import ColorVisuals
from trimesh.visual.color import vertex_to_face_color

from molop.io.types import BLOCKTYPES, QMBLOCKTYPES
from molop.logger.logger import logger
from molop.utils.consts import cutoff_epsilon_sigma_dict, k, qe

pt = Chem.GetPeriodicTable()
radius_func = {
    "vdw": pt.GetRvdw,
    "covalent": pt.GetRcovalent,
    "atomic": pt.GetRb0,
}


def clac_cutoff_epsilon_sigma(atom1: str, atom2: str) -> Tuple[float, float, float]:
    if atom1 == atom2:
        cutoff, epsilon, sigma = cutoff_epsilon_sigma_dict[atom1]
    else:
        epsilon = np.sqrt(
            cutoff_epsilon_sigma_dict[atom1][1] * cutoff_epsilon_sigma_dict[atom2][1]
        )
        sigma = (
            cutoff_epsilon_sigma_dict[atom1][2] + cutoff_epsilon_sigma_dict[atom2][2]
        ) / 2
        cutoff = 2.5 * sigma
    return cutoff, epsilon, sigma


def build_molecule_sphere_mesh(
    block: BLOCKTYPES,
    radius: Union[Literal["vdw", "covalent", "atomic"], Sequence[float]] = "vdw",
    subdivision_iterations: int = 2,
):
    """
    Build a mesh from a molecule block based on sphere model.
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


class MoleculeMeshLayer:
    def __init__(self, block: BLOCKTYPES) -> None:
        self._block: BLOCKTYPES = block
        self._mesh = None
        self._visible = True
        self.__cutoff: np.ndarray = None
        self.__epsilon: np.ndarray = None
        self.__sigma: np.ndarray = None
        self.__color_values: np.ndarray = None

    @property
    def visible(self) -> bool:
        return self._visible

    @property
    def values(self) -> np.ndarray:
        if self.__color_values is None:
            if self._mesh is None:
                raise ValueError("The mesh has not been built yet.")
            return np.zeros(self._mesh.vertices.shape[0])
        return self.__color_values

    def build_sphere_mesh(
        self,
        radius: Union[Literal["vdw", "covalent", "atomic"], Sequence[float]] = "vdw",
        subdivision_iterations: int = 2,
    ) -> trimesh.Trimesh:
        self._mesh = build_molecule_sphere_mesh(
            self._block,
            radius=radius,
            subdivision_iterations=subdivision_iterations,
        )
        return self.mesh

    @property
    def mesh(self) -> trimesh.Trimesh:
        if self._visible:
            return self._mesh
        return None

    def colorize(
        self,
        values: Union[Literal["vdw", "electric_potential"], np.ndarray] = "vdw",
        cmap: str = "coolwarm",
        scaler: Literal["minmax", "standard", "robust"] = "minmax",
        **kwargs,
    ):
        if scaler == "minmax":
            _scaler = MinMaxScaler()
        elif scaler == "standard":
            _scaler = StandardScaler()
        elif scaler == "robust":
            _scaler = RobustScaler()
        else:
            raise ValueError(f"Scaler {scaler} is not supported")
        if isinstance(values, str):
            if values == "vdw":
                _values = self.Lennard_Jones_potential(**kwargs)
            elif values == "electric_potential":
                _values = self.electric_potential(**kwargs)
        else:
            assert values.shape[0] == self.mesh.vertices.shape[0], (
                f"The number of values {values.shape[0]} must be equal to the number of vertices "
                f"{self.mesh.vertices.shape[0]}"
            )
            _values = values
        self.__color_values = _values
        colormap = mpl.colormaps[cmap]
        vertex_colors = colormap(
            _scaler.fit_transform(_values.reshape(-1, 1)).flatten()
        )[np.newaxis, :, :3][0]
        face_colors = vertex_to_face_color(vertex_colors, self.mesh.faces)
        self.mesh.visual = ColorVisuals(
            vertex_colors=vertex_colors, face_colors=face_colors
        )
        return self.mesh

    def Lennard_Jones_potential(self, probe="He"):
        if self._mesh is None:
            raise ValueError("The mesh has not been built yet.")
        self._get_cutoff_epsilon_sigma(probe)
        vertices = self._mesh.vertices
        self.__Lennard_Jones_cache = np.apply_along_axis(
            self._get_Lennard_Jones_potential, 1, vertices
        )
        return self.__Lennard_Jones_cache

    def _get_Lennard_Jones_potential(self, pos: np.ndarray):
        cutoff, epsilon, sigma = self._get_cutoff_epsilon_sigma()
        diff = self._block.dimensionless_coords - pos
        dis = np.linalg.norm(diff, axis=1)
        index_cutoff = np.where(dis < cutoff)[0]
        temp_var = sigma[index_cutoff] / dis[index_cutoff]
        return np.sum(4 * epsilon[index_cutoff] * ((temp_var) ** 12 - (temp_var) ** 6))

    def _get_cutoff_epsilon_sigma(self, probe="He"):
        if self.__cutoff is None or self.__epsilon is None or self.__sigma is None:
            cutoff_epsilon_sigma = np.array(
                [clac_cutoff_epsilon_sigma(probe, atom) for atom in self._block.atoms]
            )
            self.__cutoff, self.__epsilon, self.__sigma = (
                cutoff_epsilon_sigma[:, 0],
                cutoff_epsilon_sigma[:, 1],
                cutoff_epsilon_sigma[:, 2],
            )
        return self.__cutoff, self.__epsilon, self.__sigma

    def electric_potential(self, charge: Literal["mulliken", "npa"] = "mulliken"):
        if self._mesh is None:
            raise ValueError("The mesh has not been built yet.")
        if charge == "npa":
            if self._block.npa_charges is None:
                logger.warning(
                    "The NBO charges are not available, try Mulliken charges instead"
                )
                charges = self._block.mulliken_charges
        else:
            charges = self._block.mulliken_charges
        if charges is None:
            raise ValueError("The charges are not available")
        return np.apply_along_axis(
            self._get_ep, 1, pos=self._mesh.vertices, charges=charges
        )

    def _get_ep(self, pos: np.ndarray, charges: np.ndarray):
        diff = self._block.dimensionless_coords - pos
        dis = np.linalg.norm(diff, axis=1) * 1e-10
        charge = charges * qe
        return np.sum(k * charge / dis)

    def show(self):
        if self.mesh is not None:
            return self.mesh.show()


class MoleculeMesh:
    def __init__(self, block: BLOCKTYPES) -> None:
        self._block: BLOCKTYPES = block
        self._meshes: list[MoleculeMeshLayer] = []

    @property
    def meshes(self) -> List[MoleculeMeshLayer]:
        return self._meshes

    def add_layer(self) -> None:
        self._meshes.append(MoleculeMeshLayer(self._block))

    def __getitem__(self, key: int) -> MoleculeMeshLayer:
        return self._meshes[key]

    def __len__(self) -> int:
        return len(self._meshes)

    def __iter__(self):
        self.__index = 0
        return self

    def __delitem__(self, key: int):
        self._meshes.pop(key)

    def __next__(
        self,
    ) -> trimesh.Trimesh:
        if self.__index >= len(self):
            raise StopIteration
        else:
            self.__index += 1
            return self[self.__index - 1]

    def __repr__(self) -> str:
        return f"MoleculeMesh({self._block})"

    def show(self):
        """
        Show the all meshes of the molecule
        """
        if len(self.meshes) == 0:
            return None
        if len(self.meshes) == 1:
            return self.meshes[0].show()
        totol_mesh = self.meshes[0].mesh
        for mesh in self.meshes[1:]:
            if mesh.mesh is not None:
                totol_mesh = trimesh.util.concatenate(totol_mesh, mesh.mesh)
        return totol_mesh.show()
