"""
Author: TMJ
Date: 2024-08-31 17:52:33
LastEditors: TMJ
LastEditTime: 2024-08-31 17:54:13

Description: Re-implementation of SPMS descriptor [A Molecular Stereostructure Descriptor based on 
Spherical Projection](https://www.thieme-connect.de/products/ejournals/abstract/10.1055/s-0040-1705977).
GitHub repository: https://github.com/licheng-xu-echo/SPMS.git
"""

from typing import List, Literal, Sequence, Tuple, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, computed_field
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from rdkit.Geometry import Point3D

from molop.logger.logger import moloplogger
from molop.structure.geometry import (
    rotate_mol_anchor_to_axis,
    rotate_mol_anchor_to_plane,
    translate_mol_anchor,
    translate_mol,
)

DEBUG_TAG = "SPMS"


def check_dependencies():
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        moloplogger.error(
            f"{DEBUG_TAG}: matplotlib and seaborn are required for SPMS descriptor demonstration. "
            "Please install them by running 'pip install matplotlib seaborn'."
        )
        raise ImportError(
            f"{DEBUG_TAG}: matplotlib and seaborn are required for SPMS descriptor demonstration. "
            "Please install them by running 'pip install matplotlib seaborn'."
        )


pt = Chem.GetPeriodicTable()


# TODO no rdmol necessary


class SPMSCalculator(BaseModel):
    """
    SPMS descriptor calculator.

    Attributes:
        rdmol (Chem.rdchem.Mol):
            RDKit molecule object.
        anchor_list (Sequence[int]):
            List of anchor atom ids to keep the invariance of 3D geometry. Default is None, which means using Centroid to
            origin and the atom nearest to the centroid to X axis, and the atom farthest from the centroid to XY face.
        sphere_radius (Union[float, None]):
            Sphere radius. Default is None to use the largest possible radius for each input molecule.
            If you want to use a fixed radius, set this parameter to a float value.
        latitudinal_resolution (int):
            Number of splits on the latitudinal axis. Default is 40.
        longitude_resolution (int):
            Number of splits on the longitude axis. Default is 40.
        precision (int):
            Precision of the SPMS descriptor. Default is 8.
        custom_first_anchors (Union[Sequence[int], None]):
            List of atom ids to use as the first anchor. Default is None.
        custom_second_anchors (Union[Sequence[int], None]):
            List of atom ids to use as the second anchor. Default is None.
        custom_third_anchors (Union[Sequence[int], None]):
            List of atom ids to use as the third anchor. Default is None.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    rdmol: Chem.rdchem.Mol = Field(description="RDKit molecule object.")
    anchor_list: Union[Sequence[int], None] = Field(
        default=None,
        description="List of anchor atom ids to keep the invariance of 3D geometry. Default is None, "
        "which means using Centroid to origin and the atom nearest to the centroid to X axis, and the "
        "atom farthest from the centroid to XY face.",
    )
    sphere_radius: Union[float, None] = Field(
        default=None,
        description="Sphere radius. Default is None to use the largest possible radius for each input "
        "molecule. If you want to use a fixed radius, set this parameter to a float value.",
    )
    atom_radius: Literal["vdw", "covalent"] = Field(
        default="vdw",
        description="Atom radius type. Default is 'vdw', which means using Van der Waals radii as atom "
        "radii. If you want to use covalent radii, set this parameter to 'covalent'.",
    )
    latitudinal_resolution: int = Field(
        default=40,
        description="Number of splits on the latitudinal axis. Default is 40.",
    )
    longitude_resolution: int = Field(
        default=40,
        description="Number of splits on the longitude axis. Default is 40.",
    )
    precision: int = Field(
        default=8,
        description="Precision of the SPMS descriptor. Default is 8.",
    )
    custom_first_anchors: Union[Sequence[int], None] = Field(
        default=None,
        description="List of atom ids to use as the first anchor. Default is None.",
    )
    custom_second_anchors: Union[Sequence[int], None] = Field(
        default=None,
        description="List of atom ids to use as the second anchor. Default is None.",
    )
    custom_third_anchors: Union[Sequence[int], None] = Field(
        default=None,
        description="List of atom ids to use as the third anchor. Default is None.",
    )
    _rdmol_oriented: Chem.rdchem.Mol = PrivateAttr(default=None)

    @property
    def radius_list(self) -> np.ndarray:
        """
        Get the radius list of each atom.

        Returns:
            np.ndarray: List of atom radii.
        """
        if self.atom_radius == "vdw":
            return np.array(
                [pt.GetRvdw(atom.GetAtomicNum()) for atom in self.rdmol.GetAtoms()]
            )
        elif self.atom_radius == "covalent":
            return np.array(
                [pt.GetRcovalent(atom.GetAtomicNum()) for atom in self.rdmol.GetAtoms()]
            )
        else:
            raise ValueError(
                f"Invalid atom radius type: {self.atom_radius}. Only 'vdw' and 'covalent' are supported."
            )

    @property
    def rdmol_oriented(self) -> Chem.rdchem.Mol:
        """
        Orient the molecule to the SPMS geometry.

        Returns:
            Chem.rdchem.Mol: RDKit molecule object with SPMS geometry.
        """
        if self._rdmol_oriented is None:
            self._rdmol_oriented = geometry_initialize(
                self.rdmol,
                self.anchor_list,
                custom_first_anchors=self.custom_first_anchors,
                custom_second_anchors=self.custom_second_anchors,
                custom_third_anchors=self.custom_third_anchors,
            )
        return self._rdmol_oriented

    @property
    def sphere_radii(self) -> float:
        """
        Get the sphere radii.

        Returns:
            float: Sphere radii.
        """
        if self.sphere_radius is not None:
            return self.sphere_radius
        else:
            return get_default_sphere_radius(self.rdmol_oriented)

    @property
    def positions(self) -> np.ndarray:
        """
        Get the positions of the atoms.

        Returns:
            np.ndarray: Positions of the atoms.
        """
        return self.rdmol_oriented.GetConformer().GetPositions() + np.array(
            [0.000001, 0.000001, 0.000001]
        )

    @property
    def sphere_mesh(self) -> np.ndarray:
        """
        Get the sphere mesh.

        Returns:
            np.ndarray: Sphere mesh.
        """
        return get_sphere_mesh(
            self.sphere_radii,
            self.latitudinal_resolution,
            self.longitude_resolution,
        )

    @property
    def SPMS(self) -> np.ndarray:
        """
        Calculate SPMS descriptor.

        Returns:
            np.ndarray: SPMS descriptor.
        """
        # initialize the geometry
        positions: np.ndarray = self.positions
        radius = self.radius_list
        sphere_radius = self.sphere_radii

        # get the sphere mesh
        mesh_xyz = self.sphere_mesh

        # get the inner or outer sphere
        all_cross = np.array(
            [
                np.cross(atom_vec.reshape(-1, 3), mesh_xyz, axis=1)
                for atom_vec in positions
            ]
        ).transpose(1, 0, 2)
        mesh_xyz_h = np.linalg.norm(all_cross, axis=2) / sphere_radius
        dot = np.dot(mesh_xyz, positions.T)
        atom_vec_norm = np.linalg.norm(positions, axis=1).reshape(-1, 1)
        mesh_xyz_norm = np.linalg.norm(mesh_xyz, axis=1).reshape(-1, 1)
        orthogonal_mesh = dot / np.dot(mesh_xyz_norm, atom_vec_norm.T)
        # cross_det
        cross_det = mesh_xyz_h <= radius
        # orthogonal_det
        orthogonal_det = np.arccos(orthogonal_mesh) <= np.pi * 0.5
        double_correct = np.array([orthogonal_det, cross_det]).all(axis=0)
        double_correct_index = np.array(np.where(double_correct == True)).T
        d_1, d_2 = np.zeros_like(mesh_xyz_h), np.zeros_like(mesh_xyz_h)

        psi = np.linalg.norm(positions, axis=1)
        for item in double_correct_index:
            d_1[item[0]][item[1]] = (
                max((psi[item[1]] ** 2 - mesh_xyz_h[item[0]][item[1]] ** 2), 0) ** 0.5
            )
            d_2[item[0]][item[1]] = (
                radius[item[1]] ** 2 - mesh_xyz_h[item[0]][item[1]] ** 2
            ) ** 0.5

        sphere_descriptors = sphere_radius - d_1 - d_2
        sphere_descriptors_compact = sphere_descriptors.min(1)
        sphere_descriptors_reshaped = sphere_descriptors_compact.reshape(
            (self.latitudinal_resolution, 2 * self.longitude_resolution)
        )
        sphere_descriptors_reshaped = sphere_descriptors_reshaped.round(self.precision)
        if self.anchor_list is not None and len(self.anchor_list) == 1:
            sphere_descriptors_init = (
                np.zeros((self.latitudinal_resolution, 2 * self.longitude_resolution))
                + sphere_radius
                - radius[self.anchor_list[0]]
            )
            sphere_descriptors_final = np.min(
                np.concatenate(
                    [
                        sphere_descriptors_reshaped.reshape(
                            self.latitudinal_resolution,
                            2 * self.longitude_resolution,
                            1,
                        ),
                        sphere_descriptors_init.reshape(
                            self.latitudinal_resolution,
                            2 * self.longitude_resolution,
                            1,
                        ),
                    ],
                    axis=2,
                ),
                axis=2,
            )
        else:
            sphere_descriptors_final = sphere_descriptors_reshaped
        return sphere_descriptors_final

    @property
    def quater_descriptors(
        self,
    ) -> Tuple[float, float, float, float]:
        spms = self.SPMS

        left_top_desc_sum = spms[
            : self.latitudinal_resolution // 2, : self.latitudinal_resolution
        ].sum()
        right_top_desc_sum = spms[
            : self.latitudinal_resolution // 2, self.latitudinal_resolution :
        ].sum()
        left_bottom_desc_sum = spms[
            self.latitudinal_resolution // 2 :, : self.latitudinal_resolution
        ].sum()
        right_bottom_desc_sum = spms[
            self.latitudinal_resolution // 2 :, self.latitudinal_resolution :
        ].sum()

        _sum = spms.sum()

        left_top_desc_partial = left_top_desc_sum / _sum
        right_top_desc_partial = right_top_desc_sum / _sum
        left_bottom_desc_partial = left_bottom_desc_sum / _sum
        right_bottom_desc_partial = right_bottom_desc_sum / _sum

        return (
            left_top_desc_partial,
            right_top_desc_partial,
            left_bottom_desc_partial,
            right_bottom_desc_partial,
        )

    def draw_spms_as_heatmap(
        self,
        xy_label=True,
        cbar=True,
        scale_max: float = None,
        scale_min: float = None,
        height: float = 6,
        width: float = 12,
        cmap: str = "RdBu",
    ):
        check_dependencies()
        import seaborn as sns

        spms = self.SPMS

        if scale_max is None:
            scale_max = np.max(spms)
        if scale_min is None:
            scale_min = np.min(spms)
        g = sns.heatmap(self.SPMS, cmap=cmap, cbar=cbar, vmin=scale_min, vmax=scale_max)
        g.set_title("SPMS descriptor")
        g.figure.set_size_inches(width, height)
        if xy_label:
            g.set_xlabel("Latitudinal")
            g.set_ylabel("Longitude")
        return g

    def draw_spms_as_scatter(
        self,
        height: float = 6,
        width: float = 12,
        cmap: str = "RdBu",
        xy_label=True,
    ):
        check_dependencies()
        import matplotlib.pyplot as plt
        import seaborn as sns

        spms = self.SPMS
        rows, cols = np.where(spms != 0)
        weights = spms[rows, cols]

        g = sns.scatterplot(x=cols, y=rows, hue=weights, palette=cmap)
        g.set_title("SPMS descriptor")
        g.figure.set_size_inches(width, height)
        if xy_label:
            g.set_xlabel("Latitudinal")
            g.set_ylabel("Longitude")
        return g


def get_sphere_mesh(
    sphere_radius: float, latitudinal_resolution: int, longitude_resolution: int
) -> np.ndarray:
    theta_screenning = np.linspace(0, np.pi, latitudinal_resolution)
    phi_screenning = np.linspace(0, 2 * np.pi, 2 * longitude_resolution)
    PHI, THETA = np.meshgrid(phi_screenning, theta_screenning)
    x = sphere_radius * np.sin(THETA) * np.cos(PHI)
    y = sphere_radius * np.sin(THETA) * np.sin(PHI)
    z = sphere_radius * np.cos(THETA)
    mesh_xyz = np.concatenate(
        (
            x.flatten().reshape(1, -1),
            y.flatten().reshape(1, -1),
            z.flatten().reshape(1, -1),
        ),
        axis=0,
    ).T
    return mesh_xyz


def geometry_initialize(
    rdmol: Chem.rdchem.Mol,
    anchor_list: Union[Sequence[int], None] = None,
    *,
    custom_first_anchors: Union[Sequence[int], None] = None,
    custom_second_anchors: Union[Sequence[int], None] = None,
    custom_third_anchors: Union[Sequence[int], None] = None,
) -> Chem.rdchem.Mol:
    rwmol = Chem.RWMol(rdmol)
    # if user provides custom anchors, use them
    if all(
        anchor is not None
        for anchor in [
            custom_first_anchors,
            custom_second_anchors,
            custom_third_anchors,
        ]
    ):
        translate_mol_anchor(rwmol, custom_first_anchors)
        rotate_mol_anchor_to_axis(rwmol, custom_second_anchors, "z")
        rotate_mol_anchor_to_plane(rwmol, custom_third_anchors, "zy")
        return rwmol.GetMol()

    if anchor_list is None:
        center: Point3D = ComputeCentroid(rdmol.GetConformer(), ignoreHs=False)
    elif len(anchor_list) >= 1:
        center = np.mean(
            [rwmol.GetConformer().GetAtomPosition(i) for i in anchor_list],
            axis=0,
        )
    dis_list = np.linalg.norm(
        rwmol.GetConformer().GetPositions() - np.array(center), axis=1
    )
    farthest_atom = np.argmax(dis_list)
    dis_list[np.argmin(dis_list)] = np.max(dis_list)
    nearest_atom = np.argmin(dis_list)
    moloplogger.debug(
        f"{DEBUG_TAG}: Nearest atom: {nearest_atom}, Farthest atom: {farthest_atom}"
    )
    # put the centroid to origin
    translate_mol(rwmol, center * -1)
    # rotate the nearest atom to Z axis
    rotate_mol_anchor_to_axis(rwmol, nearest_atom, "z")
    # rotate the farthest atom to ZY plane
    rotate_mol_anchor_to_plane(rwmol, farthest_atom, "zy")
    return rwmol.GetMol()


def get_default_sphere_radius(rdmol: Chem.rdchem.Mol) -> float:
    """
    Get the default sphere radius for the input molecule.

    Parameters:
        rdmol (Chem.rdchem.Mol):
            RDKit molecule object.

    Returns:
        float: Default sphere radius.
    """
    atom_radius_vector = np.array(
        [pt.GetRvdw(atom.GetAtomicNum()) for atom in rdmol.GetAtoms()]
    )
    dis_vector = np.linalg.norm(rdmol.GetConformer().GetPositions(), axis=1, ord=2)
    return np.max(atom_radius_vector + dis_vector)
