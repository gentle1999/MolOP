'''
Author: TMJ
Date: 2024-01-13 21:26:08
LastEditors: TMJ
LastEditTime: 2024-01-13 21:58:52
Description: 请填写简介
'''
from typing import List, Dict, Literal, Union

def calc_rdkit_descs(rdmol, desc_names:List[str]):
    from rdkit.Chem import Descriptors
    from rdkit.ML.Descriptors import MoleculeDescriptors

    if desc_names is None:
        desc_names = [desc_name[0] for desc_name in Descriptors._descList]
    desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
    return {
        desc_name: value
        for desc_name, value in zip(
            desc_names, desc_calc.CalcDescriptors(rdmol)
        )
    }


def calc_dscribe_descs(
    species,
    atoms_path: str = None,
    desc_names: Literal[
        "SOAP",
        "ACSF",
        "MBTR",
        "LMBTR",
    ] = None,
):
    from ase import io
    from dscribe.descriptors import (
        SOAP,
        ACSF,
        MBTR,
        LMBTR,
    )

    descs = {
        "SOAP": SOAP(
            species=species,
            periodic=False,
            r_cut=6.0,
            n_max=8,
            l_max=6,
        ),
        "ACSF": ACSF(
            species=species,
            r_cut=6.0,
            g2_params=[[1, 1], [1, 2], [1, 3]],
            g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]],
        ),
        "MBTR": MBTR(
            species=species,
            geometry={"function": "inverse_distance"},
            grid={"min": 0, "max": 1, "n": 100, "sigma": 0.1},
            weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
            periodic=False,
            normalization="l2",
        ),
        "LMBTR": LMBTR(
            species=species,
            geometry={"function": "distance"},
            grid={"min": 0, "max": 5, "n": 100, "sigma": 0.1},
            weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
            periodic=False,
            normalization="l2",
        ),
    }
    tmp_atoms = io.read(atoms_path or ".temp.sdf", format="sdf")
    if desc_names is None:
        return {desc_name: desc.create(tmp_atoms) for desc_name, desc in descs.items()}
    else:
        return {
            desc_name: desc.create(tmp_atoms)
            for desc_name, desc in descs.items()
            if desc_name in desc_names
        }
