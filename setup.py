"""
Author: TMJ
Date: 2023-10-30 14:00:45
LastEditors: TMJ
LastEditTime: 2024-01-14 22:22:34
Description: 请填写简介
"""
"""
Author: TMJ
Date: 2023-10-30 14:00:45
LastEditors: TMJ
LastEditTime: 2024-01-11 09:42:04
Description: 请填写简介
"""
from setuptools import setup

setup(
    name="molop",
    version="0.1.1",
    description="Molcule OPerator",
    url="http://10.72.201.58:13000/TMJ/MolOP",
    author="TMJ",
    author_email="mj_t@zju.edu.cn",
    platforms=["any"],
    packages=[
        "molop",
        "molop.mol",
        "molop.utils",
        "molop.unit",
        "molop.logger",
        "molop.structure",
        "molop.descriptor",
        "molop.io",
        "molop.io.bases",
        "molop.io.qm_file",
        "molop.io.coords_file",
    ],
    install_requires=[
        # List of dependencies required by your project
        "rdkit>=2023.9.1",
        "scipy",
        "numpy",
        "pandas",
        "networkx",
        "tqdm",
        "pint",
    ],
    extras_require={"full": ["mordred", "dscribe", "ase"]},
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
