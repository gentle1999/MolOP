"""
Author: TMJ
Date: 2023-10-30 14:00:45
LastEditors: TMJ
LastEditTime: 2024-01-16 15:02:58
Description: 请填写简介
"""
from setuptools import setup

setup(
    name="molop",
    version="0.1.13",
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
        "rdkit>=2023.9.4",
        "scipy>=1.10.1",
        "numpy>=1.24.2",
        "pandas>=2.0.3",
        "networkx>=2.8.8",
        "tqdm>=4.66.1",
        "pint>=0.21.1",
    ],
    extras_require={"full": ["mordred>=1.2.0", "dscribe>=2.1.0", "ase>=3.22.1"]},
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
