'''
Author: TMJ
Date: 2023-10-30 14:00:45
LastEditors: TMJ
LastEditTime: 2024-01-10 18:56:59
Description: 请填写简介
'''
"""
Author: TMJ
Date: 2023-03-22 10:30:12
LastEditors: TMJ
LastEditTime: 2024-01-09 10:39:28
Description: 请填写简介
"""
from setuptools import setup

setup(
    name="molop",
    version="0.1.0",
    url="http://10.72.201.45:13000/TMJ/MolOP",
    author="TMJ",
    author_email="mj_t@zju.edu.cn",
    packages=[
        "molop",
        "molop.mol",
        "molop.utils",
        "molop.unit",
        "molop.logger",
        "molop.structure",
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
