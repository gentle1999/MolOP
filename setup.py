'''
Author: TMJ
Date: 2023-03-22 10:30:12
LastEditors: TMJ
LastEditTime: 2023-10-30 14:02:50
Description: 请填写简介
'''
from setuptools import setup

setup(
    name='molop',
    version='0.1.0',
    url='http://10.72.201.45:13000/TMJ/MolOP',
    author='TMJ',
    author_email='mj_t@zju.edu.cn',
    packages=['molop', 'molop.mol', 'molop.io'],
    install_requires=[
        # List of dependencies required by your project
        'rdkit',
        'scipy',
        'numpy',
        'pandas',
        'networkx',
        'tqdm',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
