<!--
 * @Author: TMJ
 * @Date: 2023-10-30 13:36:49
 * @LastEditors: TMJ
 * @LastEditTime: 2024-01-15 10:38:28
 * @Description: 请填写简介
-->
# MolOP

This repository, **Molecule OPerator**, which is the basic molecule information extraction and operation unit of the [Open TS DataBase](http://10.72.201.45:13000/tmj/OTSDB-Core) project.

## Features

- Automatically extract molecule information from the Input and Output files of the common QM calculation softwares. Offer a [molecular graph recovery algrothm](molop/structure/structure_recovery.py) from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process.
  - Coords file
    - GJF `Done`
    - XYZ `Done`
    - SDF `Done`
  - QM output file
    - G16 LOG `Done`
    - G16 IRC `Done`
    - xTB OUT `Done`
    - ORCA `TODO`
- Offer the moleculer geometry and structure edit functions. `TODO`
  - Substructure replacement `TODO`
  - Orientation change `TODO`
  - Other functions `TODO`

## Get Start

See the [Tutorial Notebook](tutorial/get_start.ipynb) for more details.

## Installation for Users

```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install --index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ molop
```

For additional descriptor calculation, you need to install the requirements below:

```bash
pip install --index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ molop[full]
```

## Installation for Developers

```bash
# clone the repository
git clone http://10.72.201.45:13000/tmj/MolOP.git
cd MolOP
conda create -n  
conda activate molop
# install the dependencies
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install .
```
