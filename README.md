<!--
 * @Author: TMJ
 * @Date: 2023-10-30 13:36:49
 * @LastEditors: TMJ
 * @LastEditTime: 2025-01-05 14:05:07
 * @Description: 请填写简介
-->
# MolOP

This repository, **Molecule OPerator**, which is a Pydantic-based molecule information extraction and operation unit. Any of the attributes of the molecule defined can be accessed with clear and unified type definition and easy-to-use API.

## Features

- Automatically extract molecule information from the Input and Output files of the common QM calculation softwares.
  - Coords file
    - GJF `Available`
    - XYZ `Available`
    - SDF `Available`
  - QM output file
    - G16 LOG `Available`
    - G16 FCHK `Available`
    - xTB OUT `Available`
    - ORCA `TODO`

- Offer a [molecular graph recovery algorithm](molop/structure/structure_recovery_alter.py) from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process. This algorithm is different from the rdDetermineBonds (Original code implemented by [Jensen group](https://github.com/jensengroup/xyz2mol) and integrated in RDKit since the 2022.09 release, which is not suitable for the free radicals and complex containing metal. See [structure_recovery_cases](docs/en/structure_recovery_cases.ipynb) to learn more about the difference).
  
  Although our algorithm overcome the free radicals and metal problem and tested on the [test_cases](tutorial/test_cases.ipynb) file, it is still not perfect. There is no denying that, rdDetermineBonds works well for normal organic molecules. Thus, we would give molecule structure recovered by rdDetermineBonds first, if error happens, we will use our algorithm to recover the molecule structure instead. We hope that this strategy can take advantage of both approaches.

- Offer the moleculer geometry and structure edit functions. `Doing`
  - Substructure replacement `Available`
  - User-friendly orientation change `Available`
  - Geometry-based TS check `Available`
  - [SMPS descriptor](https://github.com/licheng-xu-echo/SPMS.git) calculation `Available`
  - Other functions `TODO`

## Get Start

See the [Tutorial Notebook](tutorial/get_start.ipynb) for more details.

All teat cases shown in the [test_cases](tutorial/test_cases.ipynb) file.

## Installation

This package is actively developing, and this time is too early to be published to pypi. So you can install by the following.

### In ZJU intranet

You can use our self-host repository.

```bash
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop --upgrade
```

For additional descriptor calculation, you need to install the requirements below:

```bash
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop[full] --upgrade
```

### On Internet

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

or

```bash
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
pip install -e .
```

## Installation for Developers

```bash
# clone the repository
git clone https://github.com/gentle1999/MolOP.git
cd MolOP
# create a new environment (conda is not necessary)
conda create -n molop python=3.8 # The lowest python version is 3.8
conda activate molop
# install the dependencies
pip install poetry
poetry install --with dev
poetry install --all-extras
```

## Online Documentation (only avialable in ZJU intranet now)

Visit [MolOP Documentation](https://molop-gentle-7355819aa2dbcc7edc9420595fa823e6ffaebc1a874271edbf.pages.zjusct.io/).

## Start documentation server

```bash
mkdocs serve
```
