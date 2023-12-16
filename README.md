<!--
 * @Author: TMJ
 * @Date: 2023-10-30 13:36:49
 * @LastEditors: TMJ
 * @LastEditTime: 2023-11-02 16:40:51
 * @Description: 请填写简介
-->
# MolOP

This repository, **Molecule OPerator**, which is the basic molecule information extraction and operation unit of the [Open TS DataBase](http://10.72.201.45:13000/tmj/OTSDB-Core) project.

## Features
- Automatically extract molecule information from the Input and Output files of the common QM calculation softwares.
  - GJF bug fixing
  - XYZ done
  - SDF TODO
  - Other formats TODO
- Offer the moleculer geometry and structure edit functions.
  - Substructure replacement
  - Orientation change
  - Other functions
- No TS or Reaction support.

## Installation for Development

```bash
# clone the repository
git clone http://10.72.201.45:13000/tmj/MolOP.git
cd MolOP
conda create -n molop python=3.8
conda activate molop
# install the dependencies
pip install -r requirements.txt
pip install .
```

## Development

