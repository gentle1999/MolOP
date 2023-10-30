<!--
 * @Author: TMJ
 * @Date: 2023-10-30 13:36:49
 * @LastEditors: TMJ
 * @LastEditTime: 2023-10-30 17:55:58
 * @Description: 请填写简介
-->
# MolOP

This repository, **Molecule OPerator**, contains the code for the MolOP project, which is the molecule information extraction and operation unit of the [Open TS DataBase](http://10.72.201.45:13000/tmj/OTSDB-Core) project.

## Features
- Automatically extract molecule information from the Input and Output files of the common QM calculation softwares.
- Offer the moleculer geometry and structure edit functions.

## Installation for Development

```bash
# clone the repository
git clone http://10.72.201.45:13000/tmj/MolOP.git
# install the dependencies
cd MolOP
pip install -r requirements.txt
pip install .
```