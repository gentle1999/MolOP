<!--
 * @Author: TMJ
 * @Date: 2024-02-03 16:23:14
 * @LastEditors: TMJ
 * @LastEditTime: 2024-02-12 21:35:30
 * @Description: 请填写简介
-->
# Installation

This package is actively developing, and this time is too early to be published to pypi. So you can install by the following.

## In ZJU intranet
You can use our self-host repository.
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop --upgrade
```

For additional descriptor calculation, you need to install the requirements below:
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop[full] --upgrade
```

## In Internet
You can use the github repository.
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install poetry
poetry add git+https://github.com/gentle1999/MolOP.git
```