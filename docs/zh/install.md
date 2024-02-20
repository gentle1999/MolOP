<!--
 * @Author: TMJ
 * @Date: 2024-02-03 16:23:14
 * @LastEditors: cathayana populuscathayana@gmail.com
 * @LastEditTime: 2024-02-20 12:45:33
 * @Description: 请填写简介
-->
# 安装

该软件包正在积极开发中，现在发布到 pypi 上还为时过早。因此，您可以通过以下方式安装。

## 浙大内网
您可以使用我们的自托管存储库。
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop --upgrade
```

要进行额外的描述符计算，您需要安装以下环境：
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install --extra-index-url http://10.72.201.58:13000/api/packages/tmj/pypi/simple/ --trusted-host 10.72.201.58 molop[full] --upgrade
```

## 互联网
您可以使用 github 仓库。
```bash
conda install openbabel -c conda-forge # openbabel is a necessary dependence
pip install poetry
poetry add git+https://github.com/gentle1999/MolOP.git
```