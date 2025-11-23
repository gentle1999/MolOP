<!--
 * @Author: TMJ
 * @Date: 2024-01-30 16:53:29
 * @LastEditors: TMJ
 * @LastEditTime: 2024-10-11 10:29:38
 * @Description: 请填写简介
-->
# MolOP

**Molecule OPerator**，是[Open TS DataBase](http://10.72.201.58:13000/tmj/OTSDB-Core)项目的基本分子信息提取和操作单元。

**请注意，由于缺少人手维护中文文档，该部分内容是过时的，一旦冲突以英文文档为准。**

## 功能

- 自动从常用 QM 计算软件的输入和输出文件中提取分子信息。
    - 坐标文件
        - GJF `Done`
        - XYZ `Done`
        - SDF `Done`
    - QM 输出文件
        - G16 LOG `Done`
        - G16 IRC `Done`
        - G16 FCHK `Done`
        - xTB OUT `Done`
        - ORCA `TODO`

- 在[OpenBabel](https://openbabel.org/index.html)的基础工作上，从简单的原子坐标中提供一种[分子图恢复算法](structure_recovery.md)，可以方便地用于文件读取过程。该启发式分子图恢复算法的正确率已远超rdkit中提供的`rdDetermineBond`，特别是在处理自由基和金属络合物方面。请参见 [structure_recovery_cases.ipynb](structure_recovery_cases.ipynb) 了解更多差异。

- 提供分子几何和结构编辑功能。 `Doing`
    - 子结构替换 `Done`
    - 原子序号重置 `Done`
    - 分子朝向改变 `TODO`
    - 其他函数 `TODO`


## 引用
arXiv 论文即将提交。