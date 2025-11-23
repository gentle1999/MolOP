<!--
 * @Author: TMJ
 * @Date: 2024-02-16 16:20:29
 * @LastEditors: cathayana populuscathayana@gmail.com
 * @LastEditTime: 2024-02-20 14:02:18
 * @Description: 请填写简介
-->
# 结构恢复

目前MolOP中集成的启发式分子图恢复算法的正确率已远超rdkit中提供的`rdDetermineBond`，已经在molop中被移除。

`GraphReconstruction.py`中分子图恢复模块的架构可以概括如下：

```
xyz_to_rdmol (入口)
|
+--> xyz_to_separated_rwmol (处理金属)
|    |
|    +--> (如果存在金属) 遍历可能的金属价态/自由基
|    |    |
|    |    +--> xyz_to_separated_rwmol_no_metal (对于有机部分)
|    |         |
|    |         +--> xyz2omol (核心启发式算法)
|    |              |
|    |              +--> make_connections
|    |              +--> pre_clean
|    |              +--> fresh_omol_charge_radical
|    |              +--> 一系列的 `eliminate_*` 和 `clean_*` 函数
|    |              +--> get_radical_resonances
|    |              +--> process_resonance
|    |              +--> omol_score
|    |
|    +--> combine_metal_with_rwmol
|    +--> structure_score
|
+--> (如果没有金属) 直接调用 xyz_to_separated_rwmol_no_metal
|
+--> make_dative_bonds (可选)
```

::: molop.structure.GraphReconstruction