<!--
 * @Author: TMJ
 * @Date: 2024-02-16 16:20:29
 * @LastEditors: TMJ
 * @LastEditTime: 2024-10-10 18:28:26
 * @Description: 请填写简介
-->
# Structure Recovery

The accuracy of the heuristic molecular graph recovery algorithm has far surpassed that of `rdDetermineBond` provided in rdkit, and it has been removed from molop.

The architecture of the graph recovery module in `GraphReconstruction.py` can be summarized as follows:

```
xyz_to_rdmol (entry point)
|
+--> xyz_to_separated_rwmol (handles metals)
|    |
|    +--> (if metals are present) iterates through possible metal valences/radicals
|    |    |
|    |    +--> xyz_to_separated_rwmol_no_metal (for the organic part)
|    |         |
|    |         +--> xyz2omol (core heuristic algorithm)
|    |              |
|    |              +--> make_connections
|    |              +--> pre_clean
|    |              +--> fresh_omol_charge_radical
|    |              +--> a series of `eliminate_*` and `clean_*` functions
|    |              +--> get_radical_resonances
|    |              +--> process_resonance
|    |              +--> omol_score
|    |
|    +--> combine_metal_with_rwmol
|    +--> structure_score
|
+--> (if no metals) xyz_to_separated_rwmol_no_metal is called directly
|
+--> make_dative_bonds (optional)
```

::: molop.structure.GraphReconstruction
