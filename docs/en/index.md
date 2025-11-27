<!--
 * @Author: TMJ
 * @Date: 2024-01-30 16:53:29
 * @LastEditors: TMJ
 * @LastEditTime: 2024-06-28 18:51:20
 * @Description: 请填写简介
-->
# MolOP

**Molecule OPerator**, which is the basic molecule information extraction and operation library.

## Features

- Automatically extract molecule information from the Input and Output files of the common QM calculation softwares.
    - Coords file
        - GJF `Done`
        - XYZ `Done`
        - SDF `Done`
    - QM output file
        - G16 LOG `Done`
        - G16 IRC `Done`
        - G16 FCHK `Done`
        - xTB OUT `Done`
        - ORCA `TODO`

- Offer a [molecular graph recovery algorithm](structure_recovery.md) from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process. The accuracy of this heuristic molecular graph recovery algorithm has far surpassed that of `rdDetermineBond` provided in rdkit, especially for free radicals and metal complexes. See [structure_recovery_cases.ipynb](structure_recovery_cases.ipynb) to learn more about the difference).

- Offer the moleculer geometry and structure edit functions. `Doing`
    - Substructure replacement `Done`
    - Atom index reset `Done`
    - Orientation change `Done`
    - Other functions `TODO`


## Citation
arXiv paper is to be submitted soon.
