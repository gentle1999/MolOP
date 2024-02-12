<!--
 * @Author: TMJ
 * @Date: 2024-01-30 16:53:29
 * @LastEditors: TMJ
 * @LastEditTime: 2024-02-12 21:40:53
 * @Description: 请填写简介
-->
# MolOP

**Molecule OPerator**, which is the basic molecule information extraction and operation unit of the [Open TS DataBase](http://10.72.201.58:13000/tmj/OTSDB-Core) project.

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

- Offer a [molecular graph recovery algorithm](https://github.com/gentle1999/MolOP/blob/main/molop/structure/structure_recovery.py) from the simple coodinates of atoms based on the initial work by [OpenBabel](https://openbabel.org/index.html), which can be easily used in the file reading process. This algorithm is different from the rdDetermineBonds (Original code implemented by [Jensen group](https://github.com/jensengroup/xyz2mol) and integrated in RDKit from the 2022.09 release, which is not suitable for the free radicals and complex containing metal. See [rdDetermineBonds_cannot_do](https://github.com/gentle1999/MolOP/blob/main/tutorial/rdDetermineBonds_cannot_do.ipynb) to learn more about the difference).
  
  > Although our algorithm overcome the free radicals and metal problem and tested on the [test_cases](https://github.com/gentle1999/MolOP/blob/main/tutorial/test_cases.ipynb) file, it is still not perfect. There is no denying that, rdDetermineBonds works well for normal organic molecules. Thus, we would give molecule structure recovered by rdDetermineBonds first, if error happens, we will use our algorithm to recover the molecule structure instead. We hope that this strategy can take advantage of both approaches.

- Offer the moleculer geometry and structure edit functions. `Doing`
    - Substructure replacement `Done`
    - Atom index reset `Done`
    - Orientation change `TODO`
    - Other functions `TODO`


## Citation
arXiv paper is to be submitted soon.