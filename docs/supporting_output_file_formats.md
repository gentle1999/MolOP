<!--
 * @Author: TMJ
 * @Date: 2024-02-02 10:43:56
 * @LastEditors: TMJ
 * @LastEditTime: 2024-02-02 19:51:17
 * @Description: 请填写简介
-->
# Supporting Input File Formats
MolOP can easily transform input files into different formats. The following lists the supported output file formats.

Suppose we have read some g16log files:
```python
from molop import AutoParser
files = AutoParser("path/to/*.log")
```

## XYZ (.xyz) 

### Transform in batch

```python
files.to_XYZ_files(file_path="path/to/out")
```
You can determine the file path by setting the `file_path` parameter. If the `file_path` parameter is not set, the default file path will be the working directory. Then the output files will be saved as `file_path/<filename>.xyz`

### Transform in single file

```python
files[i].to_XYZ_block() # return an XYZ block with multi frames
files[i].to_XYZ_file(file_path="path/to/out.xyz") # save the XYZ block with multi frames to file
```
`file_path` there should be a file, instead of a directory.

### Transform in single frame

```python
files[i][j].to_XYZ_block() # return an XYZ block with single frame
files[i][j].to_XYZ_file(file_path="path/to/out.xyz") # save the XYZ block with single frame to file
```
`file_path` there should be a file, instead of a directory.


## SDF (.sdf, .mol)
The SDF file format (or MOL file) is a standard file format used to store chemical structures. It is a text file format that contains information about the chemical structure, such as the coordinates of atoms, the bond information, and some other chemical properties like charge and spin multiplicity. In general, we can completely recover the chemical structure from a SDF file.

There is an example of the SDF file format:
```text
charge 0 multiplicity 2
     RDKit          3D

  7  6  0  0  0  0  0  0  0  0999 V2000
   -0.6934    0.0000   -0.0019 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7936    0.0000   -0.0185 C   0  0  0  0  0  3  0  0  0  0  0  0
   -1.1062    0.8866   -0.4927 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1062   -0.8866   -0.4927 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0915   -0.0000    1.0262 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3512    0.9264    0.0407 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.3512   -0.9264    0.0407 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0
  2  1  1  0
  2  6  1  0
  2  7  1  0
  3  1  1  0
  4  1  1  0
M  RAD  1   2   2
M  END
```

## Gaussian 16 Input File (.com, .gjf, .gau, .gjc)
The Gaussian 16 input file format is a text file format used to store the input information for Gaussian 16. It contains the information about the chemical structure, such as the coordinates of atoms, the basis set, the charge and spin multiplicity, and some other chemical properties like the level of theory, the solvent, and the solvation model.

There is an example of the Gaussian 16 input file format:
```text
!Put Keywords Here, check Charge and Multiplicity.
#

 Title

0  1
C           4.40549        -0.58591        -0.00000
O           3.12691        -1.20031         0.00000
C           2.09747        -0.36319         0.00000
N           2.20839         0.94663        -0.00000
O           0.95902         1.45815        -0.00000
C           0.06922         0.45397         0.00000
C          -1.37500         0.77331         0.00000
O          -2.09755        -0.34776         0.00000
C          -3.51986        -0.17267        -0.00000
C          -4.17088        -1.54205        -0.00000
O          -1.83552         1.88339         0.00000
C           0.72761        -0.74674         0.00000
H           5.12845        -1.39862        -0.00000
H           4.53028         0.03818         0.88706
H           4.53027         0.03819        -0.88706
H          -3.80165         0.40589        -0.88655
H          -3.80165         0.40589         0.88655
H          -5.25106        -1.43016        -0.00000
H          -3.87110        -2.09999         0.88278
H          -3.87110        -2.09999        -0.88278
H           0.32154        -1.73407         0.00000

1 2 1.0 13 1.0 14 1.0 15 1.0 
2 3 1.0 
3 4 2.0 12 1.0 
4 5 1.0 
5 6 1.0 
6 7 1.0 12 2.0 
7 8 1.0 11 2.0 
8 9 1.0 
9 10 1.0 16 1.0 17 1.0 
10 18 1.0 19 1.0 20 1.0 
11 
12 21 1.0 
13 
14 
15 
16 
17 
18 
19 
20 
21 

```

