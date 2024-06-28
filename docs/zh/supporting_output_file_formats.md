<!--
 * @Author: TMJ
 * @Date: 2024-02-02 10:43:56
 * @LastEditors: TMJ
 * @LastEditTime: 2024-03-27 15:13:13
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
files.to_XYZ_file(file_path="path/to/out")
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

### Transform in batch

```python
files.to_SDF_file(file_path="path/to/out")
```
You can determine the file path by setting the `file_path` parameter. If the `file_path` parameter is not set, the default file path will be the working directory. Then the output files will be saved as `file_path/<filename>.sdf`

### Transform in single file

```python
files[i].to_SDF_block() # return a SDf block with multi frames
files[i].to_SDF_file(file_path="path/to/out.sdf") # save the SDF block with multi frames to file
```
`file_path` there should be a file, instead of a directory.

### Transform in single frame

```python
files[i][j].to_SDF_block() # return a SDF block with single frame
files[i][j].to_SDF_file(file_path="path/to/out.sdf") # save the SDF block with single frame to file
```
`file_path` there should be a file, instead of a directory.

## Gaussian 16 Input File (.com, .gjf, .gau, .gjc)

Continuing to a higher level of optimization based on the optimized structure is a very common step in quantum chemistry computational workflows. MolOP provides some simple interfaces to simplify the manpower required for this step. The specific operations are explained in detail in the section [Transform in GJF](transform_in_gjf.md).

## Chemdraw (.cdxml)

### Transform in batch

```python
files.to_chemdraw(file_path="path/to/out", frameID=-1, keep3D=True)
```
You can determine the file path by setting the `file_path` parameter. If the `file_path` parameter is not set, the default file path will be the working directory. Then the output files will be saved as `file_path/<filename>.cdxml`

You can determine the frameID by setting the `frameID` parameter. If the `frameID` parameter is not set, the default frameID will be the last frame. The `keep3D` parameter is used to determine whether to keep the 3D information of the structure. If the `keep3D` parameter is not set, the default value is `True`. If the `keep3D` parameter is set to `False`, the 3D information of the structure will be removed.

### Transform in single file

```python
files[i].to_chemdraw(file_path="path/to/out.cdxml", frameID=-1, keep3D=True) # save the specific frame to cdxml file
```
`file_path` there should be a file, instead of a directory.

### Transform in single frame

```python
files[i][j].to_chemdraw(file_path="path/to/out.cdxml", keep3D=True) # save the specific frame to cdxml file
```
`file_path` there should be a file, instead of a directory.