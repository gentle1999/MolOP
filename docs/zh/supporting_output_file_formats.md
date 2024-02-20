<!--
 * @Author: TMJ
 * @Date: 2024-02-02 10:43:56
 * @LastEditors: cathayana populuscathayana@gmail.com
 * @LastEditTime: 2024-02-20 13:59:56
 * @Description: 请填写简介
-->
# 支持的输出格式
MolOP 可以轻松地将输入文件转换成不同的格式。下面列出了支持的输出文件格式。

假设我们有一些 g16log 文件:
```python
from molop import AutoParser
files = AutoParser("path/to/*.log")
```

## XYZ (.xyz)

### 批量转换

```python
files.to_XYZ_files(file_path="path/to/out")
```

可以通过设置 `file_path` 参数来确定文件路径。如果未设置 `file_path` 参数，默认文件路径将是工作目录。那么输出文件将保存为 `file_path/<filename>.xyz` 。

### 单文件转换

```python
files[i].to_XYZ_block() # return an XYZ block with multi frames
files[i].to_XYZ_file(file_path="path/to/out.xyz") # save the XYZ block with multi frames to file
```
`file_path` 应该是一个文件，而不是一个目录。

### 单帧转换

```python
files[i][j].to_XYZ_block() # return an XYZ block with single frame
files[i][j].to_XYZ_file(file_path="path/to/out.xyz") # save the XYZ block with single frame to file
```
`file_path` 应该是一个文件，而不是一个目录。


## SDF (.sdf, .mol)

### 批量转换

```python
files.to_SDF_files(file_path="path/to/out")
```

可以通过设置 `file_path` 参数来确定文件路径。如果未设置 `file_path` 参数，默认文件路径将是工作目录。那么输出文件将保存为 `file_path/<filename>.sdf` 。

### 单文件转换

```python
files[i].to_SDF_block() # return a SDf block with multi frames
files[i].to_SDF_file(file_path="path/to/out.sdf") # save the SDF block with multi frames to file
```
`file_path` 应该是一个文件，而不是一个目录。

### 单帧转换

```python
files[i][j].to_SDF_block() # return a SDF block with single frame
files[i][j].to_SDF_file(file_path="path/to/out.sdf") # save the SDF block with single frame to file
```
`file_path` 应该是一个文件，而不是一个目录。

## Gaussian 16 Input File (.com, .gjf, .gau, .gjc)

在优化结构的基础上继续进行更高层次的优化是量子化学计算工作流程中非常常见的一步。MolOP 提供了一些简单的接口来简化这一步骤所需的人力。具体操作详见[Transform in GJF](transform_in_gjf.md)部分。

## Chemdraw (.cdxml)

### 批量转换

```python
files.to_chemdraw(file_path="path/to/out", frameID=-1, keep3D=True)
```
可以通过设置 `file_path` 参数来确定文件路径。如果未设置 `file_path` 参数，默认文件路径将是工作目录。那么输出文件将保存为 `file_path/<filename>.cdxml` 。

如果未设置 `frameID` 参数，默认帧 ID 将是最后一帧。`keep3D` 参数用于确定是否保留结构的 3D 信息。如果未设置 `keep3D` 参数，默认值为 `True`。如果将 `keep3D` 参数设置为 `False`，结构体的 3D 信息将被删除。

### 单文件转换

```python
files[i].to_chemdraw(file_path="path/to/out.cdxml", frameID=-1, keep3D=True) # save the specific frame to cdxml file
```
`file_path` 应该是一个文件，而不是一个目录。

### 单帧转换

```python
files[i][j].to_chemdraw(file_path="path/to/out.cdxml", keep3D=True) # save the specific frame to cdxml file
```
`file_path` 应该是一个文件，而不是一个目录。