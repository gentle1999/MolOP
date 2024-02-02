<!--
 * @Author: TMJ
 * @Date: 2024-02-01 20:58:56
 * @LastEditors: TMJ
 * @LastEditTime: 2024-02-02 19:37:48
 * @Description: 请填写简介
-->
# Read files

MolOP offer an easy way to read files.

```python
from molop import AutoParser
files = AutoParser("path/to/file")
```

The `path/to/file` is a wildcard of the file path. For example, you can use `path/to/*.log` to read all g16log files in the directory `path/to/`.

```python
from molop import AutoParser
files = AutoParser("path/to/*.log")
```

The file reading process can be automatically parallelized if you pass files more than your CPU cores. 

## FileParserBatch

The `files` is a `FileParserBatch` object, which contains a squence of `FileParser` objects. So you can access each file by index like `files[i]`.

## FileParser

MolOP extract all frames in a file and store them in a `FileParser` object contains a squence of `BlockParser` objects. So you can access each frame by index like `files[i][j]`.

## BlockParser

Each frame is a certain molcule with basic information like `coordinates`, `atoms` et. al. (For QM files, you can also get other information like `energy`, `orbital` et. al.)
