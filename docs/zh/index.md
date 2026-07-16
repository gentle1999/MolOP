# MolOP

MolOP 用一个 Python API 和一套命令行工具读取计算化学输出、提取结果、批量筛选并导出结构。

## 三步读取最终能量

安装当前版本：

```bash
pip install git+https://github.com/gentle1999/MolOP.git
```

下载文档使用的 [ORCA 水分子样例](../assets/examples/water_mp2.out)，将它保存为
`water_mp2.out`，然后运行：

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
frame = batch[0][-1]

print(frame.atoms)
print(frame.energies.total_energy.m_as("hartree"))
```

最后一行输出约为 `-74.999374598107`。`batch[0]` 是第一个文件，`[-1]` 是该文件的
最后一帧。

## 选择你的任务

| 要完成的任务 | 从这里开始 |
| --- | --- |
| 安装并验证环境 | [安装](getting_started/installation.md) |
| 用一个真实文件完成首次解析 | [5 分钟上手](getting_started/quickstart.md) |
| 在 Python 中读取能量、频率、布居或 NMR | [读取计算结果](guides/results.md) |
| 批量导出 CSV | [批量汇总](guides/batch.md) |
| 筛选正常结束、优化结果或过渡态 | [过滤与选择](guides/filtering.md) |
| 导出 XYZ、SDF、Gaussian 或 ORCA 输入 | [格式转换与导出](guides/conversion.md) |
| 查看某种文件能解析哪些数据 | [格式支持概览](reference/format_support.md) |

## 常见输入

MolOP 对 Gaussian log/fchk/input、ORCA output/input、xTB output、XYZ、SDF/MOL、SMILES
和 CML 等格式提供专用 reader 或 writer。具体字段和限制以
[格式支持概览](reference/format_support.md)及各格式页面为准。

MolOP 不运行 Gaussian、ORCA 或 xTB 计算。它处理已有文件，并把结果整理成可查询、可转换的
对象。

## 下一步

[安装 MolOP](getting_started/installation.md)，或直接进入
[5 分钟上手](getting_started/quickstart.md)。
