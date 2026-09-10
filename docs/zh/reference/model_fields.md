# 按科学性质查找字段

从 `frame = AutoParser(path)[0][-1]` 开始，按要提取的科学结果定位公共字段。

## 快速索引

| 结果 | 存在性检查 | 主要字段 | 常用单位 |
| --- | --- | --- | --- |
| 结构 | `bool(frame.atoms)` | `atoms`, `coords`, `rdmol` | angstrom |
| 电荷/多重度 | 直接读取 | `charge`, `multiplicity` | 无量纲 |
| 能量 | `frame.energies is not None` | `energies.total_energy`, `reference_energy`, `mp2_energy`, `ccsd_t_energy` | hartree |
| 热化学 | `frame.thermal_informations is not None` | `ZPVE`, `U_T`, `H_T`, `G_T`, `S`, `C_V` | 字段自带单位 |
| 频率 | `frame.vibrations is not None` | `frequencies`, `num_imaginary`, `vibration_modes` | cm^-1, angstrom |
| 分子轨道 | `frame.molecular_orbitals is not None` | `HOMO_energy`, `LUMO_energy`, `HOMO_LUMO_gap`, occupancies | hartree，可转 eV |
| 原子布居 | `frame.charge_spin_populations is not None` | `population_names`, `populations[name].values` | 无量纲 |
| 偶极/极化率 | `frame.polarizability is not None` | `dipole`, `polarizability_tensor`, `quadrupole` | debye, bohr^3 |
| NMR | `frame.nmr is not None` | `shielding_tensors`, `spin_spin_coupling_j/k` | ppm, Hz |
| 力/Hessian | `frame.forces/hessian is not None` | `forces`, `hessian` 及 axis metadata | hartree/bohr 等 |
| 状态 | 直接读取 | `is_normal`, `is_error`, `is_optimized`, `is_TS` | `bool | None` |

`rdmol` 是受保护的 RDKit 图副本：修改 `frame.rdmol` 返回对象不会修改 frame 的缓存。需要编辑
结构时应修改 `atoms`、`coords` 或显式拓扑字段；这些修改会使派生拓扑和 SMILES 缓存失效。

## 单位化数值

```python
energy = frame.energies.total_energy
print(energy.m_as("hartree"))

frequencies = frame.vibrations.frequencies.m_as("cm^-1")
coords = frame.coords.m_as("angstrom")
```

??? example "输出"

    随文档提供的 ORCA 水分子样例第一行结果为：

    ```text
    -74.999374598107
    ```

`.m_as(...)` 返回指定单位下的数值或 NumPy 数组。保留 quantity 本身可避免单位信息丢失。

## 布居是开放集合

```python
populations = frame.charge_spin_populations
if populations:
    for name, series in populations.population_items():
        print(name, series.scheme, series.quantity, series.values[:3])
```

??? example "输出"

    随文档提供的 ORCA 水分子样例包含两种电荷布居方案：

    ```text
    mulliken_charges mulliken charge [-0.361722  0.180857  0.180866]
    lowdin_charges lowdin charge [-0.2507    0.125348  0.125352]
    ```

`ChargeSpinPopulations` 只有 `populations` 映射，不为每种方案增加固定属性。读取 Mulliken 电荷：

```python
values = populations["mulliken_charges"].values
```

先检查 `population_names`，因为源文件可能只提供 Lowdin、Hirshfeld、CM5、NPA 或 ESP 等方案。

## 状态中的 `None`

```python
if frame.is_normal is None:
    print("源文件没有足够终止状态证据")
```

??? example "输出"

    ```text
    源文件没有足够终止状态证据
    ```

对状态字段，`None` 表示未知，不等于 `False`。

## 格式能力

字段存在不代表所有格式都能填充它。按格式核对：

- [Gaussian log](formats/g16log.md)
- [Gaussian fchk](formats/g16fchk.md)
- [ORCA output](formats/orcaout.md)
- [xTB output](formats/xtbout.md)
- [完整格式支持概览](format_support.md)

更多可运行代码见[读取计算结果](../guides/results.md)。
