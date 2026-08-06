# 读取计算结果

从最终 frame 中按科学性质查找结果。并非每种格式或计算任务都提供所有字段，始终先检查容器。

## 准备

```python
from molop import AutoParser

frame = AutoParser("calculation.log")[0][-1]
```

本页代码均打印返回值。使用共享 `water_mp2.out` 时，下面两个核心结果为：

```python
frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

??? example "输出"

    ```text
    3 (3, 3)
    -74.999374598107
    ```

## 结构、电荷与多重度

```python
print(frame.atoms)
print(frame.coords.shape)
print(frame.charge, frame.multiplicity)
```

??? example "输出"

    ```text
    [8, 1, 1]
    (3, 3)
    0 1
    ```

取得纯 NumPy 数组：

```python
coords_angstrom = frame.coords.m_as("angstrom")
print(coords_angstrom.shape)
```

??? example "输出"

    ```text
    (3, 3)
    ```

## 最终能量

```python
if frame.energies and frame.energies.total_energy is not None:
    energy_hartree = frame.energies.total_energy.m_as("hartree")
    print(energy_hartree)
```

??? example "输出"

    ```text
    -74.999374598107
    ```

`total_energy` 是 MolOP 从当前容器中选择的最高优先级可用总能量。需要区分参考能、MP2 或
CCSD(T) 时直接读取对应字段：

```python
energies = frame.energies
if energies:
    print("reference:", energies.reference_energy.m_as("hartree"))
    print("MP2:", energies.mp2_energy.m_as("hartree"))
    print("CCSD(T):", energies.ccsd_t_energy)
```

??? example "输出"

    ```text
    reference: -74.96357424008319
    MP2: -74.999374598
    CCSD(T): None
    ```

不要把 `total_energy` 当作跨方法可直接比较的统一理论水平。

## 热化学

```python
thermal = frame.thermal_informations
print("thermal available:", thermal is not None)
if thermal:
    print(thermal.ZPVE)
    print(thermal.H_T)
    print(thermal.G_T)
    print(thermal.S)
```

共享单点样例没有热化学区段：

??? example "输出"

    ```text
    thermal available: False
    ```

常用字段包括零点能 `ZPVE`、热校正 `TCE/TCH/TCG`、内能 `U_0/U_T`、焓 `H_T`、
吉布斯自由能 `G_T`、熵 `S` 和热容 `C_V`。单位随字段保存在 Pint quantity 中。

## 频率、虚频与振动模式

```python
vibrations = frame.vibrations
print("vibrations available:", vibrations is not None)
if vibrations:
    print(vibrations.frequencies.m_as("cm^-1"))
    print(vibrations.num_imaginary)
    if vibrations.vibration_modes:
        first_mode = vibrations.vibration_modes[0].m_as("angstrom")
```

共享单点样例输出：

??? example "输出"

    ```text
    vibrations available: False
    ```

有频率结果时，代码会继续输出频率数组、虚频整数和可选的 `(N, 3)` 位移数组；没有频率任务时
不会输出零频率。

`num_imaginary` 按负频率计数。过渡态判断还会结合 frame 的任务和状态信息，直接使用
`frame.is_TS` 更适合筛选。

## 分子轨道

```python
orbitals = frame.molecular_orbitals
print("molecular orbitals available:", orbitals is not None)
if orbitals:
    print(orbitals.HOMO_energy.m_as("eV"))
    print(orbitals.LUMO_energy.m_as("eV"))
    print(orbitals.HOMO_LUMO_gap.m_as("eV"))
    print(orbitals.alpha_occupancies)
```

共享水分子样例未提供结构化轨道容器：

??? example "输出"

    ```text
    molecular orbitals available: False
    ```

开壳层结果还可能提供 `beta_energies` 和 `beta_occupancies`。打印输出没有轨道系数时，MolOP
不会合成系数矩阵。

## 原子布居

```python
populations = frame.charge_spin_populations
if populations:
    print(populations.population_names)
    for name in populations.population_names:
        print(name, populations[name].values)
else:
    print("populations available: False")
```

共享 ORCA 样例的实际输出为：

??? example "输出"

    ```text
    ['mulliken_charges', 'lowdin_charges']
    mulliken_charges [-0.361722  0.180857  0.180866]
    lowdin_charges [-0.2507    0.125348  0.125352]
    ```

其他文件的名称和数值以实际打印内容为准。

布居方案是可扩展字典，名称取决于源文件实际打印内容，例如 `mulliken_charges`、
`mulliken_spins`、`lowdin_charges`、`hirshfeld_charges`、`cm5_charges`、`npa_charges` 或
`esp_charges`。不要假设某个固定名称一定存在。

## 偶极矩与极化率

```python
response = frame.polarizability
print("response available:", response is not None)
if response:
    if response.dipole is not None:
        print(response.dipole.m_as("debye"))
    if response.polarizability_tensor is not None:
        print(response.polarizability_tensor.m_as("bohr**3"))
```

共享样例只提供偶极矩：

??? example "输出"

    ```text
    response available: True
    [ 0.1502394  -0.11206605 -0.6497661 ]
    ```

根据格式和打印区段，还可能有 `isotropic_polarizability`、`anisotropic_polarizability`、
`quadrupole` 及更高多极矩。

## NMR

```python
nmr = frame.nmr
print("NMR available:", nmr is not None)
if nmr:
    for shielding in nmr.shielding_tensors:
        print(shielding.atom_index, shielding.isotropic.m_as("ppm"))

    if nmr.spin_spin_coupling_j is not None:
        coupling_hz = nmr.spin_spin_coupling_j.m_as("Hz")
```

共享水分子样例未提供 NMR 区段：

??? example "输出"

    ```text
    NMR available: False
    ```

fchk 可能只包含约化耦合 `spin_spin_coupling_k`，但缺少重建 J 所需的同位素信息。以格式页
的限制为准。

## 优化与终止状态

```python
print(frame.is_normal)
print(frame.is_error)
print(frame.is_optimized)
print(frame.is_TS)

if frame.status:
    print(frame.status.scf_converged)
if frame.geometry_optimization_status:
    print(frame.geometry_optimization_status.geometry_optimized)
```

共享单点样例输出：

??? example "输出"

    ```text
    True
    False
    False
    False
    ```

`None` 表示源文件没有足够证据，不应当强制解释为 `False`。

## 下一步

- [按科学性质查找字段](../reference/model_fields.md)
- [格式支持概览](../reference/format_support.md)
- [批量汇总](batch.md)
