# 模型字段

用于存储量子化学计算的输出，如 Gaussian 的 `.log` 文件。包含更丰富的物理化学属性。

## 文件级 (`BaseCalcFile`)

- `qm_software`: (str) 量子化学计算软件名称 (e.g., "Gaussian")。
- `qm_software_version`: (str) 软件版本号。
- `keywords`: (str) 计算任务的关键词/路由部分。
- `method`: (str) 计算方法 (e.g., "DFT", "HF")。
- `basis_set`: (str) 基组。
- `functional`: (str) 泛函。
- `charge`: (int) 分子总电荷。
- `multiplicity`: (int) 分子总多重度。
- `solvent`: (`ImplicitSolvation`) 隐式溶剂化模型信息。
  - `solvent`: (str) 溶剂名称 (e.g., "Water")。
  - `solvent_model`: (str) 溶剂模型名称 (e.g., "PCM")。
  - `atomic_radii`: (str) 溶剂模型中使用的原子半径类型。
  - `solvent_epsilon`: (float) 溶剂介电常数。
  - `solvent_epsilon_infinite`: (float) 溶剂无限介电常数。
- `temperature`: (PlainQuantity) 计算温度，单位为 `K`。
- `pressure`: (PlainQuantity) 计算压力，单位为 `atm`。
- `running_time`: (PlainQuantity) 总计算耗时，单位为 `second`。
- `status`: (`Status`) 最后一帧的计算状态。

## 帧级 (`BaseCalcFrame`)

除了包含上述 `BaseCoordsFrame` 的所有结构信息外，每一帧还可能包含以下计算属性：

- **`energies`: `Energies`** - 能量信息
  - `electronic_energy`: (PlainQuantity) 电子能量, 单位 `hartree`。
  - `scf_energy`: (PlainQuantity) SCF 能量, 单位 `hartree`。
  - `mp2_energy`: (PlainQuantity) MP2 能量, 单位 `hartree`。
  - `mp3_energy`: (PlainQuantity) MP3 能量, 单位 `hartree`。
  - `mp4_energy`: (PlainQuantity) MP4 能量, 单位 `hartree`。
  - `ccsd_energy`: (PlainQuantity) CCSD 能量, 单位 `hartree`。
  - `total_energy`: (PlainQuantity) 根据精度自动选择的总能量, 单位 `hartree`。
- **`thermal_informations`: `ThermalInformations`** - 热力学信息
  - `ZPVE`: (PlainQuantity) 零点振动能, 单位 `kcal/mol`。
  - `TCE`: (PlainQuantity) 内能热校正, 单位 `kcal/mol`。
  - `TCH`: (PlainQuantity) 焓热校正, 单位 `kcal/mol`。
  - `TCG`: (PlainQuantity) 吉布斯自由能热校正, 单位 `kcal/mol`。
  - `U_0`: (PlainQuantity) 零点能量, 单位 `kcal/mol`。
  - `U_T`: (PlainQuantity) 内能, 单位 `kcal/mol`。
  - `H_T`: (PlainQuantity) 焓, 单位 `kcal/mol`。
  - `G_T`: (PlainQuantity) 吉布斯自由能, 单位 `kcal/mol`。
  - `S`: (PlainQuantity) 熵, 单位 `cal/mol/K`。
  - `C_V`: (PlainQuantity) 定容热容, 单位 `cal/mol/K`。
- **`molecular_orbitals`: `MolecularOrbitals`** - 分子轨道信息
  - `electronic_state`: (str) 电子态。
  - `alpha_energies` / `beta_energies`: (NumpyQuantity) Alpha 和 Beta 轨道能量, 单位 `hartree`。
  - `alpha_occupancies` / `beta_occupancies`: (list[bool]) Alpha 和 Beta 轨道占据情况。
  - `alpha_symmetries` / `beta_symmetries`: (list[str]) Alpha 和 Beta 轨道对称性。
  - `HOMO_energy` / `LUMO_energy`: (PlainQuantity) 前线轨道能量。
  - `HOMO_LUMO_gap`: (PlainQuantity) 前线轨道能级差。
- **`vibrations`: `Vibrations`** - 振动分析信息
  - `frequencies`: (NumpyQuantity) 振动频率, 单位 `cm^-1`。
  - `reduced_masses`: (NumpyQuantity) 约化质量, 单位 `amu`。
  - `force_constants`: (NumpyQuantity) 力常数, 单位 `mdyne/angstrom`。
  - `IR_intensities`: (NumpyQuantity) 红外强度, 单位 `km/mol`。
  - `vibration_modes`: (list[NumpyQuantity]) 振动模式向量。
  - `num_imaginary`: (int) 虚频数量。
- **`charge_spin_populations`: `ChargeSpinPopulations`** - 布居分析
  - `mulliken_charges` / `mulliken_spins`: (list[float]) Mulliken 电荷和自旋。
  - `apt_charges`: (list[float]) 原子极化张量电荷。
  - `lowdin_charges`: (list[float]) Lowdin 电荷。
  - `hirshfeld_charges` / `hirshfeld_spins`: (list[float]) Hirshfeld 电荷和自旋。
  - `hirshfeld_q_cm5`: (list[float]) Hirshfeld CM5 电荷。
  - `npa_charges`: (list[float]) NPA 电荷。
- **`polarizability`: `Polarizability`** - 极化率和多极矩
  - `electronic_spatial_extent`: (PlainQuantity) 电子空间范围, 单位 `bohr^2`。
  - `isotropic_polarizability`: (PlainQuantity) 各向同性极化率, 单位 `bohr^3`。
  - `anisotropic_polarizability`: (PlainQuantity) 各向异性极化率, 单位 `bohr^3`。
  - `polarizability_tensor`: (PlainQuantity) 极化率张量。
  - `dipole`: (PlainQuantity) 偶极矩, 单位 `debye`。
  - `quadrupole`: (PlainQuantity) 四极矩, 单位 `debye*angstrom`。
  - `octapole`: (PlainQuantity) 八极矩, 单位 `debye*angstrom^2`。
  - `hexadecapole`: (PlainQuantity) 十六极矩, 单位 `debye*angstrom^3`。
- **`bond_orders`: `BondOrders`** - 键级信息
  - `wiberg_bond_order`: (np.ndarray) Wiberg 键级。
  - `mo_bond_order`: (np.ndarray) MO 键级。
  - `mayer_bond_order`: (np.ndarray) Mayer 键级。
- **`total_spin`: `TotalSpin`** - 总自旋信息
  - `spin_square`: (float) S\*\*2 的值。
  - `spin_quantum_number`: (float) 自旋量子数 S。
- **`single_point_properties`: `SinglePointProperties`** - 单点性质
  - `vip`: (PlainQuantity) 垂直电离势, 单位 `eV/particle`。
  - `vea`: (PlainQuantity) 垂直电子亲和能, 单位 `eV/particle`。
  - `gei`: (PlainQuantity) 全局亲电指数, 单位 `eV/particle`。
- **`geometry_optimization_status`: `GeometryOptimizationStatus`** - 几何优化状态
  - `geometry_optimized`: (bool) 几何结构是否收敛。
  - `max_force` / `rms_force`: (float) 最大和均方根力。
  - `max_displacement` / `rms_displacement`: (float) 最大和均方根位移。
  - `energy_change`: (float) 能量变化。
- **`status`: `Status`** - 当前帧的计算状态
  - `scf_converged`: (bool) SCF 是否收敛。
  - `normal_terminated`: (bool) 计算是否正常终止。
- **其他直接字段**
  - `forces`: (NumpyQuantity) 原子受力, 单位 `hartree/bohr`。
  - `hessian`: (NumpyQuantity) Hessian 矩阵, 单位 `hartree/bohr^2`。
  - `rotation_constants`: (NumpyQuantity) 转动常数, 单位 `gigahertz`。
