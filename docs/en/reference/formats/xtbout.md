# xTB Output

<!-- format-support:xtbout -->

| Item | Value |
| ---- | ----- |
| Format ID | `xtbout` |
| Extensions | `.out`, `.log`, `.xtbout` |
| Read | Yes |
| Write | No |
| Registry role | Reader |
| Data level | QM results; coordinates when printed in the output |

## Quick read

```python
from molop import AutoParser

frame = AutoParser("xtb.out", parser_detection="xtbout")[0][-1]
print(frame.qm_software_version, frame.is_normal)
print(frame.energies.total_energy.m_as("hartree"))
print(len(frame.atoms), frame.coords.shape)
```

??? example "Output contract"

    Output reports version, termination status, total energy, and geometry size. If an xTB
    single-point stdout only references an external coordinate file without printing geometry, the
    final line can be `0 (0, 3)`.

MolOP does not implicitly read adjacent coordinate files.

MolOP parses command-line standard output from xTB major version 5 and versions
`>=6`. The reader covers both the legacy xTB 5/6.1 print family and the modern
xTB `>=6` print
family. It does not define a separate xTB input format because xTB calculations
normally receive coordinates and options directly from the command line.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Version scope, setup, tasks, and status -->Versions, setup, tasks and status | Partial | Recognizes xTB major version 5 and treats versions `>=6` as the modern print family; parses version, program call, coordinate-file name, method, charge, multiplicity, OMP threads, task requests, termination, SCC status, and total wall time. | Versions before 5 are rejected. xTB 5 support is partial; compatibility with later major versions is unconfirmed. |
| <!-- feature-area:Geometry availability -->Available geometry | Partial | Parses legacy Bohr `$coord` blocks, legacy setup coordinates, modern XYZ final structures, and embedded V2000 SDF final structures. | Single-point xTB stdout often references an external coordinate file without printing coordinates. Such logs retain a property frame with empty geometry; the parser does not read adjacent files implicitly. |
| <!-- feature-area:Energies and geometry optimization -->Energies and geometry optimization | Partial | Exposes final total energy in Hartree, source-labeled energy evidence, optimization convergence, energy-change criteria, and xTB gradient norm fields. | Energy decomposition rows and full optimization trajectories are not structured; only the final geometry printed by stdout is represented. |
| <!-- feature-area:Orbitals, populations, and dipole -->Orbitals, populations, dipole and bond orders | Partial | Parses legacy and modern orbital energies/occupancies, GFN1 Mulliken/CM5 charges, GFN2 Mulliken charges, molecular dipole, and rotational constants where printed. | Thresholded Wiberg matrices are supported; unprinted pairs stay unknown. Quadrupoles, dispersion-property tables and omitted orbital ranges remain incomplete. |
| <!-- feature-area:Vibrations and thermochemistry -->Vibrations and thermochemistry | Partial | Parses projected physical frequencies, reduced masses, IR intensities, zero-point energy, total enthalpy, total free energy, heat capacity, entropy, temperature, and molecular mass where printed. | Explicit Å⁴/amu Raman activity is supported with separate frequency/activity method labels. Mode vectors, Hessian sidecars, depolarization ratios and undefined old (amu) columns remain unmapped. |
| <!-- feature-area:xTB analysis properties -->xTB analysis properties | Partial | Parses vertical IP, vertical EA, global electrophilicity index, and atom-resolved Fukui indices when requested and printed. | Additional xTB workflows such as FOD, metadynamics, MD, ONIOM, GFN-FF topology, and JSON/sidecar outputs are outside the current stdout contract. |

## xTB 5 / xTB >=6 Capability Coverage Matrix

The table below defines the field-level contract of the stdout reader. Statuses
have the following meanings:

- **Supported**: reads the listed output forms for that version family.
- **Partial**: only known print forms, final values, or a limited subset of task semantics are covered.
- **Unverified**: the parser path may apply, but no real output sample from that version family establishes a stable contract.
- **Unsupported**: no corresponding structured field is currently produced.
- **Not applicable**: the capability is outside the responsibility of the stdout reader.

“Supported” in the `xTB >=6 modern` column describes xTB 6 output layouts.
Later major versions use the same reader, but their compatibility is unconfirmed.

| xTB output capability | xTB 5 legacy | xTB >=6 modern | Structured result | Current boundary |
| --------------------- | ------------ | -------------- | ----------------- | ---------------- |
| Banner and version | Partial | Supported | `qm_software`, `qm_software_version` | Major version `5` and versions `>=6` are accepted. xTB 5 support covers selected legacy layouts; compatibility with later major versions is unconfirmed. |
| Concatenated command-line runs | Unverified | Supported | One segment and result frame per version banner | Every run must print a version banner; arbitrary shell separator text is not used for splitting. |
| Program call and coordinate-file name | Partial | Supported | `keywords`, `input_file_name`, and model-chemistry raw keywords | Only the printed `program call` and `coordinate file` lines are read; shell environment and the relative-path base are not recovered. |
| GFN Hamiltonian | Partial | Supported | `method`, `model_chemistry.method`, with method family `SEMIEMPIRICAL` | Printed GFN, GFN1, GFN2, and GFN-FF-style labels are recognized; method/version compatibility is not validated. |
| Charge, multiplicity, and OMP threads | Partial | Supported | `charge`, `multiplicity`, `request_num_cpu`, `resource_request` | Multiplicity is derived from `--uhf` or the printed unpaired-electron count plus one; common xTB defaults apply when values are absent. |
| SP, optimization, frequency, and gradient requests | Partial | Supported | `task_requests` with derivative order and requested properties | Inference covers `--opt`, `--ohess`, `--hess`, `--grad`, `--vip`, `--vea`, `--vipea`, `--fukui`, and `--vfukui`; not every xTB CLI combination is modeled. |
| Implicit solvent and temperatures | Unverified | Partial | `solvent`, `temperature`, `electron_temperature` | Only printed model, solvent, vibrational/solvent temperature, and electronic temperature values are parsed; complete ALPB/GBSA parameters are not reconstructed. |
| Legacy setup coordinates and `$coord` final structure | Partial | Partial | `atoms`, Angstrom `coords`, coordinate precision, and observed provenance | Bohr input is converted to Angstrom; original units, whitespace, and line text are not retained. |
| Modern XYZ final structure | Unverified | Supported | `atoms`, Angstrom `coords`, coordinate precision, and observed provenance | Only the final XYZ after `final structure:` is read from stdout; intermediate optimization geometries are not retained. |
| Embedded V2000 SDF final structure | Unverified | Supported | Atoms and coordinates | Only the atom block is used; bonds, formal charges, and SDF property blocks are not reconstructed by this reader. |
| External coordinates and sidecars such as `xtbopt.xyz` | Unsupported | Unsupported | Only the printed coordinate-file name is retained | Adjacent files are not read implicitly, preserving a single-source contract for stdout spans and hashes. |
| Final total energy | Partial | Supported | Hartree `energies.total_energy`; an `EnergyObservation` when source evidence is enabled | The last recognized final-energy print form is selected; individual SCC iterations are not represented as an energy trajectory. |
| SCC and energy-decomposition tables | Unsupported | Unsupported | No dedicated decomposition fields | Isotropic electrostatic, dispersion, repulsion, and similar components are not structured. |
| SCC and termination status | Partial | Supported | `status.scf_converged`, `status.normal_terminated` | SCC success may be inferred from a convergence line or final energy; normal termination depends on the `finished run` marker. |
| CPU and wall time | Partial | Supported | Independent second-valued `cpu_time` / `wall_time`; `running_time` aliases wall time | File totals aggregate selected segments separately. Printed total timing takes priority; a missing kind in a total block stays unknown. Without total timing, use the last matching row independently. Detailed module timings and scheduler time are not exposed. |
| Optimization convergence and thresholds | Partial | Supported | `geometry_optimization_status`, energy change/threshold, gradient norm/threshold | Only final printed values are retained; per-step convergence history is not structured. |
| Complete optimization trajectory | Unsupported | Unsupported | No per-step frame sequence | A command-line run currently maps to one result frame rather than one frame per optimization step. |
| Gradient vectors | Unsupported | Unsupported | Only the `gradient` task request and final gradient norm can be recorded | Per-atom stdout gradients and the `gradient` sidecar are not mapped to structured arrays. |
| Orbital energies and occupancies | Partial | Supported | Orbital energies, `total_occupancies` and available spin occupations | Reads legacy `occ./eps` and modern orbital tables. Omitted ranges cannot be recovered; total occupations split into alpha/beta channels only when spin information supports the assignment. |
| Atomic charges | Unverified | Partial | GFN1 Mulliken/CM5 and GFN2 Mulliken charges with source-labeled extensible population metadata | Spin populations, complete population tables, and every Hamiltonian-specific charge scheme are not parsed. |
| Dipole and rotational constants | Unverified | Partial | Debye dipole vector and GHz rotational constants | Dipole magnitude, quadrupole, polarizability tensors, and higher multipoles are not structured. |
| Frequencies, reduced masses, and IR intensities | Unverified | Supported | `vibrations`, including imaginary-mode count | Only printed projected physical frequencies are covered; missing columns are not synthesized. |
| Mode displacements and Hessian | Unsupported | Unsupported | No corresponding structured arrays | Displacement blocks and Hessian sidecars are not read. |
| Raman activity | Unverified | Partial | `vibrations.raman_activities` in Å⁴/amu | Reads the xTB 6.7.1 PTB activity table; no depolarization ratios are invented. |
| Wiberg bond orders | Unverified | Partial | `bond_orders["wiberg_bond_order"]` | Retains AO basis, threshold, reciprocal-rounding policy and unknown pairs. |
| Thermochemistry | Unverified | Supported | ZPVE, enthalpy, free energy, `C_V`, entropy, temperature, and molecular mass | Only summary values are retained; rotor/interpolation details and individual thermal-correction components are not structured. |
| VIP, VEA, GEI, and Fukui indices | Unverified | Supported | `single_point_properties` scalars and atom-resolved Fukui arrays | Values exist only when the corresponding analysis is requested and printed to stdout. |
| FOD, MD, metadynamics, ONIOM, and GFN-FF topology | Unsupported | Unsupported | No dedicated fields | These analyses and workflows are outside the current stdout contract. |
| JSON, Molden, Hessian, and other auxiliary output | Unsupported | Unsupported | No sidecar merge | The parser consumes only the supplied stdout text and does not scan the working directory. |
| xTB input and output writers | Not applicable | Not applicable | Reader only | MolOP does not define an xTB input format or reconstruct command lines or stdout from result models. |

## Coordinate Boundary

xTB commonly prints only the coordinate file path for single-point jobs. MolOP
does not silently read that adjacent file: doing so would mix a second source
artifact into a parser result whose source spans and hashes describe only the
stdout file. Optimization outputs that print a final structure do expose normal
coordinate fields.

## Version Boundary

The parser accepts major version `5` and routes all versions `>=6` through the
modern parsing contract. Recognized xTB outputs with a major version below `5`
raise an explicit unsupported-version error instead of falling through to
Gaussian or ORCA output readers.

PTB Raman results retain `frequency_method` and `raman_activity_method`,
including in individual modes and spectra. The stdout reader preserves its
printed precision, without silently replacing values from `vibspectrum` or
Gaussian-like `g98.out` files. See [shared data model](../data_model.md).

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [xTB CLI manual](https://github.com/grimme-lab/xtb/blob/v6.7.1/man/xtb.1.adoc); [properties](https://xtb-docs.readthedocs.io/en/latest/properties.html); [PTB](https://xtb-docs.readthedocs.io/en/latest/ptb.html).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Version scope, setup, tasks, and status -->Versions, setup, tasks and status | `Version banners`; `program call`; `Hamiltonian`; `charge/spin`; `OMP`; `solvent`; `T(el)`; `normal/converged/failure labels`; `total/component wall-time` | Map software/method/task/resource/solvation requests and conditions. Multiplicity derives from unpaired-electron setup plus one, not observed <S²>. Each banner establishes a run; Selected xTB 5 legacy layouts and xTB 6 modern layouts are supported; compatibility with later major versions is unconfirmed. `cpu_time` and `wall_time` independently read printed CPU/wall totals in seconds; `running_time` aliases wall time. If a total block exists, its missing quantities stay unknown; otherwise each kind uses its last printed row. Electron and vibrational temperatures are distinct. |
| <!-- parsing-basis:Geometry availability -->Available geometry | `final structure: XYZ / $coord / V2000 atom block`; `legacy setup coordinates` | XYZ/SDF positions use Å; legacy Bohr positions convert to Å. Only printed final geometry is observed; an external filename or an atom-population table is not coordinates. No adjacent file is loaded automatically. |
| <!-- parsing-basis:Energies and geometry optimization -->Energies and geometry optimization | `TOTAL ENERGY / total energy / total E`; `energy/gradient changes and limits`; `optimization convergence` | Use the last matching final total in Eh; optimization value/threshold booleans remain separate from normal program termination. The scalar gradient norm is not an atom-resolved force array. |
| <!-- parsing-basis:Orbitals, populations, and dipole -->Orbitals, populations, dipole and bond orders | `occ./eps and modern orbital tables`; `HOMO/LUMO labels/ellipsis`; `GFN1 Mulliken/CM5 or GFN2 q`; `dipole full`; `rotational constants`; `Wiberg/Mayer heading with Wiberg subheading` | Preserve total occupations, omitted orbital ranges and reported frontiers; spin splitting requires evidence. Charge schemes remain labeled; dipole full components use atomic units before public conversion. Rotation cm⁻¹ converts to GHz. Wiberg AO matrices retain threshold/unknown pairs and first-printed reciprocal rounding policy; complete WBO row indices can establish axes without coordinates. |
| <!-- parsing-basis:Vibrations and thermochemistry -->Vibrations and thermochemistry | `Projected eigval`; `indexed masses/IR/Raman`; `Hessian/intensity method labels`; `ZPVE/TOTAL ENTHALPY/FREE ENERGY`; `TOT CV/S`; `thermodynamic settings` | cm⁻¹/amu/km/mol and explicit Å⁴/amu Raman activities; old undefined (amu) arrays and g98 zero placeholders are omitted. A leading 5/6 near-zero external-mode prefix is removed with aligned indices as MolOP policy. Preserve GFN2-xTB/PTB method distinction, K temperature, RRHO/QRRHO conditions and HO/rotor reference. No missing displacement vectors, depolarization ratios or pressure are invented. Thermal energy quantities use printed Eh/particle, converted to kcal/mol; C_V/S use cal/mol/K. |
| <!-- parsing-basis:xTB analysis properties -->xTB analysis properties | `Vertical IP / EA / global electrophilicity index`; `f(+)/f(-)/f(0) atom rows` | `single_point_properties`: VIP/VEA/GEI use eV/particle before conversion; Fukui arrays are dimensionless in printed atom order. Only requested, printed tables are observed; FOD/MD/metadynamics are not represented by these fields. |

<!-- parsing-basis:end -->
