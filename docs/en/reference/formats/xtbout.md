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

MolOP parses command-line standard output from xTB major versions 5 and 6. The
reader covers both the legacy xTB 5/6.1 print family and the modern xTB 6 print
family. It does not define a separate xTB input format because xTB calculations
normally receive coordinates and options directly from the command line.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Version scope, setup, tasks, and status -->Version scope, setup, tasks, and status | Partial | Recognizes xTB major versions 5 and 6; parses version, program call, coordinate-file name, method, charge, multiplicity, OMP threads, task requests, termination, SCC status, and total wall time. | Versions before 5 and after 6 are intentionally rejected. The xTB 5 regression is a compact legacy-format contract rather than a complete vendor output capture. | `tests/test_xtbout_parser.py::test_xtbout_all_maintained_fixtures_parse_with_supported_major_versions`<br>`tests/test_xtbout_parser.py::test_xtbout_legacy_v5_contract_parses_metadata_geometry_and_results`<br>`tests/test_xtbout_parser.py::test_xtbout_explicitly_rejects_versions_outside_five_and_six` |
| <!-- feature-area:Geometry availability -->Geometry availability | Partial | Parses legacy Bohr `$coord` blocks, legacy setup coordinates, modern XYZ final structures, and embedded V2000 SDF final structures. | Single-point xTB stdout often references an external coordinate file without printing coordinates. Such logs retain a property frame with empty geometry; the parser does not read adjacent files implicitly. | `tests/test_xtbout_parser.py::test_xtbout_legacy_v5_contract_parses_metadata_geometry_and_results`<br>`tests/test_xtbout_parser.py::test_xtbout_modern_single_point_keeps_results_without_embedded_geometry`<br>`tests/test_xtbout_parser.py::test_xtbout_modern_optimization_parses_xyz_and_sdf_final_structures` |
| <!-- feature-area:Energies and geometry optimization -->Energies and geometry optimization | Partial | Exposes final total energy in Hartree, source-labeled energy evidence, optimization convergence, energy-change criteria, and xTB gradient norm fields. | Energy decomposition rows and full optimization trajectories are not structured; only the final geometry printed by stdout is represented. | `tests/test_xtbout_parser.py::test_xtbout_modern_optimization_parses_xyz_and_sdf_final_structures`<br>`tests/test_xtbout_parser.py::test_autoparser_detects_xtbout_and_preserves_source_evidence` |
| <!-- feature-area:Orbitals, populations, and dipole -->Orbitals, populations, and dipole | Partial | Parses legacy and modern orbital energies/occupancies, GFN1 Mulliken/CM5 charges, GFN2 Mulliken charges, molecular dipole, and rotational constants where printed. | Wiberg bond-order tables, quadrupoles, dispersion-property tables, and omitted orbital ranges are not fully reconstructed. | `tests/test_xtbout_parser.py::test_xtbout_legacy_v5_contract_parses_metadata_geometry_and_results`<br>`tests/test_xtbout_parser.py::test_xtbout_modern_single_point_keeps_results_without_embedded_geometry`<br>`tests/test_xtbout_parser.py::test_xtbout_modern_optimization_parses_xyz_and_sdf_final_structures` |
| <!-- feature-area:Vibrations and thermochemistry -->Vibrations and thermochemistry | Partial | Parses projected physical frequencies, reduced masses, IR intensities, zero-point energy, total enthalpy, total free energy, heat capacity, entropy, temperature, and molecular mass where printed. | Normal-mode displacement vectors, Hessian sidecar files, Raman intensities, and detailed rotor/interpolation tables are not structured from stdout. | `tests/test_xtbout_parser.py::test_xtbout_frequency_and_thermochemistry_are_structured` |
| <!-- feature-area:xTB analysis properties -->xTB analysis properties | Partial | Parses vertical IP, vertical EA, global electrophilicity index, and atom-resolved Fukui indices when requested and printed. | Additional xTB workflows such as FOD, metadynamics, MD, ONIOM, GFN-FF topology, and JSON/sidecar outputs are outside the current stdout contract. | `tests/test_xtbout_parser.py::test_xtbout_modern_single_point_keeps_results_without_embedded_geometry`<br>`tests/test_xtbout_parser.py::test_xtbout_fukui_indices_are_structured`<br>`tests/test_xtbout_parser.py::test_xtbout_vipea_and_gei_properties_are_structured` |

## xTB 5/6 Capability Coverage Matrix

The table below defines the field-level contract of the stdout reader. Statuses
have the following meanings:

- **Supported**: an explicit parser path and a regression sample for that version family exist.
- **Partial**: only known print forms, final values, or a limited subset of task semantics are covered.
- **Unverified**: the parser path may apply, but no real output sample from that version family establishes a stable contract.
- **Unsupported**: no corresponding structured field is currently produced.
- **Not applicable**: the capability is outside the responsibility of the stdout reader.

| xTB output capability | xTB 5 legacy | xTB 6 legacy/modern | Structured result | Current boundary |
| --------------------- | ------------ | ------------------- | ----------------- | ---------------- |
| Banner and version | Partial | Supported | `qm_software`, `qm_software_version` | Only major versions `5` and `6` are accepted. The current xTB 5 fixture is a compact legacy contract rather than a real vendor output capture. |
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
| Wall time | Partial | Supported | Second-valued `running_time`, aggregated at file level across segments | xTB wall-time lines are parsed; CPU time, module timings, and external scheduler time are not. |
| Optimization convergence and thresholds | Partial | Supported | `geometry_optimization_status`, energy change/threshold, gradient norm/threshold | Only final printed values are retained; per-step convergence history is not structured. |
| Complete optimization trajectory | Unsupported | Unsupported | No per-step frame sequence | A command-line run currently maps to one result frame rather than one frame per optimization step. |
| Gradient vectors | Unsupported | Unsupported | Only the `gradient` task request and final gradient norm can be recorded | Per-atom stdout gradients and the `gradient` sidecar are not mapped to structured arrays. |
| Orbital energies and occupancies | Partial | Supported | Alpha/beta orbital energies and occupancies | Legacy `occ./eps` and modern orbital tables are covered. Orbitals hidden by ellipses cannot be recovered, and non-spin-resolved occupancies are split into alpha/beta channels. |
| Atomic charges | Unverified | Partial | GFN1 Mulliken/CM5 and GFN2 Mulliken charges | Spin populations, complete population tables, and every Hamiltonian-specific charge scheme are not parsed. |
| Dipole and rotational constants | Unverified | Partial | Debye dipole vector and GHz rotational constants | Dipole magnitude, quadrupole, polarizability tensors, and higher multipoles are not structured. |
| Frequencies, reduced masses, and IR intensities | Unverified | Supported | `vibrations`, including imaginary-mode count | Only printed projected physical frequencies are covered; missing columns are not synthesized. |
| Mode displacements, Hessian, and Raman | Unsupported | Unsupported | No corresponding structured arrays | Stdout displacement blocks, `.hessian` sidecars, and Raman intensity tables are not read. |
| Thermochemistry | Unverified | Supported | ZPVE, enthalpy, free energy, `C_V`, entropy, temperature, and molecular mass | Only summary values are retained; rotor/interpolation details and individual thermal-correction components are not structured. |
| VIP, VEA, GEI, and Fukui indices | Unverified | Supported | `single_point_properties` scalars and atom-resolved Fukui arrays | Values exist only when the corresponding analysis is requested and printed to stdout. |
| Wiberg, FOD, MD, metadynamics, ONIOM, and GFN-FF topology | Unsupported | Unsupported | No dedicated fields | These analyses and workflows are outside the current stdout contract. |
| JSON, Molden, Hessian, and other auxiliary output | Unsupported | Unsupported | No sidecar merge | The parser consumes only the supplied stdout text and does not scan the working directory. |
| xTB input and output writers | Not applicable | Not applicable | Reader only | MolOP does not define an xTB input format or reconstruct command lines or stdout from result models. |

## Coordinate Boundary

xTB commonly prints only the coordinate file path for single-point jobs. MolOP
does not silently read that adjacent file: doing so would mix a second source
artifact into a parser result whose source spans and hashes describe only the
stdout file. Optimization outputs that print a final structure do expose normal
coordinate fields.

## Version Boundary

The parser accepts version declarations whose major version is `5` or `6`.
Recognized xTB outputs outside that range raise an explicit unsupported-version
error instead of falling through to Gaussian or ORCA output readers.
