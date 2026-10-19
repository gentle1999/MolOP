# ORCA Output

<!-- format-support:orcaout -->

| Item | Value |
| ---- | ----- |
| Format ID | `orcaout` |
| Extensions | `.out`, `.log`, `.orcaout` |
| Read | Yes |
| Write | No |
| Registry role | Reader |
| Data level | Coordinates and QM results |

## Quick read

Download the [shared ORCA example](../../../assets/examples/water_mp2.out), then run:

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", parser_detection="orcaout", n_jobs=1)[0][-1]
print(frame.qm_software, frame.method, frame.is_normal)
print(frame.energies.total_energy.m_as("hartree"))
```

??? example "Output"

    ```text
    ORCA MP2 True
    -74.999374598107
    ```

ORCA output parsing for user-facing result attributes such as energies, forces,
vibrations, populations, solvation, excited states, and response properties.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Molecular orbitals -->Molecular orbitals | Partial | Orbital energies, fractional occupations, alpha/beta channels and printed symmetry labels are structured. | Spin channels are aligned by source orbital index. Unprinted energies use `np.nan` in quantity arrays and `null` on JSON export; occupations and symmetry labels use `None`. Orbital coefficients are unavailable. |
| <!-- feature-area:Thermochemistry -->Thermochemistry | Partial | The last temperature block supplies ZPVE, thermal corrections, internal energy, enthalpy, Gibbs energy, entropy, molecular mass and rotational constants; temperature and pressure are exposed on the frame. | Entropy energy terms are divided by temperature. Heat capacities and additional partition-function diagnostics are unavailable. |
| <!-- feature-area:Software metadata, frames, and status -->Software metadata, frames and status | Partial | Reads software name, version, jobs, geometry frames, normal termination and runtime. | Only recognized banners, coordinate tables and termination markers are supported; truncated files can lack status or properties. |
| <!-- feature-area:Energies -->Energies | Partial | Preserves final total, reference and supported correlation energies in Hartree. | Method-specific corrections and decompositions may be unmapped; correlation corrections are not total energies. |
| <!-- feature-area:Forces and geometry optimization -->Forces and geometry optimization | Partial | Builds force arrays from complete Cartesian gradients and reads optimization criteria, thresholds and completion. | Not every optimizer diagnostic is mapped; missing or truncated gradients do not produce complete force arrays. |
| <!-- feature-area:Vibrations and vibrational spectra -->Vibrations and vibrational spectra | Partial | Frequencies, source mode indices, Cartesian normal modes, legacy/modern IR intensities, Raman activities and depolarization ratios are structured. | Mode-resolved arrays require complete rows aligned to stored frequencies. Raman activities omit experimental laser and temperature factors. |
| <!-- feature-area:Charge and spin populations -->Charge, spin populations and bond orders | Partial | Mulliken, Loewdin and Hirshfeld charge/spin tables, printed Mayer bond orders and the expectation value of total spin squared are structured. | Mayer matrices retain unprinted entries as NaN (JSON null) and preserve the print threshold; reduced orbital populations, Mayer valence diagnostics and other unsupported schemes remain unavailable. |
| <!-- feature-area:Solvation -->Solvation | Partial | Reads recognized CPCM/SMD models, solvent names and printed parameters. | Full solvent energy decompositions are unavailable; a solvent name or input request does not establish a free-energy correction. |
| <!-- feature-area:Excited states and multireference results -->Excited states and multireference results | Partial | Reads supported TDDFT, ADC2, EOM-CCSD, STEOM-CCSD, STEOM-DLPNO-CCSD and ROCIS state energies and electric-dipole transition properties. | SOC, magnetic-dipole and velocity forms remain separate from electric-dipole results; multireference input requests do not establish parsed output results. |
| <!-- feature-area:Response properties -->Electric response, NMR and total spin | Partial | Polarizability, NMR shielding tensors and isotropic spin-spin J matrices are structured, with shielding principal values, anisotropy and FC/SD/PSO/DSO J contributions. | Full Cartesian coupling tensors and SD/FC cross terms are not exposed. Coupling matrices retain source atom indices and use NaN (JSON null) for unprinted entries. |

Raman values are available as `frame.vibrations.raman_activities` (`angstrom^4/amu`) and
`raman_depolarization_ratios` (dimensionless). Use `to_raman_spectrum()` to obtain a spectrum
whose band intensities hold activities. The units and activity interpretation follow the
[ORCA manual](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/properties.html#raman-spectra).

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

The references below use ORCA 6.0/6.1 manuals and tutorials to explain syntax,
quantities and printed units. MolOP support is defined on this page; older
output headers and column orders can differ.

| Property and references | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Molecular orbitals --><!-- source:ORCA-OUTPUT --><!-- source:ORCA-METHOD -->Molecular orbitals ([Output controls][output], [Wavefunction types][hftype]) | `ORBITAL ENERGIES`; `SPIN UP/DOWN ORBITALS`; `occupation/E(Eh)/E(eV)/Irrep`; `omitted-orbital notices` | Keep printed occupations and energy units, reference/spin identity and source labels. Restricted sharing requires evidence; fractional open-shell totals do not imply equally split spin occupations. Missing energies use NaN; occupations use None. |
| <!-- parsing-basis:Thermochemistry --><!-- source:ORCA-THERMO -->Thermochemistry ([Thermochemistry][thermo], [Thermodynamics tutorial][thermo-tutorial]) | `Last THERMOCHEMISTRY AT block`; `Electronic/Zero point/thermal/enthalpy/Gibbs energies`; `Final entropy term` | Keep totals and corrections distinct, converting printed Eh/particle to kcal/mol. Derive U₀/corrections relative to E(el), and S from the printed entropy-energy term divided by T. Retain T/P, masses and rotation metadata; an incremental ORCA correction is not a total correction. |
| <!-- parsing-basis:Software metadata, frames, and status --><!-- source:ORCA-INPUT --><!-- source:ORCA-COORD --><!-- source:ORCA-RESOURCE --><!-- source:ORCA-UNITS -->Software metadata, frames and status ([Input structure][input], [Coordinates][coordinates], [Parallelism and memory][parallel], [Units][units]) | `Program Version`; `echoed INPUT FILE`; `coordinate/SCF/optimization boundaries`; `ORCA TERMINATED NORMALLY`; `TOTAL RUN TIME` | Read software/version, tasks/resources and status separately. Coordinates prefer the printed Å layout; AU coordinates convert from Bohr, while its MASS column remains amu. Frame/segment source spans establish membership; normal termination alone does not establish every convergence field. `wall_time` and its `running_time` alias convert printed day/hour/minute/second/millisecond fields to seconds; file totals sum selected segment durations. Missing `cpu_time` falls back to wall time × process count with `cpu_time_source="estimated"` and `cpu_time_num_processes`. Prefer the runtime MPI count, then `%pal nprocs` / `PALn`; complete serial input defaults to one. Unknown wall time or count leaves CPU time unknown. This is an estimate, following the [cclib ORCA approach][cpu-estimate], rather than a measured CPU counter. |
| <!-- parsing-basis:Energies --><!-- source:ORCA-ENERGY -->Energies ([Energies and gradients][energy], [MP2][mp2]) | `FINAL SINGLE POINT ENERGY`; `E(SCF)/MPn/CC reference/correlation/total labels` | Map final/reference/correlated totals separately in Eh, retaining source observations and distinguishing correlation correction from a total. Method identity limits which label can populate which field. |
| <!-- parsing-basis:Forces and geometry optimization --><!-- source:ORCA-OPT -->Forces and geometry optimization ([Geometry optimization][optimization], [Numerical gradients][gradient]) | `CARTESIAN GRADIENT`; `geometry-convergence values/limits`; `optimization-converged markers` | `forces` = negative gradient in Eh/Bohr; source atom order remains explicit. Criteria, thresholds and completion are observed independently; truncated/missing gradients are not padded as zeros. |
| <!-- parsing-basis:Vibrations and vibrational spectra --><!-- source:ORCA-VIB -->Vibrations and vibrational spectra ([IR and Raman][vibrations]) | `VIBRATIONAL FREQUENCIES`; `NORMAL MODES`; `IR SPECTRUM Int/T**2`; `RAMAN SPECTRUM activity/depolarization` | Align printed mode IDs with cm⁻¹ frequencies, Cartesian vectors and km/mol IR. Modern IR uses `Int` (km/mol), excluding its `T**2` atomic-unit column; legacy IR follows the older table header. Activity is Å⁴/amu; unspecified polarization uses the unsuffixed ratio field. Preserve normalization/source axes and omit incomplete optional tables. |
| <!-- parsing-basis:Charge and spin populations --><!-- source:ORCA-POP -->Charge, spin populations and bond orders ([Populations and Mayer analysis][population]) | `MULLIKEN/LOEWDIN/HIRSHFELD atomic charge/spin tables`; `Mayer bond orders larger than threshold` | Charge/spin series retain scheme and channel; Mayer records retain the print threshold. Atom-pair matrices use source atom indices and preserve NaN for unprinted pairs rather than zero. |
| <!-- parsing-basis:Solvation --><!-- source:ORCA-SOLVENT -->Solvation ([Solvation][solvation]) | `CPCM/SMD/solvation section labels and echoed input solvent options` | `solvent` reflects recognized model/name and printed metadata. Input solvation requests and output energy corrections remain distinct; a solvent name does not establish a free energy. |
| <!-- parsing-basis:Excited states and multireference results --><!-- source:ORCA-EXCITED -->Excited states and multireference results ([Excited states][excited]) | `Covered TD/ADC/EOM/STEOM/ROCIS root tables`; `electric-dipole absorption spectrum` | Populate root/state energies in eV and covered electric-dipole transitions (atomic-unit vectors converted to Debye), with explicit state identity. SOC/magnetic/velocity spectra and unsupported MRCI/CASSCF result tables are not silently read as the same electric-dipole state collection; multireference input requests do not prove parsed results. |
| <!-- parsing-basis:Response properties --><!-- source:ORCA-ELECTRIC --><!-- source:ORCA-NMR --><!-- source:ORCA-SPIN -->Electric response, NMR and total spin ([Electrical properties][electric], [NMR][nmr], [SCF][scf], [UHF spin diagnostics][spin]) | `DIPOLE MOMENT/raw polarizability`; `Total shielding tensor (ppm)`; `J[...] / isotropic coupling summary`; `<S**2>` | Electric dipole vectors convert from ea₀ to Debye; retain polarizability tensors. NMR shielding uses ppm, J uses Hz with source atom-index mapping and independent contributions. <S²> is observed; effective S is derived via sqrt(<S²>+1/4)-1/2. None of these are inferred from multiplicity alone. |

<!-- parsing-basis:end -->

[input]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html
[coordinates]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html
[parallel]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html
[output]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/output.html
[hftype]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hftype.html
[energy]: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/energygradients.html
[mp2]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/MP2.html
[population]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html
[scf]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html
[spin]: https://www.faccts.de/docs/orca/6.1/tutorials/prop/single_point.html
[vibrations]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html
[thermo]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html
[thermo-tutorial]: https://www.faccts.de/docs/orca/6.1/tutorials/prop/thermo.html
[nmr]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nmr.html
[electric]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/electric.html
[optimization]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html
[gradient]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalgradients.html
[solvation]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html
[excited]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/excitedstates.html
[units]: https://www.faccts.de/docs/orca/6.1/manual/contents/orcaarchitecture/conversionfactors.html

[cpu-estimate]: https://github.com/cclib/cclib/blob/master/cclib/parser/orcaparser.py
