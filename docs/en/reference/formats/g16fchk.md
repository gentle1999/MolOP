# Gaussian Formatted Checkpoint

<!-- format-support:g16fchk -->

| Item | Value |
| ---- | ----- |
| Format ID | `g16fchk` |
| Extensions | `.fchk`, `.fch`, `.fck` |
| Read | Yes |
| Write | No |
| Registry role | Reader |
| Data level | Single-frame Gaussian QM results and coordinates |

## Quick read

```python
from molop import AutoParser

frame = AutoParser("molecule.fchk", parser_detection="g16fchk")[0][0]
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy if frame.energies else None)
print(frame.charge_spin_populations.population_names if frame.charge_spin_populations else [])
```

??? example "Output contract"

    Output reports atom/coordinate size, available total energy, and population schemes actually
    present.

An fchk produces one current-geometry frame, not the full optimization trajectory from a log.

MolOP reads textual formatted checkpoints produced by Gaussian `formchk`. The
reader decodes fixed-width records and `N=` array lengths without depending on
adjacent field ordering. Unknown records are skipped.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:File recognition and fixture inventory -->Text records and file recognition | Partial | Reads fixed-width scalar and array records into one result frame for the current geometry. | Binary `.chk` is not accepted; unmapped records do not produce structured properties. |
| <!-- feature-area:Metadata and route semantics -->Metadata and route settings | Partial | Parses title, Gaussian version, route, charge, multiplicity, model chemistry, tasks, dispersion, solvent, and open-shell treatment; the second-line summary fills method/basis gaps in the route grammar. | Route semantics remain bounded by the shared Gaussian grammar; fchk does not contain complete Link 0 input or runtime. |
| <!-- feature-area:Geometry and Cartesian derivatives -->Geometry and Cartesian derivatives | Partial | Parses atomic numbers, Bohr coordinates, forces obtained by negating the Cartesian gradient, and a symmetric Hessian expanded from packed force constants. | An fchk represents one current geometry rather than a full optimization trajectory; coordinates and derivatives remain in source order without reconstructing the original input transform. |
| <!-- feature-area:Energies, orbitals, and spin -->Energies, orbitals and spin | Partial | Parses SCF/reference, MP2, MP3, MP4, CCSD, CCSD(T), and final total energies, alpha/beta orbital energies and occupancies, and S-squared. | MO coefficients, density matrices, natural orbitals, and method-specific energy labels outside the mapped set are not exposed. |
| <!-- feature-area:Populations and electric response -->Populations and electric response | Partial | Parses atom-count-aligned Mulliken, NPA, ESP, APT, Hirshfeld/CM5, and available spin population records, plus dipole, packed polarizability, and quadrupole values with atomic-unit conversion. | Only records present in the fchk are exposed; bond orders, hyperfine tensors, and higher electric-response arrays are not mapped. |
| <!-- feature-area:Vibrations, thermochemistry, and status -->Vibrations, thermochemistry and status | Partial | Parses frequencies, reduced masses, force constants, IR intensities, normal modes, thermal energy/enthalpy/free energy, Job Status, and optimization completion. | Raman activity and P/U ratios require complete polarizability derivatives. TCE/TCH/TCG are derived with equations; printed temperature/pressure are retained. ROA/VCD and convergence history remain unmapped. |
| <!-- feature-area:NMR shielding and spin-spin coupling -->NMR shielding and spin-spin coupling | Partial | Parses per-atom shielding tensors and isotropic, anisotropic, and principal values; FC, SD, PSO, and DSO reduced-coupling matrices; and their total K matrix in Hz. | Requires the corresponding NMR records. fchk lacks the isotope information needed to reconstruct J; referenced chemical shifts and EPR data are not mapped. |

## Capability Coverage Matrix

| fchk capability | Read | Structured result | Current boundary |
| --------------- | ---- | ----------------- | ---------------- |
| Fixed-width scalar records | Supported | `I`, `R`, `C`, `L`, and `H` values | Unknown labels are skipped. |
| `N=` array records | Supported | Numeric values or 12-character chunks read to the declared length | Truncated arrays raise a format mismatch. |
| Title, route, and Gaussian version | Supported | `title_card`, `keywords`, `qm_software_version` | Original record layout is not retained. |
| Charge and multiplicity | Supported | `charge`, `multiplicity` | Uses the fchk scalar values. |
| Method, functional, basis, and tasks | Partial | `model_chemistry`, `task_requests` | Route takes precedence; the second-line summary fills missing method/basis data, while unknown route tokens remain diagnostic. |
| File and frame count | Supported | One segment and one terminal/SP frame | fchk does not represent the multi-Link or optimization-step timeline of a Gaussian log. |
| Atoms and current Cartesian coordinates | Supported | Angstrom `coords`, source-order `atoms` | Bohr coordinates are converted to Angstrom. |
| Cartesian gradient | Supported | Hartree/Bohr `forces` | Forces are the elementwise negative of the gradient. |
| Cartesian force constants | Supported | Symmetric `(3N, 3N)` `hessian` | Only a valid packed lower-triangle length is accepted. |
| SCF, MP2-4, CCSD, CCSD(T), and total energies | Supported | `Energies` and optional source evidence | Unmapped energy labels are not assigned guessed semantics. |
| Alpha/beta orbital energies and occupancies | Supported | `MolecularOrbitals` | If a closed-shell file omits beta energies, alpha energies are copied and beta occupancy uses the beta-electron count. |
| MO coefficients and density matrices | Unsupported | None | Large raw arrays are skipped safely. |
| Atomic population records | Partial | `ChargeSpinPopulations`, including extensible ESP/NPA spin series | Mulliken, NPA, ESP, APT, Hirshfeld/CM5, and known spin labels are mapped only when arrays match the atom count. |
| S-squared | Supported | `TotalSpin` | Spin quantum number is derived from S-squared. |
| Dipole, polarizability, and quadrupole | Partial | `Polarizability` | fchk packed ordering is retained; full named tensors are not expanded. |
| Frequencies, masses, force constants, and IR | Supported | `Vibrations` | Reads the known first four per-mode groups in `Vib-E2`. |
| Normal-mode displacements | Supported | One `(N, 3)` array per mode | Marked with source-program normalization and unknown mass weighting. |
| Thermal totals | Partial | `U_T`, `H_T`, `G_T` | Full ZPVE, entropy and heat-capacity decompositions are not provided. |
| NMR magnetic shielding | Partial | `nmr.gauge`, per-atom `shielding_tensors`, isotropic, anisotropy, and principal values | Reads `NMR shielding`; tensor units are converted to ppm and Cartesian orientation remains `unknown`. Per-atom shielding values can be projected into RDKit atom properties with `qm_embedded_rdmol(embed_nmr=True)`. |
| NMR reduced spin-spin coupling | Partial | `spin_spin_coupling_k`, `spin_spin_coupling_k_components`, and `coupling_atom_indices` | Reads four packed FC/SD/PSO/DSO K contributions and sums them to total K in Hz; `spin_spin_coupling_j` and J contributions remain absent because fchk lacks isotope metadata; coupling matrices remain frame-level atom-pair data. |
| Job Status and optimization completion | Partial | `Status`, `GeometryOptimizationStatus` | Completion is inferred from an opt request plus normal Job Status; per-step convergence values are absent. |
| Source spans and hashes | Supported | One whole-file segment/frame provenance record | No `.chk`, log, or other sidecar is merged. |
| Basis, ECP, MO, and density raw arrays | Skipped | No dedicated shared fields | Being scannable does not imply structured semantic support. |
| fchk/chk writer | Unsupported | Reader only | MolOP does not rebuild fchk or process binary checkpoints. |

## Format Boundary

`.fchk` is a textual Gaussian formatted checkpoint, not a binary `.chk`. The
reader does not invoke Gaussian `formchk` or load a neighboring log/checkpoint;
all provenance describes only the supplied fchk text.

Calculated Raman data expose activity and separate P/U depolarization arrays
when complete polarizability derivatives establish their presence. FCHK thermal
corrections `TCE`, `TCH` and `TCG` retain their equations in `derivations`.
The [shared data model](../data_model.md)
describes container relationships and shared semantics.

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [Formatted checkpoint reference](https://gaussian.com/formchk/).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:File recognition and fixture inventory -->Text records and file recognition | `Fixed-width scalar/array records`; `I/R/C/L/H`; `N= count`; `header` | Decode declared record type and count, including D exponents. Only mapped labels gain semantics; binary .chk is not accepted. |
| <!-- parsing-basis:Metadata and route semantics -->Metadata and route settings | `Title and second-line job/method/basis`; `Route`; `Gaussian Version`; `Charge`; `Multiplicity` | Map metadata and input requests; header tokens fill route gaps. No runtime or complete Link0 is inferred from a formatted checkpoint. |
| <!-- parsing-basis:Geometry and Cartesian derivatives -->Geometry and Cartesian derivatives | `Atomic numbers`; `Real atomic weights`; `Current cartesian coordinates`; `Cartesian Gradient/Force Constants` | `coords`: Bohr→Å; masses: amu; `forces` = negative gradient (Eh/Bohr); packed lower-triangle Hessian expands symmetrically in Eh/Bohr². Axes remain source atom order. |
| <!-- parsing-basis:Energies, orbitals, and spin -->Energies, orbitals and spin | `SCF/MPn/Cluster/Total Energy`; `Alpha/Beta Orbital Energies`; `electron counts`; `S**2` | Energy labels map separately in Eh. Occupations derive from spin electron counts; a missing beta spectrum is shared only with an evidenced restricted reference. S = sqrt(<S²>+1/4)-1/2 is derived. |
| <!-- parsing-basis:Populations and electric response -->Populations and electric response | `Mulliken/NPA/ESP/APT/Hirshfeld/CM5 charge/spin arrays`; `Dipole Moment`; `Polarizability`; `Quadrupole Moment` | Atom-count-aligned population series retain scheme and spin labels. Dipole ea₀ and quadrupole ea₀² convert to Debye and Debye·Å; polarizability uses a₀³ and retains packed ordering. These are not bond-order matrices. |
| <!-- parsing-basis:Vibrations, thermochemistry, and status -->Vibrations, thermochemistry and status | `Number of Normal Modes`; `Vib-E2/Vib-Modes`; `Polarizability Derivatives`; `Thermal Energy/Enthalpy/Free Energy`; `Temperature/Pressure`; `Job Status` | Vib-E2 blocks 0–3 = frequency (cm⁻¹), mass (amu), force constant (mdyne/Å), IR (km/mol); 4–6 = Raman activity (Å⁴/amu), P/U, requiring complete polarizability derivatives. Mode arrays use (mode,atom,xyz) order. Thermal totals convert from printed Eh/particle to kcal/mol; TCE/TCH/TCG are totals minus Total Energy with derivations recorded. Status follows the printed label; optimization history is unavailable. |
| <!-- parsing-basis:NMR shielding and spin-spin coupling -->NMR shielding and spin-spin coupling | `NMR shielding`; `Isotropic Spin-Spin Couplings`; `NMR route gauge` | Transpose stored 3×3 shielding tensors and convert with α²×10⁶ to ppm; isotropic/anisotropic/principal values are derived. Four packed FC/SD/PSO/DSO blocks form K in Hz and their sum; J is not synthesized without isotope data. Atom indices are 0-based. |

<!-- parsing-basis:end -->
