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
| <!-- feature-area:File recognition and fixture inventory -->File recognition and fixture inventory | Fixture-covered | Reads fchk scalar and array records; every maintained fchk fixture, including large files, produces one result frame. | Binary `.chk` files are not accepted. Unknown records gain structured semantics only after an explicit shared-field mapping exists. |
| <!-- feature-area:Metadata and route semantics -->Metadata and route semantics | Partial | Parses title, Gaussian version, route, charge, multiplicity, model chemistry, tasks, dispersion, solvent, and open-shell treatment; the second-line summary fills method/basis gaps in the route grammar. | Route semantics remain bounded by the shared Gaussian grammar; fchk does not contain complete Link 0 input or runtime. |
| <!-- feature-area:Geometry and Cartesian derivatives -->Geometry and Cartesian derivatives | Partial | Parses atomic numbers, Bohr coordinates, forces obtained by negating the Cartesian gradient, and a symmetric Hessian expanded from packed force constants. | An fchk represents one current geometry rather than a full optimization trajectory; coordinates and derivatives remain in source order without reconstructing the original input transform. |
| <!-- feature-area:Energies, orbitals, and spin -->Energies, orbitals, and spin | Partial | Parses SCF/reference, MP2, MP3, MP4, CCSD, CCSD(T), and final total energies, alpha/beta orbital energies and occupancies, and S-squared. | MO coefficients, density matrices, natural orbitals, and method-specific energy labels outside the mapped set are not exposed. |
| <!-- feature-area:Populations and electric response -->Populations and electric response | Partial | Parses atom-count-aligned Mulliken, NPA, ESP, APT, Hirshfeld/CM5, and available spin population records, plus dipole, packed polarizability, and quadrupole values with atomic-unit conversion. | Only records present in the fchk are exposed; bond orders, hyperfine tensors, and higher electric-response arrays are not mapped. |
| <!-- feature-area:Vibrations, thermochemistry, and status -->Vibrations, thermochemistry, and status | Partial | Parses frequencies, reduced masses, force constants, IR intensities, normal modes, thermal energy/enthalpy/free energy, Job Status, and optimization completion. | Raman/ROA/VCD, detailed thermal corrections, temperature/pressure, and optimization convergence history are not mapped. |
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
| Thermal totals | Partial | `U_T`, `H_T`, `G_T` | Current fixtures do not contain full ZPVE, entropy, and heat-capacity decomposition. |
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
