# Gaussian Output

<!-- format-support:g16log -->

| Item | Value |
| ---- | ----- |
| Format ID | `g16log` |
| Extensions | `.log`, `.out`, `.g16`, `.gal`, `.irc`, `.gau` |
| Read | Yes |
| Write | No |
| Registry role | Reader |
| Data level | Coordinates and QM results |

## Quick read

```python
from molop import AutoParser

frame = AutoParser("calculation.log", parser_detection="g16log")[0][-1]
print(frame.is_normal)
if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))
```

??? example "Output contract"

    Output is termination status and final total energy, for example `True` and `-317.592366596`;
    the value depends on the file.

See [Find fields by scientific property](../model_fields.md) for other results.

MolOP reads Gaussian output files and extracts the data most often needed for
computational chemistry post-processing: structures, energies, thermochemistry,
vibrations, orbitals, populations, gradients, NMR and response properties, archive
records, and termination status. Coverage is described using Gaussian output
concepts familiar to computational chemists. “Partial” means common or tested
print forms are structured; it is not a claim of complete Gaussian grammar
coverage.

| Capability | Support | Extracted data | Boundaries |
| ---------- | ------- | -------------- | ---------- |
| <!-- feature-area:Gaussian job metadata -->Gaussian job metadata | Example-covered | Software/version metadata, route text, title, charge, multiplicity, and accumulated running time where present. | Fields are normalized; original log line ordering and byte-level layout are not preserved. |
| <!-- feature-area:Title card -->Title card | Example-covered | Gaussian title card text as a normalized job title. | Only the title text is advertised; surrounding raw Gaussian formatting is not preserved. |
| <!-- feature-area:Route keywords and calculation setup -->Route keywords and calculation setup | Partial | Raw route text plus recognized model chemistry, job types, dispersion, solvation, population requests, and HF handling. | Method-specific Gaussian keywords outside covered cases may remain raw text. |
| <!-- feature-area:Input and standard orientations -->Input and standard orientations | Example-covered | Atoms and coordinates from input orientation and standard orientation blocks. | Distance-matrix and stoichiometry printouts are not advertised as separate structured fields. |
| <!-- feature-area:Rotational constants -->Rotational constants | Example-covered | Rotational constants printed by Gaussian, exposed as frame-level frequency quantities. | Principal-axis diagnostic text is not preserved separately. |
| <!-- feature-area:SCF and electronic energies -->SCF and electronic energies | Partial | SCF, reference/electronic, and method-specific energies; main-log values take precedence over duplicate archive-tail values. | Only tested energy fields are advertised; not every post-HF or correction energy table is guaranteed structured. |
| <!-- feature-area:Molecular orbitals and population analysis -->Molecular orbitals and population analysis | Partial | Molecular orbital energies/occupancies; Mulliken charge/spin, APT, Lowdin, Hirshfeld/CM5, NPA, and ESP charge series; electronic spatial extent; multipoles; and early polarizability values where present. | Orbital coefficient matrices, NBO orbital/bond-order details, and schemes without a recognized atomic table are not claimed. |
| <!-- feature-area:Vibrational frequencies and IR intensities -->Vibrational frequencies and IR intensities | Partial | Frequencies, reduced masses, force constants, IR intensities, per-mode displacement vectors, and imaginary-mode flags. | Raman, VCD, and other spectrum variants are not advertised unless explicitly covered later. |
| <!-- feature-area:Thermochemistry -->Thermochemistry | Partial | Temperature, pressure, molecular mass, moments of inertia, rotational symmetry number, rotational/vibrational temperatures, rotational constants, ZPVE, thermal energy, enthalpy, Gibbs free energy, entropy, and heat capacity. | Coverage is mainly tied to frequency-style Gaussian outputs and maintained examples. |
| <!-- feature-area:Dipole and polarizability -->Dipole and polarizability | Partial | Dipole and polarizability values from covered Gaussian response-property sections. | Only example-backed response fields are advertised; not every response-property print variant is guaranteed structured. |
| <!-- feature-area:NMR shielding and spin-spin coupling -->NMR shielding and spin-spin coupling | Partial | Parses per-atom Gaussian GIAO shielding tensors, observed isotropic/anisotropic and principal values, total K/J coupling matrices, and FC, SD, PSO, and DSO contribution matrices in Hz. | Coverage targets the standard SCF GIAO print family; referenced chemical shifts, EPR, and other gauges are not structured. |
| <!-- feature-area:Cartesian gradients and forces -->Cartesian gradients and forces | Partial | Cartesian force arrays from Gaussian force sections. | Only the normalized force array is advertised; auxiliary force diagnostics are not separately recorded. |
| <!-- feature-area:Cartesian Hessian -->Cartesian Hessian | Partial | Cartesian Hessian data from Gaussian second-derivative sections. | The contract is the normalized Hessian field, not every printed second-derivative diagnostic. |
| <!-- feature-area:Geometry optimization convergence -->Geometry optimization convergence | Partial | Berny optimization summaries, including convergence thresholds, force/displacement values, energy change, and optimized-state flags. | Optimizer diagnostics outside tested examples may remain unstructured. |
| <!-- feature-area:Gaussian archive section -->Gaussian archive section | Partial | Archive-tail metadata, coordinates, energies, thermochemistry, polarizability, and Hessian fallback or augmentation fields. | Main-log fields take precedence where applicable; archive-tail data does not invent live status, temperature, or pressure. |
| <!-- feature-area:Termination status -->Termination status | Partial | Gaussian termination evidence is segment-scoped and contributes to the file-level aggregate without being copied into frames. | Frame status retains frame-local SCF evidence; detailed failure taxonomy outside structured parse diagnostics is not inferred. |
| <!-- feature-area:CPU and elapsed time -->CPU and elapsed time | Example-covered | Job CPU and elapsed-time style records accumulated into the running-time value. | Per-link timing rows are not exposed as separate timing records. |
| <!-- feature-area:Link1 multi-step jobs -->Link1 multi-step jobs | Partial | Link1 section metadata is propagated to later frames so multi-step Gaussian jobs keep the expected per-frame context. | Low-level link boundary rows are not exposed as user-facing records. |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Parsed Gaussian output can be converted to coordinate and graph formats through the registry. | Conversion quality depends on successful structure and graph recovery. |

## Gaussian Output Capability Matrix

`g16log` is a structured reader and does not provide a Gaussian log writer. In the
matrix, fakeG indicates whether the field can participate in the
[Gaussian-like renderer](fakeg.md); it does not imply reconstruction of the source log.

| Gaussian output content | Parse | Structured target | fakeG | Source fidelity | Current boundary |
| ----------------------- | ----- | ----------------- | ----- | --------------- | ---------------- |
| Gaussian software and version | Fixture-covered | `qm_software`, `qm_software_version` | Supported | Unsupported | Maintained Gaussian version banners are covered; an exact version cannot be inferred from a truncated file without a header. |
| Link0/options | Partial | File/frame `options` text | Supported | Unsupported | Printed Link0 options are stored; fakeG can derive a CPU line from `%nprocshared`, but complete structured resource projection is not guaranteed. |
| Route and title | Fixture-covered | `keywords`, `title_card`, `semantic_route` | Supported | Unsupported | Route semantics are shared with GJF; title and route are normalized without decorative lines. |
| Method, functional, basis, and tasks | Partial | `model_chemistry`, `task_requests`, and compatibility fields | Partial | Not applicable | Only shared-route keywords are structured; requested properties are never fabricated solely from the route. |
| Charge and multiplicity | Supported | Frame/file top-level fields | Supported | Unsupported | Extracted from Gaussian charge/multiplicity lines or archive metadata; main-log values take priority. |
| Link1 / multi-segment context | Partial | Segment/frame indices and propagated task metadata | Partial | Unsupported | Applicable section metadata propagates to later frames; low-level Link enter/leave rows are not public business records. |
| Input orientation | Fixture-covered | `atoms`, `coords` | Supported | Unsupported | Coordinates are normalized without preserving table whitespace or center numbers. |
| Standard orientation | Fixture-covered | `standard_coords`, transformation matrix | Supported | Unsupported | The input/standard rigid transform can be computed; standard orientation may backfill coordinates when input orientation is absent. |
| Structure-only fast parsing | Supported | Atoms, coordinates, and basic frame data | Not applicable | Not applicable | `only_extract_structure` stops after orientation and skips energies, frequencies, thermochemistry, and other expensive fields. |
| Rotational constants | Fixture-covered | `rotation_constants` | Raw/limited | Unsupported | Exposed as a GHz array; fakeG has no complete canonical rotational-constant renderer. |
| SCF/reference energy | Supported | `energies.reference_energy`, energy observations | Supported | Unsupported | Extracted from `SCF Done`; fakeG emits `E(SCF)` without the original method label or cycle history. |
| MP2, MP3, MP4, MP5 | Partial | `energies.mp*_energy` | Unsupported | Unsupported | Covered `EUMP*`/MP5 lines parse; all MP corrections and spin-component tables are not guaranteed. |
| CCSD and CCSD(T) | Partial | `energies.ccsd_energy`, `ccsd_t_energy` | Unsupported | Unsupported | Covered total-energy lines parse; other coupled-cluster variants and correction tables are not guaranteed. |
| Archive energy fallback | Partial | Merged `energies` | Indirect | Unsupported | Main-log observations take priority; archive values augment missing or duplicate fields without overriding stronger live evidence. |
| Total spin | Partial | `total_spin.spin_square`, `spin_quantum_number` | Supported | Unsupported | Only standard `S**2` / `S` print forms are covered. |
| MO energies, occupations, and symmetry | Partial | `molecular_orbitals` | Raw/unsupported | Unsupported | Alpha/beta, open-shell, and concatenated-number regressions are covered; complete MO coefficient matrices are not parsed. |
| Mulliken charge / spin populations | Partial | `charge_spin_populations` | Raw/unsupported | Unsupported | Separate and combined open-shell charge/spin tables are structured; headings, sums, and source precision are not retained. |
| APT / Lowdin populations | Partial | `charge_spin_populations` | Raw/unsupported | Unsupported | Coverage is limited to standard tables in maintained examples. |
| Hirshfeld / CM5 | Partial | `populations["hirshfeld_charges"]`, `populations["hirshfeld_spins"]`, `populations["cm5_charges"]` | Raw/unsupported | Unsupported | Hirshfeld charge/spin and CM5 charge series parse. |
| NPA / ESP atomic charges | Partial | `populations["npa_charges"]`, `populations["esp_charges"]` | Raw/unsupported | Unsupported | The last complete NPA summary and ESP atomic charge table in a frame are retained; detailed NBO orbitals and bond orders are not parsed. |
| Electronic spatial extent | Partial | `polarizability.electronic_spatial_extent` | Raw/unsupported | Unsupported | Only the standard scalar line in population output is covered. |
| Dipole and higher multipoles | Partial | Dipole, quadrupole, traceless quadrupole, octapole, and hexadecapole fields | Raw/unsupported | Unsupported | Covered field-independent/response forms parse; all units and frequency-dependent variants are not guaranteed. |
| Polarizability | Partial | Isotropic, anisotropic, and tensor fields | Raw/unsupported | Unsupported | Population, response, and archive sources merge; later explicit response data may replace an earlier approximation. |
| Harmonic frequencies | Partial | `vibrations.frequencies`, imaginary-mode flags | Supported | Unsupported | Concatenated numeric output is handled; anharmonic/VPT2 result tables are not parsed. |
| Reduced mass, force constant, and IR | Partial | Corresponding `vibrations` arrays | Supported | Unsupported | Values align by mode; Raman activity and VCD/ROA are not in the structured contract. |
| Normal-mode displacement | Partial | `vibrations.vibration_modes` and axis metadata | Supported | Unsupported | Normalized to mode/atom/Cartesian; normalization and mass weighting remain `unknown`. |
| Temperature and pressure | Partial | Frame `temperature`, `pressure` | Supported | Unsupported | Read from standard frequency/thermochemistry lines; archive data does not invent absent temperature or pressure. |
| Molecular mass, inertia, and rotational data | Partial | `thermal_informations` | Supported | Unsupported | Covers mass, moments, symmetry number, rotational temperatures, and rotational constants. |
| Vibrational temperatures | Partial | `vibrational_temperatures` and mode indices | Supported | Unsupported | Positive-frequency modes are mapped first; only demonstrable array relations are retained when counts disagree. |
| ZPVE and thermal corrections | Partial | ZPVE/TCE/TCH/TCG in `thermal_informations` | Supported | Unsupported | Values use normalized Hartree/particle units; all partition-function diagnostics are not included. |
| U0, UT, H, G, entropy, and Cv | Partial | `thermal_informations` | Supported | Unsupported | Standard Gaussian thermochemistry summaries parse; source pretty-print layout is not modeled. |
| Cartesian forces | Partial | `forces` plus axis/order/orientation metadata | Raw/unsupported | Unsupported | Shape is `(N, 3)` in source atom order; orientation is normally `unknown`. |
| Cartesian Hessian | Partial | `hessian` plus axis/order/orientation metadata | Raw/unsupported | Unsupported | Main second-derivative data take priority and archive can backfill; shape must be `(3N, 3N)`. |
| Berny convergence | Partial | `geometry_optimization_status` | Supported | Unsupported | Force/displacement values, thresholds, energy change, and optimized state are covered; other optimizers may remain unstructured. |
| SCF and termination status | Partial | Frame-local `status`, segment/file aggregate status | Synthetic | Unsupported | SCF evidence is frame-local; normal/error termination is segment/file scoped. fakeG termination is not source-status reproduction. |
| CPU / elapsed time | Fixture-covered | Accumulated `running_time` | Partial | Unsupported | File totals aggregate segments; per-Link timing records are not public, and multiframe fakeG removes frame runtime lines. |
| Archive-tail metadata/coordinates | Partial | Missing-field fallback | Indirect | Unsupported | Archive can augment metadata, coordinates, energies, thermochemistry, polarizability, and Hessian; live data take priority. |
| TD excited-state results | Unsupported | No stable public result contract | Unsupported | Unsupported | The route may identify a TD request, but excitation energies and oscillator strengths are not structured. |
| NMR magnetic shielding | Partial | `nmr.gauge`, per-atom `shielding_tensors`, isotropic, anisotropy, and principal values | Unsupported | Unsupported | Covers standard SCF GIAO `(ppm)` output; tensors align by source atom index and Cartesian orientation is currently `unknown`. Per-atom shielding values can be projected into RDKit atom properties with `qm_embedded_rdmol(embed_nmr=True)`. |
| NMR spin-spin coupling | Partial | `spin_spin_coupling_k`, `spin_spin_coupling_j`, `spin_spin_coupling_k_components`, `spin_spin_coupling_j_components`, and `coupling_atom_indices` | Unsupported | Unsupported | Expands blockwise lower-triangle total and FC/SD/PSO/DSO K/J tables into symmetric Hz matrices; coupling matrices remain frame-level atom-pair data. |
| EPR, referenced chemical shifts, NBO bond orders, and other specialized results | Unsupported | No stable public result contract | Unsupported | Unsupported | A route request does not imply result coverage; each property requires fixtures and a public semantic contract. |
| Coordinate/graph conversion | Supported | Frame structure / recovered graph | Not applicable | Not applicable | Coordinate conversion requires extracted structure; graph conversion also depends on bond-order and charge-state recovery. |
| Lossless source-log rewriting | Unsupported | Component raw snippets are inspection/fallback only | Unsupported | Unsupported | The parser is a field extractor, not a complete Gaussian-output CST, and has no byte-for-byte writer. |

`capture_source_evidence=True` retains source evidence for energies, SCF status, and
optimization convergence for audit-oriented workflows. Default parsing focuses on result
extraction. Neither mode preserves the complete original output structure.
