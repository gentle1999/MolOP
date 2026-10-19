# Data containers and cross-QM semantics

MolOP maps program inputs and outputs into shared data containers. Gaussian, ORCA
and xTB results expose common fields while retaining source identity, calculation
conditions and conventions needed to interpret their differences.

This page describes public containers and rules for comparing programs. Use the
[scientific field index](model_fields.md) to locate a property and [format support](format_support.md)
for fields each reader fills and their parsing references. Defining a field in a
container does not mean every program or calculation provides it.

## Data hierarchy

`AutoParser(path)` returns a batch of files; `batch[0]` is a file and `batch[0][-1]`
is its last retained frame. A file can contain multiple calculation steps
(segments), each with multiple geometry or result frames.

| Level | Container and access | Contents |
| ----- | -------------------- | -------- |
| Batch | `FileBatchModelDisk` | Files, selection, summaries and batch conversion. |
| File | `BaseChemFile`; `BaseCoordsFile`, `BaseQMInputFile` and `BaseCalcFile` for structures, inputs and results | Source format, provenance, steps, file-level setup and frames. |
| Frame | `BaseChemFileFrame`; `BaseCoordsFrame`, `BaseQMInputFrame` and `BaseCalcFrame` | Structure, requests, properties and status for one geometry or result. |
| Property | For example, `frame.energies` or `frame.vibrations` | Domain data, sometimes containing individual records. |
| Generic view | `PropertyTable`, `Spectrum`, `PropertyBundle` and related types | Tables, spectra and data for export to other systems. |

Choose the corresponding job and frame before comparing multistep files.
File-level summaries do not replace frame-specific conditions. `file_frame_index`
retains source-file order; `segment_index` / `segment_frame_index` identify the step
and its frame. `frame_id` belongs to the current frame collection and should not
serve as source identity after selection or reordering.

## Public data containers

The scientific types below are importable from `molop.io.base_models.DataClasses`;
generic representation types live in `molop.io.base_models.Bases`. The tables
describe purpose and access; see the [Python API](api/io.md) for full signatures.

### Input requests and calculation setup

| Container | Access or relationship | Purpose |
| --------- | ---------------------- | ------- |
| `QMModelChemistry` | `model_chemistry` | Method family, method, reference, functional, basis, dispersion, spin treatment and requested solvation. |
| `QMBasisSet` | Records in `model_chemistry.basis_sets` | Basis roles and scopes: orbital, auxiliary, ECP and element/atom-specific assignments. |
| `QMResourceRequest` | `resource_request` | CPU count, memory, program-specific options and raw directives. |
| `QMTaskRequest` | Records in `task_requests` | Single point, optimization, frequency, scan, derivative order, target state and requested properties. |
| `ExcitedStateRequest` | Records in `excited_state_requests` | Method, requested roots, target root, spin channels and root following. |
| `ActiveSpace` | `active_space` in requests or results | Active electrons, orbital count and active/inactive/frozen orbital lists. |
| `MultireferenceStateBlock` | Records in `MultireferenceRequest.state_blocks` | Root counts, reference space and active space grouped by multiplicity and symmetry. |
| `MultireferenceRequest` | Records in `multireference_requests` | Multireference method, state blocks, thresholds and requested corrections. |
| `ExplicitSolventRequest` | ORCA input `explicit_solvent_requests` | Solvent placement count, droplet/cluster settings and solvent-file requests. |

### Coordinates and comments

| Container or type | Access or relationship | Purpose |
| ----------------- | ---------------------- | ------- |
| `CoordinateContainer` | Generic ordered coordinate records | `items`, indexing, iteration and separation of real atoms from dummy/ghost atoms. |
| `AtomInInternalCoords` | One internal-coordinate atom | Reference indices, distance, angle, dihedral and dummy/ghost flags. |
| `InternalCoords` | Collection of `AtomInInternalCoords` | Internal/Cartesian coordinate conversion; lengths use Å and angles degrees. |
| `CoordinateParameter` | One coordinate parameter or scan variable | Name, raw expression, start, stop and step count. |
| `CoordinateParameters` | Collection of `CoordinateParameter` | Ordered parameters, scans and available numeric values. |
| `Comment`, `CommentContainer` | File or frame `comments` | Text, raw source, source line and ordered annotations; `CommentKind` identifies the annotation kind. |

Cartesian structure resides directly in frame `atoms`, `coords`, `charge` and
`multiplicity`; `coords` has shape `(N, 3)`. Distinguish observed coordinates,
internal-coordinate conversions and generated SMILES 2D coordinates by provenance.

### Energy, thermochemistry, solvation and status

| Container | Common access | Purpose |
| --------- | ------------- | ------- |
| `Energies` | `frame.energies` | Reference, electronic and method-resolved energies; `total_energy` selects an available total. |
| `EnergyObservation` | Records in `energies.observations` | `method`, `quantity_semantics` and `source_label` distinguish totals, correlation corrections and components. |
| `ThermalInformations` | `frame.thermal_informations` | Zero-point energy, thermal corrections, U/H/G, entropy, heat capacity and temperature/pressure/low-frequency conditions. |
| `ImplicitSolvation` | `solvent` | Recognized output solvent, model, atomic-radius scheme and dielectric constants. |
| `Status` | `status` | SCF convergence and normal termination as separate fields. |
| `GeometryOptimizationStatus` | `geometry_optimization_status` | Completion, energy/force/displacement criteria, thresholds and source-reported convergence. |

### Orbitals, populations, bond orders and total spin

| Container | Common access or relationship | Purpose |
| --------- | ----------------------------- | ------- |
| `MolecularOrbitals`, `MoleculeOrbital` | `frame.molecular_orbitals` and indexed orbitals | Alpha/beta energies, occupations, symmetry, total occupations, orbital basis and reported frontiers. |
| `NaturalAtomicOrbitals`, `NaturalAtomicOrbital` | Natural atomic orbital collection and records | Orbital/atom indices, angular momentum, type, occupancy and energy; collections can distinguish total/alpha/beta records. |
| `NaturalBondOrbitals`, `NaturalBondOrbital` | Natural bond orbital collection and records | Orbital type, bonding atoms, occupancy, energy and delocalization, with spin channels. |
| `ChargeSpinPopulations`, `AtomicPopulationSeries` | `frame.charge_spin_populations` and `populations[name]` | Extensible atom-aligned charge/spin series, retaining scheme, quantity, spin channel and source label. |
| `BondOrders`, `BondOrderMatrix` | `frame.bond_orders` and `bond_orders[name]` | Extensible atom-pair matrices retaining scheme, spin, print threshold and missing-entry reason. |
| `TotalSpin` | `frame.total_spin` | `spin_square` and `spin_quantum_number` for the spin-squared expectation and spin quantum number. |

Natural orbital types can represent independently supplied data. Consult each
format page to determine whether its reader populates them automatically.

### Vibrations, response and other properties

| Container | Common access or relationship | Purpose |
| --------- | ----------------------------- | ------- |
| `Vibrations`, `Vibration` | `frame.vibrations` and indexed modes | Frequencies, masses, force constants, IR, mode vectors, Raman activities, depolarization and mode/method identity. |
| `Polarizability` | `frame.polarizability` | Dipole, polarizability, quadrupole/higher multipoles and electronic spatial extent. |
| `NMR`, `ShieldingTensor` | `frame.nmr` and `nmr.shielding_tensors` | Per-atom shielding tensors, isotropic/principal/anisotropic values and K/J coupling matrices with components. |
| `SinglePointProperties` | `frame.single_point_properties` | Vertical IP/EA, global electrophilicity, atom-resolved Fukui indices and specialized properties such as FOD. |
| `Dispersions` | Independent dispersion-property records | `C6AA` and `C8AA`, distinct from a requested empirical dispersion correction. |

Forces and Hessians reside in `frame.forces` and `frame.hessian`. Their
`*_atom_order`, `*_axis_order` and `*_orientation` fields describe array conventions.

### Electronic states and multireference results

| Container | Common access or relationship | Purpose |
| --------- | ----------------------------- | ------- |
| `ElectronicConfiguration` | Records in `ElectronicState.configurations` | Configurations, coefficients, weights, occupations and orbital indices. |
| `ElectronicState` | Records in `ElectronicStates.states` | Total/excitation energy, root, spin/symmetry, oscillator strength, transition dipole and explicit ground-state identity. |
| `ElectronicStates` | `frame.electronic_states` | States and views as a state table, transitions or electronic spectrum. |
| `MultireferenceResult` | `frame.multireference_result` | Output method, active space, states, corrections and diagnostics. |

### Generic representations and export

| Container | Purpose | Main fields |
| --------- | ------- | ----------- |
| `BaseDataClassWithUnit` | Shared foundation for unit-aware records, numeric and summary export. | `to_unitless_dump()`, `to_unitless_dump_with_unit_keys()` |
| `PropertyPoint`, `PropertySeries` | Individual observations and ordered series. | `value`, `uncertainty`, `points`, axis labels and metadata |
| `TensorProperty` | Generic tensor with a reference frame. | `tensor`, `frame`, `metadata` |
| `PropertyTable` | Equal-length column data retaining quantity units. | `columns`, `row_labels`, `column_names`, `row()` |
| `PropertyTransition` | Transition between two states. | Initial/final indices, energy, wavelength, oscillator strength and transition dipole |
| `SpectralBand`, `Spectrum` | Bands and collections of bands/transitions. | Center, width, intensity, assignment and `bands` / `transitions` |
| `PropertyBundle` | Properties grouped by data shape. | Scalar, vector, tensor, table, spectrum, series and transition mappings |

Containers exposing `to_property_bundle()` can produce generic views; vibrations
provide `to_ir_spectrum()` and `to_raman_spectrum()`. Export retains available
metadata without completing unreported quantities. See [source evidence and
serialization](behavior/serialization.md) for export options.

### Provenance and parsing status

| Record or type | Access or relationship | Purpose |
| -------------- | ---------------------- | ------- |
| `SourceSpan` | Frame or segment `source_span` | Byte, character and line ranges for locating source text. |
| `SourceSegmentEvidence` | File `source_segments` | Step provenance, tasks, status, original frame count and retained frames. |
| `ParserProvenance` | File `parser_provenance` | Parser/dependency versions and effective configuration. |
| `ParseDiagnostic` | File or frame diagnostics | Code, severity, scope, associated field and message. |
| `ParsePresence`, `ParseCompleteness` | `parse_presence`, `parse_completeness` | Content presence and parsing completeness, distinct from calculation convergence. |

These types live in `molop.io.base_models.source`. Use `capture_source_evidence=True`
to request additional provenance. See [source evidence and serialization](behavior/serialization.md)
for source tracing and export.

## Cross-QM semantic rules

### Calculation setup and results

Prefer `model_chemistry`, `resource_request` and request containers when reading
setup across programs. Flat `method`, `functional`, `basis_set` and `keywords`
fields are convenience projections. Gaussian R/U/RO reference types reside in
`spin_treatment`; functionals retain normalized dispersion suffixes alongside
the correction identifier, for example `functional="B3LYP-GD3BJ"` and
`dispersion_correction="GD3BJ"`.

Frequency, NMR, excited-state and explicit-solvent requests describe intended
tasks. Results require corresponding output data. Basis comparisons should
include roles and scopes in `basis_sets`, not just the main basis string.
ORCA `%maxcore` is memory per process; matching MB values alone do not establish
equivalent resource requests across programs.

### Units, atom order and orientation

| Quantity | Shared default unit | Normalization |
| -------- | ------------------- | ------------- |
| Coordinates | Å | Convert Bohr to Å, preserving source atom order. |
| Electronic/orbital energies | Hartree (Eh) | Keep total, reference and correlation corrections separate. |
| Thermal corrections/totals | kcal/mol | Interpret printed Hartree per particle before conversion to molar quantities. |
| Entropy/constant-volume heat capacity | cal/(mol·K) | Temperature K and pressure atm; missing conditions remain unknown. |
| Forces/Cartesian Hessian | Eh/Bohr and Eh/Bohr² | Convert gradients through `F = -dE/dR`; printed forces are not negated again. |
| Frequency/reduced mass/force constant | cm⁻¹, amu and mdyne/Å | Preserve mode identity and normalization separately. |
| IR intensity/Raman activity | km/mol and Å⁴/amu | Preserve their distinct quantity names. |
| Dipole/polarizability | Debye and Bohr³ | Interpret tensor component ordering together with orientation. |
| Shielding/K and J coupling | ppm and Hz | Shielding differs from chemical shift; K and J remain distinct. |
| Calculation time | s | `cpu_time` and `wall_time` are independent; `running_time` aliases `wall_time`. |

Use `.to(...)` or `.m_as(...)` for another compatible unit. Equal units do not
establish equal reference frames. `forces` has shape `(N, 3)` and the Hessian
`(3N, 3N)`; compare atom order, axis order and input/standard/source orientation.
Unknown orientation stays unknown; tensors are not all rotated automatically.
Frames validate population lengths, bond-order matrix sizes and displacement atom
counts against their atoms. Each optional per-mode array is either absent or aligned
with all stored frequencies; incomplete tables must not shift values to another mode.

### Calculation time

Calculation files and frames expose `cpu_time` for CPU usage and `wall_time` for
elapsed wall-clock duration. `running_time` is an alias of `wall_time` in every
format. Missing wall time remains `None`, even when CPU time is available; CPU
time is never added to elapsed time or divided by a requested core count.

`cpu_time_source` distinguishes printed CPU durations (`observed`) from fallback
estimates (`estimated`); `None` means the source is unknown. For ORCA, missing CPU
time is estimated as wall time × process count, using the printed runtime MPI count
first, then `%pal nprocs` or `PALn`. A complete printed input without a parallel
request uses ORCA's serial default of one process. Unknown wall time or process
count leaves CPU time unknown. `cpu_time_num_processes` records the multiplier.
This estimate assumes every process uses the full elapsed duration; it does not
measure CPU utilization or exclude waiting time.

File totals sum each quantity independently across the selected source segments.
If any selected segment lacks a quantity, its file total stays unknown. With
`only_last_frame=True`, file timing follows the retained source segment, while
frame timing remains local to the stored frame. Gaussian file scopes prefer
printed job totals and frame scopes prefer local Link timings; these alternative
records are not added together. File summary columns distinguish CPU and wall time.
Any estimated segment makes the CPU total `estimated`; the process multiplier is
retained only if it is the same across all selected segments. File summaries also
include `CPUTimeSource` and `CPUTimeNumProcesses` when available.

### Energies and thermochemistry

`energies.total_energy` selects an available total in container-defined order:
`electronic_energy`, then CCSD(T), CCSD, MP5 through MP2, then `reference_energy`.
Method-specific fields remain available. Check method and source labels when
comparing protocols; this convenience field alone does not establish equal methods.

Thermochemistry uses these definitions, where `Eel` is electronic energy:

- `U_0 = Eel + ZPVE`;
- `U_T = Eel + TCE`;
- `H_T = Eel + TCH`;
- `G_T = Eel + TCG`.

Gaussian logs distinguish corrections from electronic-plus-thermal totals.
FCHK corrections derived by subtraction retain equations in `derivations`.
An incremental ORCA enthalpy correction alone is not `TCH`; its printed entropy
energy term is divided by temperature to obtain entropy. xTB retains available
thermal summaries and RRHO/QRRHO conditions.

Compare free energies alongside temperature, pressure, `entropy_model`, frequency
scaling, mode-exclusion cutoff, `rotor_reference_frequency` and imaginary-frequency
treatment. The HO/rotor interpolation reference differs from a mode-exclusion
threshold. FCHK totals alone do not establish an entropy model.

### Orbitals and spin

Alpha/beta arrays align by orbital index. Missing or shorter channels are padded
without discarding the longer channel. Padding does not establish shared spectra;
energies are shared only when a restricted shared reference is established.

`orbital_basis` distinguishes shared spatial, spin-resolved and unknown orbitals.
`total_occupancies` retains spatial-orbital totals; `occupancy_origin` identifies
observed, electron-count-derived, restricted or high-spin occupations. Closed-shell
totals can split equally; integer occupations with high-spin evidence can split
by spin. Other uncertain spin occupations remain unknown. Omitted xTB ranges are
not unoccupied orbitals; fractional occupations alone do not establish HOMO/LUMO
or SOMO identity. Reported frontiers reside in `reported_homo_index` and `reported_lumo_index`.

`spin_square` represents `<S²>`. The derived
`spin_quantum_number = sqrt(<S²> + 1/4) - 1/2` is an effective spin that can differ
from ideal spin for a contaminated wavefunction. Input multiplicity does not
replace an observed `<S²>`.

### Populations and bond orders

Enumerate populations with `population_items()` / `population_names` and bond
orders with `bond_order_items()` / `bond_order_names`. Compare `scheme`, population
`quantity` and `spin_channel`; different schemes and total/alpha/beta results
remain separate.

`populations["mulliken_charges"].values` is atom-aligned;
`bond_orders["mayer_bond_order"].values` is an atom-pair matrix. Wiberg, Mayer, MO
and other schemes use the same record type with open scheme names. Unprinted
pairs in thresholded tables remain unknown; `print_threshold` and `missing_reason`
explain omissions. Printed zero and permitted negative values remain unchanged.
Embedding bond orders in RDKit annotates existing chemical bonds without changing
the molecular graph. Values for nonbonded atom pairs remain in the frame matrix.

### Vibrations, IR and Raman

`mode_indices` are internal zero-based identifiers; `source_mode_indices` retain
printed mode numbers, with `None` for modes without a printed number. Frequencies, spectra
and displacements align by mode identity; selection preserves both identifiers.

Mode arrays use `(mode, atom, xyz)` order. Source `normalization`, `mass_weighting`,
`orientation` and atom order stay attached. Vectors are not rescaled to force
agreement across programs. The Å displacement unit is a storage convention;
normalized-vector magnitudes are not experimental absolute amplitudes.

`IR_intensities` stores IR intensity and `raman_activities` Raman activity.
Gaussian P/U ratios reside in `raman_depolarization_ratios_polarized` and
`raman_depolarization_ratios_unpolarized`; unspecified polarization uses the
unsuffixed field. Ratios are dimensionless. Activity excludes experimental laser
and temperature factors. xTB `frequency_method` and `raman_activity_method` can
differ; unprinted vectors and depolarization remain absent.

### NMR and electric response

`ShieldingTensor.atom_index` is a zero-based source-order index;
`coupling_atom_indices` identifies coupling-matrix axes. K, J and FC/SD/PSO/DSO
components remain separate. Without isotope information, FCHK K does not supply J.

Shielding tensors and printed isotropic/principal/anisotropic values are retained
separately; absent scalars can be derived from the tensor. `anisotropy_convention`
identifies the convention: shared fallback uses Haeberlen, while a format can
preserve its program's convention. A chemical shift requires a reference;
unknown tensor orientation also prevents direct component-wise comparison.
Dipole/polarizability units are converted by source section, with packed-component
ordering documented on the format page.

### Electronic states, status and missing values

`ElectronicState.state_index` is an internal unique non-negative identifier;
`root` retains the program's root number. `is_ground_state=True` / `False` / `None`
means explicit ground state, excited state or unknown. Neither list position nor
root zero alone establishes a ground state. Excited-only collections still export
transitions; an unstored initial ground state can have no internal index while
source labels remain attached.

Normal termination, SCF convergence and optimization completion are independent.
Unknown status uses `None`, distinct from `False`. A single energy or termination
line does not establish that every task succeeded.

Unknown orbital energies, bond orders and coupling cells use `np.nan`; unknown
occupations and labels use `None`. Fields allowing missing values export missing
entries as JSON `null`; `array_mode="ndarray"` can retain array NaN. Printed zeros
remain zero. An absent container, an empty array and a missing array element are
different cases when comparing or exporting results.

SMILES 2D coordinates, graphs recovered from coordinates and fakeG text are
generated data. Conversion uses existing data without performing a new quantum
calculation. Shared containers unify access and interpretation while preserving
program-specific information and physical conventions.
