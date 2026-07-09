# Model Field Map

This page maps user-visible computational chemistry data to MolOP file and frame
objects. It is written from the public API point of view: use these fields in
analysis code, notebooks, and downstream converters.

## Gaussian Output Contract

`AutoParser(..., parser_detection="g16log")` returns a batch. Each item in the
batch is a Gaussian output file object, and each file contains frame objects.
The parser stores raw Gaussian output as normalized model fields; it does not
promise byte-for-byte preservation of the original log.

### File-Level Fields

Use file-level fields for job-wide metadata and values finalized from the last
frame.

| Data | Field | Notes |
| ---- | ----- | ----- |
| Software | `file.qm_software`, `file.qm_software_version` | Usually `"Gaussian"` plus the Gaussian revision string. |
| Route text | `file.keywords` | Normalized route line text. |
| Structured route | `file.semantic_route`, `file.model_chemistry`, `file.task_requests` | Preferred structured representation for method family, basis, job types, solvation, dispersion, and related route semantics. |
| Legacy route projections | `file.method`, `file.basis_set`, `file.functional` | Compatibility fields projected from structured route data when possible. |
| Charge and multiplicity | `file.charge`, `file.multiplicity` | File-level values finalized from parsed frames. |
| Running time | `file.running_time` | Accumulated `pint` quantity when Gaussian timing records are present. |
| Status | `file.status` | File-level status finalized from the last parsed frame. |

### Frame-Level Fields

Use frame-level fields for structures and per-step QM results.

| Data | Field | Notes |
| ---- | ----- | ----- |
| Frame identity | `frame.frame_id` | Zero-based frame index inside the parsed file. |
| Structure | `frame.atoms`, `frame.atom_symbols`, `frame.coords` | Atomic numbers, symbols, and input-orientation coordinates. |
| Standard orientation | `frame.standard_coords`, `frame.standard_orientation_transformation_matrix` | Present when Gaussian prints standard orientation or when it can be reconstructed. |
| Energy data | `frame.energies` | `Energies` object with `reference_energy`, `electronic_energy`, post-HF energies, and computed `total_energy`. |
| Thermochemistry | `frame.thermal_informations` | `ThermalInformations` object with ZPVE, thermal corrections, thermodynamic energies, entropy, heat capacity, mass, and rotational/vibrational metadata. |
| Vibrations | `frame.vibrations` | Frequencies, reduced masses, force constants, IR intensities, mode vectors, and imaginary-mode counts. |
| Molecular orbitals | `frame.molecular_orbitals` | Orbital energies, occupancies, symmetries, and derived frontier-orbital quantities. |
| Populations | `frame.charge_spin_populations` | Mulliken, Lowdin, Hirshfeld, CM5, NPA, and spin population fields when printed and recognized. |
| Response properties | `frame.polarizability` | Dipole, polarizability tensor/scalars, electronic spatial extent, and multipoles where present. |
| Forces and Hessian | `frame.forces`, `frame.hessian` | Cartesian arrays with normalized units. |
| Optimization | `frame.geometry_optimization_status` | Berny convergence values and boolean convergence flags. |
| Status | `frame.status`, `frame.is_error`, `frame.is_normal`, `frame.is_TS`, `frame.is_optimized` | Public status helpers for common workflow filters. |
| Running time | `frame.running_time` | Per-frame timing when present. |

### Energy Selection

`frame.energies.total_energy` is a computed field, not an input field. It picks
the most specific available energy in this order:

`ccsd_energy -> mp5_energy -> mp4_energy -> mp3_energy -> mp2_energy -> electronic_energy -> reference_energy`

Use method-specific fields when the exact source matters. Use `total_energy` for
summary tables and filtering when “best available scalar energy” is sufficient.

### Summary Tables

`file.to_summary_df()` returns one row per frame. With the default
`brief=True`, the stable columns cover storage, charge/multiplicity, structure
summary, route metadata, environment, and status:

| Column group | Examples |
| ------------ | -------- |
| `DiskStorage` | `FilePath`, `FileFormat` |
| `General` | `Charge`, `Multiplicity`, `CanonicalSMILES`, `NumAtoms`, `FrameID` |
| `Calc Parameter` | `Software`, `Version`, `Method`, `BasisSet`, `Functional`, `Keywords` |
| `Environment` | `SolventModel`, `Solvent`, `Temperature (...)`, `Pressure (...)` |
| `Status` | `IsError`, `IsNormal`, `IsTS`, `IsOptimized` |

Use `brief=False` to add result-heavy columns such as `Energy`,
`Thermal`, `GeometryOptimizationStatus`, and `Vibration`.

For batches, `batch.to_summary_df(frameIDs="all", flatten_columns=True)` returns
one row per selected frame and uses dot-separated column names such as
`General.FrameID` and `Status.IsError`.

### fakeG Rendering

`file.format_transform("fakeg")` follows the general transform default
`frameID=-1`, so it renders the last frame unless another selector is passed.
Use `file.format_transform("fakeg", frameID="all")` or `file.render_fakeg()` for
full-file Gaussian-like output.

`fakeg` output is semantic and normalized. It is suitable for inspection,
compatibility tests, and reparse checks, but it is not a byte-for-byte Gaussian
log reproduction.

### Compatibility and Internal Data

The public API keeps some legacy flat fields for convenience:

| Compatibility field | Preferred structured field |
| ------------------- | -------------------------- |
| `method`, `basis_set`, `functional` | `model_chemistry` |
| `keywords` | `semantic_route` / `model_chemistry` / task request containers |
| `is_error`, `is_normal`, `is_TS`, `is_optimized` | `status`, `vibrations`, `geometry_optimization_status` |

The derived data used to rebuild Gaussian-like text is internal. It is not
serialized by `model_dump()`, and user code should rely on file/frame fields
instead of rendering implementation details.

## Common Field Map

| Data | Location | Type |
| ---- | -------- | ---- |
| Structure | `frame.atoms`, `frame.coords` | `list[int]`, `NumpyQuantity` |
| Bonds | `frame.bonds` | `list` |
| SMILES | `frame.to_SMILES()` | `str` |
| Total energy | `frame.energies.total_energy` | `PlainQuantity | None` |
| Reference energy | `frame.energies.reference_energy` | `PlainQuantity | None` |
| Electronic energy | `frame.energies.electronic_energy` | `PlainQuantity | None` |
| Frequencies | `frame.vibrations.frequencies` | `NumpyQuantity` |
| Imaginary frequencies | `frame.vibrations.num_imaginary` | `int` |
| Orbitals | `frame.molecular_orbitals` | `MolecularOrbitals | None` |
| Charges/spins | `frame.charge_spin_populations` | `ChargeSpinPopulations | None` |
| Optimization status | `frame.geometry_optimization_status` | `GeometryOptimizationStatus | None` |
| Calculation status | `file.status` or `frame.status` | `Status | None` |
| QM software | `file.qm_software`, `frame.qm_software` | `str` |
