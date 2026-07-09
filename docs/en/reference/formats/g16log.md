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

MolOP reads Gaussian output files and extracts the data most often needed for
computational chemistry post-processing: structures, energies, thermochemistry,
vibrations, orbitals, populations, gradients, response properties, archive
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
| <!-- feature-area:Molecular orbitals and population analysis -->Molecular orbitals and population analysis | Partial | Molecular orbital energies and occupancies, charge and spin populations, electronic spatial extent, multipole data, and early polarizability values where present. | Orbital coefficient matrices and every population-analysis flavor are not claimed. |
| <!-- feature-area:Vibrational frequencies and IR intensities -->Vibrational frequencies and IR intensities | Partial | Frequencies, reduced masses, force constants, IR intensities, per-mode displacement vectors, and imaginary-mode flags. | Raman, VCD, and other spectrum variants are not advertised unless explicitly covered later. |
| <!-- feature-area:Thermochemistry -->Thermochemistry | Partial | Temperature, pressure, molecular mass, moments of inertia, rotational symmetry number, rotational/vibrational temperatures, rotational constants, ZPVE, thermal energy, enthalpy, Gibbs free energy, entropy, and heat capacity. | Coverage is mainly tied to frequency-style Gaussian outputs and maintained examples. |
| <!-- feature-area:Dipole and polarizability -->Dipole and polarizability | Partial | Dipole and polarizability values from covered Gaussian response-property sections. | Only example-backed response fields are advertised; not every response-property print variant is guaranteed structured. |
| <!-- feature-area:Cartesian gradients and forces -->Cartesian gradients and forces | Partial | Cartesian force arrays from Gaussian force sections. | Only the normalized force array is advertised; auxiliary force diagnostics are not separately recorded. |
| <!-- feature-area:Cartesian Hessian -->Cartesian Hessian | Partial | Cartesian Hessian data from Gaussian second-derivative sections. | The contract is the normalized Hessian field, not every printed second-derivative diagnostic. |
| <!-- feature-area:Geometry optimization convergence -->Geometry optimization convergence | Partial | Berny optimization summaries, including convergence thresholds, force/displacement values, energy change, and optimized-state flags. | Optimizer diagnostics outside tested examples may remain unstructured. |
| <!-- feature-area:Gaussian archive section -->Gaussian archive section | Partial | Archive-tail metadata, coordinates, energies, thermochemistry, polarizability, and Hessian fallback or augmentation fields. | Main-log fields take precedence where applicable; archive-tail data does not invent live status, temperature, or pressure. |
| <!-- feature-area:Termination status -->Termination status | Partial | Gaussian termination information and file-level status finalization from the last parsed frame. | Status remains an aggregate signal; detailed failure taxonomy is outside the current contract. |
| <!-- feature-area:CPU and elapsed time -->CPU and elapsed time | Example-covered | Job CPU and elapsed-time style records accumulated into the running-time value. | Per-link timing rows are not exposed as separate timing records. |
| <!-- feature-area:Link1 multi-step jobs -->Link1 multi-step jobs | Partial | Link1 section metadata is propagated to later frames so multi-step Gaussian jobs keep the expected per-frame context. | Low-level link boundary rows are not exposed as user-facing records. |
| <!-- feature-area:Registry conversion -->Registry conversion | Supported | Parsed Gaussian output can be converted to coordinate and graph formats through the registry. | Conversion quality depends on successful structure and graph recovery. |
