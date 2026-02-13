# Model Field Map

The following table maps common chemical data to their locations in MolOP's data models (`File` and `Frame` objects).

| Data | Location | Type |
|------|----------|------|
| **Structure** | `frame.atoms`, `frame.coords` | `list[int]`, `NumpyQuantity` |
| **Bonds** | `frame.bonds` | `list` |
| **SMILES** | `frame.to_SMILES()` | `str` |
| **Total Energy** | `frame.energies.total_energy` | `PlainQuantity` |
| **SCF Energy** | `frame.energies.scf_energy` | `PlainQuantity` |
| **Frequencies** | `frame.vibrations.frequencies` | `NumpyQuantity` |
| **Imaginary Freqs** | `frame.vibrations.num_imaginary` | `int` |
| **Orbitals** | `frame.molecular_orbitals` | `MolecularOrbitals` |
| **Charges/Spins** | `frame.charge_spin_populations` | `ChargeSpinPopulations` |
| **Optimization Status** | `frame.geometry_optimization_status` | `GeometryOptimizationStatus` |
| **Calculation Status** | `file.status` (last frame) or `frame.status` | `Status` |
| **QM Software** | `file.qm_software` | `str` |
