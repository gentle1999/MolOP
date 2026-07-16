# Read calculation results

Find results by scientific property on the final frame. Formats and calculation types do not all
provide the same fields, so check each optional container first.

## Prepare

```python
from molop import AutoParser

frame = AutoParser("calculation.log")[0][-1]
```

Every snippet on this page prints or assigns its result. With the shared `water_mp2.out`, the two
core outputs are:

```python
frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
print(len(frame.atoms), frame.coords.shape)
print(frame.energies.total_energy.m_as("hartree"))
```

```text
3 (3, 3)
-74.999374598107
```

## Structure, charge, and multiplicity

```python
print(frame.atoms)          # atomic numbers
print(frame.coords)         # unit-aware (N, 3) coordinates
print(frame.charge)
print(frame.multiplicity)
```

The shared example has atomic numbers `[8, 1, 1]`, coordinate shape `(3, 3)`, charge `0`, and
multiplicity `1`.

Obtain a plain NumPy array with:

```python
coords_angstrom = frame.coords.m_as("angstrom")
```

## Final energy

```python
if frame.energies and frame.energies.total_energy is not None:
    energy_hartree = frame.energies.total_energy.m_as("hartree")
    print(energy_hartree)
```

The shared example prints `-74.999374598107`.

`total_energy` selects the highest-priority available total energy in the current container. Read
method-specific fields when reference, MP2, or CCSD(T) energies must remain distinct:

```python
energies = frame.energies
if energies:
    print(energies.reference_energy)
    print(energies.mp2_energy)
    print(energies.ccsd_t_energy)
```

Do not treat `total_energy` as one theory level that can be compared blindly across methods.

## Thermochemistry

```python
thermal = frame.thermal_informations
if thermal:
    print(thermal.ZPVE)
    print(thermal.H_T)
    print(thermal.G_T)
    print(thermal.S)
```

Common fields include zero-point energy `ZPVE`, corrections `TCE/TCH/TCG`, internal energies
`U_0/U_T`, enthalpy `H_T`, Gibbs free energy `G_T`, entropy `S`, and heat capacity `C_V`. Pint
quantities retain each field's unit.

## Frequencies, imaginary modes, and displacements

```python
vibrations = frame.vibrations
if vibrations:
    print(vibrations.frequencies.m_as("cm^-1"))
    print(vibrations.num_imaginary)
    if vibrations.vibration_modes:
        first_mode = vibrations.vibration_modes[0].m_as("angstrom")
```

Output is a frequency array, an integer imaginary-mode count, and an optional `(N, 3)` displacement
array. For a non-frequency job, the branch does not run; it does not print a zero frequency.

`num_imaginary` counts negative frequencies. For filtering, `frame.is_TS` is more appropriate
because it can also use task and frame status.

## Molecular orbitals

```python
orbitals = frame.molecular_orbitals
if orbitals:
    print(orbitals.HOMO_energy.m_as("eV"))
    print(orbitals.LUMO_energy.m_as("eV"))
    print(orbitals.HOMO_LUMO_gap.m_as("eV"))
    print(orbitals.alpha_occupancies)
```

Open-shell results can also provide `beta_energies` and `beta_occupancies`. MolOP does not invent a
coefficient matrix when the output prints energies but no coefficients.

## Atomic populations

```python
populations = frame.charge_spin_populations
if populations:
    print(populations.population_names)
    if "mulliken_charges" in populations:
        charges = populations["mulliken_charges"].values
        print(charges)
```

The current ORCA water population regression example outputs:

```text
['mulliken_charges', 'lowdin_charges', 'hirshfeld_charges', 'hirshfeld_spins']
[-0.642714, 0.321357, 0.321357]
```

Other files expose the names and values actually printed by their source.

Population schemes form an extensible mapping based on what the source prints. Names can include
`mulliken_charges`, `mulliken_spins`, `lowdin_charges`, `hirshfeld_charges`, `cm5_charges`,
`npa_charges`, or `esp_charges`. Do not assume one fixed name is always present.

## Dipole and polarizability

```python
response = frame.polarizability
if response:
    if response.dipole is not None:
        print(response.dipole.m_as("debye"))
    if response.polarizability_tensor is not None:
        print(response.polarizability_tensor.m_as("bohr**3"))
```

Depending on the format and printed section, the container can also include
`isotropic_polarizability`, `anisotropic_polarizability`, `quadrupole`, and higher multipoles.

## NMR

```python
nmr = frame.nmr
if nmr:
    for shielding in nmr.shielding_tensors:
        print(shielding.atom_index, shielding.isotropic.m_as("ppm"))

    if nmr.spin_spin_coupling_j is not None:
        coupling_hz = nmr.spin_spin_coupling_j.m_as("Hz")
```

An fchk file may contain only reduced `spin_spin_coupling_k` values without isotope information
needed to reconstruct J. Consult the format page for exact limitations.

## Optimization and termination status

```python
print(frame.is_normal)
print(frame.is_error)
print(frame.is_optimized)
print(frame.is_TS)

if frame.status:
    print(frame.status.scf_converged)
if frame.geometry_optimization_status:
    print(frame.geometry_optimization_status.geometry_optimized)
```

`None` means the source did not provide enough evidence; do not force it to mean `False`.

## Next steps

- [Find fields by scientific property](../reference/model_fields.md)
- [Format overview](../reference/format_support.md)
- [Batch summaries](batch.md)
