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

??? example "Output"

    ```text
    3 (3, 3)
    -74.999374598107
    ```

## Structure, charge, and multiplicity

```python
print(frame.atoms)
print(frame.coords.shape)
print(frame.charge, frame.multiplicity)
```

??? example "Output"

    ```text
    [8, 1, 1]
    (3, 3)
    0 1
    ```

Obtain a plain NumPy array with:

```python
coords_angstrom = frame.coords.m_as("angstrom")
print(coords_angstrom.shape)
```

??? example "Output"

    ```text
    (3, 3)
    ```

## Final energy

```python
if frame.energies and frame.energies.total_energy is not None:
    energy_hartree = frame.energies.total_energy.m_as("hartree")
    print(energy_hartree)
```

??? example "Output"

    ```text
    -74.999374598107
    ```

`total_energy` selects the highest-priority available total energy in the current container. Read
method-specific fields when reference, MP2, or CCSD(T) energies must remain distinct:

```python
energies = frame.energies
if energies:
    print("reference:", energies.reference_energy.m_as("hartree"))
    print("MP2:", energies.mp2_energy.m_as("hartree"))
    print("CCSD(T):", energies.ccsd_t_energy)
```

??? example "Output"

    ```text
    reference: -74.96357424008319
    MP2: -74.999374598
    CCSD(T): None
    ```

Do not treat `total_energy` as one theory level that can be compared blindly across methods.

## Thermochemistry

```python
thermal = frame.thermal_informations
print("thermal available:", thermal is not None)
if thermal:
    print(thermal.ZPVE)
    print(thermal.H_T)
    print(thermal.G_T)
    print(thermal.S)
```

The single-point sample has no thermochemistry section:

??? example "Output"

    ```text
    thermal available: False
    ```

Common fields include zero-point energy `ZPVE`, corrections `TCE/TCH/TCG`, internal energies
`U_0/U_T`, enthalpy `H_T`, Gibbs free energy `G_T`, entropy `S`, and heat capacity `C_V`. Pint
quantities retain each field's unit.

## Frequencies, imaginary modes, and displacements

```python
vibrations = frame.vibrations
print("vibrations available:", vibrations is not None)
if vibrations:
    print(vibrations.frequencies.m_as("cm^-1"))
    print(vibrations.num_imaginary)
    if vibrations.vibration_modes:
        first_mode = vibrations.vibration_modes[0].m_as("angstrom")
```

The single-point sample prints:

??? example "Output"

    ```text
    vibrations available: False
    ```

For a frequency job, the same code continues with a frequency array, an integer imaginary-mode count,
and an optional `(N, 3)` displacement array. It does not synthesize a zero frequency.

`num_imaginary` counts negative frequencies. For filtering, `frame.is_TS` is more appropriate
because it can also use task and frame status.

## Molecular orbitals

```python
orbitals = frame.molecular_orbitals
print("molecular orbitals available:", orbitals is not None)
if orbitals:
    print(orbitals.HOMO_energy.m_as("eV"))
    print(orbitals.LUMO_energy.m_as("eV"))
    print(orbitals.HOMO_LUMO_gap.m_as("eV"))
    print(orbitals.alpha_occupancies)
```

The shared water sample has no structured orbital container:

??? example "Output"

    ```text
    molecular orbitals available: False
    ```

Open-shell results can also provide `beta_energies` and `beta_occupancies`. MolOP does not invent a
coefficient matrix when the output prints energies but no coefficients.

## Atomic populations

```python
populations = frame.charge_spin_populations
if populations:
    print(populations.population_names)
    for name in populations.population_names:
        print(name, populations[name].values)
else:
    print("populations available: False")
```

The shared ORCA sample outputs:

??? example "Output"

    ```text
    ['mulliken_charges', 'lowdin_charges']
    mulliken_charges [-0.361722  0.180857  0.180866]
    lowdin_charges [-0.2507    0.125348  0.125352]
    ```

Other files expose the names and values actually printed by their source.

Population schemes form an extensible mapping based on what the source prints. Names can include
`mulliken_charges`, `mulliken_spins`, `lowdin_charges`, `hirshfeld_charges`, `cm5_charges`,
`npa_charges`, or `esp_charges`. Do not assume one fixed name is always present.

## Dipole and polarizability

```python
response = frame.polarizability
print("response available:", response is not None)
if response:
    if response.dipole is not None:
        print(response.dipole.m_as("debye"))
    if response.polarizability_tensor is not None:
        print(response.polarizability_tensor.m_as("bohr**3"))
```

The shared sample provides a dipole but no polarizability tensor:

??? example "Output"

    ```text
    response available: True
    [ 0.1502394  -0.11206605 -0.6497661 ]
    ```

Depending on the format and printed section, the container can also include
`isotropic_polarizability`, `anisotropic_polarizability`, `quadrupole`, and higher multipoles.

## NMR

```python
nmr = frame.nmr
print("NMR available:", nmr is not None)
if nmr:
    for shielding in nmr.shielding_tensors:
        print(shielding.atom_index, shielding.isotropic.m_as("ppm"))

    if nmr.spin_spin_coupling_j is not None:
        coupling_hz = nmr.spin_spin_coupling_j.m_as("Hz")
```

The shared water sample has no NMR section:

??? example "Output"

    ```text
    NMR available: False
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

The single-point sample prints:

??? example "Output"

    ```text
    True
    False
    False
    False
    ```

`None` means the source did not provide enough evidence; do not force it to mean `False`.

## Next steps

- [Find fields by scientific property](../reference/model_fields.md)
- [Format overview](../reference/format_support.md)
- [Batch summaries](batch.md)
