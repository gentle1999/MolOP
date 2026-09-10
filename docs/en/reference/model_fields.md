# Find fields by scientific property

Start with `frame = AutoParser(path)[0][-1]`, then locate a public field by the result you need.

## Quick index

| Result | Presence check | Main fields | Common units |
| --- | --- | --- | --- |
| Structure | `bool(frame.atoms)` | `atoms`, `coords`, `rdmol` | angstrom |
| Charge/multiplicity | Read directly | `charge`, `multiplicity` | dimensionless |
| Energy | `frame.energies is not None` | `energies.total_energy`, `reference_energy`, `mp2_energy`, `ccsd_t_energy` | hartree |
| Thermochemistry | `frame.thermal_informations is not None` | `ZPVE`, `U_T`, `H_T`, `G_T`, `S`, `C_V` | field-specific |
| Frequencies | `frame.vibrations is not None` | `frequencies`, `num_imaginary`, `vibration_modes` | cm^-1, angstrom |
| Molecular orbitals | `frame.molecular_orbitals is not None` | `HOMO_energy`, `LUMO_energy`, `HOMO_LUMO_gap`, occupancies | hartree, convertible to eV |
| Atomic populations | `frame.charge_spin_populations is not None` | `population_names`, `populations[name].values` | dimensionless |
| Dipole/polarizability | `frame.polarizability is not None` | `dipole`, `polarizability_tensor`, `quadrupole` | debye, bohr^3 |
| NMR | `frame.nmr is not None` | `shielding_tensors`, `spin_spin_coupling_j/k` | ppm, Hz |
| Forces/Hessian | `frame.forces/hessian is not None` | `forces`, `hessian`, and axis metadata | hartree/bohr, etc. |
| Status | Read directly | `is_normal`, `is_error`, `is_optimized`, `is_TS` | `bool | None` |

`rdmol` is a defensive RDKit graph copy: mutating the object returned by `frame.rdmol` does not
mutate the frame's cache. Edit `atoms`, `coords`, or explicit topology fields instead; those edits
invalidate derived topology and SMILES caches.

## Unit-aware values

```python
energy = frame.energies.total_energy
print(energy.m_as("hartree"))

frequencies = frame.vibrations.frequencies.m_as("cm^-1")
coords = frame.coords.m_as("angstrom")
```

??? example "Output"

    For the bundled ORCA water sample, the first line is:

    ```text
    -74.999374598107
    ```

`.m_as(...)` returns a value or NumPy array in the requested unit. Keep the quantity object when
unit metadata must remain attached.

## Populations are an open set

```python
populations = frame.charge_spin_populations
if populations:
    for name, series in populations.population_items():
        print(name, series.scheme, series.quantity, series.values[:3])
```

??? example "Output"

    The bundled ORCA water sample reports two charge-population schemes:

    ```text
    mulliken_charges mulliken charge [-0.361722  0.180857  0.180866]
    lowdin_charges lowdin charge [-0.2507    0.125348  0.125352]
    ```

`ChargeSpinPopulations` has one `populations` mapping rather than a fixed attribute per scheme. Read
Mulliken charges with:

```python
values = populations["mulliken_charges"].values
```

Check `population_names` first because a source may provide only Lowdin, Hirshfeld, CM5, NPA, or ESP
data.

## `None` status

```python
if frame.is_normal is None:
    print("the source has insufficient termination evidence")
```

??? example "Output"

    ```text
    the source has insufficient termination evidence
    ```

For status fields, `None` means unknown and is not equivalent to `False`.

## Format capability

A public field does not imply that every format fills it. Check:

- [Gaussian log](formats/g16log.md)
- [Gaussian fchk](formats/g16fchk.md)
- [ORCA output](formats/orcaout.md)
- [xTB output](formats/xtbout.md)
- [Full format overview](format_support.md)

See [Read calculation results](../guides/results.md) for runnable examples.
