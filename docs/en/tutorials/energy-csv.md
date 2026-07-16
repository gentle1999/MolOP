# Export energy and thermochemistry CSV

Extract final-frame energy, zero-point energy, enthalpy, and Gibbs free energy from each file.

## Python

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
full = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)

wanted = [
    "DiskStorage.FilePath",
    "Calc Parameter.Method",
    "Energy.total_energy.hartree",
    "Thermal.ZPVE.kilocalorie / mole",
    "Thermal.H_T.kilocalorie / mole",
    "Thermal.G_T.kilocalorie / mole",
]
energy_table = full.reindex(columns=wanted)
energy_table.to_csv("energies.csv", index=False)
print(energy_table.head().to_string(index=False))
```

## Output

```text
energies.csv
```

The header contains the six requested columns. Single-point jobs without a frequency or
thermochemistry section have empty `Thermal.*` cells rather than zero.

## Keep only files with thermochemistry

```python
thermal_batch = batch.filter_state("thermal")
thermal_table = thermal_batch.to_summary_df(
    brief=False,
    flatten_columns=True,
).reindex(columns=wanted)
print(len(batch), len(thermal_batch))
```

The two output counts are all parsed files and files with thermochemistry on at least one frame.

## Note

Flat column names include units. Temperature, pressure, standard state, and theory level still need
to be interpreted from the original calculation settings.
