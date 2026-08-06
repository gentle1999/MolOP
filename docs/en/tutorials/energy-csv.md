# Export energy and thermochemistry CSV

Extract final-frame energy, zero-point energy, enthalpy, and Gibbs free energy from each file.

## Python

```python
from molop import AutoParser

batch = AutoParser("water_mp2.out", n_jobs=1)
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
print(energy_table.drop(columns="DiskStorage.FilePath").head().to_string(index=False))
```

## Output

??? example "CSV preview"

    ```text
    Calc Parameter.Method  Energy.total_energy.hartree  Thermal.ZPVE.kilocalorie / mole  Thermal.H_T.kilocalorie / mole  Thermal.G_T.kilocalorie / mole
                      MP2                  -74.999375                              NaN                             NaN                             NaN
    ```

The file header contains the six requested columns. The values above are from the bundled
`water_mp2.out` sample; its single-point job has no frequency or thermochemistry section, so the
`Thermal.*` cells are empty in `energies.csv`, not zero.

## Keep only files with thermochemistry

```python
thermal_batch = batch.filter_state("thermal")
thermal_table = thermal_batch.to_summary_df(
    brief=False,
    flatten_columns=True,
).reindex(columns=wanted)
print(len(batch), len(thermal_batch))
```

??? example "Output"

    ```text
    1 0
    ```

The two counts are all parsed files and files with thermochemistry on at least one frame. Replace
`water_mp2.out` with `results/*.log` to process a directory of calculations.

## Note

Flat column names include units. Temperature, pressure, standard state, and theory level still need
to be interpreted from the original calculation settings.
