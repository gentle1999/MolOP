# Python API first steps

Learn the `batch -> file -> frame` access pattern used by most MolOP programs.

## Parse one file

```python
from molop import AutoParser

batch = AutoParser("calculation.log")
parsed_file = batch[0]
final_frame = parsed_file[-1]
```

`len(batch)` is the number of files added successfully. `len(parsed_file)` is the number of frames
in that file.
For the shared `water_mp2.out`, continue with:

```python
print(len(batch), len(parsed_file), parsed_file.detected_format_id)
print(len(final_frame.atoms), final_frame.coords.shape)
```

Output:

```text
1 1 orcaout
3 (3, 3)
```

## Parse multiple inputs

```python
batch = AutoParser([
    "reactants/*.log",
    "products/product.out",
    "structures/*.xyz",
])

for parsed_file in batch:
    frame = parsed_file[-1]
    print(parsed_file.filename, frame.charge, frame.multiplicity)
```

Each successfully parsed file produces one line, for example:

```text
water_mp2.out 0 1
```

Paths, glob patterns, and path lists can be mixed. Results are stably sorted by absolute path and
duplicate paths are parsed once.

## Check optional results

Different calculation types provide different fields. Check that a container exists before using
it:

```python
frame = batch[0][-1]

if frame.energies and frame.energies.total_energy is not None:
    print(frame.energies.total_energy.m_as("hartree"))

if frame.vibrations:
    print(frame.vibrations.num_imaginary)

if frame.charge_spin_populations:
    names = frame.charge_spin_populations.population_names
    print(names)
```

The shared MP2 example prints `-74.999374598107` from the energy branch. It has no frequency or
population section, so the other two branches print nothing. This is why optional containers are
checked first.

MolOP numerical results commonly carry Pint units. Use `.m_as("target unit")` to obtain a value in
a chosen unit.

## Build a table

```python
df = batch.to_summary_df(
    frame=-1,
    brief=False,
    flatten_columns=True,
)
df.to_csv("summary.csv", index=False)
```

Use `frame="all"` for every frame in every file. The default `frame=-1` selects only final frames.

## Next steps

- [Parse files](../guides/parsing.md)
- [Read calculation results](../guides/results.md)
- [Batch summaries](../guides/batch.md)
