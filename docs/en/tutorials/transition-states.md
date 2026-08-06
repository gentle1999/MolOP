# Transition-state analysis

Use this workflow for frequency outputs that contain a transition-state candidate. MolOP can select
candidate frames, inspect imaginary modes, and generate geometry-based endpoint candidates. It does
not optimize those structures or certify a reaction path.

## Select candidate frames

```python
from molop import AutoParser

batch = AutoParser("results/*.log")
normal = batch.filter_state("normal", n_jobs=-1)
ts_batch = normal.filter_state("ts", n_jobs=-1)

print("parsed:", len(batch))
print("normal:", len(normal))
print("transition-state files:", len(ts_batch))
```

??? example "Selection result (input-dependent)"

    ```text
    parsed: <count>
    normal: <count>
    transition-state files: <count>
    ```

The counts depend on the input directory. A file can contain several frames, so select the actual TS
frame rather than assuming that `parsed_file[-1]` is the one selected by `filter_state("ts")`:

```python
for parsed_file in ts_batch:
    for frame in parsed_file:
        if frame.is_TS:
            imaginary = frame.vibrations.num_imaginary if frame.vibrations else None
            print(parsed_file.filename, frame.frame_id, imaginary)
```

??? example "Candidate-frame output shape (input-dependent)"

    ```text
    <source filename> <frame_id> <imaginary-mode count or None>
    ```

`is_TS` is a programmatic candidate decision. `None` means that the source did not provide enough
evidence and is different from `False`.

## Inspect the imaginary mode

```python
ts_frame = next(
    (frame for parsed_file in ts_batch for frame in parsed_file if frame.is_TS),
    None,
)
if ts_frame is None:
    raise ValueError("No transition-state frame was found in the input")
vibrations = ts_frame.vibrations

if vibrations is None:
    raise ValueError("The TS candidate has no structured frequency data")

print("imaginary modes:", vibrations.num_imaginary)
print("frequencies (cm^-1):", vibrations.frequencies.m_as("cm^-1"))
```

??? example "Imaginary-mode output (input-dependent)"

    ```text
    imaginary modes: <count>
    frequencies (cm^-1): <frequency array>
    ```

Review the frequency sign, displacement direction, and chemical interpretation. The number of
imaginary modes alone is not a reaction-path validation.

## Generate structures along the mode

`ts_vibration` displaces coordinates along the first vibration and returns candidate `Molecule`
objects. It requires `frame.is_TS` to be true and may omit geometries that violate basic crowding
checks:

```python
candidates = ts_frame.ts_vibration(ratio=1.75, steps=7)
print("candidate geometries:", len(candidates))
```

??? example "Imaginary-mode candidates (image output)"

    ![Structures generated and reconstructed along an imaginary mode](../../assets/examples/ts_imaginary_mode.svg)

    The image is generated from a real parsed frame with one imaginary mode. Each panel is a
    `ts_vibration(...)` candidate for which graph reconstruction succeeded. The two-dimensional
    layout compares connectivity; it is not an optimized reaction path.

## Infer endpoint candidates and bond changes

```python
try:
    reactant, product = ts_frame.possible_pre_post_ts(show_3D=True)
    print(reactant.GetNumAtoms(), product.GetNumAtoms())
except ValueError as exc:
    print("endpoint inference failed:", exc)

difference = ts_frame.to_diff_rdmol()
if difference is not None:
    print("difference graph bonds:", difference.GetNumBonds())
```

??? example "Endpoint candidates and virtual-bond difference (image output)"

    ![Reactant candidate, virtual-bond difference graph, and product candidate](../../assets/examples/ts_endpoints_difference.svg)

    Red highlights identify bonds that change between the endpoint candidates. The center panel
    comes from `to_diff_rdmol()`. Its zero-order bonds are MolOP's machine-readable difference
    markers; the generator promotes them to visible lines only in the drawing copy.

`possible_pre_post_ts` returns geometry-based endpoint candidates, not optimized reactant and product
structures. `to_diff_rdmol` can return `None` when no supported bond-breaking difference is inferred.
Review atom mapping, connectivity, and the original calculation before using either result in an
automated reaction workflow.

## Related operations

- Use `filter_state("opt")` and `parsed_file.closest_optimized_frame` for optimization trajectories.
- Use `frame.vibrate(...)` when a specific vibration, rather than the first one, is required.
- Export selected structures with [Convert and export](../guides/conversion.md).
