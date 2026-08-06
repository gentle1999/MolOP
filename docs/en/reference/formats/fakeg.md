# Gaussian-like Renderer

<!-- format-support:fakeg -->

| Item | Value |
| ---- | ----- |
| Format ID | `fakeg` |
| Extensions | `.fakeg` |
| Read | No |
| Write | Yes |
| Registry role | File writer |
| Data level | Parsed Gaussian output data |

`fakeg` renders parsed Gaussian output data as Gaussian-like text for human
inspection, compatibility testing, and downstream-tool handoff. It is not a
byte-for-byte Gaussian log reproducer.

`format_transform("fakeg")` follows the general transform default `frame=-1`;
pass `frame="all"` when full-file Gaussian-like rendering is required.

```python
from molop import AutoParser

gaussian_log = AutoParser(
    "calculation.log",
    parser_detection="g16log",
    n_jobs=1,
)[0]
rendered = gaussian_log.format_transform("fakeg", frame="all")
print(type(rendered).__name__)
```

??? example "Output type"

    ```text
    str
    ```

`G16LogFile.render_fakeg()` defaults to `frame="all"`, unlike `format_transform()`.

| Capability | Support | Rendered content | Boundaries |
| ---------- | ------- | ---------------- | ---------- |
| <!-- feature-area:File-level Gaussian-like writer -->File-level Gaussian-like writer | Partial | `.fakeg` files from parsed Gaussian output data. | File-level rendering only; frame-level fakeG writing is not provided, the selected frames follow `frame`, and byte-for-byte Gaussian log reproduction is not claimed. |
| <!-- feature-area:Structure and SCF energy rendering -->Structure and SCF energy rendering | Partial | Normalized orientation and SCF-cycle text from coordinate and energy fields. | Rendered sections are semantic Gaussian-like text, not original log text. |
| <!-- feature-area:Vibrational frequency rendering -->Vibrational frequency rendering | Partial | Frequency, reduced-mass, force-constant, IR-intensity, and per-mode displacement sections. | Current coverage is limited to frequency and IR-related fields; Raman/VCD sections are not advertised. |
| <!-- feature-area:Thermochemistry rendering -->Thermochemistry rendering | Partial | Temperature, correction, energy, entropy, heat-capacity, mass, inertia, and rotational/vibrational metadata. | This is a normalized thermochemistry summary, not a complete Gaussian thermochemistry pretty-printer. |
| <!-- feature-area:Reparseable frequency and thermochemistry output -->Reparseable frequency and thermochemistry output | Supported | Rendered frequency and thermochemistry content can be parsed back into frame fields. | The round trip proves supported semantic fields, not byte-for-byte equivalence with Gaussian output. |

## fakeG Capability Matrix

The matrix distinguishes structured field-driven reconstruction from component `raw_text`
fallback. Only field-driven reconstruction can synthesize a section without the original
Gaussian log.

| fakeG capability | Data source | Structured reconstruction | Reparse | Source fidelity | Current boundary |
| ---------------- | ----------- | ------------------------- | ------- | --------------- | ---------------- |
| Registry file writer | Coordinate-bearing file model | Supported | Partial | Unsupported | Registered as a `.fakeg` file writer; rich output depends on Gaussian result fields being present on input frames. |
| Registry frame writer | Single frame | Unsupported | Not applicable | Not applicable | `frame.format_transform("fakeg")` fails because no fakeG frame writer is registered. |
| `frame.render_fakeg()` | G16Log frame | Supported | Partial | Unsupported | The model method can render one frame directly; this is distinct from registry frame writing. |
| Frame selection | `frame` selector | Supported | Not applicable | Not applicable | `format_transform()` defaults to `-1`; `G16LogFile.render_fakeg()` defaults to `"all"` and accepts int/slice/list/all. |
| Separate output | `embed_in_one_file=False` | Supported | Partial | Unsupported | `render_fakeg()` can return one string per frame; registry disk output still follows the common file-writer contract. |
| File header | Version, options, keywords, title | Supported | Partial | Unsupported | Rebuilds a Gaussian-like banner, Link0, route, and title; a missing route becomes `#p fakeg`. |
| Shared-memory CPU line | `%nprocshared` | Supported | Partial | Unsupported | Reuses the Gaussian Link0 parser; nonnumeric values produce a generic description. |
| `Symbolic Z-matrix` header | First-frame atoms/coords/charge/multiplicity | Supported | Supported | Unsupported | The label follows Gaussian convention, but current output contains Cartesian atom coordinates rather than the original Z-matrix. |
| Input/standard orientation | Atoms, coords, standard_coords | Supported | Supported | Unsupported | Standard orientation takes priority; center/type fields and precision use a fixed template. |
| SCF energy | `energies.reference_energy` | Supported | Supported | Unsupported | Emits `SCF Done: E(SCF)` with one synthetic cycle; functional labels, cycle counts, and convergence history are not retained. |
| Post-HF energy | MP2-MP5, CCSD, CCSD(T) fields | Unsupported | Unsupported | Unsupported | Corresponding Gaussian post-HF energy lines are not generated. |
| Total spin | `total_spin` | Supported | Partial | Unsupported | `S**2` and `S` can render; spin-contamination diagnostics are not reconstructed. |
| MO and populations | `molecular_orbitals`, `charge_spin_populations` | Raw/unsupported | Not guaranteed | Unsupported | No field-driven canonical renderer exists; parsed component raw text may be returned when available. |
| Dipole / polarizability / multipoles | `polarizability` | Raw/unsupported | Not guaranteed | Unsupported | The component tree may carry raw text, but pure model fields cannot synthesize the complete response section. |
| Harmonic frequencies | `vibrations.frequencies` | Supported | Supported | Unsupported | Up to three modes render per batch with synthetic `A` symmetry labels. |
| Reduced mass / force constant / IR | Corresponding `vibrations` arrays | Supported | Supported | Unsupported | Lines render only when fields exist, using fixed format and precision. |
| Normal-mode displacement | `vibration_modes`, atoms | Supported | Supported | Unsupported | Full frequency blocks render per-atom vectors; a single child-node renderer provides only a limited preview. |
| Raman / VCD / ROA | No stable fields | Unsupported | Unsupported | Unsupported | The frequency header uses Gaussian Raman wording, but no Raman activity or depolarization values are emitted. |
| Thermochemistry header | Temperature, pressure | Supported | Supported | Unsupported | Only available values render; missing conditions are not inferred. |
| Mass, inertia, rotational/vibrational temperatures | `thermal_informations` | Supported | Supported | Unsupported | Uses fixed Gaussian-like labels and precision rather than source table layout. |
| ZPVE / thermal corrections | ZPVE, TCE, TCH, TCG | Supported | Supported | Unsupported | Units and values come from normalized containers. |
| U0 / UT / H / G | Thermal summary fields | Supported | Supported | Unsupported | Generates summary lines that MolOP's Gaussian parser can recognize again. |
| Entropy / heat capacity | `S`, `C_V` | Supported | Supported | Unsupported | Emits normalized total values and a thermochemistry header, not component breakdowns. |
| Cartesian forces | `forces` | Raw/unsupported | Not guaranteed | Unsupported | No field-driven force-table renderer exists; original component raw text may survive. |
| Cartesian Hessian | `hessian` | Raw/unsupported | Not guaranteed | Unsupported | No field-driven complete second-derivative renderer exists. |
| Optimization convergence | `geometry_optimization_status` | Partial | Supported | Unsupported | File-level multiframe rendering generates step number and four convergence metrics; single-frame output mainly relies on component content. |
| Multiframe optimization selection | Optimization status + total energy + no vibrations | Supported | Supported | Unsupported | Only matching optimization frames enter the body; other selected frames may be omitted. |
| Frequency-frame selection | Last frame containing vibrations | Supported | Supported | Unsupported | Multiframe rendering selects only the last frequency frame and is not a complete frame-by-frame dump. |
| Runtime | `running_time` | Partial | Partial | Unsupported | Single-frame output can emit `Job cpu time`; multiframe file bodies remove per-frame runtime lines. |
| Termination | Renderer-generated | Supported | Supported | Unsupported | `Normal termination` is appended after multiframe optimization/frequency sections and does not prove the source calculation terminated normally. |
| Whole component-tree rendering | Parsed/synthetic components | Partial | Partial | Unsupported | Nodes with dedicated renderers use payload fields; other nodes may return `raw_text`. |
| Named component node | `component_tree.render_node(name)` | Supported | Not applicable | Unsupported | Nodes such as `l716.forceconstants`, `l716.vibration.mode[0]`, and temperature can render independently. |
| fakeG reparsing | Generated Gaussian-like text | Partial | Supported | Not applicable | Structure, optimization frames, frequency, and thermochemistry are verified; all source fields are not guaranteed to round-trip. |
| Original Gaussian-log reproduction | Source log | Unsupported | Not applicable | Unsupported | Original Link order, iterations, banners, diagnostics, timestamps, and byte layout are not retained. |

fakeG is designed to produce text that is sufficiently Gaussian-like for tests and
downstream compatibility and can be reparsed by MolOP for covered fields. It must not be
treated as evidence of source calculation success, original termination status, or complete
provenance.
