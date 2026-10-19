# Gaussian Input

<!-- format-support:gjf -->

| Item | Value |
| ---- | ----- |
| Format ID | `gjf` |
| Extensions | `.gjf`, `.gif`, `.com`, `.gau`, `.gjc` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates and QM input semantics |

Gaussian input parsing and canonical rendering cover Link0, route, title,
molecule specification, Link1, and additional sections. The writer can render an
existing GJF frame or construct Gaussian input from another coordinate-bearing format:

Download the [shared ORCA example](../../../assets/examples/water_mp2.out), then run:

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
rendered = frame.format_transform(
    "gjf",
    link0_commands={"nprocshared": "16", "mem": "32GB"},
    route_section="#p wb97xd/def2tzvp opt freq",
    title_card="optimization and frequency",
    coords_type="cartesian",
    chk=True,
)
```

`rendered` is a string beginning like:

??? example "Rendered prefix"

    ```text
    %nprocshared=16
    %mem=32GB
    %chk=optimization and frequency.chk
    #p wb97xd/def2tzvp opt freq
    ```

Provide an explicit Gaussian `route_section` for cross-format conversion. The writer
does not validate compatibility among methods, basis sets, tasks, and Gaussian versions.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Fixture parse/render inventory -->Input structure and writing | Partial | Reads and normalizes Link0, route, title, molecule specification and additional sections. | Unrecognized options can retain raw text; full Gaussian grammar and version compatibility are not guaranteed. |
| <!-- feature-area:Route semantics -->Route settings | Partial | Model chemistry, job types, dispersion, solvation, population requests, and HF handling are parsed into shared semantic containers. | `functional` is populated only for DFT/hybrid/double-hybrid semantics; unrecognized Gaussian method keywords may remain raw or unstructured. |
| <!-- feature-area:Molecule specification -->Molecule specification | Partial | Cartesian, integer Cartesian, Z-matrix variables, free delimiters, fragments, and invalid-domain diagnostics. | Only the molecule-specification forms listed below are supported; not every ONIOM/MM extension is covered. |
| <!-- feature-area:Additional sections -->Additional input sections | Partial | ModRedundant, GIC, NBO, structured rendering, and mixed-section diagnostics. | Mixed GIC/ModRedundant sections intentionally fall back to unknown with diagnostics. |
| <!-- feature-area:Link1, includes, and writer options -->Link1, includes and writer options | Partial | Link1 uses exact source spans with blank-line validation; Geom=AllCheck and checkpoint propagation during transform are supported. | `@include` crosses source artifacts and cannot be represented by the current single-artifact `SourceSpan`, so parsing rejects it with `MOL.PARSE.GJF_INCLUDE_PROVENANCE_UNSUPPORTED`; not every Gaussian writer option is covered. |

## Gaussian Input Capability Matrix

In the matrix, source fidelity means whether parsing and rendering retain original order,
whitespace, line wrapping, and token representation. It is separate from preserving the
same calculation meaning in structured fields.

| Gaussian input capability | Parse | Structured semantics | Canonical render | Source fidelity | Current boundary |
| ------------------------- | ----- | -------------------- | ---------------- | --------------- | ---------------- |
| Link0 `%key=value` | Supported | Partial | Supported | Unsupported | Stored as ordered `GJFLink0Command` items; CPU and memory project to resource requests, while other keys receive no dedicated validation. |
| `%chk` / `%oldchk` render overrides | Not applicable | Supported | Supported | Not applicable | `chk` and `old_chk` append to a temporary Link0 copy; existing directives are neither replaced nor deduplicated. |
| Single-line and multiline route sections | Supported | Partial | Supported | Unsupported | Routes are normalized to text beginning with `#`; original wrapping, whitespace, and token case are not guaranteed. |
| Method / functional / basis | Partial | Partial | Supported | Not applicable | HF, DFT/hybrid/double-hybrid, selected post-HF methods, and common basis sets are recognized; unknown text remains in the route without guaranteed semantic fields. |
| OPT / FREQ / SP / IRC / scan tasks | Partial | Partial | Supported | Not applicable | Common job types and covered parameters enter `task_requests`; all Gaussian task combinations are not enumerated. |
| `Opt(...)` options | Partial | Partial | Supported | Not applicable | Covered fields include TS, CalcFC, ReadFC, ModRedundant, MaxCycles, and coordinate-system options. |
| `Freq(...)` options | Partial | Partial | Supported | Not applicable | Covered fields include anharmonic, projected, hindered rotor, VCD/ROA/Raman, temperature/pressure, and atom selectors; task compatibility is not checked. |
| TD, SCRF, Pop, and Geom options | Partial | Partial | Supported | Not applicable | Recognized options enter dedicated semantic containers; unknown parameters remain in the raw route. |
| Title card | Supported | Supported | Supported | Unsupported | Up to five lines; model validation removes specified invalid characters, and an empty title renders as the filename or `title`. |
| Charge and multiplicity | Supported | Supported | Supported | Unsupported | Single- and multi-fragment charge/multiplicity have structured models; output uses canonical whitespace. |
| Cartesian coordinates | Supported | Supported | Supported | Unsupported | Integer and floating-point coordinates parse; the writer uses six decimal places and does not retain original delimiters or precision. |
| Atom label, atom type, atomic charge, and parameters | Partial | Partial | Partial | Unsupported | Gaussian atom labels, type, charge, `(...)` parameters, and frozen tags are stored; all ONIOM/MM extensions are not covered. |
| Dummy and ghost atoms | Partial | Partial | Partial | Unsupported | Covered `X` and `Bq` forms project to the model; implicit Cartesian/internal conversion is rejected when dummy or ghost atoms are present. |
| Numeric Z-matrix | Supported | Supported | Supported | Unsupported | Ordinary dihedrals and trailing `0/1` alternate-angle markers are handled, with bond, angle, and reference-domain validation. |
| Z-matrix variables | Partial | Partial | Partial | Unsupported | Covered labeled and unlabeled variable sections are resolved to numbers during parsing; the writer does not restore variable names. |
| Free-format coordinate delimiters | Partial | Partial | Supported | Unsupported | Covered comma, tab, and `/` forms parse; output is space-delimited. |
| Cartesian/internal coordinate conversion | Not applicable | Supported | Partial | Not applicable | `coords_type` converts ordinary real-atom fragments both ways; dummy/ghost or unsafe alignment raises `ValueError`. |
| Fragment molecule specification | Partial | Supported | Supported | Unsupported | Every atom must declare contiguous `Fragment=n`; declared fragment charge/multiplicity counts are validated. |
| `Geom=AllCheck` / `Geom=Checkpoint` | Supported | Partial | Partial | Unsupported | Parsing allows missing title or molecule blocks; canonical rendering reassembles from model fields and default-title rules. |
| Connectivity section | Partial | Unsupported | Partial | Unsupported | Original connectivity has no dedicated field; `add_gjf_connectivity=True` regenerates connectivity from the recovered molecular graph. |
| ModRedundant | Partial | Partial | Supported | Unsupported | Covered line forms are structured; complex or unknown commands may fall back to an unknown/raw additional section. |
| GIC | Partial | Partial | Supported | Unsupported | Labels, options, functions, arguments, and standalone actions are parsed; complete Gaussian GIC grammar is not claimed. |
| `$NBO ... $END` | Partial | Partial | Supported | Unsupported | Header, commands, and footer are stored; individual NBO commands are not semantically validated. |
| Gen/GenECP, custom basis, and other extra input | Partial | Unsupported | Raw | Partial | The route recognizes `Gen` / `GenECP`; detailed basis/ECP text is usually an unknown additional section. |
| Mixed or unknown additional sections | Supported | Unsupported | Raw | Partial | Mixed GIC/ModRedundant content degrades to unknown with a diagnostic; raw content can render without interpretation. |
| `--Link1--` multi-job input | Supported | Supported | Supported | Partial | Exact source spans split frames and require a blank line before the separator; the file writer uses a canonical separator. |
| Link1 checkpoint propagation | Not applicable | Partial | Supported | Not applicable | Batch/multiframe writing can generate `%chk` / `%oldchk` chains; filenames and output directories follow common transform parameters. |
| `@include` | Unsupported | Unsupported | Unsupported | Unsupported | Disk and memory parsers reject includes because multi-artifact provenance is not modeled; includes are not read recursively. |
| Build GJF from a general coordinate frame | Not applicable | Partial | Supported | Not applicable | `atoms`, `coords`, charge, and multiplicity build one Cartesian fragment; callers should supply the calculation route. |
| Gaussian version and keyword compatibility | Unsupported | Unsupported | Unsupported | Not applicable | `qm_software_version` is `Any`; acceptance by a specific Gaussian version is not checked. |

The GJF writer is a **canonical writer** for constructing and normalizing Gaussian input,
not a **source-preserving editor**. Do not use `parse -> render` as a lossless editing
workflow when all untouched source text must remain byte-for-byte identical.

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

References: [Gaussian input reference](https://gaussian.com/input/).

| Property | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Fixture parse/render inventory -->Input structure and writing | `%Link0`; `# route`; `title`; `charge/multiplicity`; `molecule and trailing sections` | Read Link0, route, title, charge, multiplicity, molecule specification and additional sections. Writing normalizes whitespace and coordinates without checking Gaussian-version compatibility. |
| <!-- parsing-basis:Route semantics -->Route settings | `HF/DFT/post-HF method/basis tokens`; `Opt/Freq/TD/SCRF/Pop` | `model_chemistry`, task/resource/excited-state/population/solvation requests derive from bounded token mappings. A requested property is not an observed result; unknown options remain raw. |
| <!-- parsing-basis:Molecule specification -->Molecule specification | `Cartesian or Z-matrix records`; `variables`; `fragment charge/spin pairs` | Retain molecule specification, fragments, labels, dummy/ghost and frozen markers; derive Cartesian coordinates only for supported internal-coordinate cases. Lengths use Å; angles use degrees. |
| <!-- parsing-basis:Additional sections -->Additional input sections | `ModRedundant/GIC records`; `$NBO ... $END`; `trailing blocks` | Classify supported additional sections into typed input objects. Mixed or unrecognized grammar is retained as unknown with diagnostics, not invented output properties. |
| <!-- parsing-basis:Link1, includes, and writer options -->Link1, includes and writer options | `--Link1--`; `blank-line boundaries`; `%chk/%oldchk`; `@include`; `Geom=AllCheck` | Preserve source spans for steps; checkpoint chains and output paths follow writer policy. external files are not expanded; @include reports an unsupported-input error. |

<!-- parsing-basis:end -->
