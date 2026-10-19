# ORCA Input

<!-- format-support:orcainp -->

| Item | Value |
| ---- | ----- |
| Format ID | `orcainp` |
| Extensions | `.inp` |
| Read | Yes |
| Write | Yes |
| Registry role | Reader, file writer, frame writer |
| Data level | Coordinates and QM input semantics |

Read ORCA coordinates, multiple jobs, model chemistry and structured input requests.

Coordinate-bearing frames can be rendered with explicit ORCA calculation settings:

Download the [shared ORCA example](../../../assets/examples/water_mp2.out), then run:

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
rendered = frame.format_transform(
    "orcainp",
    keywords="wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP",
    nprocs=16,
    maxcore=4000,
    blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
)
```

`rendered` is a string containing the ORCA simple input line, `%pal/%maxcore/%scf`, and an
`* xyz` coordinate block.

??? example "Output"
    ```text
    ! wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP

    %pal
      nprocs 16
    end

    %maxcore 4000

    %scf
      MaxIter 300
      STABPerform true
    end

    * xyz 0 1
    O      1.7849140000     1.2624220000     0.5119850000
    H      2.6482370000     1.0729290000     0.1316310000
    H      1.1831680000     1.2568160000    -0.2388350000
    *
    ```

Cross-format conversion requires `keywords`; MolOP does not reuse keyword syntax
from another QM program. Raw `%block` text is also accepted through `blocks`.

| Feature | Support | Scope | Limits |
| ------- | ------- | ----- | ------ |
| <!-- feature-area:Geometry and job splitting -->Geometry and multiple jobs | Partial | Direct coordinates, point charges, external xyzfile/pdbfile references, and $new_job splitting. | External geometry references are parsed as references; parser does not require dereferencing external files. |
| <!-- feature-area:Model chemistry and options -->Model chemistry and options | Partial | Method, functional, basis, auxiliary basis, dispersion-as-functional-suffix, mixed basis, print settings, PARAS variables, and scan coordinates. | Unsupported ORCA block options may remain raw/resource fields instead of dedicated semantic containers. |
| <!-- feature-area:ORCA 6 manual fixture families -->ORCA 6 input syntax | Partial | Single point, SCF stability, optimization, frequency, excited-state, MRCI, and Solvator input examples. | Complete coordinates require an inline geometry declaration; external geometries remain references. |
| <!-- feature-area:Excited-state and multireference requests -->Excited-state and multireference requests | Partial | Structured excited-state and multireference task semantics, including MRCI multi-job input. | ORCA excited-state and multireference keyword families outside the covered examples may still be raw. |
| <!-- feature-area:Optimization, coordinates, frequency, and Solvator structures -->Optimization, frequency and Solvator settings | Partial | Optimization constraints, internal coordinates, fragments, NEB, frequency restart, and Solvator examples. | Only the syntax and options listed below are interpreted; unrecognized content retains its raw text. |
| <!-- feature-area:Writer availability -->Input file writing | Partial | Canonical file/frame rendering from coordinate-bearing models with explicit keywords, PAL, maxcore, arbitrary percent blocks, Cartesian geometry, and `.inp` output paths. | Non-Cartesian inline geometry and arbitrary source-statement ordering are not canonicalized; cross-format conversion requires explicit ORCA keywords. |

## ORCA 6.1 Capability Matrix

The matrix distinguishes four levels so that raw-text retention is not mistaken for
structured support:

- **Supported**: the listed syntax has an explicit model and tests for reading or canonical rendering.
- **Partial**: only the listed forms are covered; all ORCA 6.1 variants are not guaranteed.
- **Raw**: original `%block` text can be passed through, but MolOP does not interpret or validate it.
- **Unsupported**: the input is rejected, information is lost, or the model cannot regenerate it.

| ORCA input capability | Parse | Structured semantics | Canonical render | Source fidelity | Current boundary |
| --------------------- | ----- | -------------------- | ---------------- | --------------- | ---------------- |
| `!` simple keyword lines | Supported | Partial | Supported | Unsupported | Keyword text is retained; method, functional, and basis recognition uses built-in vocabularies, and output is reformatted. |
| Top-level `#` comments | Partial | Not applicable | Supported | Unsupported | Standalone comment lines are stored. A second `#` or line ending closes a comment; code after a closed comment is parsed and quoted filenames retain literal `#`. Global comment position is not preserved. |
| `%maxcore` and `%pal` | Supported | Supported | Supported | Partial | Structured `maxcore` and `nprocs` overrides are supported; output uses canonical formatting. |
| Ordinary non-nested `%block ... end` | Supported | Partial | Raw | Partial | Unmodeled blocks can remain in `ORCABlock.raw_text`, but their original position relative to other statement types is not guaranteed. |
| Single-line `%` directives such as `%moinp` and `%base` | Partial | Unsupported | Raw | Partial | Generic block scanning captures them; most have no dedicated model or parameter validation beyond resource fields. |
| Nested blocks such as `%scf/SOSCF` and `%basis/NewGTO` | Partial | Partial | Partial | Unsupported | Nested recognition uses a limited allowlist; an unknown inner `end` can truncate the outer block. |
| Repeated blocks and input priority | Partial | Unsupported | Partial | Unsupported | Parsing keeps block-list order; name-based override removes old instances and appends the replacement at the end. |
| Unclassified top-level statements | Partial | Unsupported | Partial | Unsupported | They are stored in `trailing_lines` and rendered after geometry. |
| `* xyz`, `* cart`, and `* cartesian` | Supported | Supported | Supported | Unsupported | Coordinates are canonicalized as `* xyz`; whitespace, significant digits, and the original coordinate type are not retained. |
| Cartesian geometry in `%coords` | Supported | Supported | Supported | Unsupported | Parsed geometry is rendered as `* xyz`, not regenerated as `%coords`. |
| Point charges `Q q x y z` | Supported | Supported | Supported | Unsupported | Charges use six decimal places; coordinates use the writer precision setting. |
| Ghost, dummy, fragment, frozen, and atom-level basis markers | Partial | Partial | Partial | Unsupported | The listed forms are supported; original tokens, quoting, and ordering are not guaranteed. |
| Isotope `M=` and nuclear charge `Z=` | Supported | Supported | Unsupported | Unsupported | The parser stores both atom attributes, but the Cartesian writer does not emit them. |
| `* int`, `* internal`, and `* gzmt` | Partial | Partial | Unsupported | Unsupported | Covered numeric internal coordinates can be parsed and projected to Cartesian coordinates; rendering raises `NotImplementedError`. |
| `* xyzfile`, `* gzmtfile`, and `* pdbfile` | Supported | Supported | Supported | Partial | The path reference is stored and rendered without reading or validating the external file. |
| `%paras` and coordinate expressions | Partial | Partial | Partial | Unsupported | Simple parameters, ranges, and `name +/- offset` are recognized; geometry rendering replaces expressions with resolved values. |
| `$new_job` multi-job input | Supported | Supported | Supported | Unsupported | Each job maps to one frame; output uses a canonical `$new_job` separator. |
| `%Compound` workflows | Partial | Unsupported | Unsupported | Unsupported | Valid Compound-only input without top-level geometry is rejected by format probing; rendering also requires keywords and geometry. |
| Model chemistry and common tasks | Partial | Partial | Partial | Not applicable | Common SP, OPT, TS/NEB, FREQ, IRC, gradient, DFT, selected wavefunction methods, and dispersion keywords are covered. |
| Excited-state, multireference, and Solvator input | Partial | Partial | Raw | Partial | Parsing provides structured requests; rendering relies on retained raw blocks or explicit blocks and does not fully synthesize blocks from request models. |
| Full `%basis`, ECP, AuxJ/AuxJK/AuxC/CABS | Partial | Unsupported | Raw | Partial | Simple keywords and atom-level overrides are recognized; complete block contents are not modeled or validated field by field. |
| ORCA 6.1 QMMM, MD, GOAT, DOCKER, RESP, Raman, and similar features | Partial | Unsupported | Raw | Partial | Non-nested raw blocks can be supplied; nested syntax, task compatibility, and field-level construction are not guaranteed. |
| ORCA version and compatibility validation | Unsupported | Unsupported | Unsupported | Not applicable | The version is currently recorded as `Any`; keyword, method, basis, and ORCA-version compatibility are not checked. |

The current renderer is a **canonical writer** for constructing new input from coordinates
and explicit calculation settings. It is not a **source-preserving writer**. Do not rely on
round-trip behavior when editing one field in an existing complex input must leave all other
source text unchanged.

References: [ORCA 6.1 input structure](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html),
[coordinates](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html),
[basis sets](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html), and
[Compound](https://www.faccts.de/docs/orca/6.1/manual/contents/workflowsautomatization/compound.html).

## Property sources and reading conventions

<!-- parsing-basis:start -->

The table maps file contents to MolOP data. See [data containers and cross-QM
semantics](../data_model.md) for shared units, indices and missing values.

The references below use ORCA 6.0/6.1 manuals and tutorials to explain syntax,
quantities and printed units. MolOP support is defined on this page; older
output headers and column orders can differ.

| Property and references | File marker or data source | Reading or writing convention |
| -------- | -------------------------- | ----------------------------- |
| <!-- parsing-basis:Geometry and job splitting --><!-- source:ORCA-COORD --><!-- source:ORCA-INPUT --><!-- source:ORCA-UNITS -->Geometry and multiple jobs ([Coordinates][coordinates], [Input structure][input], [Units][units]) | `* xyz/internal/xyzfile`; `%coords`; `$new_job`; `charge/multiplicity`; `unit declarations` | Inline geometry maps to atoms/Å coordinates only where resolvable; internal coordinates/parameters retain typed input semantics. External paths remain references. Job delimiters define source spans. |
| <!-- parsing-basis:Model chemistry and options --><!-- source:ORCA-METHOD --><!-- source:ORCA-RESOURCE --><!-- source:ORCA-OUTPUT --><!-- source:ORCA-SOLVENT -->Model chemistry and options ([Wavefunction types][hftype], [Parallelism and memory][parallel], [Output controls][output], [Solvation][solvation]) | `! keyword lines`; `%scf/%basis/%pal`; `%maxcore`; `%output`; `PARAS variables` | Map method/spin/basis/auxiliary basis/dispersion/solvent/resource/print requests from bounded semantics. maxcore is per-process memory; unknown block options remain raw rather than invented result fields. |
| <!-- parsing-basis:ORCA 6 manual fixture families -->ORCA 6 input syntax ([Input structure][input]) | `!` keywords, `%...end` blocks and coordinate declarations | Read supported keywords, blocks and geometry declarations. Input files describe calculation requests; read energies, spectra and convergence from output files. |
| <!-- parsing-basis:Excited-state and multireference requests --><!-- source:ORCA-EXCITED -->Excited-state and multireference requests ([Excited states][excited], [Multireference input][mrci]) | `%tddft/%cis/%rocis/%casscf/%mrci and covered method keywords` | `excited_state_requests` / `multireference_requests` retain root counts, states, reference and active-space settings where mapped; roots are source request labels, not observed eigenstates. |
| <!-- parsing-basis:Optimization, coordinates, frequency, and Solvator structures --><!-- source:ORCA-OPT --><!-- source:ORCA-VIB -->Optimization, frequency and Solvator settings ([Geometry optimization][optimization], [Numerical gradients][gradient], [IR and Raman][vibrations]) | `%geom Constraints/Scan/Fragments`; `NEB`; `%freq restart`; `%solvator` | Retain supported task/constraint/scan/environment inputs, including nested block boundaries. Source coordinates and requested displacements are distinct; unrecognized options retain their raw text. |
| <!-- parsing-basis:Writer availability -->Input file writing ([Input structure][input], [Coordinates][coordinates], [Parallelism and memory][parallel]) | `Common geometry/model fields`; `explicit ! keywords`; `PAL/maxcore and %blocks` | Generate canonical Cartesian ORCA input. Cross-format keywords are explicit; a model’s Gaussian keywords do not silently become ORCA semantics. Arbitrary original ordering is not reconstructed. |

<!-- parsing-basis:end -->

[input]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html
[coordinates]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html
[parallel]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html
[output]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/output.html
[hftype]: https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hftype.html
[mrci]: https://www.faccts.de/docs/orca/6.0/manual/contents/typical/mrci.html
[vibrations]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html
[optimization]: https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html
[gradient]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalgradients.html
[solvation]: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html
[excited]: https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/excitedstates.html
[units]: https://www.faccts.de/docs/orca/6.1/manual/contents/orcaarchitecture/conversionfactors.html
