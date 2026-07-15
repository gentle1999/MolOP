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

ORCA input parsing with coordinates, multi-job splitting, model chemistry, task-family fixtures, and structured request containers.

Coordinate-bearing frames can be rendered with explicit ORCA calculation settings:

```python
rendered = frame.format_transform(
    "orcainp",
    keywords="wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP",
    nprocs=16,
    maxcore=4000,
    blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
)
```

Cross-format conversion requires `keywords`; MolOP does not reuse keyword syntax
from another QM program. Raw `%block` text is also accepted through `blocks`.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Geometry and job splitting -->Geometry and job splitting | Partial | Direct coordinates, point charges, external xyzfile/pdbfile references, and $new_job splitting. | External geometry references are parsed as references; parser does not require dereferencing external files. | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_minimal_xyz_parses_atoms_and_coords`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_point_charges_parse_into_frame_geometry`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_xyzfile_does_not_require_external_file`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_new_job_splits_frames` |
| <!-- feature-area:Model chemistry and options -->Model chemistry and options | Partial | Method, functional, basis, auxiliary basis, dispersion-as-functional-suffix, mixed basis, print settings, PARAS variables, and scan coordinates. | Unsupported ORCA block options may remain raw/resource fields instead of dedicated semantic containers. | `tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_metadata_population_recognizes_d4`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_mixed_basis_and_output_print_settings`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_paras_structures_scan_and_resolves_cartesian_variables`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_dispersion_is_functional_suffix` |
| <!-- feature-area:ORCA 6 manual fixture families -->ORCA 6 manual fixture families | Fixture-covered | Single point, SCF stability, optimization, frequency, excited-state, MRCI, and Solvator input examples. | Manual snippets without complete explicit geometry are intentionally excluded from fixture coverage. | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_single_point_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_scf_stability_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_optimization_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_frequency_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_inputs_include_orca6_manual_fixtures`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_solvator_inputs_include_orca6_manual_fixtures` |
| <!-- feature-area:Excited-state and multireference requests -->Excited-state and multireference requests | Partial | Structured excited-state and multireference task semantics, including MRCI multi-job fixtures. | ORCA excited-state and multireference keyword families outside the covered examples may still be raw. | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_manual_fixtures_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_multijob_manual_fixture_splits_frames_and_structures_last_job` |
| <!-- feature-area:Optimization, coordinates, frequency, and Solvator structures -->Optimization, coordinates, frequency, and Solvator structures | Partial | Optimization constraints, internal coordinates, fragments, NEB, frequency restart, and Solvator examples. | Coverage follows explicit ORCA manual fixtures and targeted regression examples. | `tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_constraints_block_is_not_truncated`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_internal_coords_fill_frame_atoms_and_coords`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_fragment_mixed_basis_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_neb_block_and_xtb_geometry_are_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_frequency_restart_example_is_structured`<br>`tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_solvator_fixture_is_structured` |
| <!-- feature-area:Writer availability -->Writer availability | Partial | Canonical file/frame rendering from coordinate-bearing models with explicit keywords, PAL, maxcore, arbitrary percent blocks, Cartesian geometry, and `.inp` output paths. | Non-Cartesian inline geometry and arbitrary source-statement ordering are not canonicalized; cross-format conversion requires explicit ORCA keywords. | `tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_builds_common_sp_input`<br>`tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_uses_inp_output_extension` |

## ORCA 6.1 Capability Matrix

The matrix distinguishes four levels so that raw-text retention is not mistaken for
structured support:

- **Supported**: the listed syntax has an explicit model and tests for reading or canonical rendering.
- **Partial**: only known forms or fixtures are covered; all ORCA 6.1 variants are not guaranteed.
- **Raw**: original `%block` text can be passed through, but MolOP does not interpret or validate it.
- **Unsupported**: the input is rejected, information is lost, or the model cannot regenerate it.

| ORCA input capability | Parse | Structured semantics | Canonical render | Source fidelity | Current boundary |
| --------------------- | ----- | -------------------- | ---------------- | --------------- | ---------------- |
| `!` simple keyword lines | Supported | Partial | Supported | Unsupported | Keyword text is retained; method, functional, and basis recognition uses built-in vocabularies, and output is reformatted. |
| Top-level `#` comments | Partial | Not applicable | Supported | Unsupported | Standalone comment lines are stored; global position is lost, and parsing does not resume after a closing `#` on the same line. |
| `%maxcore` and `%pal` | Supported | Supported | Supported | Partial | Structured `maxcore` and `nprocs` overrides are supported; output uses canonical formatting. |
| Ordinary non-nested `%block ... end` | Supported | Partial | Raw | Partial | Unmodeled blocks can remain in `ORCABlock.raw_text`, but their original position relative to other statement types is not guaranteed. |
| Single-line `%` directives such as `%moinp` and `%base` | Partial | Unsupported | Raw | Partial | Generic block scanning captures them; most have no dedicated model or parameter validation beyond resource fields. |
| Nested blocks such as `%scf/SOSCF` and `%basis/NewGTO` | Partial | Partial | Partial | Unsupported | Nested recognition uses a limited allowlist; an unknown inner `end` can truncate the outer block. |
| Repeated blocks and input priority | Partial | Unsupported | Partial | Unsupported | Parsing keeps block-list order; name-based override removes old instances and appends the replacement at the end. |
| Unclassified top-level statements | Partial | Unsupported | Partial | Unsupported | They are stored in `trailing_lines` and rendered after geometry. |
| `* xyz`, `* cart`, and `* cartesian` | Supported | Supported | Supported | Unsupported | Coordinates are canonicalized as `* xyz`; whitespace, significant digits, and the original coordinate type are not retained. |
| Cartesian geometry in `%coords` | Supported | Supported | Supported | Unsupported | Parsed geometry is rendered as `* xyz`, not regenerated as `%coords`. |
| Point charges `Q q x y z` | Supported | Supported | Supported | Unsupported | Charges use six decimal places; coordinates use the writer precision setting. |
| Ghost, dummy, fragment, frozen, and atom-level basis markers | Partial | Partial | Partial | Unsupported | Covered fixture forms are handled; original tokens, quoting, and ordering are not guaranteed. |
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
