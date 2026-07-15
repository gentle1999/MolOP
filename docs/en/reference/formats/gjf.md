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

```python
rendered = frame.format_transform(
    "gjf",
    link0_commands={"nprocshared": "16", "mem": "32GB"},
    route_section="#p wb97xd/def2tzvp opt freq",
    title_card="optimization and frequency",
    coords_type="cartesian",
    chk=True,
)
```

Provide an explicit Gaussian `route_section` for cross-format conversion. The writer
does not validate compatibility among methods, basis sets, tasks, and Gaussian versions.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Fixture parse/render inventory -->Fixture parse/render inventory | Fixture-covered | Maintained Gaussian input fixtures are parseable and renderable. | Coverage is fixture-anchored rather than a complete Gaussian grammar proof. | `tests/test_g16gjf_fixtures_coverage.py::test_g16gjf_all_fixtures_are_parseable_and_renderable` |
| <!-- feature-area:Route semantics -->Route semantics | Partial | Model chemistry, job types, dispersion, solvation, population requests, and HF handling are parsed into shared semantic containers. | Functional is populated only for DFT/hybrid/double-hybrid semantics; method-specific Gaussian keywords outside tests may remain raw or unstructured. | `tests/test_gaussian_route_semantics.py::test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities`<br>`tests/test_gaussian_route_semantics.py::test_hf_route_does_not_populate_functional`<br>`tests/test_gaussian_route_semantics.py::test_shared_semantic_route_exposes_structured_dispersion_and_solvation`<br>`tests/test_gaussian_route_semantics.py::test_shared_route_parser_supports_more_gaussian_job_types_and_params` |
| <!-- feature-area:Molecule specification -->Molecule specification | Partial | Cartesian, integer Cartesian, Z-matrix variables, free delimiters, fragments, and invalid-domain diagnostics. | Supported cases are the explicit molecule-spec grammar variants covered by tests. | `tests/test_gjf_spec_cases.py::test_zmat_labeled_variable_blocks_are_applied_to_atom_lines`<br>`tests/test_gjf_spec_cases.py::test_zmat_unlabeled_variable_blocks_are_applied_by_position`<br>`tests/test_gjf_spec_cases.py::test_free_format_delimiters_are_parsed_and_rendered_canonically`<br>`tests/test_gjf_spec_cases.py::test_integer_cartesian_coordinates_are_not_misparsed_as_zmatrix`<br>`tests/test_gjf_spec_cases.py::test_multifragment_assignments_match_declared_pairs`<br>`tests/test_gjf_spec_cases.py::test_invalid_zmatrix_domains_raise_clear_errors` |
| <!-- feature-area:Additional sections -->Additional sections | Partial | ModRedundant, GIC, NBO, structured rendering, and mixed-section diagnostics. | Mixed GIC/ModRedundant sections intentionally fall back to unknown with diagnostics. | `tests/test_gjf_spec_cases.py::test_modredundant_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_gic_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_nbo_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic`<br>`tests/test_gjf_spec_cases.py::test_render_prefers_structured_additional_sections` |
| <!-- feature-area:Link1, includes, and writer options -->Link1, includes, and writer options | Partial | Link1 uses exact source spans with blank-line validation; Geom=AllCheck and checkpoint propagation during transform are supported. | `@include` crosses source artifacts and cannot be represented by the current single-artifact `SourceSpan`, so parsing rejects it with `MOL.PARSE.GJF_INCLUDE_PROVENANCE_UNSUPPORTED`; not every Gaussian writer option is covered. | `tests/test_gjf_spec_strictness.py::test_link1_requires_blank_line_before_separator`<br>`tests/test_gjf_spec_strictness.py::test_geom_allcheck_allows_missing_title_and_molecule_sections`<br>`tests/test_gjf_include_handling.py::test_gjf_disk_parser_rejects_include_without_multi_artifact_spans`<br>`tests/test_format_transform_output_dir.py::test_format_transform_gjf_chk_propagation_single_file` |

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
