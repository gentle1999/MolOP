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

Gaussian input parsing/rendering with route, molecule specification, Link1, and additional-section coverage.

| Feature | Support | Scope | Limits | Test evidence |
| ------- | ------- | ----- | ------ | ------------- |
| <!-- feature-area:Fixture parse/render inventory -->Fixture parse/render inventory | Fixture-covered | Maintained Gaussian input fixtures are parseable and renderable. | Coverage is fixture-anchored rather than a complete Gaussian grammar proof. | `tests/test_g16gjf_fixtures_coverage.py::test_g16gjf_all_fixtures_are_parseable_and_renderable` |
| <!-- feature-area:Route semantics -->Route semantics | Partial | Model chemistry, job types, dispersion, solvation, population requests, and HF handling are parsed into shared semantic containers. | Functional is populated only for DFT/hybrid/double-hybrid semantics; method-specific Gaussian keywords outside tests may remain raw or unstructured. | `tests/test_gaussian_route_semantics.py::test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities`<br>`tests/test_gaussian_route_semantics.py::test_hf_route_does_not_populate_functional`<br>`tests/test_gaussian_route_semantics.py::test_shared_semantic_route_exposes_structured_dispersion_and_solvation`<br>`tests/test_gaussian_route_semantics.py::test_shared_route_parser_supports_more_gaussian_job_types_and_params` |
| <!-- feature-area:Molecule specification -->Molecule specification | Partial | Cartesian, integer Cartesian, Z-matrix variables, free delimiters, fragments, and invalid-domain diagnostics. | Supported cases are the explicit molecule-spec grammar variants covered by tests. | `tests/test_gjf_spec_cases.py::test_zmat_labeled_variable_blocks_are_applied_to_atom_lines`<br>`tests/test_gjf_spec_cases.py::test_zmat_unlabeled_variable_blocks_are_applied_by_position`<br>`tests/test_gjf_spec_cases.py::test_free_format_delimiters_are_parsed_and_rendered_canonically`<br>`tests/test_gjf_spec_cases.py::test_integer_cartesian_coordinates_are_not_misparsed_as_zmatrix`<br>`tests/test_gjf_spec_cases.py::test_multifragment_assignments_match_declared_pairs`<br>`tests/test_gjf_spec_cases.py::test_invalid_zmatrix_domains_raise_clear_errors` |
| <!-- feature-area:Additional sections -->Additional sections | Partial | ModRedundant, GIC, NBO, structured rendering, and mixed-section diagnostics. | Mixed GIC/ModRedundant sections intentionally fall back to unknown with diagnostics. | `tests/test_gjf_spec_cases.py::test_modredundant_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_gic_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_nbo_additional_section_detected_as_typed_section`<br>`tests/test_gjf_spec_cases.py::test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic`<br>`tests/test_gjf_spec_cases.py::test_render_prefers_structured_additional_sections` |
| <!-- feature-area:Link1, includes, and writer options -->Link1, includes, and writer options | Partial | Link1 splitting with blank-line validation; Geom=AllCheck; disk @include expansion; checkpoint propagation during transform. | @include expansion requires source-path context; not every Gaussian writer option is covered. | `tests/test_gjf_spec_strictness.py::test_link1_requires_blank_line_before_separator`<br>`tests/test_gjf_spec_strictness.py::test_geom_allcheck_allows_missing_title_and_molecule_sections`<br>`tests/test_gjf_include_handling.py::test_gjf_disk_parser_expands_at_include`<br>`tests/test_format_transform_output_dir.py::test_format_transform_gjf_chk_propagation_single_file` |
