from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


FormatGroup = Literal["structure", "qm-input", "qm-output", "special"]
SupportLevel = Literal["supported", "partial", "fixture-covered", "intentionally-unsupported"]


@dataclass(frozen=True)
class FeatureSupport:
    area: str
    support: SupportLevel
    scope: str
    limitations: str
    tests: tuple[str, ...]


@dataclass(frozen=True)
class FormatSupport:
    format_id: str
    group: FormatGroup
    title: str
    extensions: tuple[str, ...]
    read: bool
    write: bool
    registry_role: str
    data_level: str
    summary: str
    features: tuple[FeatureSupport, ...]


BUILTIN_READER_FORMATS = frozenset(
    {
        "g16log",
        "gjf",
        "orcainp",
        "orcaout",
        "sdf",
        "smi",
        "xyz",
    }
)
BUILTIN_WRITER_FORMATS = frozenset({"cml", "fakeg", "gjf", "sdf", "smi", "xyz"})
SPECIAL_CODEC_IDS = frozenset({"openbabel-fallback"})

FORMAT_GROUPS: dict[FormatGroup, tuple[str, ...]] = {
    "structure": ("xyz", "sdf", "smi", "cml"),
    "qm-input": ("gjf", "orcainp"),
    "qm-output": ("g16log", "orcaout", "fakeg"),
    "special": ("openbabel-fallback",),
}

FORMAT_SUPPORT: dict[str, FormatSupport] = {
    "xyz": FormatSupport(
        format_id="xyz",
        group="structure",
        title="XYZ coordinates",
        extensions=(".xyz",),
        read=True,
        write=True,
        registry_role="reader, file writer, frame writer",
        data_level="coordinates",
        summary="Standard XYZ coordinate IO with charge/multiplicity comment support.",
        features=(
            FeatureSupport(
                area="Reader",
                support="supported",
                scope="Standard multi-frame XYZ files; charge and multiplicity in comments; file-level charge/multiplicity finalized from the first frame.",
                limitations="No graph recovery; malformed atom counts or incomplete frames fail during parse-phase mismatch handling.",
                tests=(
                    "tests/test_parsers_smoke.py::test_autoparser_xyz_smoke",
                    "tests/test_io_registry_and_batch_more.py::test_base_file_parser_finalizes_file_charge_and_multiplicity_from_first_frame",
                ),
            ),
            FeatureSupport(
                area="Writer and conversion",
                support="supported",
                scope="File/frame XYZ rendering, explicit comment override, and registry conversion from coordinate-bearing parsed files.",
                limitations="Writer uses coordinate semantics; graph-only metadata is not encoded by XYZ.",
                tests=(
                    "tests/test_structure_format_features.py::test_xyz_frame_writer_uses_explicit_comment_override",
                    "tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke",
                    "tests/test_format_transform_output_dir.py::test_format_transform_batch_output_dir_writes_under_requested_dir",
                ),
            ),
            FeatureSupport(
                area="Format mismatch handling",
                support="supported",
                scope="Simple-format mismatch detection is deferred to normal parsing to avoid a second IO probe.",
                limitations="No separate file-level fingerprint is required for XYZ.",
                tests=(
                    "tests/test_io_registry_and_batch_more.py::test_simple_coordinate_readers_defer_mismatch_to_parse_phase",
                ),
            ),
        ),
    ),
    "sdf": FormatSupport(
        format_id="sdf",
        group="structure",
        title="SDF/MOL structures",
        extensions=(".sdf", ".sd", ".mol"),
        read=True,
        write=True,
        registry_role="reader, file writer, frame writer",
        data_level="coordinates on read; graph required on write",
        summary="SDF/MOL structure IO backed by RDKit graph extraction and strict graph writer semantics.",
        features=(
            FeatureSupport(
                area="Reader",
                support="supported",
                scope="SDF/MOL blocks parsed into atom, coordinate, bond, formal charge, radical, total charge, and multiplicity fields.",
                limitations="Arbitrary SD data-field round-trip is not advertised by the tested contract.",
                tests=("tests/test_autoparser_sdf_tmpfile.py::test_autoparser_sdf_tmpfile",),
            ),
            FeatureSupport(
                area="Writer graph policy",
                support="supported",
                scope="SDF file/frame writer requires graph-capable input, defaults to strict graph semantics, and supports RDKit/OpenBabel rendering engines.",
                limitations="Coords-only override is rejected when no coords-only SDF writer exists.",
                tests=(
                    "tests/test_structure_format_features.py::test_sdf_frame_writer_supports_rdkit_openbabel_engines_and_rejects_unknown_engine",
                    "tests/test_io_registry_and_batch_more.py::test_registry_builtin_default_graph_policies_match_format_semantics",
                    "tests/test_io_registry_and_batch_more.py::test_registry_raises_when_coords_override_has_no_coords_writer",
                ),
            ),
            FeatureSupport(
                area="Registry conversion",
                support="supported",
                scope="Parsed structures can be converted to SDF through the codec registry.",
                limitations="Conversion quality depends on recovered or transformed graph data.",
                tests=("tests/test_io_convert_registry_smoke.py::test_registry_convert_xyz_smoke",),
            ),
        ),
    ),
    "smi": FormatSupport(
        format_id="smi",
        group="structure",
        title="SMILES records",
        extensions=(".smi", ".txt"),
        read=True,
        write=True,
        registry_role="reader, file writer, frame writer",
        data_level="coordinates on read; graph required on write",
        summary="SMILES record IO with graph-derived charge/multiplicity and canonical SMILES writing.",
        features=(
            FeatureSupport(
                area="Reader",
                support="supported",
                scope="Non-empty SMILES records; first whitespace token parsed; 2D coordinates, graph, formal charge/radical, charge, and multiplicity derived by RDKit.",
                limitations="Trailing record names or columns are not part of the advertised parsed-data contract; input 3D coordinates cannot be preserved.",
                tests=(
                    "tests/test_autoparser_smi_tmpfile.py::test_autoparser_smi_tmpfile",
                    "tests/test_structure_format_features.py::test_smi_reader_uses_first_token_and_writer_canonicalizes",
                ),
            ),
            FeatureSupport(
                area="Writer graph policy",
                support="supported",
                scope="Canonical graph-derived SMILES writer with strict graph semantics.",
                limitations="Coords-only rendering is not the default writer contract.",
                tests=(
                    "tests/test_structure_format_features.py::test_smi_reader_uses_first_token_and_writer_canonicalizes",
                    "tests/test_io_registry_and_batch_more.py::test_registry_builtin_default_graph_policies_match_format_semantics",
                    "tests/test_io_convert_registry_smoke.py::test_registry_convert_xyz_smoke",
                ),
            ),
            FeatureSupport(
                area="Registry conversion",
                support="supported",
                scope="Gaussian output structures can be converted to SMILES through the registry.",
                limitations="Conversion depends on successful graph recovery or transformation.",
                tests=("tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke",),
            ),
        ),
    ),
    "gjf": FormatSupport(
        format_id="gjf",
        group="qm-input",
        title="Gaussian input",
        extensions=(".gjf", ".gif", ".com", ".gau", ".gjc"),
        read=True,
        write=True,
        registry_role="reader, file writer, frame writer",
        data_level="coordinates and QM input semantics",
        summary="Gaussian input parsing/rendering with route, molecule specification, Link1, and additional-section coverage.",
        features=(
            FeatureSupport(
                area="Fixture parse/render inventory",
                support="fixture-covered",
                scope="Maintained Gaussian input fixtures are parseable and renderable.",
                limitations="Coverage is fixture-anchored rather than a complete Gaussian grammar proof.",
                tests=(
                    "tests/test_g16gjf_fixtures_coverage.py::test_g16gjf_all_fixtures_are_parseable_and_renderable",
                ),
            ),
            FeatureSupport(
                area="Route semantics",
                support="partial",
                scope="Model chemistry, job types, dispersion, solvation, population requests, and HF handling are parsed into shared semantic containers.",
                limitations="Functional is populated only for DFT/hybrid/double-hybrid semantics; method-specific Gaussian keywords outside tests may remain raw or unstructured.",
                tests=(
                    "tests/test_gaussian_route_semantics.py::test_shared_gaussian_route_parser_extracts_model_chemistry_and_capabilities",
                    "tests/test_gaussian_route_semantics.py::test_hf_route_does_not_populate_functional",
                    "tests/test_gaussian_route_semantics.py::test_shared_semantic_route_exposes_structured_dispersion_and_solvation",
                    "tests/test_gaussian_route_semantics.py::test_shared_route_parser_supports_more_gaussian_job_types_and_params",
                ),
            ),
            FeatureSupport(
                area="Molecule specification",
                support="partial",
                scope="Cartesian, integer Cartesian, Z-matrix variables, free delimiters, fragments, and invalid-domain diagnostics.",
                limitations="Supported cases are the explicit molecule-spec grammar variants covered by tests.",
                tests=(
                    "tests/test_gjf_spec_cases.py::test_zmat_labeled_variable_blocks_are_applied_to_atom_lines",
                    "tests/test_gjf_spec_cases.py::test_zmat_unlabeled_variable_blocks_are_applied_by_position",
                    "tests/test_gjf_spec_cases.py::test_free_format_delimiters_are_parsed_and_rendered_canonically",
                    "tests/test_gjf_spec_cases.py::test_integer_cartesian_coordinates_are_not_misparsed_as_zmatrix",
                    "tests/test_gjf_spec_cases.py::test_multifragment_assignments_match_declared_pairs",
                    "tests/test_gjf_spec_cases.py::test_invalid_zmatrix_domains_raise_clear_errors",
                ),
            ),
            FeatureSupport(
                area="Additional sections",
                support="partial",
                scope="ModRedundant, GIC, NBO, structured rendering, and mixed-section diagnostics.",
                limitations="Mixed GIC/ModRedundant sections intentionally fall back to unknown with diagnostics.",
                tests=(
                    "tests/test_gjf_spec_cases.py::test_modredundant_additional_section_detected_as_typed_section",
                    "tests/test_gjf_spec_cases.py::test_gic_additional_section_detected_as_typed_section",
                    "tests/test_gjf_spec_cases.py::test_nbo_additional_section_detected_as_typed_section",
                    "tests/test_gjf_spec_cases.py::test_mixed_gic_and_modredundant_section_falls_back_to_unknown_with_diagnostic",
                    "tests/test_gjf_spec_cases.py::test_render_prefers_structured_additional_sections",
                ),
            ),
            FeatureSupport(
                area="Link1, includes, and writer options",
                support="partial",
                scope="Link1 splitting with blank-line validation; Geom=AllCheck; disk @include expansion; checkpoint propagation during transform.",
                limitations="@include expansion requires source-path context; not every Gaussian writer option is covered.",
                tests=(
                    "tests/test_gjf_spec_strictness.py::test_link1_requires_blank_line_before_separator",
                    "tests/test_gjf_spec_strictness.py::test_geom_allcheck_allows_missing_title_and_molecule_sections",
                    "tests/test_gjf_include_handling.py::test_gjf_disk_parser_expands_at_include",
                    "tests/test_format_transform_output_dir.py::test_format_transform_gjf_chk_propagation_single_file",
                ),
            ),
        ),
    ),
    "g16log": FormatSupport(
        format_id="g16log",
        group="qm-output",
        title="Gaussian output",
        extensions=(".log", ".out", ".g16", ".gal", ".irc", ".gau"),
        read=True,
        write=False,
        registry_role="reader",
        data_level="coordinates and QM results",
        summary="Gaussian output parsing for user-facing QM result attributes such as structure, energies, thermochemistry, vibrations, orbitals, populations, forces, and response properties.",
        features=(
            FeatureSupport(
                area="Metadata, geometry, charge, multiplicity, and status",
                support="fixture-covered",
                scope="Representative Gaussian output files expose atoms, coordinates or standard coordinates, charge, multiplicity, title/options/keywords, running time, and termination status.",
                limitations="Read-only format; output sections outside covered fixtures may remain unstructured.",
                tests=(
                    "tests/test_parsers_smoke.py::test_autoparser_g16log_smoke",
                    "tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_first_frame_parses_link1_header_payloads",
                    "tests/test_g16log_numeric_regression.py::test_g16log_file_running_time_is_accumulated_into_model_field",
                    "tests/test_io_registry_and_batch_more.py::test_g16_file_parser_finalizes_status_from_last_frame",
                ),
            ),
            FeatureSupport(
                area="Model chemistry and route-derived requests",
                support="partial",
                scope="Parsed output frames expose shared Gaussian route semantics, including model chemistry, job types, dispersion, solvation, population requests, and HF handling.",
                limitations="Route-derived semantics are limited to the shared Gaussian route parser coverage; method-specific keywords outside tests may remain raw or unstructured.",
                tests=(
                    "tests/test_gaussian_route_semantics.py::test_g16log_frame_exposes_shared_semantic_route",
                ),
            ),
            FeatureSupport(
                area="Energies",
                support="partial",
                scope="Reference/electronic and method-specific energies are populated on representative output frames; live reference energies take precedence over archive values; terminal archive-only frames do not invent energies.",
                limitations="Only tested energy fields are advertised; not every Gaussian post-HF or correction energy table is guaranteed structured.",
                tests=(
                    "tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures",
                    "tests/test_g16log_parser_v3_regression.py::test_g16log_v3_preserves_live_reference_energy_over_archive_value",
                    "tests/test_g16log_parser_v3_regression.py::test_g16log_v3_does_not_invent_archive_energies_on_terminal_frame",
                ),
            ),
            FeatureSupport(
                area="Thermochemistry",
                support="partial",
                scope="Thermal information includes temperature, molecular mass, moments of inertia, rotational symmetry number, rotational/vibrational temperatures, rotational constants, ZPVE, thermal energy, enthalpy, Gibbs free energy, entropy, and heat capacity where present.",
                limitations="Thermochemistry coverage is fixture-based and tied to frequency-style Gaussian outputs.",
                tests=(
                    "tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_l716_children_decompose_parent_payloads",
                    "tests/test_g16log_parser_v3_render.py::test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry",
                ),
            ),
            FeatureSupport(
                area="Vibrations and IR data",
                support="partial",
                scope="Vibrational frequencies, reduced masses, force constants, IR intensities, per-mode displacement vectors, and imaginary-mode flags are exposed where present.",
                limitations="Raman/VCD and other spectrum variants are not advertised unless backed by explicit tests.",
                tests=(
                    "tests/test_g16log_numeric_regression.py::test_default_g16log_parser_handles_concatenated_frequency_values",
                    "tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_l716_vibration_mode_children_expose_per_mode_semantics",
                ),
            ),
            FeatureSupport(
                area="Molecular orbitals and charge/spin populations",
                support="partial",
                scope="Molecular orbital energies/occupancies and charge/spin population data are exposed on representative closed- and open-shell outputs.",
                limitations="Orbital coefficient matrices and every population analysis flavor are not claimed by this coverage item.",
                tests=(
                    "tests/test_g16log_mo_regression.py::test_default_g16log_parser_preserves_open_shell_mo_counts",
                    "tests/test_g16log_mo_regression.py::test_default_g16log_parser_handles_concatenated_orbital_energies",
                    "tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures",
                    "tests/test_g16log_parser_v3_node_payloads.py::test_g16log_v3_major_components_parse_payloads",
                ),
            ),
            FeatureSupport(
                area="Forces, Hessian, optimization status, and polarizability",
                support="partial",
                scope="Representative outputs preserve force arrays, Hessian presence, geometry optimization status, and polarizability when present.",
                limitations="Only fixture-backed fields are advertised; dipole detail and expanded polarizability sub-blocks are not claimed.",
                tests=(
                    "tests/test_g16log_parser_v3_regression.py::test_g16log_v3_matches_v1_on_representative_fixtures",
                    "tests/test_g16log_numeric_regression.py::test_extract_labeled_float_tokens_handles_concatenated_polarizability_values",
                ),
            ),
            FeatureSupport(
                area="Link1 metadata propagation",
                support="partial",
                scope="Link1 section metadata is propagated to later frames so multi-step Gaussian jobs keep the expected per-frame context.",
                limitations="This item is limited to user-visible per-frame metadata.",
                tests=(
                    "tests/test_g16log_link1_metadata.py::test_g16log_link1_sections_propagate_section_metadata_to_later_frames",
                ),
            ),
            FeatureSupport(
                area="Registry conversion",
                support="supported",
                scope="Parsed Gaussian output can be converted to coordinate and graph formats through the registry.",
                limitations="Conversion depends on successful structure and graph recovery.",
                tests=("tests/test_io_convert_registry_smoke.py::test_registry_convert_g16log_smoke",),
            ),
        ),
    ),
    "fakeg": FormatSupport(
        format_id="fakeg",
        group="qm-output",
        title="Gaussian-like output renderer",
        extensions=(".fakeg",),
        read=False,
        write=True,
        registry_role="file writer",
        data_level="parsed Gaussian output data",
        summary="Gaussian-like text rendering from parsed Gaussian output data.",
        features=(
            FeatureSupport(
                area="Rendered Gaussian-like output sections",
                support="partial",
                scope="Renders Gaussian-like structure, charge/multiplicity, frequency, and thermochemistry sections from parsed Gaussian output data.",
                limitations="Compatibility renderer only; not a byte-for-byte Gaussian log reproduction.",
                tests=(
                    "tests/test_g16log_parser_v3_render.py::test_g16log_v3_render_contains_expected_gaussian_sections",
                    "tests/test_g16log_parser_v3_render.py::test_g16log_file_model_can_render_fakeg_for_all_frames",
                ),
            ),
            FeatureSupport(
                area="Frequency and thermochemistry reparse",
                support="supported",
                scope="Rendered frequency and thermochemistry content can be parsed back into frame fields.",
                limitations="Frame-level fakeg rendering is intentionally unsupported.",
                tests=(
                    "tests/test_g16log_parser_v3_render.py::test_g16log_file_format_transform_registers_fakeg_file_renderer_only",
                    "tests/test_g16log_parser_v3_render.py::test_g16log_rendered_fakeg_can_be_reparsed_with_frequency_and_thermochemistry",
                ),
            ),
        ),
    ),
    "orcainp": FormatSupport(
        format_id="orcainp",
        group="qm-input",
        title="ORCA input",
        extensions=(".inp",),
        read=True,
        write=False,
        registry_role="reader",
        data_level="coordinates and QM input semantics",
        summary="ORCA input parsing with coordinates, multi-job splitting, model chemistry, task-family fixtures, and structured request containers.",
        features=(
            FeatureSupport(
                area="Geometry and job splitting",
                support="partial",
                scope="Direct coordinates, point charges, external xyzfile/pdbfile references, and $new_job splitting.",
                limitations="External geometry references are parsed as references; parser does not require dereferencing external files.",
                tests=(
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_minimal_xyz_parses_atoms_and_coords",
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_point_charges_parse_into_frame_geometry",
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_xyzfile_does_not_require_external_file",
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_new_job_splits_frames",
                ),
            ),
            FeatureSupport(
                area="Model chemistry and options",
                support="partial",
                scope="Method, functional, basis, auxiliary basis, dispersion-as-functional-suffix, mixed basis, print settings, PARAS variables, and scan coordinates.",
                limitations="Unsupported ORCA block options may remain raw/resource fields instead of dedicated semantic containers.",
                tests=(
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_metadata_population_recognizes_d4",
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_mixed_basis_and_output_print_settings",
                    "tests/test_autoparser_orcainp_tmpfile.py::test_autoparser_orcainp_paras_structures_scan_and_resolves_cartesian_variables",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_dispersion_is_functional_suffix",
                ),
            ),
            FeatureSupport(
                area="ORCA 6 manual fixture families",
                support="fixture-covered",
                scope="Single point, SCF stability, optimization, frequency, excited-state, MRCI, and Solvator input examples.",
                limitations="Manual snippets without complete explicit geometry are intentionally excluded from fixture coverage.",
                tests=(
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_single_point_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_scf_stability_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_optimization_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_frequency_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_inputs_include_orca6_manual_fixtures",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_solvator_inputs_include_orca6_manual_fixtures",
                ),
            ),
            FeatureSupport(
                area="Excited-state and multireference requests",
                support="partial",
                scope="Structured excited-state and multireference task semantics, including MRCI multi-job fixtures.",
                limitations="ORCA excited-state and multireference keyword families outside the covered examples may still be raw.",
                tests=(
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_excited_state_manual_fixtures_are_structured",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_manual_fixtures_are_structured",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_mrci_multijob_manual_fixture_splits_frames_and_structures_last_job",
                ),
            ),
            FeatureSupport(
                area="Optimization, coordinates, frequency, and Solvator structures",
                support="partial",
                scope="Optimization constraints, internal coordinates, fragments, NEB, frequency restart, and Solvator examples.",
                limitations="Coverage follows explicit ORCA manual fixtures and targeted regression examples.",
                tests=(
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_constraints_block_is_not_truncated",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_manual_internal_coords_fill_frame_atoms_and_coords",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_fragment_mixed_basis_is_structured",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_optimization_neb_block_and_xtb_geometry_are_structured",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_frequency_restart_example_is_structured",
                    "tests/test_orcainp_fixtures_roundtrip.py::test_orcainp_orca6_solvator_fixture_is_structured",
                ),
            ),
            FeatureSupport(
                area="Writer availability",
                support="intentionally-unsupported",
                scope="ORCA input writer is not registered.",
                limitations="Rendering is disabled until a structured renderer can preserve ORCA input semantics.",
                tests=("tests/test_autoparser_orcainp_tmpfile.py::test_orcainp_writer_is_not_registered",),
            ),
        ),
    ),
    "orcaout": FormatSupport(
        format_id="orcaout",
        group="qm-output",
        title="ORCA output",
        extensions=(".out", ".log", ".orcaout"),
        read=True,
        write=False,
        registry_role="reader",
        data_level="coordinates and QM results",
        summary="ORCA output parsing for user-facing result attributes such as energies, forces, vibrations, populations, solvation, excited states, and response properties.",
        features=(
            FeatureSupport(
                area="Software metadata, frames, and status",
                support="fixture-covered",
                scope="Local and cclib ORCA fixtures expose software name/version, frame counts, and normal termination status across legacy and modern ORCA versions.",
                limitations="Coverage is intentionally fixture-driven; malformed or partial outputs outside this set may remain unsupported.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_versions_cover_legacy_and_modern_orca",
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                ),
            ),
            FeatureSupport(
                area="Energies",
                support="partial",
                scope="Final single-point total energies are exposed in Hartree for fixtures that declare expected final energies.",
                limitations="Only tested ORCA energy fields are advertised; method-specific correction/decomposition tables may remain unstructured.",
                tests=("tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",),
            ),
            FeatureSupport(
                area="Forces and geometry optimization",
                support="partial",
                scope="Gradient/force arrays and geometry optimization status are exposed for covered ORCA gradient and optimization outputs.",
                limitations="Coverage follows the local/cclib fixture inventory and does not imply every ORCA optimizer diagnostic is structured.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts",
                ),
            ),
            FeatureSupport(
                area="Vibrations and vibrational spectra",
                support="fixture-covered",
                scope="Vibration containers are exposed for covered frequency outputs; fixture inventory includes IR and Raman examples.",
                limitations="Feature counts prove fixture presence; only fields asserted in the structured parse contract are guaranteed structured.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad",
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts",
                ),
            ),
            FeatureSupport(
                area="Charge and spin populations",
                support="partial",
                scope="Charge/spin population containers are exposed for covered ORCA outputs that include population sections.",
                limitations="Specific population schemes outside the tested fixtures may remain unstructured.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                ),
            ),
            FeatureSupport(
                area="Solvation",
                support="partial",
                scope="Solvent or solvation model information is exposed on file or frame models for covered CPCM/SMD solvation outputs.",
                limitations="Only fixture-backed implicit-solvation fields are advertised; detailed solvent energy decompositions may remain unstructured.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad",
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts",
                ),
            ),
            FeatureSupport(
                area="Excited states and multireference results",
                support="partial",
                scope="Electronic state containers are exposed for covered TDDFT, ADC2, EOM-CCSD, STEOM-CCSD, STEOM-DLPNO-CCSD, and ROCIS outputs; multireference result containers are part of the structured contract where present.",
                limitations="Coverage is fixture-driven and does not imply every ORCA excited-state or multireference print variant is normalized.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad",
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts",
                ),
            ),
            FeatureSupport(
                area="Response properties",
                support="partial",
                scope="Polarizability is exposed as a structured field for covered outputs; fixture inventory also includes NMR and spin-spin coupling examples.",
                limitations="NMR and spin-spin coupling are currently fixture-anchored and may not yet be lifted into dedicated common fields.",
                tests=(
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_manifest_inventory_is_broad",
                    "tests/test_orca_output_fixtures.py::test_orca_output_structured_parse_contract",
                    "tests/test_orca_output_fixtures.py::test_orca_output_fixture_feature_counts",
                ),
            ),
        ),
    ),
    "cml": FormatSupport(
        format_id="cml",
        group="structure",
        title="CML writer",
        extensions=(".cml",),
        read=False,
        write=True,
        registry_role="file writer, frame writer",
        data_level="graph",
        summary="Chemical Markup Language writer backed by RDKit or OpenBabel molecule rendering.",
        features=(
            FeatureSupport(
                area="File and frame writer",
                support="supported",
                scope="Selected-frame file rendering, negative frame index selection, separate selected-frame blocks, and single-frame rendering through the registry.",
                limitations="No CML reader is registered.",
                tests=(
                    "tests/test_cml_codec.py::test_cml_writer_renders_selected_frames_in_one_file",
                    "tests/test_cml_codec.py::test_cml_writer_supports_negative_frame_id",
                    "tests/test_cml_codec.py::test_cml_writer_can_return_selected_frames_as_separate_blocks",
                    "tests/test_cml_codec.py::test_cml_frame_writer_renders_single_frame",
                ),
            ),
            FeatureSupport(
                area="Rendering engines and validation",
                support="supported",
                scope="RDKit-backed default rendering, OpenBabel-backed rendering, frame-less file rejection, and unsupported-engine diagnostics.",
                limitations="Requires a recoverable RDKit or OpenBabel molecule; arbitrary XML/CML source preservation is not covered.",
                tests=(
                    "tests/test_cml_codec.py::test_cml_writer_supports_openbabel_engine",
                    "tests/test_cml_codec.py::test_cml_writer_requires_frames",
                    "tests/test_cml_codec.py::test_cml_frame_writer_rejects_unknown_engine",
                ),
            ),
        ),
    ),
    "openbabel-fallback": FormatSupport(
        format_id="openbabel-fallback",
        group="special",
        title="OpenBabel fallback reader",
        extensions=(),
        read=True,
        write=False,
        registry_role="special fallback reader",
        data_level="coordinates",
        summary="Fallback reader for unknown extensions when OpenBabel can parse the source file.",
        features=(
            FeatureSupport(
                area="Unknown-extension fallback",
                support="partial",
                scope="Unknown extensions can be read through OpenBabel-compatible formats and converted to an XYZFile/XYZFileFrame model with detected format ID persistence and first-molecule coordinate preservation.",
                limitations="Only the first parsed molecule is converted; output is coordinate-level and OpenBabel availability controls format reach.",
                tests=(
                    "tests/test_structure_format_features.py::test_openbabel_unknown_extension_fallback_keeps_first_molecule_coordinates",
                    "tests/test_io_convert_registry_smoke.py::test_unknown_extension_fallback_smoke",
                    "tests/test_filebatchmodeldisk_codec_filter.py::test_filter_by_codec_id_uses_detected_format_id",
                ),
            ),
        ),
    ),
}

FORMAT_FEATURE_COVERAGE = {
    format_id: report.features for format_id, report in FORMAT_SUPPORT.items()
}
