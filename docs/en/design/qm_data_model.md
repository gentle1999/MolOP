# QM Common Data Model Design

This page defines the boundaries of MolOP's common QM input and result containers. The goal is to let Gaussian, ORCA, and future program parsers project into one shared model without letting duplicated fields overwrite each other.

## Layers

The common QM model has three layers:

1. Compatibility fields: `method`, `functional`, `basis_set`, `keywords`, `request_num_cpu`, `request_memory`, and `resources_raw`. They remain available for existing API users, but new parser logic should not treat them as authoritative.
2. Input request containers: `model_chemistry`, `task_requests`, `resource_request`, `excited_state_requests`, and `multireference_requests`. Parsers should write these structured fields first.
3. Result containers: `energies`, `molecular_orbitals`, `vibrations`, `single_point_properties`, `electronic_states`, and `multireference_result`. These describe calculation output and should not infer input requests in reverse.

## Source Of Truth

Structured containers are authoritative. After parsing structured semantics, parsers should call `project_common_qm_fields()` to project structured containers onto compatibility fields:

- `model_chemistry.raw_keywords` -> `keywords`
- `model_chemistry.method_family/method` -> `method`
- `model_chemistry.functional` -> `functional`
- `model_chemistry.basis_set` -> `basis_set`
- `resource_request.num_cpu/memory/raw` -> `request_num_cpu/request_memory/resources_raw`

The legacy `refresh_common_qm_containers()` method is kept for existing callers only. It backfills empty structured fields from compatibility fields and must not overwrite already parsed structured semantics.

## Input Requests

`QMModelChemistry` describes model chemistry:

- Functionals and dispersion belong to the same model-chemistry expression. Gaussian and ORCA both append the dispersion suffix to `functional`, such as `B3LYP-GD3BJ`.
- `dispersion_correction` also stores the normalized dispersion label, such as `GD3BJ` or `D4`.
- `basis_sets` stores mixed basis sets, auxiliary basis sets, ECPs, and atom-scoped overrides. `basis_set` stores only the primary global orbital basis.
- `solvation_model` and `solvent` describe the input request and are separate from output-side `BaseCalcFrame.solvent`.

`QMTaskRequest` describes only generic task types and generic task options, such as `sp`, `opt`, `freq`, `excited_state`, and `multi_reference`. Domain-specific excited-state and multi-reference parameters belong in their dedicated containers.

`ExcitedStateRequest` describes requested roots, number of states, spin channels, root following, and requested properties.

`MultireferenceRequest` describes requested multi-reference method, reference method, active space, state blocks, thresholds, and corrections.

## Results

`ElectronicStates` stores state/root-resolved output data for excited-state, MRCI, CASSCF, and related calculations. It does not replace `energies`, which remains the overall or ground-state energy container.

`MultireferenceResult` stores multi-reference output data, including method, active space, state-resolved results, corrections, and diagnostics. It does not replace input-side `MultireferenceRequest`.

## Parser Rules

New parser logic should follow this flow:

1. Parse program-specific syntax into program-specific semantic models.
2. Build common structured containers.
3. If resource or keyword data was initially parsed into legacy fields, use `backfill_common_qm_containers_from_legacy()` only to fill empty structured fields.
4. Call `project_common_qm_fields()` to update compatibility fields.
5. Test structured fields and compatibility-field projections together.

Do not store new semantics only in `keywords`, `resources_raw`, or `options`. `options` is reserved for program-specific details that do not yet have a common field.
