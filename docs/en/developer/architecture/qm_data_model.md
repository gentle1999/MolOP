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

The legacy `refresh_common_qm_containers()` method is kept for existing callers only. It backfills empty
structured fields from compatibility fields and then projects the structured values back to compatibility
fields. It never lets flat fields overwrite already parsed structured semantics. To change normalized
semantics, modify the structured container and call `project_common_qm_fields()`.

## Input Requests

`QMModelChemistry` describes model chemistry:

- `method` is the stable concrete electronic-structure method used for matching
  and protocol identity. Gaussian `R`, `U`, and `RO` prefixes are removed from
  this value and retained separately in `spin_treatment`; the original route
  token remains available in Gaussian's semantic route evidence.
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

## Result semantics across formats

These contracts apply to Gaussian logs, formatted checkpoints, ORCA logs and xTB
stdout. Tests live in `tests/test_qm_result_semantics.py` and the format-specific
numerical suites.

- **State identity:** `state_index` is an internal, unique non-negative index;
  `root` is the program's root number. Set `is_ground_state` explicitly: `True`
  identifies the ground state, `False` identifies an excited state, and `None`
  leaves identity unknown. List position and root zero do not identify a ground
  state. An excited-only collection exports every known transition. Its initial
  internal state is `None` unless a ground-state record exists; printed initial
  labels/root numbers remain in transition properties.
- **Missing observations:** unknown orbital energies, bond orders and coupling
  matrix entries remain numeric `np.nan` internally. Unknown occupations and
  symmetry labels use `None`. Only these explicitly nullable fields mask NaN
  for public JSON exports (`null`); binary ndarray exports preserve NaN. Infinities
  still fail JSON validation. A printed zero remains zero. Each `BondOrderMatrix` in `BondOrders` retains
  its print threshold and missing-entry reason alongside the matrix; RDKit
  embedding skips unobserved pairs.
- **Spin channels:** the common orbital validator aligns unequal or absent spin
  channels with unknown entries, without assuming restriction. Parsers identify
  shared restricted references explicitly. `orbital_basis`, `occupancy_origin`
  and `total_occupancies` distinguish printed spin occupations from derived
  shared-spatial occupations. Closed-shell totals split equally; evidenced
  integer high-spin totals can split into alpha/beta; other open-shell totals
  retain unknown spin occupations. xTB's omitted orbital ranges never imply zero
  occupation. Printed HOMO/LUMO labels populate `reported_homo_index` and
  `reported_lumo_index`; fractional occupations alone do not define frontiers
  or SOMOs.
- **Modes:** `mode_indices` are stable internal indices starting at zero when
  parsed. `source_mode_indices` separately retain printed numbers, or `None`
  where the source supplies only an array. Filtering retains both identities.
  Single `Vibration` records and IR/Raman bands carry source program, orientation,
  normalization and mass-weighting metadata. Source-normalized Cartesian vectors
  are not rescaled to make different programs appear equivalent.
- **Thermochemistry:** corrections and totals retain their existing definitions
  (`U_T = Eel + TCE`, `H_T = Eel + TCH`, `G_T = Eel + TCG`). The container also
  retains temperature, pressure, `entropy_model`, frequency scaling, excluded-mode
  cutoff, HO/rotor interpolation reference and imaginary-frequency cutoff when
  available. These are distinct conditions; a rotor interpolation reference is
  not a frequency-exclusion threshold. Property-bundle exports retain the conditions.
  Gaussian FCHK totals alone do not establish an entropy model.

Gaussian Raman data retain plane-polarized and unpolarized depolarization
ratios separately in `raman_depolarization_ratios_polarized` and
`raman_depolarization_ratios_unpolarized`. The unsuffixed ratio field is
reserved for source tables with unspecified polarization. xTB PTB calculations
retain distinct `frequency_method` and `raman_activity_method` labels.
FCHK thermal corrections carry their equations in `derivations`.
See [data containers and cross-QM semantics](../../reference/data_model.md) for
the public container overview and user-facing comparison rules.

Migration: callers constructing electronic states must set `is_ground_state`
explicitly. Callers needing original vibration numbers must use
`source_mode_indices` (or `Vibration.source_mode_index`), not `mode_indices`.
Do not replace missing matrix entries with zero when comparing programs.

## Extensible bond-order records

Like atomic populations, `BondOrders.bond_orders` maps stable identifiers to typed
`BondOrderMatrix` records. Each record stores its extensible `scheme` string,
`values`, optional spin channel, source label, print threshold, missing-entry reason
and program-specific metadata. No parallel dictionaries link provenance by name.

```python
import numpy as np
from molop.io.base_models.DataClasses import BondOrderMatrix, BondOrders

orders = BondOrders(bond_orders={
    "mayer_bond_order": BondOrderMatrix(
        scheme="mayer",
        values=np.array([[np.nan, 0.9], [0.9, np.nan]]),
        source_label="Mayer bond orders larger than 0.1",
        print_threshold=0.1,
        missing_reason="unprinted atom pairs",
    ),
})
matrix = orders["mayer_bond_order"].values
```

New schemes need only another dictionary entry. Serialization, property bundles
and RDKit embedding do not enumerate supported scheme names. Matrices must be
non-empty square arrays of the same shape; a frame with geometry also checks them
against its source atom count. Access records through `bond_order_items()`,
`bond_order_names`, `get_bond_order()` or indexing. RDKit skips missing observations
and preserves observed zero values on existing bonds.

Migrate `BondOrders(mayer_bond_order=matrix)` to the dictionary constructor above,
and replace `frame.bond_orders.mayer_bond_order` with
`frame.bond_orders["mayer_bond_order"].values`. Old flat fields are rejected rather
than silently ignored. Serialized matrices live under `bond_orders[name]["values"]`;
NaN entries still export as JSON `null`.

## Parser Rules

New parser logic should follow this flow:

1. Declare a canonical `format_id`, implement `_quick_check_file_format()` and `_locate_segments()` on the file parser, and locate original source ranges directly.
2. Implement `_parse_artifact_metadata()` or `_parse_segment_metadata()` only when that scope has metadata.
3. Parse program-specific syntax into semantic models and build common structured containers.
4. If resource or keyword data was initially parsed into legacy fields, use `backfill_common_qm_containers_from_legacy()` only to fill empty structured fields.
5. Call `project_common_qm_fields()` and test structured fields and compatibility projections together.

Do not reintroduce `_split_file()` or `_parse_metadata()`. The locator is the
only source of truth for segment/frame boundaries and the source slices sent to
the frame parser.

Do not store new semantics only in `keywords`, `resources_raw`, or `options`. `options` is reserved for program-specific details that do not yet have a common field.
