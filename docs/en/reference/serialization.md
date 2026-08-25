# Source Evidence and Serialization

Export file metadata, frame results, and source evidence to a database or audit system. For ordinary
analysis, use `AutoParser`, `to_summary_df`, or frame fields; use this page when parsed facts must be
preserved across systems.

## Export payloads

```python
from molop import AutoParser

chem_file = AutoParser(
    "water_mp2.out",
    capture_source_evidence=True,
    release_file_content=True,
    n_jobs=1,
)[0]
file_payload = chem_file.to_unitless_dump_with_unit_keys(exclude_none=True)
frame_payloads = [
    frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    for frame in chem_file
]

print(file_payload["schema_version"], file_payload["source_format"])
print(len(frame_payloads), frame_payloads[0]["file_frame_index"])
```

??? example "Real output"

    ```text
    molop-calculation-export-v1 orcaout
    1 0
    ```

Keep file and frame payloads separate. The file payload carries provenance and file-level metadata;
the frame payload carries coordinates, topology, scientific results, and frame-level source locators.

## Values and units

`to_unitless_dump_with_unit_keys()` converts Pint quantities to magnitudes and puts canonical units in
keys such as `Energy.total_energy.hartree` or `coords (angstrom)`. Arrays become lists by default;
pass `array_mode="ndarray"` when a database stores numeric sidecars.

| Payload | Purpose | Stable fields |
| --- | --- | --- |
| File payload | File identity, provenance, and file metadata | `schema_version`, `source_format`, `parser_provenance` |
| Frame payload | Structure and calculation results | `file_frame_index`, `source_span`, `parse_presence` |
| Source evidence | Source locations and content digests | `source_segments`, `source_block_sha256` |

## Indexing and boundaries

- Use `file_frame_index` for stable source-file frame order.
- `frame_id` belongs to the retained frame collection; it can change after selection or reordering and
  is not a global identity.
- `capture_source_evidence=True` requests additional source locators and parser provenance.
- Source-evidence capture does not read `frame.rdmol` or reconstruct a molecular graph; topology
  remains lazy until a graph-dependent operation requests it.
- MolOP supplies parsed facts; the receiving system owns database identity, payload validation, array
  encoding, artifact verification, and admission/QC policy.

See the [parser contract](../developer/parser-contract.md) for source spans and parser lifecycle rules.

## Related pages

- [API Contracts](api_contracts.md)
- [Find fields by scientific property](model_fields.md)
- [Format overview](format_support.md)
