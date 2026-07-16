# Developer documentation

This section targets contributors changing MolOP source, extending readers/writers, or maintaining
data models.

## Responsibility boundaries

```text
source -> codec/reader -> common File/Frame models -> summary/transform/DTO
                       ↘ format-specific models
```

- A reader extracts source evidence and populates common or format-specific models.
- Common models define cross-format scientific results and units.
- A writer promises only information expressible by its target format.
- Database identity, admission policy, and reaction-path construction do not belong in a parser.

## Development entry points

| Task | Documentation |
| --- | --- |
| Add a reader or writer | [Plugin development](plugins.md) |
| Read the complete parser/model constraints | [Parser contract](parser-contract.md) |
| Change documentation | [Documentation contributions](documentation.md) |
| Run tests and quality gates | [Development environment and quality gates](quality.md) |
| Understand the common QM model | [Common QM data model](../design/qm_data_model.md) |
| Understand the ORCA input model | [ORCA input model](../design/orca_input_model.md) |
| Understand source lifecycle | [API contracts](../reference/api_contracts.md) |

User workflows must not depend on this section. Public API changes need corresponding user-guide
and executable-example updates.
