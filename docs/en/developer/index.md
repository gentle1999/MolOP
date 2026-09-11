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
| Add a reader or writer | [Plugin development](extensions/plugins.md) |
| Read the complete parser/model constraints | [Parser contract](contracts/parser.md) |
| Change documentation | [Documentation contributions](contributing/documentation.md) |
| Run tests and quality gates | [Development environment and quality gates](contributing/quality.md) |
| Understand the common QM model | [Common QM data model](architecture/qm_data_model.md) |
| Understand the ORCA input model | [ORCA input model](architecture/orca_input_model.md) |
| Understand source lifecycle | [API contracts](contracts/api.md) |

User workflows must not depend on this section. Public API changes need corresponding user-guide
and executable-example updates.
