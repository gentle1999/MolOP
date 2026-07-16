# Reader/writer plugin development

Before adding a format, confirm that an existing dedicated codec or OpenBabel fallback cannot meet
the requirement.

## Minimum reader responsibilities

1. Declare a stable `format_id`, extensions, priority, and content probe.
2. Read source through the shared parsing lifecycle.
3. Produce at least one valid File/Frame model.
4. Record presence, diagnostics, and supported source evidence for optional scientific fields.
5. Add parsing, error, and detection tests with a real or minimal fixture.

## Minimum writer responsibilities

1. Declare target format ID, extension, and capability level.
2. Implement frame selection, embedding/splitting, and render/write behavior.
3. State coordinate-, graph-, and format-specific information boundaries.
4. Declare dynamic options and CLI help in the registry.
5. Test string rendering, disk output, naming, and failure behavior.

## Registration location

Built-in codecs are registered through `src/molop/io/codecs/catalog.py` and their format packages.
Do not maintain a second format list in prose. Capability status comes from
`tests/format_feature_coverage/support_matrix.py` and generated sections on format pages.

## Required contract

The [Parser contract](parser-contract.md) preserves full class hierarchy, registration order, probe,
source span, diagnostics, transform, and test constraints.

Run the [quality gates](quality.md) and update paired format pages and user examples with the
implementation.
