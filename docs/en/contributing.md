# Contributing to MolOP

Choose the developer guide that matches the task:

- [Developer overview](developer/overview.md)
- [Reader/writer plugin development](developer/plugins.md)
- [Complete parser contract](developer/parser-contract.md)
- [Documentation contributions](developer/documentation.md)
- [Development environment and quality gates](developer/quality.md)

For issues, provide a minimal reproducer, exact command, version, and expected result. For code
changes, keep scope clear and add tests plus paired user documentation for new or changed public
behavior.

## Before opening a pull request

1. Run the focused test for the changed behavior.
2. Run the relevant format-coverage or documentation checks.
3. Run `make check` when source changes cross module boundaries.
4. Review `git diff --check` and confirm that generated files, local logs, and unrelated work are not included.

For a new reader or writer, start with [the complete parser contract](developer/parser-contract.md)
and [the plugin guide](developer/plugins.md). For a documentation change, maintain the Chinese and
English pages together and follow [the documentation contribution rules](developer/documentation.md).
