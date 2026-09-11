# Contributing to MolOP

Choose the developer guide that matches the task:

- [Developer overview](../index.md)
- [Reader/writer plugin development](../extensions/plugins.md)
- [Complete parser contract](../contracts/parser.md)
- [Documentation contributions](documentation.md)
- [Development environment and quality gates](quality.md)

For issues, provide a minimal reproducer, exact command, version, and expected result. For code
changes, keep scope clear and add tests plus paired user documentation for new or changed public
behavior.

## Before opening a pull request

1. Run the focused test for the changed behavior.
2. Run the relevant format-coverage or documentation checks.
3. Run `make check` when source changes cross module boundaries.
4. Review `git diff --check` and confirm that generated files, local logs, and unrelated work are not included.

For a new reader or writer, start with [the complete parser contract](../contracts/parser.md)
and [the plugin guide](../extensions/plugins.md). For a documentation change, maintain the Chinese and
English pages together and follow [the documentation contribution rules](documentation.md).
