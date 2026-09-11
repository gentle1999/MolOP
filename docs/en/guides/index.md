# Workflows

This section is organized by the task to complete. Each page starts with a small runnable example,
then explains common variants, output shapes, and information boundaries.

## Everyday operations

| Task | Page | Use it for |
| --- | --- | --- |
| Parse files | [Parse files](parsing.md) | Single files, globs, mixed inputs, and parse options |
| Read results | [Read calculation results](results.md) | Energies, thermochemistry, vibrations, orbitals, populations, and NMR |
| Build summaries | [Batch summaries](batch.md) | DataFrames, CSV/JSON output, and parallel processing |
| Filter and select | [Filter and select](filtering.md) | Normal, optimized, transition-state, and custom conditions |
| Convert and export | [Convert and export](conversion.md) | XYZ, SDF, Gaussian/ORCA input, and writing files |
| Recover structures | [Structure recovery](structure-recovery.md) | Coordinate topology, dative bonds, and TS candidates |
| CLI recipes | [CLI task recipes](cli-recipes.md) | Reproducible shell command chains |
| Diagnose issues | [Troubleshooting](troubleshooting.md) | Paths, formats, fields, conversion, and native libraries |

## Complete task tutorials

When a task combines several operations, open the [complete task tutorials](../tutorials/index.md).
Each tutorial starts from input files and finishes with a selection, summary, or export workflow that
can be adapted to a real project.

## Choose a starting page

- Reading one file: start with [Parse files](parsing.md).
- Already have a batch: read [calculation results](results.md), [batch summaries](batch.md), or [filter and select](filtering.md).
- Writing a new file: use [Convert and export](conversion.md), then check the target [format reference](../reference/index.md).
- Tuning logging, parallelism, or topology recovery: see the [configuration reference](../reference/config.md).
