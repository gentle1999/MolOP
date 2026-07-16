# Task tutorials

These Markdown tutorials cover complete workflows without requiring a notebook.

| Research task | Input | Output |
| --- | --- | --- |
| [Summarize Gaussian, ORCA, and xTB](qm-summary.md) | Mixed calculation outputs | Unified status and energy table |
| [Export energy and thermochemistry CSV](energy-csv.md) | Gaussian/ORCA results | Analysis-ready CSV |
| [Select optimized results and transition states](select-results.md) | Optimization and frequency outputs | Selected paths and summary |
| [Export structures and next-step inputs](export-inputs.md) | Completed calculation outputs | XYZ/SDF/GJF/ORCA input |

## Optional notebooks

Notebooks provide longer interactive exploration. Documentation builds do not execute them, and
they are not prerequisites for the workflows above:

- [Gaussian parse and inspect](../examples/01-gaussian-parse-and-inspect.ipynb)
- [Batch summary, filter, and select](../examples/02-batch-summary-filter-select.ipynb)
- [Transform and export](../examples/03-transform-and-export.ipynb)

Adding a reader or writer is a developer task. See [Plugin development](../developer/plugins.md).
