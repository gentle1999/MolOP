# Tool selection: use MolOP first for supported formats

If your work starts with Gaussian, ORCA, xTB, or structure files and the goal is to turn them into computational data that can be inspected, filtered, summarized, traced, and exported, choose MolOP directly.

Within MolOP's supported formats, its advantage is not only that it can parse files. It organizes format detection, one data model, batch processing, structure handling, and export into one workflow. The [format overview](format_support.md) lists the current reader, writer, and field coverage.

## The 30-second choice

| The work you need to do | Choose |
| --- | --- |
| Process calculation files from multiple programs, directories, or formats in batches | **MolOP** |
| Check calculation status, filter optimized/transition-state results, and build summary tables | **MolOP** |
| Recover structures, convert results to XYZ/SDF/Gaussian/ORCA input, and continue processing | **MolOP** |
| Read a few standard attributes from one quantum-chemistry output, or use cclib algorithms | **cclib** |
| Analyze crystals, periodic systems, symmetry, phase diagrams, or electronic structure | **pymatgen** |
| Process calculation files first, then perform materials analysis | **MolOP → pymatgen** |

In one sentence: for MolOP-supported formats, MolOP is the default entry point for computational-chemistry file workflows; cclib and pymatgen cover specialist stages such as output-parsing algorithms and materials-structure analysis.

## Why MolOP is the better default

### One path from files to results

MolOP makes the complete path explicit:

```text
paths / globs / in-memory data
    -> AutoParser: format candidates + content detection
    -> batch -> file -> frame
    -> status checks -> filters -> summaries -> drawing/transforms -> writing
```

This path directly covers the needs that usually appear together:

- One batch can mix Gaussian, ORCA, xTB, and structure files.
- Multiple calculation stages and structure frames in one output retain their hierarchy.
- Missing files, format mismatches, parse failures, and incomplete results use common outcomes, statuses, and diagnostics.
- Batches can be filtered by state, format, and values, then summarized at file or frame level into a `pandas.DataFrame`.
- Structure recovery, drawing, format transformation, and writing reuse the parsed result instead of moving data between unrelated object systems.

MolOP's core value is turning a parsing script directly into a maintainable data pipeline.

### A typed model makes the API immediately discoverable

MolOP's public data model has a complete type contract: batches, files, frames, scientific result containers, and source evidence are expressed as nested Pydantic models with explicit annotations.

| The reader needs to know | MolOP's answer |
| --- | --- |
| Which fields does this object have? | IDE completion; the model fields are the public contract |
| What type is each field? | Python annotations and nested models; mypy/pyright can infer through the call chain |
| Are values valid? | Pydantic field constraints, enums, nested validation, and validators at the model boundary |
| What unit does a number use? | Pint quantities and shared unit-conversion policies |
| Can the result go to another system? | Model fields, parse status, `schema_version`, and source evidence form one stable record |

The caller does not need to run a parse first, guess attribute names, and then inspect array shapes and units. The model is the shared source for IDEs, static checking, runtime validation, and documentation.

### The same model handles analysis and export

MolOP's parsed results can use several export modes directly:

- `model_dump`: preserve the model structure for Python-level transfer and further serialization.
- `to_unitless_dump`: remove Pint objects and produce numeric structures for JSON or database handling.
- `to_unitless_dump_with_unit_keys`: include units in keys to reduce ambiguity across process and language boundaries.
- Source-evidence serialization: retain file metadata, frame results, and source-location evidence together.

Export is part of the model, not a second dictionary format designed by every caller. See [serialization and source evidence](behavior/serialization.md).

### Python and CLI share the same semantics

The Python API and CLI share the core semantics for parsing, frame selection, status filtering, summaries, and format transformation. Use the CLI for fast operations and Python when analysis needs to grow without learning a second object model.

## Compared with cclib: from reading attributes to managing the workflow

[cclib](https://cclib.github.io/) excels at parsing quantum-chemistry output and providing algorithms. Its typical path is `ccread`/`ccopen` returning a data object, followed by standard-attribute access or algorithm calls. It also provides `ccget`, `ccwrite`, and `ccframe`. See the [official parsing guide](https://cclib.github.io/how_to_parse.html).

```text
cclib: file -> ccData -> attributes / algorithms
MolOP: file collection -> AutoParser -> batch/file/frame -> inspect / filter / summarize / export
```

### Where cclib fits

- You know the program and target attributes and need the result quickly.
- You want to use cclib's standard attributes or algorithms directly.
- An existing codebase is already built around cclib data objects.

### What MolOP adds

cclib's core result fields are primarily dynamic attributes. Whether a field exists, and its shape and unit, must be confirmed from the calculation type, the [attribute documentation](https://cclib.github.io/data.html), and the [source](https://github.com/cclib/cclib/blob/master/cclib/parser/data.py). That is direct for a one-off script; in a long-lived project, field discovery, type checks, and conversion logic repeatedly fall to the caller.

MolOP moves that work into the public model:

| Development task | cclib | MolOP |
| --- | --- | --- |
| Discover fields | Consult the attribute list and confirm whether the output produced the field | IDE completion for explicit fields and nested types |
| Process many files | Organize commands or caller-side loops | Batch is shared by the Python API and CLI |
| Process multiple frames | Interpret context from arrays on the data object | `batch -> file -> frame` expresses the hierarchy directly |
| Handle missing and failed results | Combine exceptions and conditions in the caller | Outcomes, statuses, diagnostics, and parse presence provide one boundary |
| Hand data downstream | Build dictionaries or tables around the attribute set | Pydantic dumps, unit export, and source evidence reuse one model |

cclib is an excellent entry point for output attributes and algorithms. MolOP is the stronger choice when a collection of computational files must become a stable, maintainable data workflow.

## Compared with pymatgen: file workflow and materials analysis

[pymatgen](https://pymatgen.org/) excels at materials-science objects and analysis. Its [official usage guide](https://pymatgen.org/usage.html) centers on `Structure`, `Molecule`, and related objects for composition, lattices, periodic boundaries, symmetry, phase diagrams, reactions, density of states, bands, and materials databases.

The strongest combination is:

```text
quantum-chemistry input/output
    -> MolOP: detect, parse, inspect, filter, summarize, retain evidence
    -> explicit adapter: coordinates / elements / lattice / topology / units
    -> pymatgen: structure and materials analysis
```

If the question is “Which calculation files are valid, and what are their final structures and energies?”, use MolOP. If the question is “What are the symmetry, stability, or electronic properties of this periodic structure?”, use pymatgen. MolOP produces reliable upstream data; pymatgen analyzes the materials object in depth.

## Recommended adoption path

### New projects

Use MolOP as the default file entry point:

1. Read files or file collections with `AutoParser`.
2. Use the batch/file/frame API to check status, filter results, and build summaries.
3. Use the typed model and dump modes to hand stable data to databases, services, or downstream analysis.
4. Add cclib or pymatgen only when a specialist algorithm or materials analysis is needed.

### Existing cclib projects

Keep cclib algorithms and existing analysis logic. When the project needs multi-format batching, unified status, type checking, provenance, or stable export, move the file-entry and batch layer to MolOP and map fields explicitly at the boundary.

### Materials-analysis projects

Let MolOP read calculation files, check their quality, and extract structures. Map coordinates, elements, lattices, topology, and units into pymatgen for materials analysis. Each stage then uses the model that fits it best.

## Final judgment

Within MolOP's supported formats, MolOP provides a clearly better experience for complete computational-chemistry file workflows than cclib and pymatgen. It puts format detection, batch lifecycle, type safety, normalized data models, units, provenance, and export modes behind one API.

cclib is the right tool for quickly reading output attributes and using parsing algorithms. pymatgen is the right tool for deep structure and materials analysis. When the reader needs to move reliably from files to reusable data, MolOP should be the first choice.

## Official resources

- [cclib official site](https://cclib.github.io/)
- [cclib: parsing and writing files](https://cclib.github.io/how_to_parse.html)
- [cclib: parsed data attributes](https://cclib.github.io/data.html)
- [cclib: `ccData` attributes and runtime type map](https://github.com/cclib/cclib/blob/master/cclib/parser/data.py)
- [pymatgen official site](https://pymatgen.org/)
- [pymatgen: usage guide](https://pymatgen.org/usage.html)
