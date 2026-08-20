# Structure recovery

Recover bonds, bond orders, formal charges, and radical states from elements, three-dimensional
coordinates, charge, and multiplicity.

## Metal-complex example

The example is a Gaussian 16 single-point calculation:
[download `mn_complex_sp.log`](../../assets/examples/mn_complex_sp.log). The source has no molecular
graph; MolOP invokes MolGR when `frame.rdmol` is first accessed.

```python
from rdkit import Chem

from molop import AutoParser

frame = AutoParser("mn_complex_sp.log", n_jobs=1)[0][-1]
mol = frame.rdmol

if mol is None:
    raise RuntimeError("structure recovery failed")

mn = next(atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Mn")
dative_count = sum(bond.GetBondType() == Chem.BondType.DATIVE for bond in mol.GetBonds())

print(frame.formula)
print(f"{mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")
print(f"Mn charge {mn.GetFormalCharge():+d}, degree {mn.GetDegree()}")
print(f"{dative_count} dative bonds")
print(frame.topology_reconstruction_backend, frame.topology_reconstruction_status)
```

??? example "Real output and molecular graph"

    ```text
    C12H15MnO3P+
    32 atoms, 36 bonds
    Mn charge +1, degree 8
    8 dative bonds
    cpp succeeded
    ```

    ![Raw coordinates, RDKit distance connectivity, and the MolGR reconstruction for a metal complex](../../assets/examples/mn_complex_graph_reconstruction.svg)

    | Treatment | Bonds | Mn coordination | Result |
    | --- | ---: | --- | --- |
    | Raw XYZ | 0 | None | Elements and coordinates only |
    | RDKit `DetermineConnectivity` | 36 | 8 `SINGLE` bonds | Distance connectivity without bond-order semantics |
    | RDKit `DetermineBonds(charge=1)` | - | - | Raises `ValueError` because Mn has no predefined valence |
    | MolGR | 36 | 8 `DATIVE` bonds | Also recovers ligand bond orders, formal charges, and coordination bonds |

    The SVG is generated from the three real molecular objects for the same calculation with
    MolOP's default `rdkit-dof` drawer. The script does not edit SVG paths by hand.

The distinction is not whether 36 atom pairs can be connected, but which chemical semantics are
recovered. RDKit can connect atoms by distance but cannot assign a valence to Mn. MolGR represents
all eight Mn-C bonds as coordination bonds while assigning ligand bond orders and charges.

## Configure recovery globally

Recovery settings belong to the process-wide `molopconfig`, not to an individual `Molecule`:

```python
from molop import molopconfig

molopconfig.graph_reconstruction_backend = "python"  # default: "cpp"
molopconfig.make_dative_bonds = False  # default: True
molopconfig.make_stereochemistry = False  # default: True
molopconfig.prewarm_topologies = True  # default: False
```

Set these values before the first access to `rdmol` for a frame. Lazy reconstruction reads the
current global values and records the effective settings on the frame. Restore the defaults or use a
separate process when different workflows require different policies.

## Optional native parallel prewarming

Graph-dependent operations use lazy reconstruction by default. MolOP's process-parallel entry points
use spawn-like `loky` workers, so the worker that consumes a graph can build it without a separate
parent-process warmup step. Set `molopconfig.prewarm_topologies = True` when a workflow benefits from
an eagerly populated parent-side cache. Existing graphs and already-attempted frames are skipped, and
frames with different backend, dative-bond, or stereochemistry settings are submitted as separate
homogeneous batches so provenance remains correct.

When enabled, prewarming is used by:

- batch `format_transform()` for graph writers (`sdf`, `smi`, `cml`);
- Gaussian `format_transform("gjf", add_gjf_connectivity=True)`;
- frame-level `to_summary_df()` for a file or a batch;
- trajectory, vibration, and TS-vibration animation;
- TS pre/post candidate inference and its summary/export paths remain lazy: these endpoints are
  generated dynamically from each TS frame, so a fresh loky worker may reconstruct temporary
  candidates instead of receiving a parent-built candidate map.
- `filter_custom()` and `groupby()` callbacks. Because these callbacks are arbitrary, MolOP warms
  every eligible frame in the input snapshot before dispatching them to loky, even when a particular
  callback happens to inspect metadata only.
- file-batch parsing with `capture_source_evidence=True` deliberately stays in one process. Those
  parsers create and inspect frames while source spans are being attached, so their graphs cannot be
  completely prewarmed before dispatch.

Operations with an enumerable source-frame set prewarm those frames with MolGR's native batch API in
the parent before outer joblib/loky processes start only when the option is enabled. Operations that
generate temporary graph candidates while running, such as TS endpoints, keep the normal public call
path and may perform single-molecule lazy reconstruction inside a fresh loky worker. The worker process
must be created by spawn/loky; a fork child is rejected by MolGR's PID guard before native code runs.

This boundary also applies to explicit `threading` backends and nested parallel calls: if a callback
tries to trigger an unprewarmed reconstruction while MolOP parallel work is active, MolOP fails fast
with a concurrency error instead of waiting on a self-deadlocking pool. The molecule remains
retryable after this scheduling error.

The native guard covers both loky and joblib's multiprocessing child processes. It also refuses to
start while an unrelated external loky executor still has pending work; consume or close that result
before starting a graph-dependent operation. An idle reusable pool is drained before MolGR starts.

The guard also rejects live Python child processes that are not managed by MolOP, including external
`multiprocessing` and `loky.ProcessPoolExecutor` workers. Join or close those processes first; without
that boundary MolGR's native pool cannot be proven independent from them.

Do not call `os.fork()` after MolGR native prewarming, or switch a graph-dependent task to a
fork-based `multiprocessing` backend. Forking can copy native runtime state that Python cannot
fully guard. MolOP-managed process-parallel entry points explicitly select joblib's `loky` backend;
loky uses an independent-interpreter, spawn-like boundary and does not inherit the parent's MolGR
native runtime through POSIX `fork`. MolOP deliberately does not mutate Python's global
`multiprocessing` start method or wrap joblib's legacy `multiprocessing` backend with standard
`spawn`, which would disrupt notebooks, interactive entry points, and caller-owned process
configuration. MolOP entry points normalize only joblib's legacy `multiprocessing` configuration to
loky; outer threading and third-party backends are preserved, and an explicit `backend` argument on
the MolOP parallel API remains the compatibility override. Keep graph-dependent work on loky and
workers read prewarmed source-frame caches where available and may lazily reconstruct dynamic
single-molecule candidates. The native batch adapter itself remains a parent-process operation and
rejects invocation from a worker, preventing nested MolGR native pools.

The native batch adapter validates request identity and completion. Duplicate or unknown results are
treated as a batch error; if the iterator stops early, every missing frame is marked `failed` and will
not retry native reconstruction from a worker. Closing or exhausting an outer joblib iterator releases
its process-pool guard. MolGR does not currently expose a Python timeout or hard cancellation hook, so
a hung native call remains a process-level failure boundary rather than something MolOP can interrupt
in the same process.

`suspicious_fallback` keeps the private native molecule only as raw evidence and never invents a
disconnected graph from empty ordinary topology fields. If a user-defined serializer drops that
private cache, a worker fails closed with an empty graph instead of fabricating an incorrect topology.

When enabled, the operation's existing `n_jobs` value controls the native prewarm limit for enumerable
batch summaries and format conversion. Coordinate-only `xyz`, ORCA input, and Gaussian coordinate
renders keep lazy graph behavior and do not force reconstruction merely because they are converted.

## Status values

| Status | Meaning |
| --- | --- |
| `provided` | The source format provided bonds, formal charges, and radicals |
| `succeeded` | MolGR obtained a normal coordinate-derived candidate |
| `suspicious_fallback` | A usable fallback candidate needs review |
| `failed` | No RDKit molecule could be built; `rdmol` is `None` |
| `None` | Topology access has not run yet, or no status applies |

When MolGR reports a non-fatal item failure or suspicious fallback,
`frame.topology_reconstruction_diagnostics` retains its stable `code`, `stage`, `backend`,
`counts`, `details`, and cause fields. This ordinary dictionary is preserved when a prewarmed frame
is serialized to a loky worker.

## Why review matters

Recovering topology from geometry is not unique. Close contacts, ion pairs, radicals, metal
complexes, and distorted geometries can admit several plausible bonding assignments. A
`suspicious_fallback` may still be useful, but should not enter a high-confidence database without
inspection.

## Provide charge and multiplicity

Override missing source values at the parser entry point:

```python
batch = AutoParser(
    "radical.xyz",
    total_charge=0,
    total_multiplicity=2,
)
mol = batch[0][0].rdmol
```

These values affect recovery. Do not guess them merely to make reconstruction succeed.

## Check before exporting

```python
frame = batch[0][-1]
mol = frame.rdmol

if mol is not None and frame.topology_reconstruction_status != "suspicious_fallback":
    print(frame.to_canonical_SMILES())
```

??? example "Output shape (input-dependent)"

    ```text
    <canonical SMILES>
    ```

Graph-level SDF, SMILES, and CML writers need a usable molecular graph. XYZ and coordinate-mode
Gaussian/ORCA inputs can still export coordinates without preserving full topology semantics.

## Next steps

- [Convert and export](conversion.md)
- [SDF/MOL format](../reference/formats/sdf.md)
- [Structure API](../reference/api/structure.md)
