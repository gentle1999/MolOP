# Structure recovery

Infer bonds, bond orders, formal charges, and radical states from a frame containing only elements
and three-dimensional coordinates.

## Shortest example

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
mol = frame.rdmol

if mol is None:
    raise RuntimeError("structure recovery failed")

print(mol.GetNumAtoms(), mol.GetNumBonds())
print(frame.smiles)
print(frame.topology_reconstruction_status)
```

Shared example output:

```text
3 2
[H]O[H]
succeeded
```

When the source does not provide a molecular graph, accessing `frame.rdmol` asks MolGR to recover
one from the coordinates.

## Status values

| Status | Meaning |
| --- | --- |
| `provided` | The source format provided bonds, formal charges, and radicals |
| `succeeded` | MolGR obtained a normal coordinate-derived candidate |
| `suspicious_fallback` | A usable fallback candidate needs review |
| `failed` | No RDKit molecule could be built; `rdmol` is `None` |
| `None` | Topology access has not run yet, or no status applies |

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

Graph-level SDF, SMILES, and CML writers need a usable molecular graph. XYZ and coordinate-mode
Gaussian/ORCA inputs can still export coordinates without preserving full topology semantics.

## Next steps

- [Convert and export](conversion.md)
- [SDF/MOL format](../reference/formats/sdf.md)
- [Structure API](../reference/api/structure.md)
