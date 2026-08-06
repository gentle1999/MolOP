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
```

Set these values before the first access to `rdmol` for a frame. Lazy reconstruction reads the
current global values and records the effective settings on the frame. Restore the defaults or use a
separate process when different workflows require different policies.

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
