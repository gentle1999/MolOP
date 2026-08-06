"""Generate the checked-in SVGs used by the bilingual documentation."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit_dof import MolsToGridDofImage

from molop import AutoParser


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "docs" / "assets" / "examples"
METAL_COMPLEX_LOG = ASSETS / "mn_complex_sp.log"
TS_LOG = (
    ROOT
    / "tests"
    / "test_files"
    / "g16log"
    / ("000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log")
)


def write_grid(
    path: Path,
    molecules: list[Chem.Mol],
    legends: list[str],
    *,
    highlight_atoms: list[list[int]] | None = None,
    highlight_bonds: list[list[int]] | None = None,
    mols_per_row: int = 3,
    sub_image_size: tuple[int, int] = (280, 240),
) -> None:
    """Render RDKit molecules as SVG and write the generated result."""

    svg = MolsToGridDofImage(
        molecules,
        molsPerRow=mols_per_row,
        subImgSize=sub_image_size,
        legends=legends,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
        use_svg=True,
        return_image=False,
    )
    if not isinstance(svg, str) or "<svg" not in svg:
        raise RuntimeError(f"RDKit did not return SVG output for {path}")
    path.write_text(svg, encoding="utf-8")


def bond_index_by_pair(mol: Chem.Mol, pair: tuple[int, int]) -> int | None:
    start, end = pair
    bond = mol.GetBondBetweenAtoms(start, end)
    return None if bond is None else bond.GetIdx()


def generate_reconstruction_image() -> None:
    frame = AutoParser(METAL_COMPLEX_LOG, n_jobs=1)[0][-1]
    molecule = frame.rdmol
    if molecule is None:
        raise RuntimeError("The documentation metal complex did not reconstruct a graph")

    raw_xyz = Chem.MolFromXYZBlock(frame.to_XYZ())
    if raw_xyz is None or raw_xyz.GetNumBonds() != 0:
        raise RuntimeError("The raw XYZ baseline is unavailable or already contains bonds")

    distance_connectivity = Chem.Mol(raw_xyz)
    rdDetermineBonds.DetermineConnectivity(distance_connectivity)

    mn_indices = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() == "Mn"]
    if len(mn_indices) != 1:
        raise RuntimeError("The documentation metal complex must contain exactly one Mn atom")
    mn_index = mn_indices[0]
    mn_bonds = [
        bond
        for bond in molecule.GetBonds()
        if mn_index in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    ]
    if not mn_bonds or any(bond.GetBondType() != Chem.BondType.DATIVE for bond in mn_bonds):
        raise RuntimeError("MolGR did not assign the expected Mn dative bonds")

    write_grid(
        ASSETS / "mn_complex_graph_reconstruction.svg",
        [raw_xyz, distance_connectivity, molecule],
        [
            f"raw XYZ | {raw_xyz.GetNumAtoms()} atoms, 0 bonds",
            f"RDKit connectivity | {distance_connectivity.GetNumBonds()} single bonds",
            f"MolGR | {molecule.GetNumBonds()} bonds, {len(mn_bonds)} Mn-C dative",
        ],
        highlight_atoms=[[mn_index], [mn_index], [mn_index]],
        sub_image_size=(380, 320),
    )


def generate_transition_state_images() -> None:
    parsed = AutoParser(TS_LOG, n_jobs=1)
    frame = next((candidate for candidate in parsed[0] if candidate.is_TS), None)
    if frame is None:
        raise RuntimeError("The documentation TS fixture has no TS frame")

    vibration_molecules = [candidate.rdmol for candidate in frame.ts_vibration()]
    vibration_molecules = [molecule for molecule in vibration_molecules if molecule is not None]
    if len(vibration_molecules) < 2:
        raise RuntimeError("The documentation TS fixture produced too few vibration candidates")
    write_grid(
        ASSETS / "ts_imaginary_mode.svg",
        vibration_molecules,
        [f"imaginary-mode candidate {index}" for index in range(1, len(vibration_molecules) + 1)],
    )

    reactant, product = frame.possible_pre_post_ts(show_3D=True)
    difference = frame.to_diff_rdmol()
    if difference is None:
        raise RuntimeError("The documentation TS fixture did not produce a difference graph")

    changed_pairs: list[tuple[int, int]] = []
    display_difference = Chem.RWMol(difference)
    for bond in display_difference.GetBonds():
        if bond.GetBondType() != Chem.BondType.ZERO:
            continue
        changed_pairs.append(
            (
                min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            )
        )
        # ZERO is MolOP's machine-readable difference marker. Promote it only
        # in this display copy so RDKit can draw the virtual bond in red.
        bond.SetBondType(Chem.BondType.SINGLE)
    display_difference = display_difference.GetMol()

    reactant_highlights = [
        bond_index
        for pair in changed_pairs
        if (bond_index := bond_index_by_pair(reactant, pair)) is not None
    ]
    product_highlights = [
        bond_index
        for pair in changed_pairs
        if (bond_index := bond_index_by_pair(product, pair)) is not None
    ]
    difference_highlights = [
        bond.GetIdx()
        for bond in display_difference.GetBonds()
        if tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) in changed_pairs
    ]
    write_grid(
        ASSETS / "ts_endpoints_difference.svg",
        [reactant, display_difference, product],
        ["reactant candidate", "virtual-bond difference graph", "product candidate"],
        highlight_bonds=[reactant_highlights, difference_highlights, product_highlights],
    )


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    generate_reconstruction_image()
    generate_transition_state_images()


if __name__ == "__main__":
    main()
