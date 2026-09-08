"""Substituent-site discovery and input normalization."""

from __future__ import annotations

from collections.abc import Iterable

from rdkit import Chem

from molop.utils.types import RdMol

from .substitution_models import ReplacementAnchor, ReplacementSite, ReplacementSiteSet


GENERATED_ATOM_PROP = "_molop_generated_substituent"
SUPPORTED_BOND_TYPES = (
    Chem.rdchem.BondType.SINGLE,
    Chem.rdchem.BondType.DOUBLE,
    Chem.rdchem.BondType.TRIPLE,
    Chem.rdchem.BondType.DATIVE,
    Chem.rdchem.BondType.DATIVEL,
    Chem.rdchem.BondType.DATIVER,
    Chem.rdchem.BondType.DATIVEONE,
)


def normalize_query(query: str | RdMol, *, require_connected: bool = True) -> RdMol:
    """Return a query molecule with a predictable SMARTS interpretation."""
    if isinstance(query, str):
        query_mol = Chem.MolFromSmarts(query)
    elif isinstance(query, (Chem.Mol, Chem.RWMol)):
        query_mol = Chem.MolFromSmarts(Chem.MolToSmarts(query))
    else:
        query_mol = None
    if query_mol is None or query_mol.GetNumAtoms() == 0:
        raise ValueError("query must be a valid SMARTS string or RDKit molecule")
    if require_connected and len(Chem.GetMolFrags(query_mol)) != 1:
        raise ValueError("query must describe one connected substituent")
    return query_mol


def normalize_replacement(replacement: str | RdMol) -> RdMol:
    """Return a copy of a replacement molecule parsed as SMILES when needed."""
    if isinstance(replacement, str):
        replacement_mol = Chem.MolFromSmiles(replacement)
    elif isinstance(replacement, (Chem.Mol, Chem.RWMol)):
        replacement_mol = Chem.Mol(replacement)
    else:
        replacement_mol = None
    if replacement_mol is None or replacement_mol.GetNumAtoms() == 0:
        raise ValueError("replacement must be a valid SMILES string or RDKit molecule")
    return replacement_mol


def _validate_atom_idx(mol: RdMol, atom_idx: int, name: str) -> None:
    if not 0 <= atom_idx < mol.GetNumAtoms():
        raise ValueError(f"{name} is outside the molecule atom range")


def _site_from_bond(
    mol: RdMol,
    scaffold_atom: int,
    substituent_atom: int,
    *,
    query_match: Iterable[int] = (),
    allow_ring: bool = False,
) -> ReplacementSite:
    _validate_atom_idx(mol, scaffold_atom, "scaffold_atom")
    _validate_atom_idx(mol, substituent_atom, "substituent_atom")
    bond = mol.GetBondBetweenAtoms(scaffold_atom, substituent_atom)
    if bond is None:
        raise ValueError("the replacement site atoms must be bonded")
    if bond.IsInRing() and not allow_ring:
        raise ValueError("substituent replacement does not support ring bonds")
    if bond.GetBondType() not in SUPPORTED_BOND_TYPES:
        raise ValueError(f"unsupported replacement bond type: {bond.GetBondType()}")
    return ReplacementSite(
        scaffold_atom=scaffold_atom,
        substituent_atom=substituent_atom,
        bond_type=bond.GetBondType(),
        bond_stereo=bond.GetStereo(),
        query_match=tuple(query_match),
    )


def explicit_site(mol: RdMol, scaffold_atom: int, substituent_atom: int) -> ReplacementSite:
    """Resolve an explicitly supplied scaffold-to-substituent bond."""
    return _site_from_bond(mol, scaffold_atom, substituent_atom)


def _contains_generated_atom(mol: RdMol, atom_indices: Iterable[int]) -> bool:
    return any(
        mol.GetAtomWithIdx(atom_idx).HasProp(GENERATED_ATOM_PROP) for atom_idx in atom_indices
    )


def find_sites(
    mol: RdMol,
    query_mol: RdMol,
    *,
    bind_idx: int | None = None,
    exclude_generated: bool = False,
) -> list[ReplacementSite]:
    """Find matches that describe exactly one external, non-ring attachment edge.

    A match is accepted only when exactly one bond connects it to the rest of
    the molecule.  This is the defining restriction of substituent replacement:
    the query must represent a pendant fragment, not an arbitrary subgraph.
    """
    sites: dict[tuple[int, int], ReplacementSite] = {}
    for match in mol.GetSubstructMatches(query_mol, uniquify=True):
        match_set = set(match)
        if exclude_generated and _contains_generated_atom(mol, match):
            continue
        boundary: list[tuple[int, int]] = []
        for match_atom_idx in match:
            for neighbor in mol.GetAtomWithIdx(match_atom_idx).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Explicit hydrogens are valence bookkeeping, not an external
                # molecular attachment.  They are commonly present after a
                # coordinate-bearing molecule has been passed through AddHs.
                if neighbor_idx not in match_set and neighbor.GetAtomicNum() != 1:
                    boundary.append((match_atom_idx, neighbor_idx))
        if len(boundary) != 1:
            continue
        substituent_atom, scaffold_atom = boundary[0]
        if exclude_generated and mol.GetAtomWithIdx(scaffold_atom).HasProp(GENERATED_ATOM_PROP):
            continue
        site = _site_from_bond(
            mol,
            scaffold_atom=scaffold_atom,
            substituent_atom=substituent_atom,
            query_match=match,
        )
        if bind_idx is not None and site.scaffold_atom != bind_idx:
            continue
        sites[(site.scaffold_atom, site.substituent_atom)] = site
    return sorted(sites.values(), key=lambda site: (site.scaffold_atom, site.substituent_atom))


def find_multisite_sites(
    mol: RdMol,
    query_mol: RdMol,
    anchor_count: int,
    *,
    exclude_generated: bool = False,
) -> list[ReplacementSiteSet]:
    """Find query matches with an exact number of external edges.

    Unlike :func:`find_sites`, the query may be disconnected.  This allows a
    multi-anchor edit to select several independent pendant fragments in the
    parent molecule and replace them as one geometric assembly.
    """
    if anchor_count < 2:
        raise ValueError("anchor_count must be >= 2")

    sites: dict[tuple[tuple[int, ...], tuple[tuple[int, int], ...]], ReplacementSiteSet] = {}
    for match in mol.GetSubstructMatches(query_mol, uniquify=True):
        match_set = set(match)
        if exclude_generated and _contains_generated_atom(mol, match):
            continue
        boundary: set[tuple[int, int]] = set()
        for match_atom_idx in match:
            for neighbor in mol.GetAtomWithIdx(match_atom_idx).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Explicit hydrogens are not molecular attachment edges.
                if neighbor_idx not in match_set and neighbor.GetAtomicNum() != 1:
                    boundary.add((match_atom_idx, neighbor_idx))
        if len(boundary) != anchor_count:
            continue

        anchors: list[ReplacementAnchor] = []
        try:
            for substituent_atom, scaffold_atom in sorted(
                boundary, key=lambda item: (item[1], item[0])
            ):
                if exclude_generated and mol.GetAtomWithIdx(scaffold_atom).HasProp(
                    GENERATED_ATOM_PROP
                ):
                    raise ValueError("generated atom is part of the boundary")
                site = _site_from_bond(
                    mol,
                    scaffold_atom=scaffold_atom,
                    substituent_atom=substituent_atom,
                    query_match=match,
                    allow_ring=True,
                )
                anchors.append(
                    ReplacementAnchor(
                        scaffold_atom=site.scaffold_atom,
                        substituent_atom=site.substituent_atom,
                        bond_type=site.bond_type,
                        bond_stereo=site.bond_stereo,
                    )
                )
        except ValueError:
            continue
        site_set = ReplacementSiteSet(anchors=tuple(anchors), query_match=tuple(match))
        key = (
            tuple(sorted(match)),
            tuple((item.scaffold_atom, item.substituent_atom) for item in anchors),
        )
        sites[key] = site_set
    return sorted(
        sites.values(),
        key=lambda site: tuple(
            (anchor.scaffold_atom, anchor.substituent_atom) for anchor in site.anchors
        ),
    )


def mark_generated_atoms(mol: Chem.RWMol) -> None:
    """Mark replacement atoms so batch replacement cannot re-match its output."""
    for atom in mol.GetAtoms():
        atom.SetBoolProp(GENERATED_ATOM_PROP, True)


def clear_generated_atom_marks(mol: Chem.Mol) -> Chem.Mol:
    """Remove internal batch-processing markers from a returned molecule."""
    result = Chem.RWMol(mol)
    for atom in result.GetAtoms():
        if atom.HasProp(GENERATED_ATOM_PROP):
            atom.ClearProp(GENERATED_ATOM_PROP)
    return result.GetMol()
