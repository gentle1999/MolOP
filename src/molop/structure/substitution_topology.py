"""Pure graph operations for substituent replacement."""

from __future__ import annotations

from rdkit import Chem

from molop.utils.types import RdMol

from .substitution_models import (
    MultiScaffoldFragment,
    ReplacementSite,
    ReplacementSiteSet,
    ScaffoldFragment,
)
from .substitution_sites import SUPPORTED_BOND_TYPES
from .utils import bond_type_mapping


def cut_scaffold(mol: RdMol, site: ReplacementSite) -> ScaffoldFragment:
    """Cut one replacement edge and retain every component except the substituent side."""
    working = Chem.RWMol(mol)
    bond = working.GetBondBetweenAtoms(site.scaffold_atom, site.substituent_atom)
    if bond is None:
        raise ValueError("replacement site bond is not present")
    if bond.IsInRing():
        raise ValueError("substituent replacement does not support ring bonds")
    if bond.GetBondType() not in SUPPORTED_BOND_TYPES:
        raise ValueError(f"unsupported replacement bond type: {bond.GetBondType()}")

    radical_delta = int(round(bond_type_mapping[bond.GetBondType()]))
    working.RemoveBond(site.scaffold_atom, site.substituent_atom)
    for atom_idx in (site.scaffold_atom, site.substituent_atom):
        atom = working.GetAtomWithIdx(atom_idx)
        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + radical_delta)

    fragment_indices = Chem.GetMolFrags(working, asMols=False, sanitizeFrags=False)
    fragments = Chem.GetMolFrags(working, asMols=True, sanitizeFrags=False)
    target_fragment_idx = next(
        (
            fragment_idx
            for fragment_idx, atom_indices in enumerate(fragment_indices)
            if site.substituent_atom in atom_indices
        ),
        None,
    )
    if target_fragment_idx is None:
        raise RuntimeError("could not identify the substituent fragment after cutting")

    retained = [
        (fragment, tuple(atom_indices))
        for fragment_idx, (fragment, atom_indices) in enumerate(
            zip(fragments, fragment_indices, strict=True)
        )
        if fragment_idx != target_fragment_idx
    ]
    if not retained:
        raise RuntimeError("the replacement site does not leave a scaffold")

    scaffold = retained[0][0]
    original_atom_indices = [int(atom_idx) for atom_idx in retained[0][1]]
    for fragment, atom_indices in retained[1:]:
        scaffold = Chem.CombineMols(scaffold, fragment)
        original_atom_indices.extend(int(atom_idx) for atom_idx in atom_indices)
    try:
        attachment_atom = original_atom_indices.index(site.scaffold_atom)
    except ValueError as exc:
        raise RuntimeError("the scaffold no longer contains its attachment atom") from exc
    return ScaffoldFragment(
        mol=Chem.Mol(scaffold),
        attachment_atom=attachment_atom,
        original_atom_indices=tuple(original_atom_indices),
    )


def combine_scaffold_and_replacement(
    scaffold: ScaffoldFragment,
    replacement: RdMol,
    bond_type: Chem.rdchem.BondType,
) -> tuple[Chem.RWMol, int]:
    """Join a prepared replacement (attachment atom zero) to a scaffold."""
    if replacement.GetNumAtoms() == 0:
        raise ValueError("replacement must contain at least one atom")
    combined = Chem.CombineMols(replacement, scaffold.mol)
    new_scaffold_atom = replacement.GetNumAtoms() + scaffold.attachment_atom
    result = Chem.RWMol(combined)
    if result.GetBondBetweenAtoms(0, new_scaffold_atom) is not None:
        raise ValueError("replacement attachment bond already exists")
    result.AddBond(0, new_scaffold_atom, bond_type)

    radical_delta = int(round(bond_type_mapping[bond_type]))
    scaffold_atom = result.GetAtomWithIdx(new_scaffold_atom)
    available = scaffold_atom.GetNumRadicalElectrons()
    if available < radical_delta:
        raise ValueError(
            f"scaffold attachment atom has {available} radical electrons; "
            f"{radical_delta} are required for a {bond_type} bond"
        )
    scaffold_atom.SetNumRadicalElectrons(available - radical_delta)
    # The replacement attachment radical was consumed while preparing its
    # open valence.  It is intentionally zero here, rather than subtracted a
    # second time during graph assembly.
    result.GetAtomWithIdx(0).SetNumRadicalElectrons(0)
    return result, new_scaffold_atom


def cut_multisite_scaffold(mol: RdMol, site: ReplacementSiteSet) -> MultiScaffoldFragment:
    """Cut all boundary edges of an experimental multi-anchor target set."""
    if len(site.anchors) < 2:
        raise ValueError("a multi-anchor site must contain at least two anchors")
    if len({anchor.scaffold_atom for anchor in site.anchors}) != len(site.anchors):
        raise ValueError("multi-anchor sites require distinct scaffold atoms")
    if len({anchor.substituent_atom for anchor in site.anchors}) != len(site.anchors):
        raise ValueError("multi-anchor sites require distinct substituent atoms")

    working = Chem.RWMol(mol)
    for anchor in site.anchors:
        bond = working.GetBondBetweenAtoms(anchor.scaffold_atom, anchor.substituent_atom)
        if bond is None:
            raise ValueError("a multi-anchor replacement edge is not present")
        if bond.GetBondType() not in SUPPORTED_BOND_TYPES:
            raise ValueError(f"unsupported replacement bond type: {bond.GetBondType()}")
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            raise ValueError("the experimental multi-anchor prototype supports single bonds only")

    for anchor in site.anchors:
        working.RemoveBond(anchor.scaffold_atom, anchor.substituent_atom)
        for atom_idx in (anchor.scaffold_atom, anchor.substituent_atom):
            atom = working.GetAtomWithIdx(atom_idx)
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + 1)

    fragment_indices = Chem.GetMolFrags(working, asMols=False, sanitizeFrags=False)
    fragments = Chem.GetMolFrags(working, asMols=True, sanitizeFrags=False)
    target_fragment_indices = {
        fragment_idx
        for fragment_idx, atom_indices in enumerate(fragment_indices)
        if any(anchor.substituent_atom in atom_indices for anchor in site.anchors)
    }
    if not target_fragment_indices:
        raise RuntimeError("could not identify the multi-anchor target fragments")
    target_indices = set().union(
        *(set(fragment_indices[fragment_idx]) for fragment_idx in target_fragment_indices)
    )
    if not all(anchor.substituent_atom in target_indices for anchor in site.anchors):
        raise ValueError("could not identify all multi-anchor target fragments")
    query_indices = set(site.query_match)
    if not query_indices.issubset(target_indices):
        raise ValueError("every multi-anchor query component must have at least one external edge")
    unexpected_heavy_atoms = [
        atom_idx
        for atom_idx in target_indices - query_indices
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 1
    ]
    if unexpected_heavy_atoms:
        raise ValueError(
            "the query must cover every heavy atom in the multi-anchor target fragments"
        )

    retained = [
        (fragment, tuple(atom_indices))
        for fragment_idx, (fragment, atom_indices) in enumerate(
            zip(fragments, fragment_indices, strict=True)
        )
        if fragment_idx not in target_fragment_indices
    ]
    if not retained:
        raise RuntimeError("the multi-anchor replacement does not leave a scaffold")
    scaffold = retained[0][0]
    original_atom_indices = list(retained[0][1])
    for fragment, atom_indices in retained[1:]:
        scaffold = Chem.CombineMols(scaffold, fragment)
        original_atom_indices.extend(atom_indices)
    attachment_atoms: list[int] = []
    for anchor in site.anchors:
        try:
            attachment_atoms.append(original_atom_indices.index(anchor.scaffold_atom))
        except ValueError as exc:
            raise RuntimeError("the scaffold lost a multi-anchor atom") from exc
    return MultiScaffoldFragment(
        mol=Chem.Mol(scaffold),
        attachment_atoms=tuple(attachment_atoms),
        original_atom_indices=tuple(original_atom_indices),
    )


def combine_scaffold_and_multisite_replacement(
    scaffold: MultiScaffoldFragment,
    replacement: RdMol,
    site: ReplacementSiteSet,
    replacement_indices: tuple[int, ...],
) -> Chem.RWMol:
    """Join an experimental replacement through all corresponding anchors."""
    if len(site.anchors) != len(replacement_indices):
        raise ValueError("replacement indices must match the number of anchors")
    combined = Chem.CombineMols(replacement, scaffold.mol)
    result = Chem.RWMol(combined)
    replacement_count = replacement.GetNumAtoms()
    for anchor, replacement_idx, scaffold_idx in zip(
        site.anchors, replacement_indices, scaffold.attachment_atoms, strict=True
    ):
        new_replacement_atom = replacement_idx
        new_scaffold_atom = replacement_count + scaffold_idx
        if result.GetBondBetweenAtoms(new_replacement_atom, new_scaffold_atom) is not None:
            raise ValueError("multi-anchor replacement bond already exists")
        result.AddBond(new_replacement_atom, new_scaffold_atom, anchor.bond_type)

        scaffold_atom = result.GetAtomWithIdx(new_scaffold_atom)
        available = scaffold_atom.GetNumRadicalElectrons()
        if available < 1:
            raise ValueError("scaffold anchor has no open valence for the replacement")
        scaffold_atom.SetNumRadicalElectrons(available - 1)

        replacement_atom = result.GetAtomWithIdx(new_replacement_atom)
        if replacement_atom.GetNumRadicalElectrons() > 0:
            replacement_atom.SetNumRadicalElectrons(replacement_atom.GetNumRadicalElectrons() - 1)
    return result
