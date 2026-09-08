"""Geometry-aware replacement of molecular fragments.

Replacement molecules use RDKit dummy atoms (``[*]``) to mark the atom-side
of each new attachment bond.  The dummy atoms are retained while the
replacement conformer is prepared, then removed immediately before graph
assembly.  ``replace_multisite_substituent`` is an experimental prototype
for replacing one or more parent-side fragments through two or more anchors
while freezing the retained parent coordinates.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdDistGeom import EmbedMolecule

from molop.config import moloplogger
from molop.utils.decorators import experimental
from molop.utils.functions import find_rigid_transform, invert_transform_coords, transform_coords
from molop.utils.types import RdMol

from . import GeometryTransformation
from .substitution_models import (
    ReplacementAnchor,
    ReplacementCandidate,
    ReplacementOptions,
    ReplacementSite,
    ReplacementSiteSet,
    StereoPolicy,
)
from .substitution_sites import (
    SUPPORTED_BOND_TYPES,
    clear_generated_atom_marks,
    explicit_site,
    find_multisite_sites,
    find_sites,
    mark_generated_atoms,
    normalize_query,
    normalize_replacement,
)
from .substitution_topology import (
    combine_scaffold_and_multisite_replacement,
    combine_scaffold_and_replacement,
    cut_multisite_scaffold,
    cut_scaffold,
)
from .topology import reset_atom_index
from .utils import bond_stereo_mapping, estimate_bond_length
from .validation import check_crowding, get_crowding_score


DEBUG_TAG = "[STRUCTURE]"
_ROTATABLE_BOND_TYPES = (
    Chem.rdchem.BondType.SINGLE,
    Chem.rdchem.BondType.DATIVE,
    Chem.rdchem.BondType.DATIVEL,
    Chem.rdchem.BondType.DATIVER,
    Chem.rdchem.BondType.DATIVEONE,
)
_DATIVE_BOND_TYPES = (
    Chem.rdchem.BondType.DATIVE,
    Chem.rdchem.BondType.DATIVEL,
    Chem.rdchem.BondType.DATIVER,
    Chem.rdchem.BondType.DATIVEONE,
)


class SubstituentReplacementError(RuntimeError):
    """Raised when no geometry candidate satisfies the replacement contract."""


def _copy_mol(mol: RdMol) -> Chem.Mol:
    return Chem.Mol(mol)


def _ensure_conformer(mol: Chem.Mol, random_seed: int) -> None:
    if mol.GetNumConformers() == 0:
        status = EmbedMolecule(mol, randomSeed=random_seed)
        if status < 0:
            raise SubstituentReplacementError("could not generate a replacement conformer")
    if mol.GetNumConformers() == 0:
        raise SubstituentReplacementError("replacement has no conformer")


def _geometry_bond_type(bond_type: Chem.rdchem.BondType) -> Chem.rdchem.BondType:
    return Chem.rdchem.BondType.DATIVE if bond_type in _DATIVE_BOND_TYPES else bond_type


def _resolve_dummy_indices(
    replacement_mol: RdMol,
    *,
    expected_count: int | None = None,
    requested_indices: Sequence[int] | None = None,
) -> tuple[int, ...]:
    """Resolve and validate all ``[*]`` atoms used as attachment markers."""
    dummy_indices = tuple(
        atom.GetIdx() for atom in replacement_mol.GetAtoms() if atom.GetAtomicNum() == 0
    )
    if requested_indices is None:
        selected = dummy_indices
    else:
        selected = tuple(int(index) for index in requested_indices)
        if len(set(selected)) != len(selected):
            raise ValueError("replacement dummy indices must be distinct")
        if set(selected) != set(dummy_indices):
            raise ValueError(
                "replacement dummy indices must identify every [*] atom and no other atom"
            )
    if expected_count is not None and len(selected) != expected_count:
        marker_name = "marker" if expected_count == 1 else "markers"
        raise ValueError(
            f"replacement must contain exactly {expected_count} [*] attachment {marker_name}"
        )
    if not selected:
        raise ValueError("replacement must contain at least one [*] attachment marker")
    for dummy_idx in selected:
        if dummy_idx < 0 or dummy_idx >= replacement_mol.GetNumAtoms():
            raise ValueError("replacement dummy index is outside the molecule")
        if replacement_mol.GetAtomWithIdx(dummy_idx).GetAtomicNum() != 0:
            raise ValueError("replacement dummy indices must identify [*] atoms")
    return selected


def _dummy_attachment_neighbors(
    replacement_mol: RdMol, dummy_indices: Sequence[int]
) -> tuple[int, ...]:
    """Return the real atom attached to each dummy marker."""
    attachment_indices: list[int] = []
    for dummy_idx in dummy_indices:
        dummy = replacement_mol.GetAtomWithIdx(dummy_idx)
        neighbors = list(dummy.GetNeighbors())
        if len(neighbors) != 1:
            raise ValueError("each replacement [*] marker must have exactly one molecular neighbor")
        neighbor = neighbors[0]
        if neighbor.GetAtomicNum() in (0, 1):
            raise ValueError("each replacement [*] marker must be attached to a non-hydrogen atom")
        if neighbor.GetIdx() in attachment_indices:
            raise ValueError("replacement [*] markers must attach to distinct atoms")
        attachment_indices.append(neighbor.GetIdx())
    return tuple(attachment_indices)


def _validate_dummy_bonds(
    replacement_mol: RdMol,
    dummy_indices: Sequence[int],
    bond_type: Chem.rdchem.BondType,
) -> None:
    """Ensure dummy bonds encode the requested parent attachment bond."""
    # RDKit SMILES does not provide a useful dative ``[*]`` bond notation for
    # this marker.  A single marker bond is therefore the input convention for
    # all dative parent bonds.
    expected_type = (
        Chem.rdchem.BondType.SINGLE
        if bond_type in _DATIVE_BOND_TYPES
        else _geometry_bond_type(bond_type)
    )
    for dummy_idx in dummy_indices:
        neighbor = next(iter(replacement_mol.GetAtomWithIdx(dummy_idx).GetNeighbors()))
        bond = replacement_mol.GetBondBetweenAtoms(dummy_idx, neighbor.GetIdx())
        if bond is None or bond.GetBondType() != expected_type:
            raise ValueError(
                "the bond from each replacement [*] marker must match the parent "
                f"attachment bond ({expected_type})"
            )


def _remove_dummy_markers(
    replacement: RdMol,
    dummy_indices: Sequence[int],
    attachment_indices: Sequence[int],
) -> tuple[Chem.Mol, tuple[int, ...]]:
    """Remove markers and return the reindexed attachment atoms."""
    if len(dummy_indices) != len(attachment_indices):
        raise ValueError("replacement dummy and attachment index counts must match")
    working = Chem.RWMol(replacement)
    for dummy_idx in sorted(dummy_indices, reverse=True):
        working.RemoveAtom(dummy_idx)
    adjusted_indices = tuple(
        attachment_idx - sum(dummy_idx < attachment_idx for dummy_idx in dummy_indices)
        for attachment_idx in attachment_indices
    )
    for attachment_idx in adjusted_indices:
        # The marker bond supplied the open valence.  The scaffold bond will
        # consume it during graph assembly, so no replacement radical remains.
        working.GetAtomWithIdx(attachment_idx).SetNumRadicalElectrons(0)
    return working.GetMol(), adjusted_indices


def _prepare_marked_replacement(
    replacement_mol: RdMol,
    bond_type: Chem.rdchem.BondType,
    *,
    random_seed: int,
    start_atomic_num: int,
    marker_indices: tuple[int, ...],
) -> Chem.Mol:
    """Prepare one dummy-marked replacement and move its attachment atom to zero."""
    if bond_type not in SUPPORTED_BOND_TYPES:
        raise ValueError(f"unsupported replacement bond type: {bond_type}")
    attachment_indices = _dummy_attachment_neighbors(replacement_mol, marker_indices)
    _validate_dummy_bonds(replacement_mol, marker_indices, bond_type)

    replacement = Chem.AddHs(_copy_mol(replacement_mol), addCoords=True)
    _ensure_conformer(replacement, random_seed)
    working = Chem.RWMol(replacement)
    marker_idx = marker_indices[0]
    attachment_idx = attachment_indices[0]
    GeometryTransformation.standard_orient(working, [marker_idx, attachment_idx])
    rdMolTransforms.SetBondLength(
        working.GetConformer(),
        marker_idx,
        attachment_idx,
        estimate_bond_length(
            working.GetAtomWithIdx(attachment_idx).GetAtomicNum(),
            start_atomic_num,
            bond_type,
        ),
    )
    prepared, (prepared_attachment_idx,) = _remove_dummy_markers(
        working.GetMol(), marker_indices, attachment_indices
    )
    return reset_atom_index(prepared, [prepared_attachment_idx])


def prepare_replacement(
    replacement_mol: RdMol,
    bond_type: Chem.rdchem.BondType,
    *,
    random_seed: int,
    start_atomic_num: int,
) -> Chem.Mol:
    """Normalize, embed, and orient a replacement with one ``[*]`` marker."""
    marker_indices = _resolve_dummy_indices(
        replacement_mol,
        expected_count=1,
    )
    return _prepare_marked_replacement(
        replacement_mol,
        bond_type,
        random_seed=random_seed,
        start_atomic_num=start_atomic_num,
        marker_indices=marker_indices,
    )


def _prepare_multisite_replacement(
    replacement_mol: RdMol,
    anchors: tuple[ReplacementAnchor, ...],
    *,
    random_seed: int,
    marker_indices: tuple[int, ...] | None = None,
) -> tuple[Chem.Mol, tuple[int, ...]]:
    """Prepare a multi-anchor replacement and return its attachment indices."""
    if any(anchor.bond_type != Chem.rdchem.BondType.SINGLE for anchor in anchors):
        raise ValueError("the experimental multi-anchor prototype supports single bonds only")
    marker_indices = _resolve_dummy_indices(
        replacement_mol,
        expected_count=len(anchors),
        requested_indices=marker_indices,
    )
    attachment_indices = _dummy_attachment_neighbors(replacement_mol, marker_indices)
    _validate_dummy_bonds(
        replacement_mol,
        marker_indices,
        Chem.rdchem.BondType.SINGLE,
    )
    # Multi-anchor placement only needs the heavy-atom conformer.  Removing
    # the marker before embedding also avoids RDKit distance-geometry issues
    # for replacement graphs with several dummy atoms.
    replacement, replacement_indices = _remove_dummy_markers(
        _copy_mol(replacement_mol), marker_indices, attachment_indices
    )
    Chem.SanitizeMol(replacement)
    replacement = Chem.AddHs(replacement, addCoords=True)
    _ensure_conformer(replacement, random_seed)
    working = Chem.RWMol(replacement)
    hydrogen_indices: list[int] = []
    for replacement_idx in replacement_indices:
        hydrogen = next(
            (
                neighbor.GetIdx()
                for neighbor in working.GetAtomWithIdx(replacement_idx).GetNeighbors()
                if neighbor.GetAtomicNum() == 1
            ),
            None,
        )
        if hydrogen is None:
            raise ValueError("each replacement [*] attachment atom must have an available valence")
        hydrogen_indices.append(hydrogen)
    for hydrogen_idx in sorted(set(hydrogen_indices), reverse=True):
        working.RemoveAtom(hydrogen_idx)
    Chem.SanitizeMol(working)
    return working.GetMol(), replacement_indices


def _multisite_target_positions(
    mol: RdMol,
    replacement: RdMol,
    site: ReplacementSiteSet,
    replacement_indices: tuple[int, ...],
) -> np.ndarray:
    """Estimate replacement attachment positions from the old bond directions."""
    source_positions = np.asarray(mol.GetConformer().GetPositions())
    targets: list[np.ndarray] = []
    for anchor, replacement_idx in zip(site.anchors, replacement_indices, strict=True):
        direction = (
            source_positions[anchor.substituent_atom] - source_positions[anchor.scaffold_atom]
        )
        direction_length = float(np.linalg.norm(direction))
        if direction_length <= 1e-8:
            raise SubstituentReplacementError("multi-anchor source bond has zero length")
        bond_length = estimate_bond_length(
            mol.GetAtomWithIdx(anchor.scaffold_atom).GetAtomicNum(),
            replacement.GetAtomWithIdx(replacement_idx).GetAtomicNum(),
            anchor.bond_type,
        )
        targets.append(
            source_positions[anchor.scaffold_atom] + direction / direction_length * bond_length
        )
    return np.asarray(targets)


def _rotate_positions_about_axis(
    positions: np.ndarray,
    origin: np.ndarray,
    axis: np.ndarray,
    angle: float,
) -> np.ndarray:
    """Rotate points around an axis using Rodrigues' formula."""
    axis_length = float(np.linalg.norm(axis))
    if axis_length <= 1e-8:
        raise SubstituentReplacementError("multi-anchor target points are coincident")
    unit_axis = axis / axis_length
    shifted = positions - origin
    cosine = np.cos(angle)
    sine = np.sin(angle)
    return (
        shifted * cosine
        + np.cross(unit_axis, shifted) * sine
        + np.outer(shifted @ unit_axis, unit_axis) * (1.0 - cosine)
        + origin
    )


def _set_positions(mol: Chem.Mol, positions: np.ndarray) -> Chem.Mol:
    conformer = mol.GetConformer()
    for atom_idx, position in enumerate(positions):
        conformer.SetAtomPosition(atom_idx, tuple(float(value) for value in position))
    return mol


def _multisite_replacement_poses(
    replacement: RdMol,
    replacement_indices: tuple[int, ...],
    targets: np.ndarray,
    *,
    angle_split: int,
    anchor_tolerance: float,
) -> list[Chem.Mol]:
    """Generate rigidly fitted poses for a multi-anchor replacement.

    Two anchors leave one rotational degree of freedom around their axis, so
    that mode is sampled.  Three or more anchors must define a non-collinear
    geometry and determine one proper rigid pose.
    """
    source_positions = np.asarray(replacement.GetConformer().GetPositions())
    source_anchor_positions = source_positions[list(replacement_indices)]
    if len(replacement_indices) >= 3:
        source_singular_values = np.linalg.svd(
            source_anchor_positions - source_anchor_positions.mean(axis=0),
            compute_uv=False,
        )
        target_singular_values = np.linalg.svd(
            targets - targets.mean(axis=0),
            compute_uv=False,
        )
        source_non_collinear = source_singular_values[1] > max(
            1e-6, source_singular_values[0] * 1e-6
        )
        target_non_collinear = target_singular_values[1] > max(
            1e-6, target_singular_values[0] * 1e-6
        )
        if not source_non_collinear or not target_non_collinear:
            raise SubstituentReplacementError(
                "three or more anchors must define non-collinear geometries"
            )
    transform = find_rigid_transform(source_anchor_positions, targets)
    fitted_positions = transform_coords(source_positions, transform)
    fit_error = float(
        np.sqrt(
            np.mean(np.sum((fitted_positions[list(replacement_indices)] - targets) ** 2, axis=1))
        )
    )
    if fit_error > anchor_tolerance:
        raise SubstituentReplacementError(
            f"replacement anchor fit RMSD {fit_error:.3f} exceeds {anchor_tolerance:.3f} Å"
        )

    poses: list[Chem.Mol] = []
    if len(replacement_indices) >= 3:
        candidate = Chem.Mol(replacement)
        return [_set_positions(candidate, fitted_positions)]

    axis = targets[1] - targets[0]
    for angle in np.linspace(0.0, 2.0 * np.pi, angle_split, endpoint=False):
        candidate = Chem.Mol(replacement)
        candidate_positions = _rotate_positions_about_axis(
            fitted_positions, targets[0], axis, float(angle)
        )
        poses.append(_set_positions(candidate, candidate_positions))
    return poses


def _validate_multisite_candidate(
    source: RdMol,
    candidate: Chem.Mol,
    scaffold_atom_count: int,
    scaffold_original_indices: tuple[int, ...],
    scaffold_attachment_atoms: tuple[int, ...],
    site: ReplacementSiteSet,
    replacement_indices: tuple[int, ...],
    *,
    anchor_tolerance: float,
) -> None:
    """Check frozen scaffold coordinates and the newly formed anchor bonds."""
    source_positions = np.asarray(source.GetConformer().GetPositions())
    candidate_positions = np.asarray(candidate.GetConformer().GetPositions())
    scaffold_positions = candidate_positions[:scaffold_atom_count]
    expected_positions = source_positions[list(scaffold_original_indices)]
    scaffold_drift = float(np.max(np.abs(scaffold_positions - expected_positions)))
    if scaffold_drift > 1e-6:
        raise SubstituentReplacementError(
            f"frozen scaffold coordinates drifted by {scaffold_drift:.3e} Å"
        )
    for anchor, replacement_idx, scaffold_idx in zip(
        site.anchors, replacement_indices, scaffold_attachment_atoms, strict=True
    ):
        result_replacement_idx = scaffold_atom_count + replacement_idx
        result_scaffold_idx = scaffold_idx
        bond_length = float(
            np.linalg.norm(
                candidate_positions[result_replacement_idx]
                - candidate_positions[result_scaffold_idx]
            )
        )
        expected_length = estimate_bond_length(
            source.GetAtomWithIdx(anchor.scaffold_atom).GetAtomicNum(),
            candidate.GetAtomWithIdx(result_replacement_idx).GetAtomicNum(),
            anchor.bond_type,
        )
        if abs(bond_length - expected_length) > anchor_tolerance:
            raise SubstituentReplacementError(
                f"multi-anchor bond length error {abs(bond_length - expected_length):.3f} "
                f"exceeds {anchor_tolerance:.3f} Å"
            )


def _assemble_multisite_candidate(
    mol: RdMol,
    site: ReplacementSiteSet,
    replacement_mol: RdMol,
    replacement_marker_indices: tuple[int, ...],
    options: ReplacementOptions,
    seed: int,
    *,
    anchor_tolerance: float,
) -> Chem.Mol:
    """Experimental multi-anchor replacement with a frozen scaffold."""
    if len(site.anchors) < 2:
        raise ValueError("the experimental prototype requires at least two anchors")
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        raise ValueError("the experimental multi-anchor replacement requires a 3D parent conformer")

    scaffold = cut_multisite_scaffold(mol, site)
    replacement, replacement_indices = _prepare_multisite_replacement(
        replacement_mol,
        site.anchors,
        random_seed=seed,
        marker_indices=replacement_marker_indices,
    )
    targets = _multisite_target_positions(mol, replacement, site, replacement_indices)
    poses = _multisite_replacement_poses(
        replacement,
        replacement_indices,
        targets,
        angle_split=options.angle_split,
        anchor_tolerance=anchor_tolerance,
    )

    candidates: list[ReplacementCandidate] = []
    for pose in poses:
        combined = combine_scaffold_and_multisite_replacement(
            scaffold,
            pose,
            site,
            replacement_indices,
        )
        final = _finalize_atom_order(combined, pose.GetNumAtoms())
        Chem.SanitizeMol(final)
        _validate_multisite_candidate(
            mol,
            final,
            scaffold.mol.GetNumAtoms(),
            scaffold.original_atom_indices,
            scaffold.attachment_atoms,
            site,
            replacement_indices,
            anchor_tolerance=anchor_tolerance,
        )
        if not check_crowding(final, options.crowding_threshold):
            continue
        candidates.append(
            ReplacementCandidate(
                mol=final,
                seed=seed,
                crowding_score=_safe_crowding_score(final),
            )
        )
    if not candidates:
        raise SubstituentReplacementError("no valid multi-anchor geometry candidate was found")
    return max(candidates, key=lambda item: item.crowding_score).mol


def _orient_host(
    mol: RdMol, site: ReplacementSite, random_seed: int
) -> tuple[Chem.Mol, np.ndarray]:
    """Orient the host temporarily and return the transform needed to undo it."""
    oriented = Chem.RWMol(mol)
    _ensure_conformer(oriented, random_seed)
    source_positions = np.array(oriented.GetConformer().GetPositions(), copy=True)
    if site.bond_type == Chem.rdchem.BondType.DOUBLE:
        neighbor = next(
            (
                atom.GetIdx()
                for atom in oriented.GetAtomWithIdx(site.scaffold_atom).GetNeighbors()
                if atom.GetIdx() != site.substituent_atom
            ),
            None,
        )
        indices = (
            [site.scaffold_atom, site.substituent_atom, neighbor]
            if neighbor is not None
            else [site.scaffold_atom, site.substituent_atom]
        )
    else:
        indices = [site.scaffold_atom, site.substituent_atom]
    GeometryTransformation.standard_orient(oriented, indices)
    oriented_mol = oriented.GetMol()
    transform = find_rigid_transform(
        source_positions, np.asarray(oriented_mol.GetConformer().GetPositions())
    )
    return oriented_mol, transform


def _restore_host_orientation(mol: Chem.Mol, transform: np.ndarray) -> Chem.Mol:
    """Restore the source host frame after geometry assembly in a temporary frame."""
    for conformer in mol.GetConformers():
        positions = invert_transform_coords(np.asarray(conformer.GetPositions()), transform)
        for atom_idx, position in enumerate(positions):
            conformer.SetAtomPosition(atom_idx, tuple(float(value) for value in position))
    return mol


def _orient_double_replacement(replacement: RdMol) -> list[Chem.Mol]:
    """Generate orientations that expose each substituent at the double bond."""
    anchors = [atom.GetIdx() for atom in replacement.GetAtomWithIdx(0).GetNeighbors()]
    if not anchors:
        anchors = [None]
    oriented: list[Chem.Mol] = []
    for anchor in anchors:
        candidate = Chem.RWMol(replacement)
        dummy_idx = candidate.AddAtom(Chem.Atom("H"))
        indices = [dummy_idx, 0] if anchor is None else [dummy_idx, 0, anchor]
        GeometryTransformation.standard_orient(candidate, indices)
        candidate.RemoveAtom(dummy_idx)
        oriented.append(candidate.GetMol())
    return oriented


def _stereo_from_geometry(
    mol: RdMol, begin: int, end: int
) -> tuple[Chem.rdchem.BondStereo, tuple[int, ...]]:
    """Infer the connecting double-bond stereo and its defining atom indices."""
    checked = _copy_mol(mol)
    Chem.SanitizeMol(checked)
    Chem.AssignStereochemistryFrom3D(checked, replaceExistingTags=True)
    Chem.DetectBondStereochemistry(checked)
    bond = checked.GetBondBetweenAtoms(begin, end)
    if bond is None:
        raise SubstituentReplacementError("connecting bond disappeared during stereo detection")
    return bond.GetStereo(), tuple(bond.GetStereoAtoms())


def _target_stereo(site: ReplacementSite, policy: str) -> Chem.rdchem.BondStereo | None:
    if policy == "E":
        return bond_stereo_mapping["E"]
    if policy == "Z":
        return bond_stereo_mapping["Z"]
    if policy == "preserve" and site.bond_stereo in (
        bond_stereo_mapping["E"],
        bond_stereo_mapping["Z"],
    ):
        return site.bond_stereo
    return None


def set_best_dihedral(
    rmol: Chem.RWMol,
    threshold: float,
    angle_split: int,
    start: int,
    end: int,
) -> float | None:
    """Choose the lowest-UFF-energy non-crowded torsion around one bond."""
    bond = rmol.GetBondBetweenAtoms(start, end)
    if bond is None:
        raise ValueError("the rotatable bond is not found")
    if bond.GetBondType() not in _ROTATABLE_BOND_TYPES:
        raise ValueError("the selected bond type should be rotatable")
    first_atom_idx = next(
        (
            atom.GetIdx()
            for atom in rmol.GetAtomWithIdx(start).GetNeighbors()
            if atom.GetIdx() != end
        ),
        None,
    )
    fourth_atom_idx = next(
        (
            atom.GetIdx()
            for atom in rmol.GetAtomWithIdx(end).GetNeighbors()
            if atom.GetIdx() != start
        ),
        None,
    )
    if first_atom_idx is None or fourth_atom_idx is None:
        return None

    Chem.SanitizeMol(rmol)
    best_angle: float | None = None
    best_score = float("-inf")
    for angle in np.linspace(0.0, 2.0 * np.pi, angle_split, endpoint=False):
        rdMolTransforms.SetDihedralRad(
            rmol.GetConformer(),
            first_atom_idx,
            start,
            end,
            fourth_atom_idx,
            float(angle),
        )
        if not check_crowding(rmol, threshold):
            continue
        try:
            score = get_crowding_score(rmol)
        except Exception:
            score = 0.0
        if best_angle is None or score > best_score:
            best_angle = float(angle)
            best_score = score
    if best_angle is None:
        raise SubstituentReplacementError("no non-crowded dihedral angle was found")
    rdMolTransforms.SetDihedralRad(
        rmol.GetConformer(),
        first_atom_idx,
        start,
        end,
        fourth_atom_idx,
        best_angle,
    )
    return best_angle


def _assemble_single_candidate(
    mol: RdMol,
    site: ReplacementSite,
    replacement_mol: RdMol,
    options: ReplacementOptions,
    seed: int,
    *,
    mark_replacement: bool,
) -> Chem.Mol:
    oriented_host, host_transform = _orient_host(mol, site, seed)
    scaffold = cut_scaffold(oriented_host, site)
    replacement = prepare_replacement(
        replacement_mol,
        site.bond_type,
        random_seed=seed,
        start_atomic_num=mol.GetAtomWithIdx(site.scaffold_atom).GetAtomicNum(),
    )
    if mark_replacement:
        marked = Chem.RWMol(replacement)
        mark_generated_atoms(marked)
        replacement = marked.GetMol()
    combined, scaffold_atom = combine_scaffold_and_replacement(
        scaffold, replacement, site.bond_type
    )

    if site.bond_type in _ROTATABLE_BOND_TYPES:
        set_best_dihedral(
            combined,
            threshold=options.crowding_threshold,
            angle_split=options.angle_split,
            start=scaffold_atom,
            end=0,
        )
    elif site.bond_type == Chem.rdchem.BondType.TRIPLE:
        for anchor in combined.GetAtomWithIdx(0).GetNeighbors():
            if anchor.GetIdx() != scaffold_atom:
                set_best_dihedral(
                    combined,
                    threshold=options.crowding_threshold,
                    angle_split=options.angle_split,
                    start=0,
                    end=anchor.GetIdx(),
                )
    final = _finalize_atom_order(combined, replacement.GetNumAtoms())
    return _restore_host_orientation(final, host_transform)


def _assemble_double_candidates(
    mol: RdMol,
    site: ReplacementSite,
    replacement_mol: RdMol,
    options: ReplacementOptions,
    seed: int,
    *,
    mark_replacement: bool,
) -> list[ReplacementCandidate]:
    oriented_host, host_transform = _orient_host(mol, site, seed)
    scaffold = cut_scaffold(oriented_host, site)
    replacement = prepare_replacement(
        replacement_mol,
        site.bond_type,
        random_seed=seed,
        start_atomic_num=mol.GetAtomWithIdx(site.scaffold_atom).GetAtomicNum(),
    )
    target_stereo = _target_stereo(site, options.stereo_policy)
    candidates: list[ReplacementCandidate] = []
    for oriented_replacement in _orient_double_replacement(replacement):
        if mark_replacement:
            marked = Chem.RWMol(oriented_replacement)
            mark_generated_atoms(marked)
            oriented_replacement = marked.GetMol()
        combined, scaffold_atom = combine_scaffold_and_replacement(
            scaffold, oriented_replacement, site.bond_type
        )
        Chem.SanitizeMol(combined)
        stereo, stereo_atoms = _stereo_from_geometry(combined, scaffold_atom, 0)
        if target_stereo is not None and stereo != target_stereo:
            continue
        connecting_bond = combined.GetBondBetweenAtoms(scaffold_atom, 0)
        if connecting_bond is None:
            raise SubstituentReplacementError("connecting double bond disappeared")
        connecting_bond.SetStereo(stereo)
        if len(stereo_atoms) == 2:
            connecting_bond.SetStereoAtoms(stereo_atoms[0], stereo_atoms[1])
        for anchor in combined.GetAtomWithIdx(0).GetNeighbors():
            if anchor.GetIdx() != scaffold_atom:
                set_best_dihedral(
                    combined,
                    threshold=options.crowding_threshold,
                    angle_split=options.angle_split,
                    start=0,
                    end=anchor.GetIdx(),
                )
        final = _finalize_atom_order(combined, oriented_replacement.GetNumAtoms())
        final = _restore_host_orientation(final, host_transform)
        candidates.append(
            ReplacementCandidate(
                mol=final,
                seed=seed,
                stereo=stereo,
                crowding_score=_safe_crowding_score(final),
            )
        )
    if not candidates:
        wanted = options.stereo_policy
        raise SubstituentReplacementError(f"could not generate a double-bond {wanted} candidate")
    return candidates


def _safe_crowding_score(mol: RdMol) -> float:
    try:
        return get_crowding_score(mol)
    except Exception:
        return 0.0


def _finalize_atom_order(mol: Chem.RWMol | RdMol, replacement_count: int) -> Chem.Mol:
    """Return scaffold atoms first while preserving the generated coordinates."""
    atom_count = mol.GetNumAtoms()
    mapping = list(range(replacement_count, atom_count)) + list(range(replacement_count))
    return reset_atom_index(mol, mapping)


def replace_substituent_once(
    mol: RdMol,
    replacement_mol: RdMol,
    site: ReplacementSite,
    options: ReplacementOptions,
    *,
    seed: int,
    mark_replacement: bool = False,
) -> Chem.Mol:
    """Apply one already-resolved replacement site."""
    if site.bond_type == Chem.rdchem.BondType.DOUBLE:
        candidates = _assemble_double_candidates(
            mol,
            site,
            replacement_mol,
            options,
            seed,
            mark_replacement=mark_replacement,
        )
        candidate = max(candidates, key=lambda item: item.crowding_score)
        result = candidate.mol
    else:
        result = _assemble_single_candidate(
            mol,
            site,
            replacement_mol,
            options,
            seed,
            mark_replacement=mark_replacement,
        )
    Chem.SanitizeMol(result)
    if not check_crowding(result, options.crowding_threshold):
        raise SubstituentReplacementError("assembled substituent is too crowded")
    return result


def _seed_order(options: ReplacementOptions) -> list[int]:
    seeds = list(range(1, options.attempt_num + 1))
    np.random.RandomState(options.random_seed).shuffle(seeds)
    return seeds


def _replace_all_once(
    mol: RdMol,
    query_mol: RdMol,
    replacement_mol: RdMol,
    options: ReplacementOptions,
    *,
    seed: int,
) -> Chem.Mol:
    current = _copy_mol(mol)
    initial_sites = find_sites(current, query_mol)
    if not initial_sites:
        raise ValueError("query does not identify a pendant substituent")
    maximum_replacements = len(initial_sites)
    replaced = 0
    while True:
        sites = find_sites(current, query_mol, exclude_generated=True)
        if not sites:
            break
        if replaced >= maximum_replacements:
            raise SubstituentReplacementError(
                "replacement generated more eligible sites than the original molecule"
            )
        current = replace_substituent_once(
            current,
            replacement_mol,
            sites[0],
            options,
            seed=seed,
            mark_replacement=True,
        )
        replaced += 1
    result = clear_generated_atom_marks(current)
    Chem.SanitizeMol(result)
    if not check_crowding(result, options.crowding_threshold):
        raise SubstituentReplacementError("assembled substituents are too crowded")
    return result


@experimental
def replace_multisite_substituent(
    mol: RdMol,
    query: str | RdMol,
    replacement: str | RdMol,
    replacement_absolute_indices: Sequence[int] | None = None,
    *,
    attempt_num: int = 10,
    crowding_threshold: float = 0.75,
    angle_split: int = 12,
    randomSeed: int = 114514,
    anchor_tolerance: float = 0.25,
) -> RdMol:
    """Experimentally replace one or more parent fragments through two or more anchors.

    The replacement must contain two or more ``[*]`` markers, one for each
    new attachment bond. Marker order follows input atom order by default;
    ``replacement_absolute_indices`` may explicitly provide that order using
    the input dummy-atom indices. The query may be connected or disconnected,
    and must identify exactly as many external single-bond attachments as
    markers. This permits several independent parent-side fragments to be
    replaced by one coordinate-bearing replacement assembly.
    Two anchors sample axial rotation; three or more anchors require
    non-collinear anchor geometries and use one rigid fit. The parent scaffold
    coordinates are frozen.
    """
    if anchor_tolerance < 0:
        raise ValueError("anchor_tolerance must be >= 0")
    replacement_mol = normalize_replacement(replacement)
    replacement_marker_indices = _resolve_dummy_indices(
        replacement_mol,
        requested_indices=(
            tuple(int(index) for index in replacement_absolute_indices)
            if replacement_absolute_indices is not None
            else None
        ),
    )
    if len(replacement_marker_indices) < 2:
        raise ValueError("the experimental prototype requires at least two [*] markers")
    query_mol = normalize_query(query, require_connected=False)
    sites = find_multisite_sites(mol, query_mol, len(replacement_marker_indices))
    if not sites:
        raise ValueError("query does not identify a unique multi-anchor replacement site")
    if len(sites) > 1:
        raise ValueError(
            "query identifies multiple multi-anchor substituents; use a more specific query"
        )
    options = ReplacementOptions(
        attempt_num=attempt_num,
        crowding_threshold=crowding_threshold,
        angle_split=angle_split,
        random_seed=randomSeed,
    )

    failures: list[str] = []
    for seed in _seed_order(options):
        try:
            result = _assemble_multisite_candidate(
                mol,
                sites[0],
                replacement_mol,
                replacement_marker_indices,
                options,
                seed,
                anchor_tolerance=anchor_tolerance,
            )
        except (SubstituentReplacementError, ValueError, RuntimeError) as exc:
            failures.append(f"seed {seed}: {exc}")
            continue
        moloplogger.debug(
            f"{DEBUG_TAG} | Experimental multi-anchor replacement succeeded at seed {seed}."
        )
        return result
    detail = "; ".join(failures[-3:])
    raise SubstituentReplacementError(
        f"no valid multi-anchor replacement was generated after {options.attempt_num} attempts"
        + (f": {detail}" if detail else "")
    )


def replace_substituent(
    mol: RdMol,
    query: str | RdMol | None,
    replacement: str | RdMol,
    bind_idx: int | None = None,
    replace_all: bool = False,
    attempt_num: int = 10,
    crowding_threshold: float = 0.75,
    angle_split: int = 10,
    randomSeed: int = 114514,
    start_idx: int | None = None,
    end_idx: int | None = None,
    *,
    stereo_policy: StereoPolicy = "preserve",
) -> RdMol:
    """Replace one or more pendant substituents with a dummy-marked replacement.

    Query-based replacement requires a connected query match with exactly one
    external non-ring bond. ``bind_idx`` is the atom index on the retained
    scaffold. Alternatively, ``start_idx`` and ``end_idx`` explicitly name
    the scaffold and substituent sides of one bond. The replacement SMILES
    must contain exactly one ``[*]`` marker; the atom next to that marker is
    connected to the retained scaffold.
    """
    if (start_idx is None) != (end_idx is None):
        raise ValueError("start_idx and end_idx must be provided together")
    if replace_all and (start_idx is not None or bind_idx is not None):
        raise ValueError("replace_all cannot be combined with explicit site selection")
    replacement_mol = normalize_replacement(replacement)
    query_mol: RdMol | None = None
    options = ReplacementOptions(
        attempt_num=attempt_num,
        crowding_threshold=crowding_threshold,
        angle_split=angle_split,
        random_seed=randomSeed,
        stereo_policy=stereo_policy,
    )

    if start_idx is not None and end_idx is not None:
        sites = [explicit_site(mol, start_idx, end_idx)]
    else:
        if query is None:
            raise ValueError("query is required when the replacement site is not explicit")
        query_mol = normalize_query(query)
        sites = find_sites(mol, query_mol, bind_idx=bind_idx)
        if not sites:
            raise ValueError("query does not identify a unique pendant substituent")
        if not replace_all and len(sites) > 1:
            raise ValueError(
                "query identifies multiple pendant substituents; use bind_idx or replace_all"
            )

    failures: list[str] = []
    for seed in _seed_order(options):
        try:
            if replace_all:
                if query_mol is None:
                    raise RuntimeError("replace_all requires a query-based replacement")
                result = _replace_all_once(
                    mol,
                    query_mol,
                    replacement_mol,
                    options,
                    seed=seed,
                )
            else:
                result = replace_substituent_once(
                    mol,
                    replacement_mol,
                    sites[0],
                    options,
                    seed=seed,
                )
        except (SubstituentReplacementError, ValueError, RuntimeError) as exc:
            failures.append(f"seed {seed}: {exc}")
            continue
        moloplogger.debug(f"{DEBUG_TAG} | Substituent replacement succeeded at seed {seed}.")
        return result
    detail = "; ".join(failures[-3:])
    raise SubstituentReplacementError(
        f"no valid substituent replacement was generated after {options.attempt_num} attempts"
        + (f": {detail}" if detail else "")
    )


__all__ = [
    "SubstituentReplacementError",
    "replace_multisite_substituent",
    "replace_substituent",
]
