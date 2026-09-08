"""Value objects used by substituent replacement."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from rdkit import Chem

from molop.utils.types import RdMol


StereoPolicy = Literal["preserve", "E", "Z", "any"]


@dataclass(frozen=True)
class ReplacementAnchor:
    """One scaffold-to-fragment boundary edge in an experimental multisite edit."""

    scaffold_atom: int
    substituent_atom: int
    bond_type: Chem.rdchem.BondType
    bond_stereo: Chem.rdchem.BondStereo


@dataclass(frozen=True)
class ReplacementSiteSet:
    """All external edges selected by one experimental multi-anchor match.

    ``query_match`` may contain several disconnected query components.
    """

    anchors: tuple[ReplacementAnchor, ...]
    query_match: tuple[int, ...] = ()


@dataclass(frozen=True)
class ReplacementSite:
    """A single, acyclic bond through which a substituent is attached.

    ``scaffold_atom`` belongs to the part that is kept and ``substituent_atom``
    belongs to the part that is removed.  Atom indices refer to the molecule
    from which the site was discovered.
    """

    scaffold_atom: int
    substituent_atom: int
    bond_type: Chem.rdchem.BondType
    bond_stereo: Chem.rdchem.BondStereo
    query_match: tuple[int, ...] = ()


@dataclass(frozen=True)
class ReplacementOptions:
    """Validated options shared by one or more replacement attempts."""

    attempt_num: int = 10
    crowding_threshold: float = 0.75
    angle_split: int = 10
    random_seed: int = 114514
    stereo_policy: StereoPolicy = "preserve"

    def __post_init__(self) -> None:
        if self.attempt_num < 1:
            raise ValueError("attempt_num must be >= 1")
        if self.angle_split < 1:
            raise ValueError("angle_split must be >= 1")
        if self.crowding_threshold < 0:
            raise ValueError("crowding_threshold must be >= 0")
        if self.stereo_policy not in ("preserve", "E", "Z", "any"):
            raise ValueError("stereo_policy must be one of 'preserve', 'E', 'Z', or 'any'")


@dataclass(frozen=True)
class ScaffoldFragment:
    """The retained molecular graph after a replacement site is cut."""

    mol: RdMol
    attachment_atom: int
    original_atom_indices: tuple[int, ...]


@dataclass(frozen=True)
class MultiScaffoldFragment:
    """The retained scaffold after cutting all edges of a multi-anchor site.

    The retained graph may initially contain several connected components;
    the replacement assembly can bridge them through its anchor bonds.
    """

    mol: RdMol
    attachment_atoms: tuple[int, ...]
    original_atom_indices: tuple[int, ...]


@dataclass(frozen=True)
class ReplacementCandidate:
    """A fully assembled and validated candidate geometry."""

    mol: RdMol
    seed: int
    stereo: Chem.rdchem.BondStereo = Chem.rdchem.BondStereo.STEREONONE
    crowding_score: float = 0.0
