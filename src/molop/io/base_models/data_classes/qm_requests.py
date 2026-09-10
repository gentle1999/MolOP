"""Structured quantum-chemistry request and input-semantic data classes."""

from __future__ import annotations

from typing import Any, ClassVar

from pint._typing import UnitLike
from pint.facets.plain import PlainQuantity
from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.unit import atom_ureg


class QMResourceRequest(BaseDataClassWithUnit):
    default_units: ClassVar[dict[str, UnitLike]] = {"memory": atom_ureg.megabyte}
    set_default_units: ClassVar[bool] = True

    num_cpu: int | None = Field(default=None, description="Requested number of CPU cores")
    memory: PlainQuantity | None = Field(
        default=None, description="Requested memory, normalized to megabyte"
    )
    raw: str = Field(default="", description="Raw resource directives")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific resource options"
    )


class QMBasisSet(BaseDataClassWithUnit):
    name: str = Field(default="", description="Basis-set name or expression")
    role: str = Field(
        default="orbital", description="Basis role, such as orbital, auxiliary, or ecp"
    )
    scope: str | None = Field(
        default=None, description="Scope for this basis assignment, such as global or atom"
    )
    atom_indices: list[int] = Field(
        default_factory=list, description="0-based atom indices covered by this basis"
    )
    element_symbols: list[str] = Field(
        default_factory=list, description="Element symbols covered by this basis"
    )
    raw: str = Field(default="", description="Raw basis-set directive")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific basis-set options"
    )


class QMModelChemistry(BaseDataClassWithUnit):
    method_family: str | None = Field(default=None, description="Canonical method family")
    method: str | None = Field(default=None, description="Concrete method name")
    reference_method: str | None = Field(
        default=None, description="Reference wavefunction or method"
    )
    functional: str | None = Field(default=None, description="DFT or double-hybrid functional")
    basis_set: str | None = Field(default=None, description="Primary orbital basis set")
    auxiliary_basis_set: str | None = Field(
        default=None, description="Primary auxiliary or fitting basis set"
    )
    basis_sets: list[QMBasisSet] = Field(
        default_factory=list, description="Structured global or local basis-set assignments"
    )
    dispersion_correction: str | None = Field(
        default=None, description="Empirical or nonlocal dispersion correction"
    )
    solvation_model: str | None = Field(default=None, description="Solvation model")
    solvent: str | None = Field(default=None, description="Solvent name")
    relativistic: str | None = Field(default=None, description="Relativistic treatment")
    spin_treatment: str | None = Field(default=None, description="Restricted/open-shell treatment")
    raw_keywords: str = Field(default="", description="Raw model chemistry keyword source")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific model chemistry options"
    )


class QMTaskRequest(BaseDataClassWithUnit):
    task_type: str = Field(description="Canonical task type, such as sp, opt, freq, or scan")
    enabled: bool = Field(default=True, description="Whether this task is requested")
    derivative_order: int | None = Field(
        default=None, description="Requested derivative order when applicable"
    )
    target_state: int | None = Field(default=None, description="Target electronic state/root")
    transition_state: bool = Field(
        default=False, description="Whether TS optimization is requested"
    )
    scan: bool = Field(default=False, description="Whether this task includes a coordinate scan")
    properties: list[str] = Field(default_factory=list, description="Requested properties")
    source_keywords: list[str] = Field(default_factory=list, description="Source keyword tokens")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific task options"
    )


class ExcitedStateRequest(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether excited-state treatment is requested")
    family: str | None = Field(
        default=None, description="Excited-state method family, such as TDDFT, CIS, or EOM-CCSD"
    )
    reference_method: str | None = Field(default=None, description="Reference method")
    nroots: int | None = Field(default=None, description="Requested number of roots")
    root: int | None = Field(default=None, description="Primary target root")
    secondary_root: int | None = Field(default=None, description="Secondary target root")
    roots: list[int] = Field(default_factory=list, description="Explicit requested roots")
    singlets: bool = Field(default=False, description="Whether singlet states are requested")
    triplets: bool = Field(default=False, description="Whether triplet states are requested")
    spin_flip: bool = Field(default=False, description="Whether spin-flip treatment is requested")
    follow_root: bool = Field(default=False, description="Whether root following is requested")
    properties: list[str] = Field(
        default_factory=list, description="Requested excited-state properties"
    )
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific excited-state options"
    )


class ActiveSpace(BaseDataClassWithUnit):
    electrons: int | None = Field(default=None, description="Number of active electrons")
    orbitals: int | None = Field(default=None, description="Number of active orbitals")
    roots: int | None = Field(default=None, description="Number of averaged or targeted roots")
    active_orbitals: list[int] = Field(
        default_factory=list, description="0-based active orbital indices"
    )
    inactive_orbitals: list[int] = Field(
        default_factory=list, description="0-based inactive orbital indices"
    )
    frozen_orbitals: list[int] = Field(
        default_factory=list, description="0-based frozen orbital indices"
    )
    raw: str = Field(default="", description="Raw active-space expression")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific active-space options"
    )


class MultireferenceStateBlock(BaseDataClassWithUnit):
    multiplicity: int | None = Field(default=None, description="Spin multiplicity")
    irrep: str | None = Field(default=None, description="Irreducible representation label")
    nroots: int | None = Field(default=None, description="Number of roots in this block")
    excitations: str | None = Field(default=None, description="Excitation selection")
    refs: str | None = Field(default=None, description="Reference-space definition")
    active_space: ActiveSpace | None = Field(default=None, description="Active-space summary")
    raw_lines: list[str] = Field(default_factory=list, description="Raw state-block lines")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific state-block options"
    )


class MultireferenceRequest(BaseDataClassWithUnit):
    enabled: bool = Field(
        default=False, description="Whether multi-reference treatment is requested"
    )
    method: str | None = Field(default=None, description="Multi-reference method")
    reference_method: str | None = Field(default=None, description="Reference method")
    ci_type: str | None = Field(default=None, description="Configuration interaction type")
    active_space: ActiveSpace | None = Field(
        default=None, description="Global active-space summary"
    )
    state_blocks: list[MultireferenceStateBlock] = Field(
        default_factory=list, description="Requested spin/symmetry/root blocks"
    )
    thresholds: dict[str, float] = Field(default_factory=dict, description="Numeric thresholds")
    corrections: list[str] = Field(default_factory=list, description="Requested corrections")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific multi-reference options"
    )


class ExplicitSolventRequest(BaseDataClassWithUnit):
    enabled: bool = Field(
        default=False, description="Whether explicit-solvent placement is requested"
    )
    solvent_model: str | None = Field(default=None, description="Underlying implicit solvent model")
    solvent: str | None = Field(
        default=None, description="Solvent name chosen for explicit placement"
    )
    solvent_file: str | None = Field(default=None, description="Custom solvent file, if provided")
    nsolv: int | None = Field(default=None, description="Number of explicit solvent molecules")
    cluster_mode: str | None = Field(default=None, description="SOLVATOR cluster mode")
    droplet: bool = Field(default=False, description="Whether droplet mode is enabled")
    radius: float | None = Field(default=None, description="Target droplet radius")
    fixsolute: bool = Field(default=True, description="Whether the solute is kept frozen")
    vacuumsearch: bool = Field(default=False, description="Whether vacuum search is enabled")
    randomsolv: bool = Field(default=False, description="Whether random placement is enabled")
    printlevel: str | None = Field(default=None, description="SOLVATOR print level")
    source_blocks: list[str] = Field(default_factory=list, description="Source input block names")
    options: dict[str, Any] = Field(
        default_factory=dict, description="Program-specific explicit-solvent options"
    )


__all__ = [
    "ActiveSpace",
    "ExplicitSolventRequest",
    "ExcitedStateRequest",
    "MultireferenceRequest",
    "MultireferenceStateBlock",
    "QMBasisSet",
    "QMModelChemistry",
    "QMResourceRequest",
    "QMTaskRequest",
]
