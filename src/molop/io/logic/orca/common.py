from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Literal

import numpy as np
from pydantic import BaseModel, Field, model_validator
from rdkit import Chem
from typing_extensions import Self

from molop.io.base_models.DataClasses import (
    AtomInInternalCoords,
    CoordinateContainer,
    CoordinateParameters,
    ImplicitSolvation,
    InternalCoords,
)
from molop.unit import atom_ureg


class ORCAKeywordLine(BaseModel):
    text: str = Field(default="", description="ORCA keyword line without leading !")


class ORCACommentLine(BaseModel):
    text: str = Field(default="", description="ORCA comment line without leading #")


class ORCABlockLine(BaseModel):
    text: str = Field(default="", description="Line inside an ORCA % block")


class ORCABlock(BaseModel):
    name: str = Field(default="", description="Block name without leading %")
    lines: list[ORCABlockLine] = Field(default_factory=list, description="Block body lines")
    raw_header: str = Field(default="", description="Original block header line")
    raw_text: str = Field(default="", description="Original block text")

    def body_text(self) -> str:
        return "\n".join(line.text for line in self.lines)


class ORCAAtomBasisOverride(BaseModel):
    kind: Literal["newgto", "newauxgto"] = Field(description="Atom-level basis directive")
    basis_set: str | None = Field(default=None, description="Named basis set, if present")
    tokens: list[str] = Field(default_factory=list, description="Directive tokens")


class ORCAOutputPrintSetting(BaseModel):
    target: str = Field(description="ORCA output print target inside Print[...]")
    value: str = Field(description="Assigned output print value")


class ORCAMultiReferenceNewBlock(BaseModel):
    multiplicity: int | None = Field(default=None, description="Block multiplicity")
    irrep: str | None = Field(default=None, description="Block irrep label")
    nroots: int | None = Field(default=None, description="Number of roots for this block")
    excitations: str | None = Field(default=None, description="Excitation selection mode")
    refs: str | None = Field(default=None, description="Reference space definition")
    raw_lines: list[str] = Field(default_factory=list, description="Raw lines for this NewBlock")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled NewBlock options"
    )


class ORCAMultiReferenceSemantic(BaseModel):
    enabled: bool = Field(default=False, description="Whether a multi-reference task was found")
    method: str | None = Field(default=None, description="Canonical multi-reference method")
    ci_type: str | None = Field(default=None, description="Raw CIType value from the %mrci block")
    reference_method: str | None = Field(
        default=None, description="Underlying reference method when it can be inferred"
    )
    source_blocks: list[str] = Field(
        default_factory=list, description="Multi-reference related block names used for parsing"
    )
    block_options: dict[str, dict[str, str | None]] = Field(
        default_factory=dict, description="Normalized option maps per multi-reference block"
    )
    new_blocks: list[ORCAMultiReferenceNewBlock] = Field(
        default_factory=list, description="Structured CI blocks from %mrci"
    )
    ewin: tuple[float, float] | None = Field(
        default=None, description="Selected orbital energy window, if available"
    )
    tsel: float | None = Field(default=None, description="Selection threshold")
    tpre: float | None = Field(default=None, description="Pre-diagonalization threshold")
    tnat: float | None = Field(default=None, description="Natural orbital threshold")
    etol: float | None = Field(default=None, description="Energy convergence threshold")
    rtol: float | None = Field(default=None, description="Residual convergence threshold")
    solver: str | None = Field(default=None, description="Solver selection")
    int_mode: str | None = Field(default=None, description="Integral transformation mode")
    use_ivos: bool = Field(default=False, description="Whether IVOs are enabled")
    all_singles: bool = Field(default=False, description="Whether all singles are included")
    do_ddcimp2: bool = Field(default=False, description="Whether the DDCI-MP2 correction is on")
    do_nat_orbs: bool = Field(default=False, description="Whether natural orbitals are requested")
    eunsel_opt: str | None = Field(default=None, description="Unselected energy correction mode")
    davidson_opt: str | None = Field(default=None, description="Davidson correction mode")
    partitioning: str | None = Field(default=None, description="Partitioning choice")
    fopt: str | None = Field(default=None, description="Fock operator choice")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled multi-reference options"
    )


class ORCAExcitedStateSemantic(BaseModel):
    enabled: bool = Field(default=False, description="Whether an excited-state task was found")
    family: str | None = Field(
        default=None, description="Canonical excited-state family such as TDDFT or EOM-CCSD"
    )
    reference_method: str | None = Field(
        default=None, description="Underlying reference method when it is spelled out"
    )
    source_blocks: list[str] = Field(
        default_factory=list, description="Excited-state related block names used for parsing"
    )
    block_options: dict[str, dict[str, str | None]] = Field(
        default_factory=dict, description="Normalized option maps per excited-state block"
    )
    nroots: int | None = Field(default=None, description="Number of excited roots")
    iroot: int | None = Field(default=None, description="Target root index")
    jroot: int | None = Field(default=None, description="Secondary root index")
    followiroot: bool = Field(default=False, description="Whether IROOT following is enabled")
    triplets: bool = Field(default=False, description="Whether triplets are requested")
    sf: bool = Field(default=False, description="Whether spin-flip is requested")
    nacme: bool = Field(default=False, description="Whether NACME evaluation is requested")
    etf: bool = Field(default=False, description="Whether ETF is requested")
    dosoc: bool = Field(default=False, description="Whether SOC is requested in CIS")
    doalpha: bool = Field(
        default=False, description="Whether alpha-channel only calculation is requested"
    )
    rootwise: bool = Field(default=False, description="Whether rootwise solving is requested")
    do_dbfilter: bool = Field(
        default=False, description="Whether STEOM doubly excited filtering is enabled"
    )
    do_store_steom: bool = Field(
        default=False, description="Whether STEOM intermediates are stored"
    )
    do_simple_dens: bool = Field(
        default=False, description="Whether STEOM simple density is disabled"
    )
    add_l2_term: bool = Field(default=False, description="Whether DLPNO STEOM L2 term is enabled")
    do_full_semiclassical: bool = Field(
        default=False, description="Whether full semiclassical treatment is enabled"
    )
    do_higher_moments: bool = Field(default=False, description="Whether higher moments are enabled")
    firkeepfirstref: bool = Field(
        default=False, description="Whether first reference is kept in follow-root mode"
    )
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled excited-state block options"
    )


class ORCAExplicitSolventSemantic(BaseModel):
    enabled: bool = Field(default=False, description="Whether a SOLVATOR block was found")
    solvent_model: str | None = Field(default=None, description="Implicit model used by SOLVATOR")
    solvent: str | None = Field(default=None, description="Explicit solvent name")
    solvent_file: str | None = Field(default=None, description="Custom solvent file")
    nsolv: int | None = Field(default=None, description="Number of solvent molecules")
    cluster_mode: str | None = Field(default=None, description="Cluster mode")
    droplet: bool = Field(default=False, description="Droplet mode")
    radius: float | None = Field(default=None, description="Droplet radius")
    fixsolute: bool = Field(default=True, description="Fix solute during placement")
    vacuumsearch: bool = Field(default=False, description="Vacuum search")
    randomsolv: bool = Field(default=False, description="Random solvent placement")
    printlevel: str | None = Field(default=None, description="Print level")
    source_blocks: list[str] = Field(default_factory=list, description="Source block names")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled explicit-solvent options"
    )


class ORCAGeometryAtom(BaseModel):
    symbol: str = Field(description="Atom symbol")
    atomic_number: int | None = Field(default=None, description="Atomic number")
    x: float | None = Field(default=None, description="X coordinate")
    y: float | None = Field(default=None, description="Y coordinate")
    z: float | None = Field(default=None, description="Z coordinate")
    charge: float | None = Field(default=None, description="Point charge, if present")
    is_dummy: bool = Field(default=False, description="Whether the atom is a dummy atom")
    is_ghost: bool = Field(default=False, description="Whether the atom is a ghost atom")
    fragment_id: int | None = Field(default=None, description="Fragment id, if present")
    frozen: bool = Field(default=False, description="Whether the atom is frozen")
    isotope: str | None = Field(default=None, description="Isotope token")
    nuclear_charge: str | None = Field(default=None, description="Nuclear charge token")
    basis_set: str | None = Field(default=None, description="Atom-level orbital basis override")
    auxiliary_basis_set: str | None = Field(
        default=None, description="Atom-level auxiliary basis override"
    )
    basis_overrides: list[ORCAAtomBasisOverride] = Field(default_factory=list)
    internal_coord: AtomInInternalCoords | None = Field(
        default=None,
        description="Internal-coordinate row for non-Cartesian ORCA geometry",
    )


class ORCAGeometry(CoordinateContainer[ORCAGeometryAtom]):
    ctype: Literal[
        "xyz", "cart", "cartesian", "int", "internal", "gzmt", "xyzfile", "gzmtfile", "pdbfile"
    ] = "xyz"
    charge: int = 0
    multiplicity: int = 1
    units: str | None = None
    external_path: str | None = None
    internal_coords: InternalCoords | None = None
    coordinate_parameters: CoordinateParameters | None = None
    point_charges: list[dict[str, float]] = Field(default_factory=list)
    source: Literal["star", "percent_coords", "unknown"] = "unknown"


def populate_common_orca_qm_fields(
    target: Any,
    *,
    default_version: str | None = None,
) -> None:
    """Populate ORCA common QM fields from existing legacy/container state."""
    target.qm_software = "ORCA"
    if default_version is not None and not target.qm_software_version:
        target.qm_software_version = default_version
    target.backfill_common_qm_containers_from_legacy()
    target.project_common_qm_fields()


class ORCACommonQMFieldsMixin:
    qm_software: str = Field(default="ORCA")
    auxiliary_basis_set: str = Field(
        default="",
        description="ORCA auxiliary basis keyword derived from structured input metadata",
    )
    dispersion_correction: str = Field(
        default="",
        description="ORCA dispersion correction keyword derived from structured input metadata",
    )

    @model_validator(mode="after")
    def _normalize_common_orca_fields(self) -> Self:
        populate_common_orca_qm_fields(self)
        return self


class ORCAOutputQMFieldsMixin(ORCACommonQMFieldsMixin):
    input_file_name: str = Field(default="", description="ORCA input file name printed in output")


def project_orca_geometry_to_qm_frame(
    target: Any,
    geometry: ORCAGeometry | None,
    *,
    default_version: str | None = None,
) -> bool:
    """Project ORCA input geometry fields onto a QM frame-like target."""
    if geometry is None:
        populate_common_orca_qm_fields(target, default_version=default_version)
        return False

    target.charge = geometry.charge
    target.multiplicity = geometry.multiplicity
    has_mixed_basis = any(atom.basis_overrides for atom in geometry)
    if geometry.internal_coords is not None:
        pt = Chem.GetPeriodicTable()
        target.atoms = [pt.GetAtomicNumber(symbol) for symbol in geometry.get_symbols()]
        target.coords = geometry.internal_coords.to_cartesian_coords()
    elif geometry:
        pt = Chem.GetPeriodicTable()
        atoms: list[int] = []
        coords: list[list[float]] = []
        for atom in geometry.real_atoms():
            if atom.x is None or atom.y is None or atom.z is None:
                continue
            atomic_number = atom.atomic_number or pt.GetAtomicNumber(atom.symbol)
            if atomic_number <= 0:
                continue
            atoms.append(atomic_number)
            coords.append([atom.x, atom.y, atom.z])
        if atoms:
            target.atoms = atoms
            target.coords = np.asarray(coords, dtype=float) * atom_ureg.angstrom
    populate_common_orca_qm_fields(target, default_version=default_version)
    return has_mixed_basis


ORCA_PRINTED_INPUT_METADATA_FIELDS = (
    "keywords",
    "method",
    "basis_set",
    "functional",
    "model_chemistry",
    "task_requests",
    "excited_state_requests",
    "multireference_requests",
    "resources_raw",
    "request_num_cpu",
    "request_memory",
    "resource_request",
    "charge",
    "multiplicity",
    "auxiliary_basis_set",
    "dispersion_correction",
)


def project_orca_printed_input_metadata(frame: Any | Mapping[str, Any]) -> dict[str, Any]:
    """Project an ORCA input frame onto ORCA output file metadata fields."""
    metadata: dict[str, Any] = {}
    for field in ORCA_PRINTED_INPUT_METADATA_FIELDS:
        if isinstance(frame, Mapping):
            if field in frame:
                metadata[field] = frame[field]
            continue
        if hasattr(frame, field):
            metadata[field] = getattr(frame, field)

    geometry = (
        frame.get("geometry") if isinstance(frame, Mapping) else getattr(frame, "geometry", None)
    )
    if geometry is not None:
        metadata.setdefault("charge", getattr(geometry, "charge", 0))
        metadata.setdefault("multiplicity", getattr(geometry, "multiplicity", 1))

    model_chemistry = (
        frame.get("model_chemistry")
        if isinstance(frame, Mapping)
        else getattr(frame, "model_chemistry", None)
    )
    if solvent_model := getattr(model_chemistry, "solvation_model", None):
        metadata["solvent"] = ImplicitSolvation(
            solvent_model=solvent_model,
            solvent=getattr(model_chemistry, "solvent", None),
        )
    return metadata
