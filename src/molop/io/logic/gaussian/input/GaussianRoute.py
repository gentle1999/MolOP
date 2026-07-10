from __future__ import annotations

from typing import Any, Protocol

from pydantic import Field, model_validator
from typing_extensions import Self

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.base_models.DataClasses import (
    ExcitedStateRequest,
    QMModelChemistry,
    QMTaskRequest,
)


SEMI_EMPIRICAL_METHODS = ("am1", "pm3", "pm6", "pm7", "pddg", "indo", "cndo", "pm3mm")


class GaussianRouteToken(BaseDataClassWithUnit):
    raw: str = Field(description="Original route token")
    normalized: str = Field(description="Normalized route token")
    kind: str = Field(default="unknown", description="Token classification")
    key: str | None = Field(default=None, description="Normalized token key")
    scalar_value: str | None = Field(default=None, description="Scalar token value")
    param_map: dict[str, str | None] | None = Field(
        default=None, description="Structured token parameters"
    )


class GaussianModelChemistry(BaseDataClassWithUnit):
    method_token: str | None = Field(default=None, description="Raw method token")
    method_family: str | None = Field(default=None, description="Method family")
    functional: str | None = Field(default=None, description="DFT or double-hybrid functional")
    basis_set: str | None = Field(default=None, description="Basis-set token")
    auxiliary_basis_set: str | None = Field(default=None, description="Auxiliary fitting basis")
    spin_qualifier: str | None = Field(default=None, description="Restricted/unrestricted prefix")
    basis_family: str | None = Field(default=None, description="Basis-set family")
    basis_has_diffuse: bool = Field(
        default=False, description="Whether basis includes diffuse augmentation"
    )
    basis_polarization: list[str] = Field(
        default_factory=list, description="Basis polarization markers"
    )
    low_level: GaussianModelChemistry | None = Field(
        default=None, description="Lower-level model chemistry for composite `high//low` routes"
    )


class GaussianRouteOptionMap(BaseDataClassWithUnit):
    keyword: str = Field(description="Keyword name")
    scalar_value: str | None = Field(default=None, description="Scalar value if present")
    params: dict[str, str | None] = Field(default_factory=dict, description="Structured parameters")


class GaussianRouteDiagnostics(BaseDataClassWithUnit):
    confidence: float = Field(default=0.0, description="0-1 confidence for semantic extraction")
    messages: list[str] = Field(default_factory=list, description="Diagnostic messages")


class GaussianOptOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether Opt was requested")
    restart: bool = Field(default=False, description="Whether Restart was requested")
    transition_state: bool = Field(
        default=False, description="Whether TS optimization is requested"
    )
    saddle_order: int | None = Field(default=None, description="Requested saddle-point order")
    very_tight: bool = Field(
        default=False, description="Whether VeryTight convergence is requested"
    )
    tight: bool = Field(default=False, description="Whether Tight convergence is requested")
    loose: bool = Field(default=False, description="Whether Loose convergence is requested")
    max_cycles: int | None = Field(default=None, description="Maximum optimization cycles")
    max_step: int | None = Field(default=None, description="Maximum optimization step")
    recalc_fc: int | None = Field(default=None, description="Force constant recalculation interval")
    read_fc: bool = Field(default=False, description="Whether ReadFC was requested")
    read_cartesian_fc: bool = Field(
        default=False, description="Whether RCFC/ReadCartesianFC was requested"
    )
    calc_fc: bool = Field(default=False, description="Whether CalcFC was requested")
    calc_all: bool = Field(default=False, description="Whether CalcAll was requested")
    calc_hf_fc: bool = Field(default=False, description="Whether CalcHFFC was requested")
    has_modredundant: bool = Field(
        default=False, description="Whether ModRedundant-style constraints are implied"
    )
    expert: bool | None = Field(
        default=None, description="Expert/NoExpert setting when explicitly controlled"
    )
    eigen_test: bool | None = Field(
        default=None, description="EigenTest setting when explicitly controlled"
    )
    coordinate_system: str | None = Field(
        default=None, description="Optimization coordinate system"
    )
    qst_mode: str | None = Field(default=None, description="QST2/QST3 mode if present")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled opt options"
    )


class GaussianFreqOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether Freq was requested")
    anharmonic: bool = Field(default=False, description="Whether anharmonic analysis is requested")
    read_anharm: bool = Field(default=False, description="Whether ReadAnharm was requested")
    projected: bool = Field(
        default=False, description="Whether projected frequencies are requested"
    )
    tprojected: bool = Field(default=False, description="Whether TProjected was requested")
    hindered_rotor: bool = Field(default=False, description="Whether HinderedRotor was requested")
    vibrot: bool = Field(default=False, description="Whether VibRot was requested")
    polar: bool = Field(default=False, description="Whether Polar was requested")
    hpmodes: bool = Field(
        default=False, description="Whether high-precision mode vectors are requested"
    )
    read_fcht: bool = Field(default=False, description="Whether ReadFCHT was requested")
    read_isotopes: bool = Field(default=False, description="Whether ReadIsotopes was requested")
    select_normal_modes: bool = Field(
        default=False, description="Whether SelectNormalModes was requested"
    )
    save_normal_modes: bool = Field(
        default=False, description="Whether SaveNormalModes was requested"
    )
    vcd: bool = Field(default=False, description="Whether VCD is requested")
    roa: bool = Field(default=False, description="Whether ROA is requested")
    raman: bool = Field(default=False, description="Whether Raman analysis is requested")
    no_raman: bool = Field(
        default=False, description="Whether Raman analysis is explicitly disabled"
    )
    cphf_rd_freq: bool = Field(default=False, description="Whether CPHF=RdFreq was requested")
    layer: str | None = Field(default=None, description="Layer selector")
    atoms: str | None = Field(default=None, description="Included atoms selector")
    not_atoms: str | None = Field(default=None, description="Excluded atoms selector")
    temperature: float | None = Field(default=None, description="Requested temperature")
    pressure: float | None = Field(default=None, description="Requested pressure")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled freq options"
    )


class GaussianTDOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether TD was requested")
    nstates: int | None = Field(default=None, description="Number of excited states")
    root: int | None = Field(default=None, description="Target root")
    singlets: bool = Field(default=False, description="Whether singlets were requested")
    triplets: bool = Field(default=False, description="Whether triplets were requested")
    tda: bool = Field(default=False, description="Whether TDA was requested")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled TD options"
    )


class GaussianSCRFOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether SCRF was requested")
    model: str | None = Field(default=None, description="SCRF model, e.g. smd")
    model_family: str | None = Field(default=None, description="Canonical SCRF model family")
    solvent: str | None = Field(default=None, description="Solvent name if present")
    read: bool = Field(default=False, description="Whether SCRF=Read is requested")
    iefpcm: bool = Field(default=False, description="Whether IEFPCM is requested")
    cpcm: bool = Field(default=False, description="Whether CPCM is requested")
    smd: bool = Field(default=False, description="Whether SMD is requested")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled SCRF options"
    )


class GaussianPopOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether population analysis was requested")
    none: bool = Field(default=False, description="Whether Pop=None was requested")
    full: bool = Field(default=False, description="Whether full population output is requested")
    nbo: bool = Field(default=False, description="Whether NBO output is requested")
    nbo_read: bool = Field(default=False, description="Whether NBORead was requested")
    nbo6_read: bool = Field(default=False, description="Whether NBO6Read was requested")
    nbo7_read: bool = Field(default=False, description="Whether NBO7Read was requested")
    hirshfeld: bool = Field(
        default=False, description="Whether Hirshfeld population analysis is requested"
    )
    cm5: bool = Field(default=False, description="Whether CM5 charges are requested")
    mk: bool = Field(default=False, description="Whether MK charges are requested")
    chelpg: bool = Field(default=False, description="Whether CHelpG charges are requested")
    orbitals: int | None = Field(default=None, description="Requested number of orbitals to print")
    read_radii: bool = Field(default=False, description="Whether ReadRadii was requested")
    read_at_radii: bool = Field(default=False, description="Whether ReadAtRadii was requested")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled Pop options"
    )


class GaussianGeomOptions(BaseDataClassWithUnit):
    enabled: bool = Field(default=False, description="Whether Geom was requested")
    mode: str | None = Field(
        default=None, description="Primary Geom mode, e.g. allcheck/checkpoint/check"
    )
    checkpoint: bool = Field(default=False, description="Whether Checkpoint was requested")
    allcheck: bool = Field(default=False, description="Whether AllCheck was requested")
    check: bool = Field(default=False, description="Whether Check was requested")
    huge: bool = Field(default=False, description="Whether Huge was requested")
    modify: bool = Field(default=False, description="Whether Modify was requested")
    new_definition: bool = Field(default=False, description="Whether NewDefinition was requested")
    new_redundant: bool = Field(default=False, description="Whether NewRedundant was requested")
    no_test: bool = Field(default=False, description="Whether NoTest was requested")
    gic: bool = Field(default=False, description="Whether GIC was requested")
    add_gic: bool = Field(default=False, description="Whether AddGIC was requested")
    read_all_gic: bool = Field(default=False, description="Whether ReadAllGIC was requested")
    no_gic: bool = Field(default=False, description="Whether NoGIC was requested")
    connectivity: bool = Field(default=False, description="Whether Connectivity was requested")
    mod_connectivity: bool = Field(
        default=False, description="Whether ModConnectivity was requested"
    )
    gen_connectivity: bool = Field(
        default=False, description="Whether GenConnectivity was requested"
    )
    zm_connectivity: bool = Field(default=False, description="Whether ZMConnectivity was requested")
    distance: bool = Field(default=False, description="Whether Distance was requested")
    no_distance: bool = Field(default=False, description="Whether NoDistance was requested")
    cangle: bool = Field(default=False, description="Whether CAngle was requested")
    angle: bool = Field(default=False, description="Whether Angle was requested")
    no_angle: bool = Field(default=False, description="Whether NoAngle was requested")
    cdihedral: bool = Field(default=False, description="Whether CDihedral was requested")
    dihedral: bool = Field(default=False, description="Whether Dihedral was requested")
    no_dihedral: bool = Field(default=False, description="Whether NoDihedral was requested")
    print_input_orient: bool = Field(
        default=False, description="Whether PrintInputOrient was requested"
    )
    print: bool = Field(default=False, description="Whether Print was requested")
    step: int | None = Field(default=None, description="Step index if requested")
    ngeom: int | None = Field(default=None, description="NGeom index if requested")
    extra_options: dict[str, str | None] = Field(
        default_factory=dict, description="Unmodeled Geom options"
    )


class GaussianRouteSemantic(BaseDataClassWithUnit):
    raw_route: str = Field(default="", description="Raw Gaussian route section")
    normalized_route: str = Field(default="", description="Whitespace-normalized route string")
    dieze_tag: str | None = Field(default=None, description="Normalized route tag")
    tokens: list[GaussianRouteToken] = Field(
        default_factory=list, description="Normalized route tokens"
    )
    model_chemistry: GaussianModelChemistry = Field(
        default_factory=GaussianModelChemistry, description="Model chemistry summary"
    )
    job_types: list[str] = Field(default_factory=list, description="Detected job types")
    route_modifiers: list[str] = Field(default_factory=list, description="Detected route modifiers")
    capabilities: list[str] = Field(default_factory=list, description="Capability tags")
    option_maps: dict[str, GaussianRouteOptionMap] = Field(
        default_factory=dict, description="Structured keyword option maps"
    )
    opt_options: GaussianOptOptions = Field(default_factory=GaussianOptOptions)
    freq_options: GaussianFreqOptions = Field(default_factory=GaussianFreqOptions)
    td_options: GaussianTDOptions = Field(default_factory=GaussianTDOptions)
    scrf_options: GaussianSCRFOptions = Field(default_factory=GaussianSCRFOptions)
    pop_options: GaussianPopOptions = Field(default_factory=GaussianPopOptions)
    geom_options: GaussianGeomOptions = Field(default_factory=GaussianGeomOptions)
    solvation_model: str | None = Field(default=None, description="Detected SCRF model")
    empirical_dispersion: str | None = Field(default=None, description="Detected dispersion model")
    checkpoint_geometry_mode: str | None = Field(
        default=None, description="Detected Geom=Check/AllCheck mode"
    )
    external_program: str | None = Field(
        default=None, description="Detected external program reference"
    )
    unknown_tokens: list[str] = Field(
        default_factory=list, description="Tokens not semantically classified"
    )
    diagnostics: GaussianRouteDiagnostics = Field(
        default_factory=GaussianRouteDiagnostics, description="Best-effort parser diagnostics"
    )

    def to_route_dict(self) -> dict[str, Any]:
        projected: dict[str, Any] = {}
        for token in self.tokens:
            if token.kind == "dieze-tag":
                continue

            if token.kind == "model-chemistry":
                if self.model_chemistry.method_token:
                    projected.setdefault(self.model_chemistry.method_token.lower(), None)
                if self.model_chemistry.basis_set:
                    projected.setdefault(self.model_chemistry.basis_set.lower(), None)
                if self.model_chemistry.auxiliary_basis_set:
                    projected.setdefault(self.model_chemistry.auxiliary_basis_set.lower(), None)
                continue

            if token.param_map is not None and token.key is not None:
                if token.kind == "basis-set" and all(k in {"d", "p"} for k in token.param_map):
                    projected[token.normalized] = None
                else:
                    projected[token.key] = token.param_map
                continue

            if token.scalar_value is not None and token.key is not None:
                projected[token.key] = token.scalar_value
                continue

            if token.key is not None:
                projected[token.key] = None
            else:
                projected[token.normalized] = None

        return projected


class _CommonQMInputTarget(Protocol):
    keywords: str
    method: str
    basis_set: str
    functional: str
    model_chemistry: QMModelChemistry
    task_requests: list[QMTaskRequest]
    excited_state_requests: list[ExcitedStateRequest]

    def backfill_common_qm_containers_from_legacy(self) -> None: ...

    def project_common_qm_fields(self) -> None: ...


def as_common_qm_input_target(target: Any) -> _CommonQMInputTarget:
    return target


_GAUSSIAN_TASK_TYPES = {
    "sp": "sp",
    "opt": "opt",
    "freq": "freq",
    "force": "force",
    "irc": "irc",
    "ircmax": "irc",
    "scan": "scan",
    "polar": "property",
    "td": "excited_state",
    "nmr": "property",
    "pop": "population_analysis",
    "stable": "wavefunction_stability",
    "volume": "property",
    "admp": "dynamics",
    "bomd": "dynamics",
    "oniom": "composite",
}


def _with_dispersion_suffix(functional: str | None, dispersion: str | None) -> str | None:
    if not functional:
        return functional
    if not dispersion:
        return functional
    suffix = dispersion.upper()
    if functional.upper().endswith(f"-{suffix}"):
        return functional
    return f"{functional}-{suffix}"


def _method_family_allows_functional(method_family: str | None) -> bool:
    if method_family is None:
        return False
    normalized = method_family.upper().replace("_", "-").replace(" ", "-")
    return normalized in {"DFT", "DOUBLE-HYBRID", "DOUBLE-HYBRID-DFT"}


def _semantic_or_legacy_functional(
    semantic_model: Any,
    legacy_functional: str,
    legacy_method: str,
) -> str | None:
    if semantic_model.functional:
        return semantic_model.functional
    method_family = semantic_model.method_family or legacy_method or None
    if legacy_functional and _method_family_allows_functional(method_family):
        return legacy_functional
    return None


def _source_keywords_for_job(semantic_route: GaussianRouteSemantic, job_type: str) -> list[str]:
    return [
        token.raw
        for token in semantic_route.tokens
        if token.normalized.split("=", 1)[0].split("(", 1)[0] == job_type
    ]


def build_gaussian_model_chemistry(
    semantic_route: GaussianRouteSemantic,
    *,
    keywords: str,
    legacy_method: str = "",
    legacy_basis_set: str = "",
    legacy_functional: str = "",
) -> QMModelChemistry:
    semantic_model = semantic_route.model_chemistry
    functional = _with_dispersion_suffix(
        _semantic_or_legacy_functional(semantic_model, legacy_functional, legacy_method),
        semantic_route.empirical_dispersion,
    )
    basis_set = semantic_model.basis_set or legacy_basis_set or None
    if basis_set and basis_set.lower() == "genecp" and legacy_basis_set == "pseudopotential":
        basis_set = legacy_basis_set
    return QMModelChemistry(
        method_family=semantic_model.method_family or legacy_method or None,
        method=semantic_model.method_token or legacy_method or None,
        functional=functional,
        basis_set=basis_set,
        auxiliary_basis_set=semantic_model.auxiliary_basis_set,
        dispersion_correction=(
            semantic_route.empirical_dispersion.upper()
            if semantic_route.empirical_dispersion
            else None
        ),
        solvation_model=semantic_route.scrf_options.model_family,
        solvent=semantic_route.scrf_options.solvent,
        spin_treatment=semantic_model.spin_qualifier,
        raw_keywords=keywords,
        options={
            "dieze_tag": semantic_route.dieze_tag,
            "basis_family": semantic_model.basis_family,
            "basis_has_diffuse": semantic_model.basis_has_diffuse,
            "basis_polarization": semantic_model.basis_polarization,
            "route_modifiers": semantic_route.route_modifiers,
            "capabilities": semantic_route.capabilities,
            "unknown_tokens": semantic_route.unknown_tokens,
        },
    )


def build_gaussian_task_requests(semantic_route: GaussianRouteSemantic) -> list[QMTaskRequest]:
    job_types = list(semantic_route.job_types)
    if semantic_route.opt_options.enabled and "opt" not in job_types:
        job_types.append("opt")
    if semantic_route.freq_options.enabled and "freq" not in job_types:
        job_types.append("freq")
    if semantic_route.td_options.enabled and "td" not in job_types:
        job_types.append("td")
    if semantic_route.pop_options.enabled and "pop" not in job_types:
        job_types.append("pop")
    if not job_types and semantic_route.model_chemistry.method_token:
        job_types.append("sp")

    tasks: list[QMTaskRequest] = []
    for job_type in dict.fromkeys(job_types):
        task_type = _GAUSSIAN_TASK_TYPES.get(job_type, job_type)
        options: dict[str, Any] = {}
        derivative_order: int | None = None
        transition_state = False
        scan = False
        properties: list[str] = []
        target_state: int | None = None

        if job_type == "opt":
            opt = semantic_route.opt_options
            derivative_order = 1
            transition_state = opt.transition_state
            scan = opt.has_modredundant or opt.qst_mode is not None
            options = opt.model_dump()
        elif job_type == "freq":
            freq = semantic_route.freq_options
            derivative_order = 2
            properties = [
                name
                for name in (
                    "anharmonic",
                    "projected",
                    "polar",
                    "vcd",
                    "roa",
                    "raman",
                )
                if getattr(freq, name)
            ]
            options = freq.model_dump()
        elif job_type == "force":
            derivative_order = 1
        elif job_type == "td":
            target_state = semantic_route.td_options.root
            properties = ["tda"] if semantic_route.td_options.tda else []
            options = semantic_route.td_options.model_dump()
        elif job_type == "pop":
            options = semantic_route.pop_options.model_dump()
        elif job_type == "stable":
            properties = ["stability"]

        tasks.append(
            QMTaskRequest(
                task_type=task_type,
                derivative_order=derivative_order,
                target_state=target_state,
                transition_state=transition_state,
                scan=scan,
                properties=properties,
                source_keywords=_source_keywords_for_job(semantic_route, job_type),
                options=options,
            )
        )
    return tasks


def build_gaussian_excited_state_requests(
    semantic_route: GaussianRouteSemantic,
) -> list[ExcitedStateRequest]:
    td_options = semantic_route.td_options
    if not td_options.enabled:
        return []
    family = "TDDFT" if semantic_route.model_chemistry.method_family == "DFT" else "TD"
    return [
        ExcitedStateRequest(
            enabled=True,
            family=family,
            nroots=td_options.nstates,
            root=td_options.root,
            roots=[td_options.root] if td_options.root is not None else [],
            singlets=td_options.singlets,
            triplets=td_options.triplets,
            properties=["tda"] if td_options.tda else [],
            source_blocks=["route"],
            options=td_options.model_dump(),
        )
    ]


def populate_common_gaussian_qm_containers(
    target: _CommonQMInputTarget,
    semantic_route: GaussianRouteSemantic,
) -> None:
    target.model_chemistry = build_gaussian_model_chemistry(
        semantic_route,
        keywords=target.keywords,
        legacy_method=target.method,
        legacy_basis_set=target.basis_set,
        legacy_functional=target.functional,
    )
    target.task_requests = build_gaussian_task_requests(semantic_route)
    target.excited_state_requests = build_gaussian_excited_state_requests(semantic_route)
    target.backfill_common_qm_containers_from_legacy()
    target.project_common_qm_fields()


def _infer_gaussian_legacy_method(
    semantic_route: GaussianRouteSemantic,
    *,
    keywords: str,
    method: str,
    functional: str,
) -> str:
    if method:
        return method
    lowered_keywords = keywords.lower()
    lowered_functional = functional.lower()
    if any(semi in lowered_keywords for semi in SEMI_EMPIRICAL_METHODS):
        return "SEMI-EMPIRICAL"
    if semantic_route.model_chemistry.method_family is not None:
        return semantic_route.model_chemistry.method_family
    if lowered_functional in {"hf", "rhf", "uhf", "rohf"} or lowered_functional.endswith("hf"):
        return "HF"
    if lowered_functional.endswith("fc"):
        return "FC"
    if functional:
        return "DFT"
    return ""


def populate_gaussian_legacy_qm_fields_from_semantic(
    target: Any,
    semantic_route: GaussianRouteSemantic,
) -> None:
    """Populate legacy flat QM fields from Gaussian route semantics, then sync containers."""
    qm_target = as_common_qm_input_target(target)
    semantic_model = semantic_route.model_chemistry
    if qm_target.basis_set.lower() == "genecp":
        qm_target.basis_set = "pseudopotential"
    elif not qm_target.basis_set and semantic_model.basis_set is not None:
        qm_target.basis_set = semantic_model.basis_set

    if not qm_target.functional and semantic_model.functional is not None:
        qm_target.functional = semantic_model.functional

    inferred_method = _infer_gaussian_legacy_method(
        semantic_route,
        keywords=qm_target.keywords,
        method=qm_target.method,
        functional=qm_target.functional,
    )
    if inferred_method:
        qm_target.method = inferred_method

    if semantic_route.empirical_dispersion:
        qm_target.functional = (
            _with_dispersion_suffix(
                qm_target.functional,
                semantic_route.empirical_dispersion,
            )
            or qm_target.functional
        )

    populate_common_gaussian_qm_containers(qm_target, semantic_route)


class GaussianRouteSemanticFieldsMixin:
    qm_software: str = Field(default="Gaussian")
    options: str = Field(default="", description="options comment")
    title_card: str = Field(default="", description="title card")
    job_type: str = Field(default="", description="Job type")
    keywords: str = Field(
        default="",
        description="Keywords for the QM parameters",
    )
    method: str = Field(
        default="",
        description="QM method used to perform the calculation. e.g. DFT or GFN2-xTB",
    )
    basis_set: str = Field(
        default="",
        description="Basis set used in the QM calculation, only for DFT calculations",
    )
    functional: str = Field(
        default="",
        description="Functional used in the QM calculation, only for DFT calculations",
    )
    semantic_route: GaussianRouteSemantic = Field(
        default_factory=GaussianRouteSemantic,
        description="Structured Gaussian route semantics",
    )

    @model_validator(mode="after")
    def _normalize_gaussian_route_semantic_fields(self) -> Self:
        populate_gaussian_legacy_qm_fields_from_semantic(self, self.semantic_route)
        return self

    @property
    def dieze_tag(self) -> str | None:
        return self.semantic_route.dieze_tag
