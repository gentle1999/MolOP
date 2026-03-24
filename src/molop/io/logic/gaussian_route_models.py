from __future__ import annotations

import re
from typing import Any

from pydantic import Field

from molop.io.base_models.Bases import BaseDataClassWithUnit
from molop.io.patterns.G16Patterns import SEMI_EMPIRICAL_METHODS


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
    functional: str | None = Field(default=None, description="DFT functional or method name")
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


_DFT_FUNCTIONALS = {
    "b3lyp",
    "pbe",
    "pbe0",
    "bp86",
    "tpssh",
    "m06",
    "m06-2x",
    "wb97x-d",
    "r2scan",
    "dsdpbep86",
    "revdsdpbep86",
    "cam-b3lyp",
    "b2plyp",
}

_METHOD_FAMILY_MAP = {
    "hf": "HF",
    "rhf": "HF",
    "uhf": "HF",
    "rohf": "HF",
    "mp2": "MP2",
    "ump2": "MP2",
    "rmp2": "MP2",
    "ccsd": "CCSD",
    "ccsd(t)": "CCSD(T)",
    "uccsd": "CCSD",
    "rccsd": "CCSD",
    "uccsd(t)": "CCSD(T)",
    "rccsd(t)": "CCSD(T)",
    "cis": "CIS",
    "td": "TD",
}

_JOB_TYPE_KEYS = {
    "sp": "SinglePointEnergy",
    "opt": "GeometryOptimization",
    "freq": "Frequencies",
    "force": "ForceConstants",
    "irc": "ReactionPath",
    "ircmax": "ReactionPathMaximum",
    "scan": "PotentialEnergyScan",
    "polar": "Polarizability",
    "td": "ExcitedState",
    "nmr": "NMR",
    "pop": "PopulationAnalysis",
    "stable": "WavefunctionStability",
    "volume": "MolecularVolume",
    "admp": "DirectDynamics",
    "bomd": "BornOppenheimerDynamics",
    "oniom": "ONIOM",
}

_ROUTE_MODIFIER_KEYS = {"ts", "readfc", "calcfc", "modredundant", "qst2", "qst3", "gic", "addgic"}


def _classify_basis_family(token: str) -> str | None:
    lowered = token.lower()
    if lowered.startswith(("def2", "ma-def2")):
        return "def2"
    if "cc-pv" in lowered:
        return "dunning"
    if re.match(r"^\d-\d+g", lowered):
        return "pople"
    if lowered in {"gen", "genecp"}:
        return "general"
    if lowered in {"lanl2dz", "lanl2mb", "sdd", "sddall"}:
        return "ecp"
    if lowered == "ugbs":
        return "ugbs"
    return None


def _basis_polarization_markers(token: str) -> list[str]:
    lowered = token.lower()
    markers: list[str] = []
    if "**" in lowered:
        markers.append("**")
    elif "*" in lowered:
        markers.append("*")
    if m := re.search(r"\(([^\)]+)\)", lowered):
        markers.extend(part.strip() for part in m.group(1).split(",") if part.strip())
    return list(dict.fromkeys(markers))


def _normalize_route_for_semantics(route: str) -> str:
    stripped_lines = [line.strip() for line in route.splitlines() if line.strip()]
    normalized = re.sub(r"\s+", " ", " ".join(stripped_lines)).strip()
    return re.sub(r"^#\s+([nNpPtT])\b", lambda m: f"#{m.group(1)}", normalized)


def _top_level_tokens(route: str) -> list[str]:
    tokens: list[str] = []
    current: list[str] = []
    depth = 0
    for ch in route:
        if ch.isspace() and depth == 0:
            if current:
                tokens.append("".join(current))
                current = []
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(depth - 1, 0)
        current.append(ch)
    if current:
        tokens.append("".join(current))
    return tokens


def _split_top_level_csv(raw: str) -> list[str]:
    if not raw.strip():
        return []
    parts: list[str] = []
    current: list[str] = []
    depth = 0
    for ch in raw:
        if ch == "," and depth == 0:
            token = "".join(current).strip()
            if token:
                parts.append(token)
            current = []
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(depth - 1, 0)
        current.append(ch)
    token = "".join(current).strip()
    if token:
        parts.append(token)
    return parts


def _split_top_level_assignment(token: str) -> tuple[str, str] | None:
    depth = 0
    for idx, ch in enumerate(token):
        if ch in "([{":
            depth += 1
            continue
        if ch in ")]}":
            depth = max(depth - 1, 0)
            continue
        if ch == "=" and depth == 0:
            return token[:idx].strip(), token[idx + 1 :].strip()
    return None


def _parse_param_map(raw: str) -> dict[str, str | None]:
    params: dict[str, str | None] = {}
    for part in _split_top_level_csv(raw):
        assignment = _split_top_level_assignment(part)
        if assignment is None:
            params[part.lower()] = None
        else:
            key, value = assignment
            params[key.lower()] = value.lower()
    return params


def _parse_route_token(token: str) -> tuple[str | None, str | None, dict[str, str | None] | None]:
    stripped = token.strip().lstrip("#")
    assignment = _split_top_level_assignment(stripped)
    if assignment is not None:
        left, right = assignment
        if right.startswith("(") and right.endswith(")"):
            return left.lower(), None, _parse_param_map(right[1:-1])
        return left.lower(), right.lower(), None

    if "(" in stripped and stripped.endswith(")"):
        head, tail = stripped.split("(", 1)
        return head.lower(), None, _parse_param_map(tail[:-1])

    return stripped.lower(), None, None


def _is_basis_like(token: str) -> bool:
    lowered = token.lower()
    return any(
        re.match(pattern, lowered)
        for pattern in (
            r"^(?:ma-)?def2[a-z0-9+\-]*$",
            r"^(?:ma-)?def2-[a-z0-9+\-]+$",
            r"^(?:aug-|jun-|jul-|may-|apr-)?cc-pv[a-z0-9]+z(?:-pp)?(?:fit|jkfit|ri)?$",
            r"^\d-\d+g(?:\([a-z0-9,+'\-]+\)|\*\*?|\+\+?g?(?:\([a-z0-9,+'\-]+\))?)?$",
            r"^(?:gen|genecp|sto-3g|3-21g|4-31g|6-21g|6-31g|6-311g|lanl2dz|lanl2mb|sdd|sddall|ugbs)$",
        )
    )


def _normalize_spin_and_method(token: str) -> tuple[str, str | None, str | None]:
    lowered = token.lower()
    spin_qualifier = None
    for prefix in ("ro", "ru", "rh", "uh", "u", "r"):
        if (
            lowered.startswith(prefix)
            and lowered not in _DFT_FUNCTIONALS
            and lowered not in _METHOD_FAMILY_MAP
        ):
            continue
    for prefix, tag in (("ro", "RO"), ("u", "U"), ("r", "R")):
        if (
            lowered.startswith(prefix)
            and lowered not in _DFT_FUNCTIONALS
            and lowered[len(prefix) :] in _METHOD_FAMILY_MAP
        ):
            spin_qualifier = tag
            lowered = lowered[len(prefix) :]
            break

    if lowered in _DFT_FUNCTIONALS:
        return token, spin_qualifier, "DFT"
    if lowered in _METHOD_FAMILY_MAP:
        return token, spin_qualifier, _METHOD_FAMILY_MAP[lowered]
    if lowered in {semi.lower() for semi in SEMI_EMPIRICAL_METHODS}:
        return token, spin_qualifier, "SEMI-EMPIRICAL"
    return token, spin_qualifier, None


def _dedupe(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))


def _build_model_chemistry(
    method_token: str, basis_token: str, auxiliary_basis: str | None = None
) -> GaussianModelChemistry:
    method_token_norm, spin_qualifier, method_family = _normalize_spin_and_method(method_token)
    return GaussianModelChemistry(
        method_token=method_token_norm,
        spin_qualifier=spin_qualifier,
        method_family=method_family,
        functional=method_token_norm if method_family in {"DFT", "HF", "FC"} else None,
        basis_set=basis_token,
        auxiliary_basis_set=auxiliary_basis,
        basis_family=_classify_basis_family(basis_token),
        basis_has_diffuse=(
            "+" in basis_token.lower()
            or basis_token.lower().startswith(("aug-", "jun-", "jul-", "may-", "apr-"))
        ),
        basis_polarization=_basis_polarization_markers(basis_token),
    )


def _int_or_none(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except Exception:
        return None


def _float_or_none(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


def _build_opt_options(
    option_map: GaussianRouteOptionMap | None, route_modifiers: list[str], job_types: list[str]
) -> GaussianOptOptions:
    if (
        option_map is None
        and "opt" not in job_types
        and not any(mod in route_modifiers for mod in ["ts", "readfc", "calcfc", "qst2", "qst3"])
    ):
        return GaussianOptOptions()
    params = dict(option_map.params) if option_map else {}
    keys = set(params)
    keys.update(route_modifiers)
    return GaussianOptOptions(
        enabled=option_map is not None
        or "opt" in job_types
        or any(mod in route_modifiers for mod in ["ts", "readfc", "calcfc", "qst2", "qst3"]),
        restart=("restart" in keys),
        transition_state=("ts" in keys),
        saddle_order=_int_or_none(params.get("saddle")),
        very_tight=("verytight" in keys),
        tight=("tight" in keys),
        loose=("loose" in keys),
        max_cycles=_int_or_none(params.get("maxcycles") or params.get("maxcycle")),
        max_step=_int_or_none(params.get("maxstep")),
        recalc_fc=_int_or_none(params.get("recalcfc")),
        read_fc=("readfc" in keys),
        read_cartesian_fc=("rcfc" in keys or "readcartesianfc" in keys),
        calc_fc=("calcfc" in keys),
        calc_all=("calcall" in keys),
        calc_hf_fc=("calchffc" in keys),
        has_modredundant=("modredundant" in route_modifiers),
        expert=(False if "noexpert" in keys else (True if "expert" in keys else None)),
        eigen_test=(False if "noeigentest" in keys else (True if "eigentest" in keys else None)),
        coordinate_system=(
            "cartesian"
            if "cartesian" in keys
            else (
                "internal" if "internal" in keys else ("redundant" if "redundant" in keys else None)
            )
        ),
        qst_mode="qst3" if "qst3" in keys else ("qst2" if "qst2" in keys else None),
        extra_options={
            key: value
            for key, value in params.items()
            if key
            not in {
                "ts",
                "restart",
                "saddle",
                "verytight",
                "tight",
                "loose",
                "readfc",
                "rcfc",
                "readcartesianfc",
                "calcfc",
                "calcall",
                "calchffc",
                "qst2",
                "qst3",
                "maxcycles",
                "maxcycle",
                "maxstep",
                "recalcfc",
                "expert",
                "noexpert",
                "eigentest",
                "noeigentest",
                "cartesian",
                "internal",
                "redundant",
            }
        },
    )


def _build_freq_options(option_map: GaussianRouteOptionMap | None) -> GaussianFreqOptions:
    if option_map is None:
        return GaussianFreqOptions()
    params = option_map.params
    keys = set(params)
    cphf_scalar = params.get("cphf")
    return GaussianFreqOptions(
        enabled=True,
        anharmonic=("anharmonic" in keys),
        read_anharm=("readanharm" in keys),
        projected=("projected" in keys),
        tprojected=("tprojected" in keys),
        hindered_rotor=("hinderedrotor" in keys),
        vibrot=("vibrot" in keys),
        polar=("polar" in keys),
        hpmodes=("hpmodes" in keys),
        read_fcht=("readfcht" in keys),
        read_isotopes=("readisotopes" in keys),
        select_normal_modes=("selectnormalmodes" in keys),
        save_normal_modes=("savenormalmodes" in keys),
        vcd=("vcd" in keys),
        roa=("roa" in keys),
        raman=("raman" in keys),
        no_raman=("noraman" in keys),
        cphf_rd_freq=(cphf_scalar == "rdfreq"),
        layer=params.get("layer"),
        atoms=params.get("atoms"),
        not_atoms=params.get("notatoms"),
        temperature=_float_or_none(params.get("temperature") or params.get("temp")),
        pressure=_float_or_none(params.get("pressure") or params.get("press")),
        extra_options={
            key: value
            for key, value in params.items()
            if key
            not in {
                "anharmonic",
                "readanharm",
                "projected",
                "tprojected",
                "hinderedrotor",
                "vibrot",
                "polar",
                "hpmodes",
                "readfcht",
                "readisotopes",
                "selectnormalmodes",
                "savenormalmodes",
                "vcd",
                "roa",
                "raman",
                "noraman",
                "cphf",
                "layer",
                "atoms",
                "notatoms",
                "temperature",
                "temp",
                "pressure",
                "press",
            }
        },
    )


def _build_td_options(option_map: GaussianRouteOptionMap | None) -> GaussianTDOptions:
    if option_map is None:
        return GaussianTDOptions()
    params = option_map.params
    keys = set(params)
    return GaussianTDOptions(
        enabled=True,
        nstates=_int_or_none(params.get("nstates")),
        root=_int_or_none(params.get("root")),
        singlets=("singlets" in keys),
        triplets=("triplets" in keys),
        tda=("tda" in keys),
        extra_options={
            key: value
            for key, value in params.items()
            if key not in {"nstates", "root", "singlets", "triplets", "tda"}
        },
    )


def _build_scrf_options(option_map: GaussianRouteOptionMap | None) -> GaussianSCRFOptions:
    if option_map is None:
        return GaussianSCRFOptions()
    params = option_map.params
    model = next(
        (key for key, value in params.items() if value is None and key not in {"solvent"}), None
    )
    option_keys = set(params)
    model_family = (
        "smd"
        if "smd" in option_keys
        else ("cpcm" if "cpcm" in option_keys else ("iefpcm" if "iefpcm" in option_keys else model))
    )
    return GaussianSCRFOptions(
        enabled=True,
        model=model,
        model_family=model_family,
        solvent=params.get("solvent"),
        read=("read" in option_keys),
        iefpcm=("iefpcm" in option_keys),
        cpcm=("cpcm" in option_keys),
        smd=("smd" in option_keys),
        extra_options={
            key: value
            for key, value in params.items()
            if key not in {"solvent", model, "read", "iefpcm", "cpcm", "smd"}
        },
    )


def _build_pop_options(option_map: GaussianRouteOptionMap | None) -> GaussianPopOptions:
    if option_map is None:
        return GaussianPopOptions()
    params = option_map.params
    keys = set(params)
    return GaussianPopOptions(
        enabled=True,
        none=("none" in keys),
        full=("full" in keys),
        nbo=("nbo" in keys),
        nbo_read=("nboread" in keys),
        nbo6_read=("nbo6read" in keys),
        nbo7_read=("nbo7read" in keys),
        hirshfeld=("hirshfeld" in keys),
        cm5=("cm5" in keys),
        mk=("mk" in keys),
        chelpg=("chelpg" in keys),
        orbitals=_int_or_none(params.get("orbitals")),
        read_radii=("readradii" in keys),
        read_at_radii=("readatradii" in keys),
        extra_options={
            key: value
            for key, value in params.items()
            if key
            not in {
                "none",
                "full",
                "nbo",
                "nboread",
                "nbo6read",
                "nbo7read",
                "hirshfeld",
                "cm5",
                "mk",
                "chelpg",
                "orbitals",
                "readradii",
                "readatradii",
            }
        },
    )


def _build_geom_options(option_map: GaussianRouteOptionMap | None) -> GaussianGeomOptions:
    if option_map is None:
        return GaussianGeomOptions()
    params = option_map.params
    keys = set(params)
    mode = (
        option_map.scalar_value
        if option_map.scalar_value in {"allcheck", "checkpoint", "check", "modify"}
        else (
            "allcheck"
            if "allcheck" in keys
            else (
                "checkpoint"
                if "checkpoint" in keys
                else ("check" if "check" in keys else ("modify" if "modify" in keys else None))
            )
        )
    )
    step = _int_or_none(params.get("step"))
    ngeom = _int_or_none(params.get("ngeom"))
    if step is not None and ngeom is None:
        ngeom = step + 1
    return GaussianGeomOptions(
        enabled=True,
        mode=mode,
        checkpoint=("checkpoint" in keys or option_map.scalar_value == "checkpoint"),
        allcheck=("allcheck" in keys or option_map.scalar_value == "allcheck"),
        check=("check" in keys or option_map.scalar_value == "check"),
        huge=("huge" in keys or option_map.scalar_value == "huge"),
        modify=("modify" in keys or option_map.scalar_value == "modify"),
        new_definition=("newdefinition" in keys),
        new_redundant=("newredundant" in keys),
        no_test=("notest" in keys),
        gic=("gic" in keys),
        add_gic=("addgic" in keys),
        read_all_gic=("readallgic" in keys),
        no_gic=("nogic" in keys),
        connectivity=("connectivity" in keys),
        mod_connectivity=("modconnectivity" in keys),
        gen_connectivity=("genconnectivity" in keys),
        zm_connectivity=("zmconnectivity" in keys),
        distance=("distance" in keys),
        no_distance=("nodistance" in keys),
        cangle=("cangle" in keys),
        angle=("angle" in keys),
        no_angle=("noangle" in keys),
        cdihedral=("cdihedral" in keys),
        dihedral=("dihedral" in keys),
        no_dihedral=("nodihedral" in keys),
        print_input_orient=("printinputorient" in keys),
        print=("print" in keys),
        step=step,
        ngeom=ngeom,
        extra_options={
            key: value
            for key, value in params.items()
            if key
            not in {
                "allcheck",
                "checkpoint",
                "check",
                "huge",
                "modify",
                "newdefinition",
                "newredundant",
                "notest",
                "gic",
                "addgic",
                "readallgic",
                "nogic",
                "connectivity",
                "modconnectivity",
                "genconnectivity",
                "zmconnectivity",
                "distance",
                "nodistance",
                "cangle",
                "angle",
                "noangle",
                "cdihedral",
                "dihedral",
                "nodihedral",
                "printinputorient",
                "print",
                "step",
                "ngeom",
            }
        },
    )


def parse_gaussian_route_semantic(route: str) -> GaussianRouteSemantic:
    normalized_route = _normalize_route_for_semantics(route)
    raw_tokens = _top_level_tokens(normalized_route)
    dieze_tag = None
    semantic = GaussianRouteSemantic(
        raw_route=route,
        normalized_route=normalized_route,
        dieze_tag=None,
    )

    tokens: list[GaussianRouteToken] = []
    job_types: list[str] = []
    route_modifiers: list[str] = []
    capabilities: list[str] = []
    unknown_tokens: list[str] = []
    solvation_model: str | None = None
    empirical_dispersion: str | None = None
    checkpoint_geometry_mode: str | None = None
    external_program: str | None = None
    option_maps: dict[str, GaussianRouteOptionMap] = {}
    model_chemistry = GaussianModelChemistry()

    for token_index, token in enumerate(raw_tokens):
        stripped = token.strip()
        upper = stripped.upper()
        normalized = stripped.lower().lstrip("#")
        if token_index == 0 and upper in {"#", "#N", "#P", "#T"}:
            dieze_tag = "#N" if upper == "#" else upper
            tokens.append(
                GaussianRouteToken(
                    raw=stripped,
                    normalized=normalized,
                    kind="dieze-tag",
                    key=normalized,
                )
            )
            continue

        token_kind = "unknown"
        token_key, token_scalar_value, token_param_map = _parse_route_token(stripped)

        if "//" in normalized and model_chemistry.method_token is None:
            high_raw, low_raw = stripped.split("//", 1)
            high_parts = [part for part in high_raw.split("/") if part]
            low_parts = [part for part in low_raw.split("/") if part]
            if len(high_parts) >= 2 and len(low_parts) >= 2:
                high_mc = _build_model_chemistry(
                    high_parts[0], high_parts[1], high_parts[2] if len(high_parts) >= 3 else None
                )
                low_mc = _build_model_chemistry(
                    low_parts[0], low_parts[1], low_parts[2] if len(low_parts) >= 3 else None
                )
                high_mc.low_level = low_mc
                model_chemistry = high_mc
                token_kind = "model-chemistry"
                job_types.extend(["opt", "sp"])
                capabilities.extend(["GeometryOptimization", "SinglePointEnergy"])

        if (
            "/" in normalized
            and not normalized.startswith("scrf")
            and not normalized.startswith("oniom")
            and "//" not in normalized
        ):
            parts = [part for part in normalized.split("/") if part]
            if len(parts) >= 2 and model_chemistry.method_token is None:
                left, right = parts[0], parts[1]
                if _is_basis_like(right) or right in {"gen", "genecp"}:
                    model_chemistry = _build_model_chemistry(
                        left, parts[1], parts[2] if len(parts) >= 3 else None
                    )
                    token_kind = "model-chemistry"

        if token_kind == "unknown":
            method_token, spin_qualifier, method_family = _normalize_spin_and_method(normalized)
            if model_chemistry.method_token is None and method_family is not None:
                model_chemistry.method_token = method_token
                model_chemistry.spin_qualifier = spin_qualifier
                model_chemistry.method_family = method_family
                model_chemistry.functional = (
                    method_token if method_family in {"DFT", "HF", "FC"} else None
                )
                token_kind = "method"
            elif model_chemistry.basis_set is None and _is_basis_like(normalized):
                model_chemistry.basis_set = stripped
                model_chemistry.basis_family = _classify_basis_family(stripped)
                model_chemistry.basis_has_diffuse = "+" in normalized or normalized.startswith(
                    ("aug-", "jun-", "jul-", "may-", "apr-")
                )
                model_chemistry.basis_polarization = _basis_polarization_markers(stripped)
                token_kind = "basis-set"

        key = normalized.split("=", 1)[0].split("(", 1)[0]
        if key in _JOB_TYPE_KEYS:
            job_types.append(key)
            capabilities.append(_JOB_TYPE_KEYS[key])
            token_kind = "job-type"

        if token_key is not None and (
            token_param_map is not None or token_scalar_value is not None
        ):
            option_maps[token_key] = GaussianRouteOptionMap(
                keyword=token_key,
                scalar_value=token_scalar_value,
                params=token_param_map or {},
            )

        if key in _ROUTE_MODIFIER_KEYS:
            route_modifiers.append(key)
            token_kind = "route-modifier"

        if key == "scrf":
            solvation_model = stripped
            capabilities.append("Solvation")
            token_kind = "solvation"

        if key in {"em", "empiricaldispersion"}:
            capabilities.append("Dispersion")
            route_modifiers.append(key)
            if token_scalar_value is not None:
                empirical_dispersion = token_scalar_value
            token_kind = "route-modifier"

        if key in {"geom"} and any(
            flag in normalized for flag in ("allcheck", "checkpoint", "check")
        ):
            capabilities.append("CheckpointGeometry")
            if token_scalar_value is not None:
                checkpoint_geometry_mode = token_scalar_value
            token_kind = "route-modifier"

        if key == "geom" and token_param_map is not None:
            if any(flag in token_param_map for flag in ("allcheck", "checkpoint", "check")):
                capabilities.append("CheckpointGeometry")
                if "allcheck" in token_param_map:
                    checkpoint_geometry_mode = "allcheck"
                elif "checkpoint" in token_param_map:
                    checkpoint_geometry_mode = "checkpoint"
                elif "check" in token_param_map:
                    checkpoint_geometry_mode = "check"
            token_kind = "route-modifier"

        if key == "external":
            if token_scalar_value is not None:
                external_program = token_scalar_value
            token_kind = "route-modifier"

        if token_kind == "unknown":
            unknown_tokens.append(stripped)

        tokens.append(
            GaussianRouteToken(
                raw=stripped,
                normalized=normalized,
                kind=token_kind,
                key=token_key,
                scalar_value=token_scalar_value,
                param_map=token_param_map,
            )
        )

    if model_chemistry.method_token is None and model_chemistry.basis_set is not None:
        semantic.diagnostics.messages.append("Basis set detected without a confident method token")

    if model_chemistry.method_family is None and model_chemistry.method_token is not None:
        semantic.diagnostics.messages.append("Method token detected but method family is unknown")

    semantic.dieze_tag = dieze_tag
    semantic.tokens = tokens
    semantic.model_chemistry = model_chemistry
    semantic.job_types = _dedupe(job_types)
    semantic.route_modifiers = _dedupe(route_modifiers)
    semantic.capabilities = _dedupe(capabilities)
    semantic.option_maps = option_maps
    semantic.opt_options = _build_opt_options(
        option_maps.get("opt"), semantic.route_modifiers, semantic.job_types
    )
    semantic.freq_options = _build_freq_options(option_maps.get("freq"))
    semantic.td_options = _build_td_options(option_maps.get("td"))
    semantic.scrf_options = _build_scrf_options(option_maps.get("scrf"))
    semantic.pop_options = _build_pop_options(option_maps.get("pop"))
    semantic.geom_options = _build_geom_options(option_maps.get("geom"))
    semantic.solvation_model = solvation_model
    semantic.empirical_dispersion = empirical_dispersion
    semantic.checkpoint_geometry_mode = checkpoint_geometry_mode
    semantic.external_program = external_program
    semantic.unknown_tokens = _dedupe(unknown_tokens)
    confidence = 0.2
    if semantic.model_chemistry.method_token:
        confidence += 0.3
    if semantic.model_chemistry.basis_set:
        confidence += 0.2
    if semantic.job_types:
        confidence += 0.1
    if semantic.capabilities:
        confidence += 0.1
    if semantic.unknown_tokens:
        confidence -= min(0.2, 0.05 * len(semantic.unknown_tokens))
        semantic.diagnostics.messages.append(
            f"Unclassified route tokens: {', '.join(semantic.unknown_tokens)}"
        )
    semantic.diagnostics.confidence = max(0.0, min(1.0, confidence))
    return semantic
