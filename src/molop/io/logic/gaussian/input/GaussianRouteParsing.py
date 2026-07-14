from __future__ import annotations

from molop.io.base_models.SearchPattern import MolOPPattern
from molop.io.logic.gaussian.input.GaussianRoute import (
    SEMI_EMPIRICAL_METHODS,
    GaussianFreqOptions,
    GaussianGeomOptions,
    GaussianModelChemistry,
    GaussianOptOptions,
    GaussianPopOptions,
    GaussianRouteOptionMap,
    GaussianRouteSemantic,
    GaussianRouteToken,
    GaussianSCRFOptions,
    GaussianTDOptions,
)


class GaussianRoutePatterns:
    """Gaussian route grammar fragments."""

    POPLE_BASIS_PREFIX = MolOPPattern(content_pattern=r"^\d-\d+g")
    BASIS_POLARIZATION_GROUP = MolOPPattern(content_pattern=r"\((?P<markers>[^\)]+)\)")
    DIEZE_TAG_WITH_SPACING = MolOPPattern(content_pattern=r"^#\s+(?P<tag>[nNpPtT])\b")
    BASIS_LIKE = (
        MolOPPattern(content_pattern=r"^(?P<basis>(?:ma-)?def2[a-z0-9+\-]*)$"),
        MolOPPattern(content_pattern=r"^(?P<basis>(?:ma-)?def2-[a-z0-9+\-]+)$"),
        MolOPPattern(
            content_pattern=r"^(?P<basis>(?:aug-|jun-|jul-|may-|apr-)?"
            r"cc-pv[a-z0-9]+z(?:-pp)?(?:fit|jkfit|ri)?)$"
        ),
        MolOPPattern(
            content_pattern=r"^(?P<basis>\d-\d+g(?:\([a-z0-9,+'\-]+\)|"
            r"\*\*?|\+\+?g?(?:\([a-z0-9,+'\-]+\))?)?)$"
        ),
        MolOPPattern(
            content_pattern=r"^(?P<basis>gen|genecp|sto-3g|3-21g|4-31g|6-21g|"
            r"6-31g|6-311g|lanl2dz|lanl2mb|sdd|sddall|ugbs)$"
        ),
    )


gaussian_route_patterns = GaussianRoutePatterns()


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
    if gaussian_route_patterns.POPLE_BASIS_PREFIX.match(lowered):
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
    if matched := gaussian_route_patterns.BASIS_POLARIZATION_GROUP.search(lowered):
        markers.extend(part.strip() for part in matched.group("markers").split(",") if part.strip())
    return list(dict.fromkeys(markers))


def _normalize_route_for_semantics(route: str) -> str:
    stripped_lines = [line.strip() for line in route.splitlines() if line.strip()]
    normalized = " ".join(" ".join(stripped_lines).split()).strip()
    if matched := gaussian_route_patterns.DIEZE_TAG_WITH_SPACING.match(normalized):
        return f"#{matched.group('tag')}{normalized[matched.end() :]}"
    return normalized


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
    return any(pattern.match(lowered) is not None for pattern in gaussian_route_patterns.BASIS_LIKE)


def _normalize_spin_and_method(token: str) -> tuple[str, str | None, str | None]:
    lowered = token.lower()
    spin_qualifier = None
    semi_empirical_methods = {semi.lower() for semi in SEMI_EMPIRICAL_METHODS}
    for prefix, tag in (("ro", "RO"), ("u", "U"), ("r", "R")):
        unprefixed = lowered[len(prefix) :]
        if (
            lowered.startswith(prefix)
            and lowered not in _DFT_FUNCTIONALS
            and (
                unprefixed in _DFT_FUNCTIONALS
                or unprefixed in _METHOD_FAMILY_MAP
                or unprefixed in semi_empirical_methods
            )
        ):
            spin_qualifier = tag
            lowered = unprefixed
            break

    if lowered in _DFT_FUNCTIONALS:
        return token, spin_qualifier, "DFT"
    if lowered in _METHOD_FAMILY_MAP:
        return token, spin_qualifier, _METHOD_FAMILY_MAP[lowered]
    if lowered in semi_empirical_methods:
        return token, spin_qualifier, "SEMI-EMPIRICAL"
    return token, spin_qualifier, None


def _method_family_allows_functional(method_family: str | None) -> bool:
    if method_family is None:
        return False
    normalized = method_family.upper().replace("_", "-").replace(" ", "-")
    return normalized in {"DFT", "DOUBLE-HYBRID", "DOUBLE-HYBRID-DFT"}


def _without_spin_prefix(method_token: str, spin_qualifier: str | None) -> str:
    if spin_qualifier == "RO":
        return method_token[2:]
    if spin_qualifier in {"R", "U"}:
        return method_token[1:]
    return method_token


def _functional_from_method_token(
    method_token: str,
    method_family: str | None,
    spin_qualifier: str | None,
) -> str | None:
    if _method_family_allows_functional(method_family):
        return _without_spin_prefix(method_token, spin_qualifier)
    return None


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
        functional=_functional_from_method_token(method_token_norm, method_family, spin_qualifier),
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
                model_chemistry.functional = _functional_from_method_token(
                    method_token, method_family, spin_qualifier
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
