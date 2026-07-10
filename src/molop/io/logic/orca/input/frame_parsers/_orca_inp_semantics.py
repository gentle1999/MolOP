from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, cast

from molop.io.base_models.DataClasses import (
    ActiveSpace,
    ExcitedStateRequest,
    ExplicitSolventRequest,
    MultireferenceRequest,
    MultireferenceStateBlock,
    QMBasisSet,
    QMModelChemistry,
    QMTaskRequest,
)
from molop.io.logic.orca.common import (
    ORCABlock,
    ORCAExcitedStateSemantic,
    ORCAExplicitSolventSemantic,
    ORCAGeometry,
    ORCAMultiReferenceNewBlock,
    ORCAMultiReferenceSemantic,
)
from molop.io.logic.orca.input.frame_parsers._orca_inp_patterns import orca_inp_patterns
from molop.io.logic.orca.input.frame_parsers._orca_inp_tokens import (
    parse_orca_float,
    parse_orca_int,
    strip_orca_inline_comment,
)


_KNOWN_BASIS_PREFIXES = ("def2-", "ma-", "cc-", "aug-", "sto-", "pc-", "ano-", "minix")
_KNOWN_BASIS_TOKENS = {"sv", "svp", "sv(p)", "tzvp", "tzvpp", "qzvp", "qzvpp"}
_KNOWN_FUNCTIONALS = {
    "b3lyp",
    "bp",
    "bhlyp",
    "b2plyp",
    "dsd-blyp",
    "dsd-pbep86",
    "dsdpbep86",
    "pbe0",
    "pbe",
    "opbe",
    "bp86",
    "blyp",
    "tpss",
    "b97",
    "wb97x",
    "wb97m-v",
    "wb97m(2)",
    "m06",
}
_KNOWN_WAVEFUNCTION_METHODS = {
    "hf",
    "rhf",
    "uhf",
    "mp2",
    "ri-mp2",
    "dlnpo-mp2",
    "dlpno-mp2",
    "ccsd",
    "ccsd(t)",
    "qcisd",
    "qcisd(t)",
    "mracpf",
    "sorci",
}
_KNOWN_MULTI_REFERENCE_METHODS = {
    "cepa1": "CEPA1",
    "cepa2": "CEPA2",
    "cepa3": "CEPA3",
    "mracpf": "MRACPF",
    "mracpf2": "MRACPF2",
    "mracpf2a": "MRACPF2a",
    "mraqcc": "MRAQCC",
    "mrcepa0": "MRCEPA_0",
    "mrcepa_0": "MRCEPA_0",
    "mrcepar": "MRCEPA_R",
    "mrcepa_r": "MRCEPA_R",
    "mrci": "MRCI",
    "mrci+q": "MRCI+Q",
    "mrddci1": "MRDDCI1",
    "mrddci2": "MRDDCI2",
    "mrddci3": "MRDDCI3",
    "mrmp2": "MRMP2",
    "mrmp3": "MRMP3",
    "mrre2": "MRRE2",
    "mrre3": "MRRE3",
    "mrre4": "MRRE4",
    "sorci": "SORCI",
    "sorcp": "SORCP",
}
_KNOWN_REFERENCE_METHODS = {
    "casscf",
    "hf",
    "rhf",
    "uhf",
    "rohf",
    "rks",
    "roks",
}
_EXCITED_STATE_FAMILIES = {
    "adc2": "ADC2",
    "bt-pno-eom-ccsd": "BT-PNO-EOM-CCSD",
    "cis": "CIS",
    "eom-ccsd": "EOM-CCSD",
    "ih-fsmr-ccsd": "IH-FSMR-CCSD",
    "ip-eom-ccsd": "IP-EOM-CCSD",
    "mcrpa": "MCRPA",
    "rocis": "ROCIS",
    "steom-ccsd": "STEOM-CCSD",
    "steom-dlpno-ccsd": "STEOM-DLPNO-CCSD",
    "tddft": "TDDFT",
    "td-dft": "TDDFT",
}
_CC_EXCITED_STATE_FAMILIES = {
    "ADC2",
    "BT-PNO-EOM-CCSD",
    "EOM-CCSD",
    "IH-FSMR-CCSD",
    "IP-EOM-CCSD",
    "STEOM-CCSD",
    "STEOM-DLPNO-CCSD",
}
_KNOWN_DISPERSION_CORRECTIONS = {
    "d3": "D3",
    "d3bj": "D3BJ",
    "d3bjm": "D3BJM",
    "d3bjabc": "D3BJABC",
    "d3zero": "D3ZERO",
    "d3zerom": "D3ZEROM",
    "d4": "D4",
    "nl": "NL",
    "vv10": "VV10",
}
_KNOWN_AUXILIARY_BASIS_SUFFIXES = ("/c", "/j", "/jk")
_KNOWN_SOLVATION_MODELS = {
    "cpcm": "CPCM",
    "cpcmc": "CPCM",
    "smd": "SMD",
    "cosmors": "COSMORS",
    "alpb": "ALPB",
    "ddcosmo": "ddCOSMO",
    "cpcm-x": "CPCM-X",
}


def _is_auxiliary_basis_token(lowered: str) -> bool:
    return any(lowered.endswith(suffix) for suffix in _KNOWN_AUXILIARY_BASIS_SUFFIXES)


def _normalize_multi_reference_method(token: str) -> str | None:
    lowered = token.strip().lower()
    if not lowered:
        return None
    candidates: list[str] = [lowered, lowered.replace("-", "").replace("_", "")]
    for prefix in ("ri-", "f12-ri-", "f12-", "cim-"):
        if not lowered.startswith(prefix):
            continue
        stripped = lowered[len(prefix) :]
        candidates.extend((stripped, stripped.replace("-", "").replace("_", "")))

    for candidate in candidates:
        if normalized := _KNOWN_MULTI_REFERENCE_METHODS.get(candidate):
            return normalized
    return None


def _derive_keyword_semantics(keyword_text: str) -> tuple[str, str, str, str, str]:
    tokens = keyword_text.split()
    method = ""
    functional = ""
    basis_set = ""
    auxiliary_basis_set = ""
    dispersion_correction = ""
    for token in tokens:
        cleaned = token.strip()
        lowered = cleaned.lower()
        if not cleaned:
            continue
        if not dispersion_correction:
            dispersion_correction = _KNOWN_DISPERSION_CORRECTIONS.get(lowered, "")
        if not auxiliary_basis_set and _is_auxiliary_basis_token(lowered):
            auxiliary_basis_set = cleaned
            continue
        if not basis_set and (
            lowered.startswith(_KNOWN_BASIS_PREFIXES)
            or lowered in _KNOWN_BASIS_TOKENS
            or lowered in {"6-31g", "6-31g*", "6-31g(d)"}
        ):
            basis_set = cleaned
        functional_candidate = lowered.removeprefix("dlpno-")
        if not functional and (
            functional_candidate in _KNOWN_FUNCTIONALS
            or any(functional_candidate.startswith(prefix) for prefix in _KNOWN_FUNCTIONALS)
        ):
            functional = cleaned.upper() if cleaned.islower() else cleaned
            continue
        multi_reference_method = _normalize_multi_reference_method(cleaned)
        if multi_reference_method:
            method = multi_reference_method
    if functional:
        method = "DFT"
        if dispersion_correction:
            functional = f"{functional}-{dispersion_correction}"
    elif not method:
        for token in tokens:
            token_lower = token.lower()
            compact_lower = token_lower.replace("-", "").replace("_", "")
            if token_lower in _KNOWN_WAVEFUNCTION_METHODS:
                method = token.upper()
                break
            if "dlpnoccsd" in compact_lower:
                method = "DLPNO-CCSD(T)" if "(t)" in token_lower else "DLPNO-CCSD"
                break
            if "mp2" in compact_lower:
                method = "MP2"
                break
    return method, functional, basis_set, auxiliary_basis_set, dispersion_correction


def _normalize_orca_option_key(key: str) -> str:
    return "".join(key.split()).lower()


def _parse_orca_block_option_line(line: str) -> tuple[str, str | None] | None:
    content = strip_orca_inline_comment(line).strip()
    if not content:
        return None
    lowered = content.lower()
    if lowered in {"end", "step_end"}:
        return None

    if lowered.startswith("%"):
        content = content[1:].strip()
        if not content:
            return None
        parts = content.split(maxsplit=1)
        if len(parts) <= 1:
            return None
        content = parts[1].strip()
        if not content:
            return None

    if content.lower().endswith(" end"):
        content = content[: -len(" end")].rstrip()
    if not content:
        return None

    key: str
    value: str | None
    if "=" in content:
        key, value = content.split("=", 1)
    else:
        parts = content.split(maxsplit=1)
        key = parts[0]
        value = parts[1] if len(parts) > 1 else None

    normalized_key = _normalize_orca_option_key(key)
    if not normalized_key:
        return None

    if value is None:
        return normalized_key, None
    cleaned_value = value.strip().strip('"').strip("'")
    return normalized_key, cleaned_value or None


def _parse_orca_block_option_map(block: ORCABlock) -> dict[str, str | None]:
    option_map: dict[str, str | None] = {}
    header_content = strip_orca_inline_comment(block.raw_header).strip()
    if header_content.startswith("%"):
        header_content = header_content[1:].strip()
        if header_content:
            parts = header_content.split(maxsplit=1)
            header_content = parts[1].strip() if len(parts) > 1 else ""

    for line in [header_content, *(line.text for line in block.lines)]:
        if not line:
            continue
        parsed = _parse_orca_block_option_line(line)
        if parsed is None:
            continue
        key, value = parsed
        option_map[key] = value
    return option_map


def _parse_bool_option(options: Mapping[str, str | None], *keys: str) -> bool:
    for key in keys:
        if key not in options:
            continue
        value = options[key]
        if value is None:
            return True
        lowered = value.strip().lower()
        if lowered in {"1", "true", "t", "yes", "y", "on"}:
            return True
        return lowered not in {"0", "false", "f", "no", "n", "off"}
    return False


def _parse_int_option(options: Mapping[str, str | None], *keys: str) -> int | None:
    for key in keys:
        if key not in options or options[key] is None:
            continue
        try:
            return int(float(cast(str, options[key])))
        except ValueError:
            continue
    return None


def _parse_float_option(options: Mapping[str, str | None], *keys: str) -> float | None:
    for key in keys:
        if key not in options or options[key] is None:
            continue
        parsed = parse_orca_float(cast(str, options[key]))
        if parsed is not None:
            return parsed
    return None


def _parse_float_range_option(
    options: Mapping[str, str | None], *keys: str
) -> tuple[float, float] | None:
    for key in keys:
        value = options.get(key)
        if value is None:
            continue
        parts = [part.strip() for part in value.split(",")]
        if len(parts) < 2:
            continue
        start = parse_orca_float(parts[0])
        end = parse_orca_float(parts[1])
        if start is None or end is None:
            continue
        return float(start), float(end)
    return None


def _parse_orca_explicit_solvent_semantics(
    keyword_text: str, blocks: Sequence[ORCABlock]
) -> ORCAExplicitSolventSemantic:
    block = next((candidate for candidate in blocks if candidate.name.lower() == "solvator"), None)
    if block is None:
        return ORCAExplicitSolventSemantic()

    options = _parse_orca_block_option_map(block)
    solvent_model = None
    solvent = None
    for token in _keyword_tokens(keyword_text):
        matched = orca_inp_patterns.SOLVATION_TOKEN.match(token)
        if matched is None:
            continue
        if matched.group("model").lower() in _KNOWN_SOLVATION_MODELS:
            solvent_model = _KNOWN_SOLVATION_MODELS[matched.group("model").lower()]
            solvent = _normalize_solvation_solvent_name(matched.group("solvent"))
            break

    if options.get("solvent"):
        solvent = _normalize_solvation_solvent_name(cast(str, options["solvent"]))
    solvent_file = cast(str, options["solventfile"]) if options.get("solventfile") else None

    return ORCAExplicitSolventSemantic(
        enabled=True,
        solvent_model=solvent_model,
        solvent=solvent,
        solvent_file=solvent_file,
        nsolv=_parse_int_option(options, "nsolv"),
        cluster_mode=options.get("clustermode"),
        droplet=_parse_bool_option(options, "droplet"),
        radius=_parse_float_option(options, "radius"),
        fixsolute=not _parse_bool_option(options, "nofixsolute")
        if "nofixsolute" in options
        else _parse_bool_option(options, "fixsolute") or True,
        vacuumsearch=_parse_bool_option(options, "vacuumsearch"),
        randomsolv=_parse_bool_option(options, "randomsolv"),
        printlevel=options.get("printlevel"),
        source_blocks=[block.name.lower()],
        extra_options={
            key: value
            for key, value in options.items()
            if key
            not in {
                "solvent",
                "solventfile",
                "nsolv",
                "clustermode",
                "droplet",
                "radius",
                "fixsolute",
                "nofixsolute",
                "vacuumsearch",
                "randomsolv",
                "printlevel",
            }
        },
    )


def _normalize_solvation_solvent_name(value: str) -> str:
    return value.strip().strip('"').strip("'").rstrip(",").lower()


def _parse_orca_solvation_semantics(
    keyword_text: str, blocks: Sequence[ORCABlock]
) -> tuple[str | None, str | None, dict[str, Any]]:
    solvation_model: str | None = None
    solvent: str | None = None
    options: dict[str, Any] = {}

    for token in _keyword_tokens(keyword_text):
        token_upper = token.strip().upper()
        if token_upper == "NOCPCM":
            return None, None, {"enabled": False}

        matched = orca_inp_patterns.SOLVATION_TOKEN.match(token)
        if matched is None:
            continue

        model_key = matched.group("model").lower()
        solvent_name = _normalize_solvation_solvent_name(matched.group("solvent"))
        if model_key not in _KNOWN_SOLVATION_MODELS:
            continue

        solvation_model = solvation_model or _KNOWN_SOLVATION_MODELS[model_key]
        solvent = solvent or solvent_name
        if model_key == "cpcmc":
            options["epsilon_function"] = "COSMO"
            options["solvation_variant"] = "CPCMC"

    for block in blocks:
        if block.name.lower() != "cpcm":
            continue

        block_options = _parse_orca_block_option_map(block)
        if _parse_bool_option(block_options, "smd"):
            solvation_model = solvation_model or "SMD"
            options["smd"] = True
        if block_options.get("smdsolvent"):
            solvent = solvent or _normalize_solvation_solvent_name(
                cast(str, block_options["smdsolvent"])
            )
        if any(
            key in block_options
            for key in ("epsilon", "refrac", "rsolv", "rmin", "pmin", "fepstype", "xfeps")
        ):
            solvation_model = solvation_model or "CPCM"
        for key in (
            "epsilon",
            "refrac",
            "rsolv",
            "rmin",
            "pmin",
            "fepstype",
            "xfeps",
            "surfacetype",
            "scale_gauss",
            "cpcmccm",
            "draco",
            "draco_charges",
        ):
            if key in block_options and block_options[key] is not None:
                options[key] = block_options[key]
        fepstype = block_options.get("fepstype")
        if isinstance(fepstype, str) and fepstype.strip().lower() == "cosmo":
            options["epsilon_function"] = "COSMO"
            solvation_model = solvation_model or "CPCM"
        if _parse_bool_option(block_options, "draco"):
            options["draco"] = True
            if solvation_model is None:
                solvation_model = "CPCM"

    return solvation_model, solvent, options


def _collect_excited_state_blocks(blocks: Sequence[ORCABlock]) -> list[ORCABlock]:
    relevant_names = {"tddft", "cis", "rocis", "casscf", "mcrpa", "mdci", "scf"}
    return [block for block in blocks if block.name.lower() in relevant_names]


def _detect_excited_state_family(
    keyword_text: str, blocks: Sequence[ORCABlock]
) -> tuple[str | None, str | None]:
    tokens = keyword_text.split()
    family: str | None = None
    reference_method: str | None = None

    for idx, token in enumerate(tokens):
        lowered = token.lower()
        if lowered in _EXCITED_STATE_FAMILIES:
            family = _EXCITED_STATE_FAMILIES[lowered]
            if idx > 0:
                previous = tokens[idx - 1].strip().lower()
                if previous in _KNOWN_REFERENCE_METHODS:
                    reference_method = tokens[idx - 1].strip().upper()
            break

    if family is None:
        block_names = {block.name.lower() for block in blocks}
        if "mcrpa" in block_names:
            family = "MCRPA"
        elif "tddft" in block_names:
            family = "TDDFT"
        elif "cis" in block_names:
            family = "CIS"
        elif "rocis" in block_names:
            family = "ROCIS"

    if reference_method is None:
        for block in blocks:
            if block.name.lower() != "scf":
                continue
            options = _parse_orca_block_option_map(block)
            candidate = options.get("hftyp")
            if candidate:
                reference_method = candidate.strip().upper()
                break
            for key in ("hf", "rhf", "uhf", "rohf", "rks", "roks"):
                if key in options:
                    reference_method = key.upper()
                    break
            if reference_method:
                break

    if reference_method is None and family == "MCRPA":
        reference_method = "CASSCF"

    return family, reference_method


def _resolve_method_for_excited_state(
    base_method: str, family: str | None, reference_method: str | None, blocks: Sequence[ORCABlock]
) -> str:
    if family in _CC_EXCITED_STATE_FAMILIES or family in {"CIS", "ROCIS"}:
        return family
    if family == "MCRPA":
        if any(block.name.lower() == "casscf" for block in blocks):
            return "CASSCF"
        return base_method or reference_method or "MCRPA"
    if family == "TDDFT" and not base_method:
        return reference_method or "TDDFT"
    return base_method


def _resolve_multi_reference_method(
    keyword_method: str, ci_type: str | None, blocks: Sequence[ORCABlock]
) -> str:
    ci_type_method = _normalize_multi_reference_method(ci_type or "")
    if ci_type_method:
        return ci_type_method

    keyword_method_normalized = _normalize_multi_reference_method(keyword_method)
    if keyword_method_normalized:
        return keyword_method_normalized

    if any(block.name.lower() == "mrci" for block in blocks):
        return "MRCI"

    return keyword_method


def _parse_multi_reference_new_blocks(block: ORCABlock) -> list[ORCAMultiReferenceNewBlock]:
    parsed: list[ORCAMultiReferenceNewBlock] = []
    lines = [line.text for line in block.lines]
    idx = 0

    while idx < len(lines):
        stripped = lines[idx].strip()
        if not stripped or stripped.startswith("#"):
            idx += 1
            continue
        tokens = stripped.split()
        if tokens[0].lower() != "newblock":
            idx += 1
            continue

        multiplicity = parse_orca_int(tokens[1]) if len(tokens) > 1 else None
        irrep = tokens[2] if len(tokens) > 2 else None
        nroots = (
            parse_orca_int(tokens[4]) if len(tokens) > 4 and tokens[3].lower() == "nroots" else None
        )
        raw_lines = [stripped]
        refs: str | None = None
        excitations: str | None = None
        extra_options: dict[str, str | None] = {}

        inline_tail = tokens[3:] if len(tokens) > 3 else []
        if inline_tail:
            inline_text = " ".join(inline_tail)
            matched_exc = orca_inp_patterns.MRCI_INLINE_EXCITATIONS.search(inline_text)
            if matched_exc:
                excitations = matched_exc.group("excitations")
            matched_refs = orca_inp_patterns.MRCI_INLINE_REFS.search(inline_text)
            if matched_refs:
                refs = matched_refs.group("refs").strip()

        idx += 1
        refs_parts: list[str] = [refs] if refs else []
        in_refs = False

        while idx < len(lines):
            current = lines[idx]
            current_stripped = current.strip()
            current_lower = current_stripped.lower()
            if not current_stripped or current_stripped.startswith("#"):
                raw_lines.append(current_stripped)
                idx += 1
                continue
            if current_lower.startswith("newblock"):
                break

            raw_lines.append(current_stripped)
            if current_lower == "end":
                idx += 1
                break
            if current_lower.startswith("refs"):
                refs_value = current_stripped[4:].strip()
                refs_end = orca_inp_patterns.TRAILING_END.search(refs_value)
                refs_inline = (
                    refs_value[: refs_end.start()] if refs_end is not None else refs_value
                ).strip()
                if refs_inline:
                    refs_parts.append(refs_inline)
                if orca_inp_patterns.END_TOKEN_AT_LINE_END.search(current_stripped) is None:
                    in_refs = True
                idx += 1
                continue
            if in_refs:
                if current_lower == "end":
                    in_refs = False
                else:
                    refs_parts.append(current_stripped)
                idx += 1
                continue

            parsed_option = _parse_orca_block_option_line(current_stripped)
            if parsed_option is not None:
                key, value = parsed_option
                if key == "nroots":
                    nroots = parse_orca_int(value or "")
                elif key == "excitations":
                    excitations = value
                else:
                    extra_options[key] = value
            idx += 1
        else:
            idx = len(lines)

        parsed.append(
            ORCAMultiReferenceNewBlock(
                multiplicity=multiplicity,
                irrep=irrep,
                nroots=nroots,
                excitations=excitations,
                refs=" ".join(part for part in refs_parts if part).strip() or None,
                raw_lines=raw_lines,
                extra_options=extra_options,
            )
        )

    return parsed


def _build_multi_reference_semantic(
    keyword_text: str, blocks: Sequence[ORCABlock]
) -> ORCAMultiReferenceSemantic:
    keyword_tokens = keyword_text.split()
    has_keyword_multi_reference_method = any(
        _normalize_multi_reference_method(token) for token in keyword_tokens
    )
    has_mrci_block = any(block.name.lower() == "mrci" for block in blocks)
    if not has_mrci_block and not has_keyword_multi_reference_method:
        return ORCAMultiReferenceSemantic()

    relevant_names = {"casscf", "mrci", "method", "base", "paras"}
    relevant_blocks = [block for block in blocks if block.name.lower() in relevant_names]
    block_options = {
        block.name.lower(): _parse_orca_block_option_map(block) for block in relevant_blocks
    }
    mrci_block = next((block for block in relevant_blocks if block.name.lower() == "mrci"), None)
    mrci_options = block_options.get("mrci", {})
    canonical_ci_type = _normalize_multi_reference_method(
        mrci_options.get("citype") or ""
    ) or mrci_options.get("citype")

    reference_method: str | None = None
    if "casscf" in block_options:
        reference_method = "CASSCF"
    else:
        for token in keyword_tokens:
            lowered = token.strip().lower()
            if lowered in _KNOWN_REFERENCE_METHODS:
                reference_method = token.strip().upper()
                break

    keyword_multi_reference_method = next(
        (
            token.strip()
            for token in keyword_tokens
            if _normalize_multi_reference_method(token.strip())
        ),
        "",
    )
    method = _resolve_multi_reference_method(
        keyword_multi_reference_method, canonical_ci_type, relevant_blocks
    )

    known_keys = {
        "citype",
        "ewin",
        "tsel",
        "tpre",
        "tnat",
        "etol",
        "rtol",
        "solver",
        "intmode",
        "useivos",
        "allsingles",
        "doddcimp2",
        "donatorbs",
        "eunselopt",
        "davidsonopt",
        "partitioning",
        "fopt",
    }
    return ORCAMultiReferenceSemantic(
        enabled=True,
        method=method or None,
        ci_type=canonical_ci_type,
        reference_method=reference_method,
        source_blocks=[block.name.lower() for block in relevant_blocks],
        block_options=block_options,
        new_blocks=_parse_multi_reference_new_blocks(mrci_block) if mrci_block is not None else [],
        ewin=_parse_float_range_option(mrci_options, "ewin"),
        tsel=_parse_float_option(mrci_options, "tsel"),
        tpre=_parse_float_option(mrci_options, "tpre"),
        tnat=_parse_float_option(mrci_options, "tnat"),
        etol=_parse_float_option(mrci_options, "etol"),
        rtol=_parse_float_option(mrci_options, "rtol"),
        solver=mrci_options.get("solver"),
        int_mode=mrci_options.get("intmode"),
        use_ivos=_parse_bool_option(mrci_options, "useivos"),
        all_singles=_parse_bool_option(mrci_options, "allsingles"),
        do_ddcimp2=_parse_bool_option(mrci_options, "doddcimp2"),
        do_nat_orbs=_parse_bool_option(mrci_options, "donatorbs"),
        eunsel_opt=mrci_options.get("eunselopt"),
        davidson_opt=mrci_options.get("davidsonopt"),
        partitioning=mrci_options.get("partitioning"),
        fopt=mrci_options.get("fopt"),
        extra_options={key: value for key, value in mrci_options.items() if key not in known_keys},
    )


def _build_excited_state_semantic(
    keyword_text: str, blocks: Sequence[ORCABlock]
) -> ORCAExcitedStateSemantic:
    relevant_blocks = _collect_excited_state_blocks(blocks)
    block_options = {
        block.name.lower(): _parse_orca_block_option_map(block) for block in relevant_blocks
    }
    family, reference_method = _detect_excited_state_family(keyword_text, relevant_blocks)

    source_block_name: str | None = None
    for candidate in (
        family.lower() if family else None,
        "mdci" if family in _CC_EXCITED_STATE_FAMILIES else None,
        "mcrpa" if family == "MCRPA" else None,
        "tddft" if family == "TDDFT" else None,
        "cis" if family == "CIS" else None,
        "rocis" if family == "ROCIS" else None,
        "casscf" if family == "CASSCF" else None,
    ):
        if candidate and candidate in block_options:
            source_block_name = candidate
            break

    if source_block_name is None and "mdci" in block_options:
        source_block_name = "mdci"

    if family is None and source_block_name is None:
        return ORCAExcitedStateSemantic()

    source_options = block_options.get(source_block_name, {}) if source_block_name else {}
    family = family or ("MDCI" if source_block_name == "mdci" else None)
    semantic = ORCAExcitedStateSemantic(
        enabled=family is not None or source_block_name is not None,
        family=family,
        reference_method=reference_method,
        source_blocks=[block.name.lower() for block in relevant_blocks],
        block_options=block_options,
        nroots=_parse_int_option(source_options, "nroots"),
        iroot=_parse_int_option(source_options, "iroot"),
        jroot=_parse_int_option(source_options, "jroot"),
        followiroot=_parse_bool_option(source_options, "followiroot"),
        triplets=_parse_bool_option(source_options, "triplets"),
        sf=_parse_bool_option(source_options, "sf"),
        nacme=_parse_bool_option(source_options, "nacme"),
        etf=_parse_bool_option(source_options, "etf"),
        dosoc=_parse_bool_option(source_options, "dosoc"),
        doalpha=_parse_bool_option(source_options, "doalpha"),
        rootwise=_parse_bool_option(source_options, "rootwise", "dorootwise"),
        do_dbfilter=_parse_bool_option(source_options, "dodbfilter", "dodbfilter"),
        do_store_steom=_parse_bool_option(source_options, "dostoresteom"),
        do_simple_dens=_parse_bool_option(source_options, "dosimpledens"),
        add_l2_term=_parse_bool_option(source_options, "addl2term"),
        do_full_semiclassical=_parse_bool_option(source_options, "dofullsemiclassical"),
        do_higher_moments=_parse_bool_option(source_options, "dohighermoments"),
        firkeepfirstref=_parse_bool_option(source_options, "firkeepfirstref"),
        extra_options={
            key: value
            for key, value in source_options.items()
            if key
            not in {
                "nroots",
                "iroot",
                "jroot",
                "followiroot",
                "triplets",
                "sf",
                "nacme",
                "etf",
                "dosoc",
                "doalpha",
                "rootwise",
                "dorootwise",
                "dodbfilter",
                "dostoresteom",
                "dosimpledens",
                "addl2term",
                "dofullsemiclassical",
                "dohighermoments",
                "firkeepfirstref",
            }
        },
    )
    return semantic


def _keyword_tokens(keyword_text: str) -> list[str]:
    return keyword_text.split()


def _orca_method_family(method: str, functional: str, multi_reference_enabled: bool) -> str | None:
    if functional:
        return "DFT"
    normalized = method.upper()
    if not normalized:
        return None
    if multi_reference_enabled:
        return "multi-reference"
    if normalized in {"HF", "RHF", "UHF", "ROHF"}:
        return "HF"
    if "MP2" in normalized:
        return "MP2"
    if "CCSD" in normalized:
        return "CCSD"
    return method


def _collect_orca_basis_sets(
    basis_set: str,
    auxiliary_basis_set: str,
    geometry: ORCAGeometry | None,
) -> list[QMBasisSet]:
    basis_sets: list[QMBasisSet] = []
    if basis_set:
        basis_sets.append(QMBasisSet(name=basis_set, role="orbital", scope="global", raw=basis_set))
    if auxiliary_basis_set:
        basis_sets.append(
            QMBasisSet(
                name=auxiliary_basis_set,
                role="auxiliary",
                scope="global",
                raw=auxiliary_basis_set,
            )
        )
    if geometry is None:
        return basis_sets

    real_atom_index = 0
    for atom in geometry:
        if atom.is_dummy or atom.is_ghost:
            continue
        for override in atom.basis_overrides:
            role = "auxiliary" if override.kind == "newauxgto" else "orbital"
            basis_sets.append(
                QMBasisSet(
                    name=override.basis_set or "",
                    role=role,
                    scope="atom",
                    atom_indices=[real_atom_index],
                    element_symbols=[atom.symbol],
                    raw=" ".join(override.tokens),
                    options={"kind": override.kind},
                )
            )
        real_atom_index += 1
    return basis_sets


def _build_orca_model_chemistry(
    *,
    keywords: str,
    method: str,
    functional: str,
    basis_set: str,
    auxiliary_basis_set: str,
    dispersion_correction: str,
    blocks: Sequence[ORCABlock],
    geometry: ORCAGeometry | None,
    multi_reference_enabled: bool,
) -> QMModelChemistry:
    solvation_model, solvent, solvation_options = _parse_orca_solvation_semantics(keywords, blocks)
    options: dict[str, Any] = {
        "has_mixed_basis": any(atom.basis_overrides for atom in geometry or [])
    }
    if solvation_options:
        options["solvation"] = solvation_options
    return QMModelChemistry(
        method_family=_orca_method_family(method, functional, multi_reference_enabled),
        method=method or None,
        functional=functional or None,
        basis_set=basis_set or None,
        auxiliary_basis_set=auxiliary_basis_set or None,
        basis_sets=_collect_orca_basis_sets(basis_set, auxiliary_basis_set, geometry),
        dispersion_correction=dispersion_correction or None,
        solvation_model=solvation_model,
        solvent=solvent,
        raw_keywords=keywords,
        options=options,
    )


def _build_orca_task_requests(
    keyword_text: str,
    blocks: Sequence[ORCABlock],
    excited_state_semantic: ORCAExcitedStateSemantic,
    multi_reference_semantic: ORCAMultiReferenceSemantic,
) -> list[QMTaskRequest]:
    tokens = _keyword_tokens(keyword_text)
    lowered_tokens = [token.lower() for token in tokens]
    block_names = [block.name.lower() for block in blocks]
    tasks: list[QMTaskRequest] = []

    def matching_tokens(*needles: str) -> list[str]:
        return [
            token
            for token, lowered in zip(tokens, lowered_tokens, strict=True)
            if any(needle in lowered for needle in needles)
        ]

    scan_requested = (
        "geom" in block_names
        and any("scan" in line.text.lower() for block in blocks for line in block.lines)
    ) or any("scan" in token for token in lowered_tokens)
    ts_requested = any("ts" in token for token in lowered_tokens) or any(
        name in block_names for name in ("neb", "irc")
    )

    opt_tokens = matching_tokens("opt", "neb-ts", "scants", "surfcrossopt")
    if opt_tokens:
        tasks.append(
            QMTaskRequest(
                task_type="opt",
                derivative_order=1,
                transition_state=ts_requested,
                scan=scan_requested,
                source_keywords=opt_tokens,
                source_blocks=[name for name in block_names if name in {"geom", "neb"}],
            )
        )

    freq_tokens = matching_tokens("freq", "anfreq", "numfreq", "surfcrossnumfreq")
    if freq_tokens:
        tasks.append(
            QMTaskRequest(
                task_type="freq",
                derivative_order=2,
                source_keywords=freq_tokens,
                source_blocks=[name for name in block_names if name == "freq"],
                options={
                    "analytic": any(token.lower() == "anfreq" for token in freq_tokens),
                    "numeric": any("numfreq" in token.lower() for token in freq_tokens),
                },
            )
        )

    if any(token == "irc" for token in lowered_tokens):
        tasks.append(QMTaskRequest(task_type="irc", source_keywords=matching_tokens("irc")))

    if "neb" in block_names and not any(task.task_type == "neb" for task in tasks):
        tasks.append(
            QMTaskRequest(
                task_type="neb",
                transition_state=ts_requested,
                source_keywords=matching_tokens("neb"),
                source_blocks=["neb"],
            )
        )

    if excited_state_semantic.enabled:
        tasks.append(
            QMTaskRequest(
                task_type="excited_state",
                target_state=excited_state_semantic.iroot,
                properties=[
                    name
                    for name, enabled in {
                        "nacme": excited_state_semantic.nacme,
                        "etf": excited_state_semantic.etf,
                        "soc": excited_state_semantic.dosoc,
                        "spin_flip": excited_state_semantic.sf,
                    }.items()
                    if enabled
                ],
                source_keywords=matching_tokens(
                    "tddft",
                    "cis",
                    "rocis",
                    "mcrpa",
                    "eom",
                    "adc",
                    "steom",
                ),
                source_blocks=excited_state_semantic.source_blocks,
            )
        )

    if multi_reference_semantic.enabled:
        tasks.append(
            QMTaskRequest(
                task_type="multi_reference",
                source_keywords=matching_tokens("mr", "sorci", "cepa"),
                source_blocks=multi_reference_semantic.source_blocks,
            )
        )

    gradient_tokens = matching_tokens("engrad", "grad")
    if gradient_tokens and not tasks:
        tasks.append(
            QMTaskRequest(
                task_type="sp",
                derivative_order=1,
                properties=["gradient"],
                source_keywords=gradient_tokens,
            )
        )

    if not tasks:
        tasks.append(QMTaskRequest(task_type="sp", derivative_order=0))
    return tasks


def _build_orca_excited_state_requests(
    semantic: ORCAExcitedStateSemantic,
) -> list[ExcitedStateRequest]:
    if not semantic.enabled:
        return []
    properties = [
        name
        for name, enabled in {
            "nacme": semantic.nacme,
            "etf": semantic.etf,
            "soc": semantic.dosoc,
            "alpha_only": semantic.doalpha,
            "rootwise": semantic.rootwise,
            "dbfilter": semantic.do_dbfilter,
            "store_steom": semantic.do_store_steom,
            "simple_density": semantic.do_simple_dens,
            "l2_term": semantic.add_l2_term,
            "full_semiclassical": semantic.do_full_semiclassical,
            "higher_moments": semantic.do_higher_moments,
        }.items()
        if enabled
    ]
    return [
        ExcitedStateRequest(
            enabled=True,
            family=semantic.family,
            reference_method=semantic.reference_method,
            nroots=semantic.nroots,
            root=semantic.iroot,
            secondary_root=semantic.jroot,
            roots=[root for root in (semantic.iroot, semantic.jroot) if root is not None],
            triplets=semantic.triplets,
            spin_flip=semantic.sf,
            follow_root=semantic.followiroot,
            properties=properties,
            source_blocks=semantic.source_blocks,
            options={
                key: value
                for key, value in semantic.model_dump().items()
                if key
                not in {
                    "enabled",
                    "family",
                    "reference_method",
                    "source_blocks",
                    "nroots",
                    "iroot",
                    "jroot",
                    "triplets",
                    "sf",
                    "followiroot",
                }
            },
        )
    ]


def _active_space_from_refs(refs: str | None) -> ActiveSpace | None:
    if refs is None:
        return None
    match = orca_inp_patterns.CAS_REFS.search(refs)
    if match is None:
        return ActiveSpace(raw=refs)
    return ActiveSpace(
        electrons=int(match.group("electrons")),
        orbitals=int(match.group("orbitals")),
        raw=refs,
    )


def _active_space_from_options(options: Mapping[str, str | None]) -> ActiveSpace | None:
    electrons = parse_orca_int(options.get("nel") or "")
    orbitals = parse_orca_int(options.get("norb") or "")
    roots = parse_orca_int(options.get("nroots") or "")
    if electrons is None and orbitals is None and roots is None:
        return None
    return ActiveSpace(
        electrons=electrons,
        orbitals=orbitals,
        roots=roots,
        options={
            key: value
            for key, value in options.items()
            if key in {"mult", "irootmult", "actorbs", "closedorbs"}
        },
    )


def _build_orca_multireference_requests(
    semantic: ORCAMultiReferenceSemantic,
) -> list[MultireferenceRequest]:
    if not semantic.enabled:
        return []

    casscf_active_space = _active_space_from_options(semantic.block_options.get("casscf", {}))
    state_blocks = [
        MultireferenceStateBlock(
            multiplicity=block.multiplicity,
            irrep=block.irrep,
            nroots=block.nroots,
            excitations=block.excitations,
            refs=block.refs,
            active_space=_active_space_from_refs(block.refs),
            raw_lines=block.raw_lines,
            options=block.extra_options,
        )
        for block in semantic.new_blocks
    ]
    thresholds = {
        key: value
        for key, value in {
            "tsel": semantic.tsel,
            "tpre": semantic.tpre,
            "tnat": semantic.tnat,
            "etol": semantic.etol,
            "rtol": semantic.rtol,
        }.items()
        if value is not None
    }
    corrections = [
        correction
        for correction, enabled in {
            "DDCI-MP2": semantic.do_ddcimp2,
            f"Davidson:{semantic.davidson_opt}": semantic.davidson_opt is not None,
            f"EUnsel:{semantic.eunsel_opt}": semantic.eunsel_opt is not None,
        }.items()
        if enabled
    ]
    options: dict[str, Any] = {
        "block_options": semantic.block_options,
        "extra_options": semantic.extra_options,
    }
    if semantic.ewin is not None:
        options["ewin"] = semantic.ewin
    for key in ("solver", "int_mode", "partitioning", "fopt"):
        if value := getattr(semantic, key):
            options[key] = value

    return [
        MultireferenceRequest(
            enabled=True,
            method=semantic.method,
            reference_method=semantic.reference_method,
            ci_type=semantic.ci_type,
            active_space=casscf_active_space,
            state_blocks=state_blocks,
            thresholds=thresholds,
            corrections=corrections,
            source_blocks=semantic.source_blocks,
            options=options,
        )
    ]


def _build_orca_explicit_solvent_requests(
    semantic: ORCAExplicitSolventSemantic,
) -> list[ExplicitSolventRequest]:
    if not semantic.enabled:
        return []
    return [
        ExplicitSolventRequest(
            enabled=True,
            solvent_model=semantic.solvent_model,
            solvent=semantic.solvent,
            solvent_file=semantic.solvent_file,
            nsolv=semantic.nsolv,
            cluster_mode=semantic.cluster_mode,
            droplet=semantic.droplet,
            radius=semantic.radius,
            fixsolute=semantic.fixsolute,
            vacuumsearch=semantic.vacuumsearch,
            randomsolv=semantic.randomsolv,
            printlevel=semantic.printlevel,
            source_blocks=semantic.source_blocks,
            options={"extra_options": semantic.extra_options},
        )
    ]


def build_orca_input_semantic_payload(
    keyword_text: str,
    blocks: Sequence[ORCABlock],
    geometry: ORCAGeometry | None,
) -> dict[str, Any]:
    method, functional, basis_set, auxiliary_basis_set, dispersion_correction = (
        _derive_keyword_semantics(keyword_text)
    )
    excited_state_semantic = _build_excited_state_semantic(keyword_text, blocks)
    if (
        excited_state_semantic.enabled
        and excited_state_semantic.reference_method is None
        and method
        and excited_state_semantic.family in {"CIS", "ROCIS"}
    ):
        excited_state_semantic.reference_method = method
    if excited_state_semantic.enabled:
        method = _resolve_method_for_excited_state(
            method,
            excited_state_semantic.family,
            excited_state_semantic.reference_method,
            blocks,
        )

    multi_reference_semantic = _build_multi_reference_semantic(keyword_text, blocks)
    if multi_reference_semantic.enabled and multi_reference_semantic.method:
        method = multi_reference_semantic.method

    explicit_solvent_semantic = _parse_orca_explicit_solvent_semantics(keyword_text, blocks)
    model_chemistry = _build_orca_model_chemistry(
        keywords=keyword_text,
        method=method,
        functional=functional,
        basis_set=basis_set,
        auxiliary_basis_set=auxiliary_basis_set,
        dispersion_correction=dispersion_correction,
        blocks=blocks,
        geometry=geometry,
        multi_reference_enabled=multi_reference_semantic.enabled,
    )

    return {
        "qm_software": "ORCA",
        "qm_software_version": "Any",
        "keywords": keyword_text,
        "method": method,
        "functional": functional,
        "basis_set": basis_set,
        "auxiliary_basis_set": auxiliary_basis_set,
        "dispersion_correction": dispersion_correction,
        "excited_state_semantic": excited_state_semantic,
        "multi_reference_semantic": multi_reference_semantic,
        "explicit_solvent_semantic": explicit_solvent_semantic,
        "model_chemistry": model_chemistry,
        "task_requests": _build_orca_task_requests(
            keyword_text,
            blocks,
            excited_state_semantic,
            multi_reference_semantic,
        ),
        "excited_state_requests": _build_orca_excited_state_requests(excited_state_semantic),
        "multireference_requests": _build_orca_multireference_requests(multi_reference_semantic),
        "explicit_solvent_requests": _build_orca_explicit_solvent_requests(
            explicit_solvent_semantic
        ),
    }
