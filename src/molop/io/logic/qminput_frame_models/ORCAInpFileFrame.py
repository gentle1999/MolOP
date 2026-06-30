from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from typing import Any, Literal, Protocol, cast

import numpy as np
from pydantic import BaseModel, Field, model_validator
from rdkit import Chem
from typing_extensions import Self

from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.DataClasses import (
    ActiveSpace,
    AtomInInternalCoords,
    CoordinateContainer,
    CoordinateParameters,
    ExcitedStateRequest,
    ExplicitSolventRequest,
    InternalCoords,
    MultireferenceRequest,
    MultireferenceStateBlock,
    QMBasisSet,
    QMModelChemistry,
    QMTaskRequest,
)
from molop.io.base_models.Mixins import DiskStorageMixin, MemoryStorageMixin
from molop.unit import atom_ureg


class _HasExplicitSolventRequests(Protocol):
    explicit_solvent_requests: list[ExplicitSolventRequest]


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
    doalpha: bool = Field(default=False, description="Whether alpha-channel only calculation is requested")
    rootwise: bool = Field(default=False, description="Whether rootwise solving is requested")
    do_dbfilter: bool = Field(default=False, description="Whether STEOM doubly excited filtering is enabled")
    do_store_steom: bool = Field(default=False, description="Whether STEOM intermediates are stored")
    do_simple_dens: bool = Field(default=False, description="Whether STEOM simple density is disabled")
    add_l2_term: bool = Field(default=False, description="Whether DLPNO STEOM L2 term is enabled")
    do_full_semiclassical: bool = Field(
        default=False, description="Whether full semiclassical treatment is enabled"
    )
    do_higher_moments: bool = Field(
        default=False, description="Whether higher moments are enabled"
    )
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


def _parse_float(token: str) -> float | None:
    cleaned = token.strip().strip('"').strip("'").rstrip(",")
    if "," in cleaned:
        cleaned = cleaned.split(",", 1)[0]
    if not cleaned:
        return None
    try:
        return float(cleaned.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _parse_int(token: str) -> int | None:
    cleaned = token.strip().strip('"').strip("'").rstrip(",")
    if not cleaned:
        return None
    try:
        return int(cleaned)
    except ValueError:
        try:
            as_float = float(cleaned)
        except ValueError:
            return None
        if not as_float.is_integer():
            return None
        return int(as_float)


def _derive_keyword_semantics(keyword_text: str) -> tuple[str, str, str, str, str]:
    tokens = re.split(r"\s+", keyword_text.strip())
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


def _parse_output_print_settings(blocks: list[ORCABlock]) -> list[ORCAOutputPrintSetting]:
    settings: list[ORCAOutputPrintSetting] = []
    pattern = re.compile(
        r"^Print\s*\[\s*(?P<target>[^\]]+?)\s*\]\s*(?:=|\s+)\s*(?P<value>.+)$",
        re.I,
    )
    for block in blocks:
        if block.name.lower() != "output":
            continue
        for line in block.lines:
            matched = pattern.match(line.text.strip())
            if matched is None:
                continue
            settings.append(
                ORCAOutputPrintSetting(
                    target=matched.group("target").strip(),
                    value=matched.group("value").strip(),
                )
            )
    return settings


def _parse_request_num_cpu(blocks: list[ORCABlock]) -> int | None:
    for block in blocks:
        if block.name.lower() != "pal":
            continue
        for line in block.lines:
            tokens = line.text.split()
            if len(tokens) >= 2 and tokens[0].lower() == "nprocs":
                try:
                    return int(float(tokens[1]))
                except ValueError:
                    continue
    return None


def _parse_request_memory(blocks: list[ORCABlock]):
    for block in blocks:
        if block.name.lower() != "maxcore":
            continue
        values: list[str] = []
        header_tokens = block.raw_header.split()
        if len(header_tokens) > 1:
            values.extend(header_tokens[1:])
        for line in block.lines:
            values.extend(line.text.split())
        for value in values:
            try:
                return float(value) * atom_ureg.megabyte
            except ValueError:
                continue
    return None


def _strip_inline_comment(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return ""
    if "#" in line:
        return line.split("#", 1)[0]
    return line


def _normalize_orca_option_key(key: str) -> str:
    return re.sub(r"\s+", "", key).strip().lower()


def _parse_orca_block_option_line(line: str) -> tuple[str, str | None] | None:
    content = _strip_inline_comment(line).strip()
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
    header_content = _strip_inline_comment(block.raw_header).strip()
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
        parsed = _parse_float(cast(str, options[key]))
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
        start = _parse_float(parts[0])
        end = _parse_float(parts[1])
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
        matched = re.match(r"^(?P<model>[A-Za-z][A-Za-z0-9-]*)\((?P<solvent>[^()]+)\)$", token)
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
        fixsolute=not _parse_bool_option(options, "nofixsolute") if "nofixsolute" in options else _parse_bool_option(options, "fixsolute") or True,
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

        matched = re.match(r"^(?P<model>[A-Za-z][A-Za-z0-9-]*)\((?P<solvent>[^()]+)\)$", token)
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
        for key in ("epsilon", "refrac", "rsolv", "rmin", "pmin", "fepstype", "xfeps", "surfacetype", "scale_gauss", "cpcmccm", "draco", "draco_charges"):
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


def _detect_excited_state_family(keyword_text: str, blocks: Sequence[ORCABlock]) -> tuple[str | None, str | None]:
    tokens = [token for token in re.split(r"\s+", keyword_text.strip()) if token]
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

        multiplicity = _parse_int(tokens[1]) if len(tokens) > 1 else None
        irrep = tokens[2] if len(tokens) > 2 else None
        nroots = _parse_int(tokens[4]) if len(tokens) > 4 and tokens[3].lower() == "nroots" else None
        raw_lines = [stripped]
        refs: str | None = None
        excitations: str | None = None
        extra_options: dict[str, str | None] = {}

        inline_tail = tokens[3:] if len(tokens) > 3 else []
        if inline_tail:
            inline_text = " ".join(inline_tail)
            matched_exc = re.search(r"\bexcitations\s+(\S+)", inline_text, flags=re.I)
            if matched_exc:
                excitations = matched_exc.group(1)
            matched_refs = re.search(r"\brefs\s+(.+?)\s+end\b", inline_text, flags=re.I)
            if matched_refs:
                refs = matched_refs.group(1).strip()

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
                refs_inline = re.sub(r"\s+end\s*$", "", refs_value, flags=re.I).strip()
                if refs_inline:
                    refs_parts.append(refs_inline)
                if not re.search(r"\bend\s*$", current_stripped, flags=re.I):
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
                    nroots = _parse_int(value or "")
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
    keyword_tokens = [token for token in re.split(r"\s+", keyword_text.strip()) if token]
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
    canonical_ci_type = _normalize_multi_reference_method(mrci_options.get("citype") or "") or mrci_options.get(
        "citype"
    )

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
    return [token for token in re.split(r"\s+", keyword_text.strip()) if token]


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
    options: dict[str, Any] = {"has_mixed_basis": any(atom.basis_overrides for atom in geometry or [])}
    if solvation_options:
        options["solvation"] = solvation_options
    return QMModelChemistry(
        method_family=_orca_method_family(
            method, functional, multi_reference_enabled
        ),
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
    match = re.search(r"\bcas\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)", refs, flags=re.I)
    if match is None:
        return ActiveSpace(raw=refs)
    return ActiveSpace(electrons=int(match.group(1)), orbitals=int(match.group(2)), raw=refs)


def _active_space_from_options(options: Mapping[str, str | None]) -> ActiveSpace | None:
    electrons = _parse_int(options.get("nel") or "")
    orbitals = _parse_int(options.get("norb") or "")
    roots = _parse_int(options.get("nroots") or "")
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


def _populate_common_orca_qm_containers(frame: BaseQMInputFrame, mixin: ORCAInpFileFrameMixin) -> None:
    frame.model_chemistry = _build_orca_model_chemistry(
        keywords=frame.keywords,
        method=frame.method,
        functional=frame.functional,
        basis_set=frame.basis_set,
        auxiliary_basis_set=mixin.auxiliary_basis_set,
        dispersion_correction=mixin.dispersion_correction,
        blocks=mixin.blocks,
        geometry=mixin.geometry,
        multi_reference_enabled=mixin.multi_reference_semantic.enabled,
    )
    frame.task_requests = _build_orca_task_requests(
        frame.keywords,
        mixin.blocks,
        mixin.excited_state_semantic,
        mixin.multi_reference_semantic,
    )
    frame.excited_state_requests = _build_orca_excited_state_requests(
        mixin.excited_state_semantic
    )
    frame.multireference_requests = _build_orca_multireference_requests(
        mixin.multi_reference_semantic
    )
    if hasattr(frame, "explicit_solvent_requests"):
        cast(_HasExplicitSolventRequests, frame).explicit_solvent_requests = (
            _build_orca_explicit_solvent_requests(mixin.explicit_solvent_semantic)
        )
    frame.backfill_common_qm_containers_from_legacy()
    frame.project_common_qm_fields()


class ORCAInpFileFrameMixin:
    comment_lines: list[ORCACommentLine] = Field(default_factory=list)
    keyword_lines: list[ORCAKeywordLine] = Field(default_factory=list)
    blocks: list[ORCABlock] = Field(default_factory=list)
    geometry: ORCAGeometry | None = Field(default=None)
    excited_state_semantic: ORCAExcitedStateSemantic = Field(
        default_factory=ORCAExcitedStateSemantic,
        description="Structured excited-state task semantics",
    )
    multi_reference_semantic: ORCAMultiReferenceSemantic = Field(
        default_factory=ORCAMultiReferenceSemantic,
        description="Structured multi-reference task semantics",
    )
    explicit_solvent_semantic: ORCAExplicitSolventSemantic = Field(
        default_factory=ORCAExplicitSolventSemantic,
        description="Structured explicit-solvent semantics",
    )
    explicit_solvent_requests: list[ExplicitSolventRequest] = Field(default_factory=list)
    trailing_lines: list[str] = Field(default_factory=list)
    auxiliary_basis_set: str = Field(
        default="",
        description="ORCA auxiliary basis keyword derived from the input keywords",
    )
    dispersion_correction: str = Field(
        default="",
        description="ORCA dispersion correction keyword derived from the input keywords",
    )
    has_mixed_basis: bool = Field(
        default=False,
        description="Whether atom-level basis overrides were found",
    )
    output_print_settings: list[ORCAOutputPrintSetting] = Field(default_factory=list)

    def _render(self, **kwargs) -> str:
        _ = kwargs
        raise NotImplementedError(
            f"{self.__class__.__name__} does not support ORCA input rendering yet."
        )

    @model_validator(mode="after")
    def set_orca_properties(self) -> Self:
        typed_self = cast(BaseQMInputFrame, self)
        typed_self.qm_software = "ORCA"
        typed_self.qm_software_version = "Any"

        keyword_text = "\n".join(line.text for line in self.keyword_lines if line.text.strip())
        typed_self.keywords = keyword_text
        (
            method,
            functional,
            basis_set,
            auxiliary_basis_set,
            dispersion_correction,
        ) = _derive_keyword_semantics(keyword_text)
        typed_self.method = method
        typed_self.functional = functional
        typed_self.basis_set = basis_set
        self.auxiliary_basis_set = auxiliary_basis_set
        self.dispersion_correction = dispersion_correction
        self.excited_state_semantic = _build_excited_state_semantic(keyword_text, self.blocks)
        if (
            self.excited_state_semantic.enabled
            and self.excited_state_semantic.reference_method is None
            and method
            and self.excited_state_semantic.family in {"CIS", "ROCIS"}
        ):
            self.excited_state_semantic.reference_method = method
        if self.excited_state_semantic.enabled:
            typed_self.method = _resolve_method_for_excited_state(
                typed_self.method,
                self.excited_state_semantic.family,
                self.excited_state_semantic.reference_method,
                self.blocks,
            )
        self.multi_reference_semantic = _build_multi_reference_semantic(keyword_text, self.blocks)
        if self.multi_reference_semantic.enabled and self.multi_reference_semantic.method:
            typed_self.method = self.multi_reference_semantic.method
        self.explicit_solvent_semantic = _parse_orca_explicit_solvent_semantics(
            keyword_text, self.blocks
        )

        resources = []
        for block in self.blocks:
            if block.raw_text:
                resources.append(block.raw_text)
                continue
            header = block.raw_header or f"%{block.name}"
            body = block.body_text()
            resources.append(header if not body else f"{header}\n{body}")
        typed_self.resources_raw = "\n".join(resources)
        typed_self.request_num_cpu = _parse_request_num_cpu(self.blocks)
        typed_self.request_memory = _parse_request_memory(self.blocks)
        self.output_print_settings = _parse_output_print_settings(self.blocks)

        if self.geometry is None:
            _populate_common_orca_qm_containers(typed_self, self)
            return self

        typed_self.charge = self.geometry.charge
        typed_self.multiplicity = self.geometry.multiplicity
        self.has_mixed_basis = any(atom.basis_overrides for atom in self.geometry)
        if self.geometry.internal_coords is not None:
            typed_self.atoms = [
                Chem.GetPeriodicTable().GetAtomicNumber(symbol)
                for symbol in self.geometry.get_symbols()
            ]
            typed_self.coords = self.geometry.internal_coords.to_cartesian_coords()
        elif self.geometry:
            pt = Chem.GetPeriodicTable()
            atoms: list[int] = []
            coords: list[list[float]] = []
            for atom in self.geometry.real_atoms():
                if atom.x is None or atom.y is None or atom.z is None:
                    continue
                atomic_number = atom.atomic_number or pt.GetAtomicNumber(atom.symbol)
                if atomic_number <= 0:
                    continue
                atoms.append(atomic_number)
                coords.append([atom.x, atom.y, atom.z])
            if atoms:
                typed_self.atoms = atoms
                typed_self.coords = np.asarray(coords, dtype=float) * atom_ureg.angstrom
        _populate_common_orca_qm_containers(typed_self, self)
        return self


class ORCAInpFileFrameMemory(
    MemoryStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameMemory"]
): ...


class ORCAInpFileFrameDisk(
    DiskStorageMixin, ORCAInpFileFrameMixin, BaseQMInputFrame["ORCAInpFileFrameDisk"]
): ...
