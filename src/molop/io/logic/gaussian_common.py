from __future__ import annotations

from typing import Any, Protocol

from molop.io.base_models.DataClasses import (
    ExcitedStateRequest,
    QMModelChemistry,
    QMTaskRequest,
)
from molop.io.logic.gaussian_route_models import GaussianRouteSemantic


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
