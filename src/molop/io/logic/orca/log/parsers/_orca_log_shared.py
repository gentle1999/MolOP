from __future__ import annotations

from typing import Any

from molop.io.base_models.DataClasses import Status
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns
from molop.unit import atom_ureg


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def extract_orca_status(text: str, *, include_termination: bool = True) -> Status | None:
    status_dict: dict[str, bool] = {}
    if include_termination:
        normal_pos = text.rfind("****ORCA TERMINATED NORMALLY****")
        abnormal_pos = max(
            text.rfind("ORCA finished by error termination"),
            text.rfind("ORCA TERMINATED ABNORMALLY"),
        )
        if normal_pos >= 0 or abnormal_pos >= 0:
            status_dict["normal_terminated"] = normal_pos > abnormal_pos

    scf_evidence = [
        (matched.start(), True) for matched in orca_log_patterns.SCF_CONVERGED.find_matches(text)
    ]
    scf_evidence.extend(
        (matched.start(), False) for matched in orca_log_patterns.SCF_FAILED.find_matches(text)
    )
    if scf_evidence:
        status_dict["scf_converged"] = max(scf_evidence, key=lambda evidence: evidence[0])[1]

    return Status.model_validate(status_dict) if status_dict else None


def extract_orca_running_time(text: str) -> Any | None:
    if matches := orca_log_patterns.RUN_TIME.find_matches(text):
        matched = matches[-1]
        seconds = (
            int(matched.group("days")) * 86400
            + int(matched.group("hours")) * 3600
            + int(matched.group("minutes")) * 60
            + int(matched.group("seconds"))
            + int(matched.group("msec")) / 1000
        )
        return seconds * atom_ureg.second
    return None
