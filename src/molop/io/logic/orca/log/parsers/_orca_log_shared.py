from __future__ import annotations

from typing import Any

from molop.io.base_models.DataClasses import Status
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns
from molop.unit import atom_ureg


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def extract_orca_status(text: str) -> Status | None:
    if "****ORCA TERMINATED NORMALLY****" in text:
        return Status(normal_terminated=True, scf_converged=True)
    if "ORCA finished by error termination" in text or "ORCA TERMINATED ABNORMALLY" in text:
        return Status(normal_terminated=False, scf_converged=False)
    return None


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
