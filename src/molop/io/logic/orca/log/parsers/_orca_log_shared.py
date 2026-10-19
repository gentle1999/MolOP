"""ORCA explicit status and elapsed-time extraction.

Basis: output fixture termination, SCF and TOTAL RUN TIME anchors. See
``docs/en/reference/formats/orcaout.md`` for verified
official manuals and tutorials. For field mappings and MolOP conventions, see
``docs/en/reference/data_model.md``.
"""

from __future__ import annotations

from typing import Any

from pint.facets.plain import PlainQuantity

from molop.io.base_models.DataClasses import Status
from molop.io.logic.orca.log.parsers._orca_log_patterns import orca_log_patterns
from molop.unit import atom_ureg


def _as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def extract_orca_status(text: str, *, include_termination: bool = True) -> Status | None:
    """Use explicit output markers; the last observed status wins within scope.

    Evidence: corpus termination anchors and synthetic contradictory-SCF tests;
    see the reference index. Normal termination does not imply SCF convergence.
    """
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
    """Convert fixture TOTAL RUN TIME day/hour/minute/second/msec fields to seconds."""
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


def extract_orca_cpu_timing(
    text: str,
    wall_time: PlainQuantity | None,
    *,
    cpu_time: PlainQuantity | None = None,
    requested_num_cpu: int | None = None,
    inherited_num_cpu: int | None = None,
    keywords: str = "",
    serial_input: bool = False,
) -> dict[str, Any]:
    """Estimate missing CPU usage from elapsed time and the process count.

    This is a fallback estimate, not a measured utilization counter. Runtime
    MPI counts precede input requests; complete serial input defaults to one
    process. The multiplier remains available with the timing provenance.
    """
    num_cpu = inherited_num_cpu
    if matches := orca_log_patterns.RUNTIME_PROCESS_COUNT.find_matches(text):
        num_cpu = int(matches[-1].group("count"))
    elif num_cpu is None:
        if requested_num_cpu is not None:
            num_cpu = requested_num_cpu
        elif matches := orca_log_patterns.PAL_PROCESS_COUNT.find_matches(keywords):
            num_cpu = int(matches[-1].group("count"))
        elif serial_input:
            num_cpu = 1
    if num_cpu is not None and num_cpu <= 0:
        num_cpu = None
    source = "observed" if cpu_time is not None else None
    if cpu_time is None and wall_time is not None and num_cpu is not None:
        cpu_time = wall_time * num_cpu
        source = "estimated"
    return {"cpu_time": cpu_time, "cpu_time_source": source, "cpu_time_num_processes": num_cpu}
