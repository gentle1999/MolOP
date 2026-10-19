"""Independent CPU and elapsed durations from Gaussian timing records."""

from typing import Any

from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns
from molop.unit import atom_ureg


def extract_g16_timings(text: str, *, prefer_job_totals: bool = True) -> dict[str, Any]:
    """Keep CPU and wall time separate, without adding summaries to Link rows.

    Basis: Job cpu time / Elapsed time and Leave Link cpu / elap labels;
    see docs/en/reference/formats/g16log.md. File scopes prefer job totals;
    frame scopes prefer their local Link rows. Each quantity falls back
    independently when its preferred record is absent.
    """
    job_seconds: dict[str, list[float]] = {"cpu_time": [], "wall_time": []}
    for matched in g16_log_patterns.JOB_TIME.find_matches(text):
        field = "cpu_time" if matched.group("label") == "Job cpu time" else "wall_time"
        job_seconds[field].append(
            int(matched.group("days")) * 86400
            + int(matched.group("hours")) * 3600
            + int(matched.group("minutes")) * 60
            + float(matched.group("seconds"))
        )
    link_seconds: dict[str, list[float]] = {"cpu_time": [], "wall_time": []}
    for matched in g16_log_patterns.PROCEDURE_TIME.find_matches(text):
        link_seconds["cpu_time"].append(float(matched.group("cpu")))
        link_seconds["wall_time"].append(float(matched.group("elapsed")))
    primary, fallback = (
        (job_seconds, link_seconds) if prefer_job_totals else (link_seconds, job_seconds)
    )
    timings = {
        field: sum(values) * atom_ureg.second
        if (values := primary[field] or fallback[field])
        else None
        for field in ("cpu_time", "wall_time")
    }
    source: str | None = "observed" if timings["cpu_time"] is not None else None
    num_processes = None
    if source is not None:
        provenance = g16_log_patterns.CPU_TIME_SOURCE.find_matches(text)
        if any(row.group("source") == "estimated" for row in provenance):
            source = "estimated"
        counts = {int(row.group("count")) if row.group("count") else None for row in provenance}
        if (
            len(counts) == 1
            and None not in counts
            and len(provenance) == len(job_seconds["cpu_time"])
        ):
            num_processes = next(iter(counts))
    return timings | {"cpu_time_source": source, "cpu_time_num_processes": num_processes}
