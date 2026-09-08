from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from molop.config import available_cpu_count


ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = ROOT / "tests" / "test_files" / "g16log"
# Discover the complete Gaussian fixture corpus instead of selecting a few
# representative files.  This keeps the stress test aligned with additions
# to the test data without requiring another hard-coded list.
LOG_PATHS = tuple(path.relative_to(ROOT).as_posix() for path in sorted(LOG_DIR.glob("*.log")))
EXPECTED_FILE_COUNT = len(LOG_PATHS)


_STRESS_WORKER = textwrap.dedent(
    """
    import json
    import sys

    from molop import AutoParser, molopconfig

    molopconfig.show_progress_bar = False
    assert molopconfig.prewarm_topologies is False
    n_jobs = int(sys.argv[1])
    paths = %r
    batch = AutoParser(paths, n_jobs=n_jobs, only_last_frame=True)
    summary = batch.to_summary_df(frame="all", n_jobs=n_jobs)
    rendered = batch.format_transform("smi", frame="all", n_jobs=n_jobs)
    frame_count = sum(len(file_model) for file_model in batch)
    assert frame_count == len(paths)
    assert len(summary) == frame_count
    assert len(rendered) == len(paths)
    # A source log can be valid while its coordinates do not yield a chemical
    # graph (for example, a metal-containing or otherwise incomplete result).
    # Keep those empty strings in the result and only require complete output
    # coverage here; the subprocess exit status is the crash sentinel.
    assert all(value is not None for value in rendered.values())
    nonempty_rendered = sum(bool(value) for value in rendered.values())
    print(json.dumps({
        "files": len(batch),
        "frames": frame_count,
        "summary_rows": len(summary),
        "rendered_files": len(rendered),
        "nonempty_rendered": nonempty_rendered,
        "empty_rendered": len(rendered) - nonempty_rendered,
        "n_jobs": n_jobs,
    }))
    """
) % (list(LOG_PATHS),)


@pytest.mark.stress
def test_gaussian_log_lazy_reconstruction_stress() -> None:
    """Repeat the spawn/loky log workflow in isolated processes."""

    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[variable] = "1"
    source_path = str(ROOT / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        path for path in (source_path, environment.get("PYTHONPATH")) if path
    )
    assert EXPECTED_FILE_COUNT > 0
    # Keep the stress workload bounded while using the available runner
    # capacity up to four processes. A one-CPU runner falls back to serial mode.
    n_jobs = max(1, min(4, available_cpu_count(), EXPECTED_FILE_COUNT))

    observed_counts: list[tuple[int, int]] = []
    for iteration in range(3):
        completed = subprocess.run(
            [sys.executable, "-c", _STRESS_WORKER, str(n_jobs)],
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        assert completed.returncode == 0, (
            f"stress iteration {iteration + 1} failed with exit code "
            f"{completed.returncode}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
        output_lines = [line for line in completed.stdout.splitlines() if line.strip()]
        assert output_lines, f"stress iteration {iteration + 1} produced no result"
        result = json.loads(output_lines[-1])
        assert result["files"] == EXPECTED_FILE_COUNT
        assert result["frames"] == EXPECTED_FILE_COUNT
        assert result["summary_rows"] == EXPECTED_FILE_COUNT
        assert result["rendered_files"] == EXPECTED_FILE_COUNT
        assert result["nonempty_rendered"] + result["empty_rendered"] == EXPECTED_FILE_COUNT
        assert result["n_jobs"] == n_jobs
        observed_counts.append((result["nonempty_rendered"], result["empty_rendered"]))

    # Repeated isolated workers must produce the same success/failure split;
    # otherwise the stress run is exposing nondeterministic native behavior.
    assert len(set(observed_counts)) == 1
