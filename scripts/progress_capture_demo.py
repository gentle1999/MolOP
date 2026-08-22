"""Run the MolOP progress-capture examples without chemistry fixtures.

Run from the repository root with::

    uv run --frozen python scripts/progress_capture_demo.py
"""

from __future__ import annotations

import argparse
import json
import time
from collections.abc import Sequence
from pathlib import Path

from molop.utils.progressbar import (
    AdaptiveProgress,
    ProgressEvent,
    ProgressRecorder,
    parallel_map,
    progress_jsonl_sink,
    progress_listener,
)


def _work(value: int) -> int:
    time.sleep(0.02)
    return value * 2


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--items", type=int, default=5, help="Number of demo items.")
    parser.add_argument(
        "--jsonl",
        type=Path,
        default=Path("progress_capture_demo.jsonl"),
        help="JSONL output path for the sink example.",
    )
    parser.add_argument(
        "--show-progress",
        action="store_true",
        help="Render terminal progress bars in addition to emitting events.",
    )
    args = parser.parse_args(argv)
    if args.items < 1:
        parser.error("--items must be positive")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    values = list(range(args.items))
    disable = not args.show_progress

    callback_events: list[ProgressEvent] = []
    list(
        AdaptiveProgress(
            values,
            desc="Direct callback",
            total=len(values),
            disable=disable,
            progress_callback=callback_events.append,
        )
    )
    print(f"direct callback: {len(callback_events)} events")

    listener_events: list[ProgressEvent] = []
    with progress_listener(listener_events.append):
        list(
            AdaptiveProgress(
                values,
                desc="Global listener",
                total=len(values),
                disable=disable,
            )
        )
    print(f"global listener: {len(listener_events)} events")

    recorder = ProgressRecorder()
    with recorder:
        results = list(
            parallel_map(
                _work,
                values,
                n_jobs=1,
                desc="Progress recorder",
                total=len(values),
                disable=disable,
                return_results=True,
            )
            or []
        )
    print(
        f"recorder: results={results}, events={len(recorder.history())}, "
        f"active={recorder.active_tasks()}"
    )

    args.jsonl.parent.mkdir(parents=True, exist_ok=True)
    before = 0
    if args.jsonl.exists():
        before = len(args.jsonl.read_text(encoding="utf-8").splitlines())
    with progress_jsonl_sink(args.jsonl):
        list(
            AdaptiveProgress(
                values,
                desc="JSONL sink",
                total=len(values),
                disable=disable,
            )
        )
    records = [
        json.loads(line) for line in args.jsonl.read_text(encoding="utf-8").splitlines()[before:]
    ]
    print(f"JSONL sink: wrote {len(records)} events to {args.jsonl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
