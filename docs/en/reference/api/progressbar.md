# molop.utils.progressbar

Progress bars and parallel execution helpers.

A directly runnable demonstration of every capture mode is available as
[`scripts/progress_capture_demo.py`](https://github.com/gentle1999/MolOP/blob/main/scripts/progress_capture_demo.py):

```bash
uv run --frozen python scripts/progress_capture_demo.py
```

MolOP renders progress with [tqdm](https://github.com/tqdm/tqdm), but every bar
also broadcasts machine-readable
[`ProgressEvent`][molop.utils.progressbar.ProgressEvent] notifications.
Third-party tools can capture MolOP progress without parsing terminal output:

- **In-process** — register a listener with
  [`register_progress_listener`][molop.utils.progressbar.register_progress_listener]
  (or use the
  [`progress_listener`][molop.utils.progressbar.progress_listener] context
  manager) and receive every event synchronously.
- **In-process, live state** — use the thread-safe
  [`ProgressRecorder`][molop.utils.progressbar.ProgressRecorder] to keep an
  up-to-date snapshot of every running bar, readable from any thread at any
  time (even while the bars are still running).
- **Cross-process** — use
  [`progress_jsonl_sink`][molop.utils.progressbar.progress_jsonl_sink] to append
  JSON Lines to a file that other tools can tail or read.

Events are emitted even when the terminal bar is suppressed
(`molopconfig.show_progress_bar = False`), because event emission is
independent of rendering. A single bar emits `start` → (`update` /
`description`)* → `close` events carrying a stable `task_id`, the description,
`done`/`total` counters and a timestamp.

```python
import molop.utils.progressbar as pb


def on_progress(event: pb.ProgressEvent) -> None:
    print(f"{event.description}: {event.done}/{event.total}")


with pb.progress_listener(on_progress):
    list(pb.AdaptiveProgress(range(5), desc="parse", total=5, disable=True))

# Cross-process capture: tail -f progress.jsonl from another tool
with pb.progress_jsonl_sink("progress.jsonl"):
    list(pb.AdaptiveProgress(range(5), desc="parse", total=5, disable=True))
```

??? example "Example callback output"

    ```text
    parse: 0/5
    parse: 1/5
    parse: 2/5
    parse: 3/5
    parse: 4/5
    parse: 5/5
    ```

## Keeping a live, thread-safe snapshot (`ProgressRecorder`)

Every event is a *snapshot*: it carries the full `done` / `total` /
`description` state of its bar, never a delta. A consumer therefore only needs
to overwrite its record with each event — no accumulation, no drift, and a
missed event self-heals on the next one.

[`ProgressRecorder`][molop.utils.progressbar.ProgressRecorder] is the official
implementation of that pattern. It is a listener itself, so it works as a
context manager; from any other thread you can query
[`snapshot()`][molop.utils.progressbar.ProgressRecorder.snapshot] — for example
a GUI poller or a dashboard — and it reflects the bars' current state even
while they are still running.

```python
from molop.utils.progressbar import ProgressRecorder, parallel_map

recorder = ProgressRecorder()

with recorder:                      # registers + unregisters automatically
    results = parallel_map(
        lambda value: value * 2,
        range(5),
        n_jobs=1,
        total=5,
        disable=True,
        return_results=True,
    )

# After the run: full event timeline
print(recorder.history())

# Mid-run, from another thread (e.g. a GUI update loop):
#   for task_id, event in recorder.snapshot().items():
#       widget.set_value(event.done / event.total if event.total else 0)
#   percent = recorder.percent(task_id)   # 0.0-100.0, or None if unknown
```

??? example "Example snapshot output"

    ```text
    (ProgressEvent(task_id=1, event='start', done=0, total=5), ...)
    ```

The same recorder works as a plain listener without the context manager:

```python
from molop.utils.progressbar import parallel_map, register_progress_listener

recorder = ProgressRecorder()
unregister = register_progress_listener(recorder)
try:
    parallel_map(lambda value: value * 2, range(5), n_jobs=1, disable=True)
finally:
    unregister()
```

Key points for third-party consumers:

- Key the recorded state by `task_id` — concurrent bars never collide.
- `close` does **not** imply 100%: an interrupted run closes with `done < total`.
  Always use the event's own `done` / `total` values.
- The listener is invoked synchronously from the thread driving the bar. A
  serial bar runs on the caller's thread; `parallel_map` may emit from
  joblib's completion callback thread. Marshal UI updates back to your
  UI/event-loop thread, or simply poll `snapshot()`.
- Events may arrive faster than the terminal redraws (tqdm throttles
  rendering, not emission). Compute percentages from `done` / `total`, not from
  event counts, and throttle your own rendering if needed.

::: molop.utils.progressbar
