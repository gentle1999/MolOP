# molop.utils.progressbar

进度条与并行执行辅助函数。native 重建边界和 loky 执行器生命周期由独立执行层管理；本模块保留
旧的 guard 导入路径以兼容已有调用方。

所有捕获方式的可直接运行演示位于
[`scripts/progress_capture_demo.py`](https://github.com/gentle1999/MolOP/blob/main/scripts/progress_capture_demo.py)：

```bash
uv run --frozen python scripts/progress_capture_demo.py
```

MolOP 使用 [tqdm](https://github.com/tqdm/tqdm) 渲染进度条，但每个进度条同时还会广播
机器可读的 [`ProgressEvent`][molop.utils.progressbar.ProgressEvent] 事件。
第三方工具无需解析终端输出即可捕获 MolOP 的处理进度：

- **进程内** — 通过
  [`register_progress_listener`][molop.utils.progressbar.register_progress_listener]
  注册监听器（或使用
  [`progress_listener`][molop.utils.progressbar.progress_listener] 上下文管理器），
  同步接收每一个事件。
- **进程内、实时状态** — 使用线程安全的
  [`ProgressRecorder`][molop.utils.progressbar.ProgressRecorder] 维护所有运行中
  进度条的实时快照，任何线程、任何时刻（包括进度条尚未跑完时）都可以查询。
- **跨进程** — 使用
  [`progress_jsonl_sink`][molop.utils.progressbar.progress_jsonl_sink]
  将 JSON Lines 追加写入文件，其他工具可以 tail 或读取该文件。

MolOP 只会关闭自己创建且已经空闲的 loky 执行器。宿主程序已有的空闲执行器不会被关闭；如果
宿主执行器仍有待处理任务，native 重建会 fail-closed 并要求先消费或关闭该结果流。

即使终端进度条被关闭（`molopconfig.show_progress_bar = False`），事件仍然会发出，
因为事件发送与渲染相互独立。每个进度条发出 `start` → (`update` / `description`)*
→ `close` 事件，携带稳定的 `task_id`、描述、`done`/`total` 计数和时间戳。

```python
import molop.utils.progressbar as pb


def on_progress(event: pb.ProgressEvent) -> None:
    print(f"{event.description}: {event.done}/{event.total}")


with pb.progress_listener(on_progress):
    list(pb.AdaptiveProgress(range(5), desc="parse", total=5, disable=True))

# 跨进程捕获：其他工具执行 tail -f progress.jsonl
with pb.progress_jsonl_sink("progress.jsonl"):
    list(pb.AdaptiveProgress(range(5), desc="parse", total=5, disable=True))
```

??? example "回调输出示例"

    ```text
    parse: 0/5
    parse: 1/5
    parse: 2/5
    parse: 3/5
    parse: 4/5
    parse: 5/5
    ```

## 实时、线程安全的快照（`ProgressRecorder`）

每个事件都是**快照**：它携带该进度条完整的 `done` / `total` / `description`
状态，而不是增量。因此消费者只需要用事件覆盖自己的记录即可——无需累加、不会漂移，
即使漏掉某个事件，下一个事件到来时也会自动重新对齐。

[`ProgressRecorder`][molop.utils.progressbar.ProgressRecorder] 是该模式的官方实现。
它本身就是一个监听器，可以直接用作上下文管理器；任何其他线程都可以通过
[`snapshot()`][molop.utils.progressbar.ProgressRecorder.snapshot] 查询——例如 GUI
轮询或仪表盘——即使在进度条尚未跑完时，也能拿到与条完全一致的实时状态。

```python
from molop.utils.progressbar import ProgressRecorder, parallel_map

recorder = ProgressRecorder()

with recorder:                      # 自动注册 + 注销监听器
    results = parallel_map(
        lambda value: value * 2,
        range(5),
        n_jobs=1,
        total=5,
        disable=True,
        return_results=True,
    )

# 运行结束后：完整事件时间线
print(recorder.history())

# 运行中，从其他线程（例如 GUI 更新循环）：
#   for task_id, event in recorder.snapshot().items():
#       widget.set_value(event.done / event.total if event.total else 0)
#   percent = recorder.percent(task_id)   # 0.0-100.0；未知时为 None
```

??? example "快照输出示例"

    ```text
    (ProgressEvent(task_id=1, event='start', done=0, total=5), ...)
    ```

同一个 recorder 也可以不用上下文管理器，作为普通监听器使用：

```python
from molop.utils.progressbar import parallel_map, register_progress_listener

recorder = ProgressRecorder()
unregister = register_progress_listener(recorder)
try:
    parallel_map(lambda value: value * 2, range(5), n_jobs=1, disable=True)
finally:
    unregister()
```

第三方消费者需要注意的要点：

- 按 `task_id` 区分记录——并发的多个进度条不会互相覆盖。
- `close` 并不代表 100%：被中断的运行会以 `done < total` 关闭，请始终使用事件
  自带的 `done` / `total` 数值。
- 监听器在驱动进度条的线程里同步调用。串行进度条使用调用方线程；`parallel_map`
  可能在 Joblib 的完成回调线程中发出事件。如果消费者需要更新 UI，请把更新
  marshal 回 UI/事件循环线程，或者直接轮询 `snapshot()`。
- 事件频率可能高于终端重绘频率（tqdm 只节流渲染、不节流事件）。计算百分比请用
  `done` / `total`，不要数事件条数；必要时自己节流渲染。

::: molop.utils.progressbar
