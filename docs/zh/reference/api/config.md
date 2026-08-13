# molop.config

本页面提供了 `molop.config` 模块的 API 参考。

模块暴露进程级 `molopconfig` 对象，用于控制进度显示、并行任务上限、日志和结构恢复选项。
单次调用的行为优先使用显式函数参数；只有需要影响后续调用的全局策略时才修改全局配置。

## 并行默认值

`max_jobs` 默认是 `None`，表示自动检测当前进程可用的 CPU。MolOP 会在 joblib 的调度器/容器
配额、CPU affinity 和 `os.cpu_count()` 中取最严格的限制。因此 `n_jobs=-1` 会使用当前进程可用
的全部 CPU，但不会超过调度器分配的额度。

需要设置进程级上限时传入正整数：

```python
from molop import molopconfig

molopconfig.max_jobs = 32
```

由环境创建默认配置时，可设置 `MOLOP_MAX_JOBS=32`。即使环境变量存在，显式构造
`MolOPConfig(max_jobs=None)` 仍会保持自动检测。正数 `n_jobs` 仍表示单次调用上限，并受
`effective_max_jobs` 限制；`n_jobs=1` 继续表示串行排查模式。CLI 可通过
`molop --max-jobs 32 parse ...` 为单次命令设置相同的进程级上限。

::: molop.config
