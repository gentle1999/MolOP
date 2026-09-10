# molop.config

本页面提供了 `molop.config` 模块的 API 参考。

模块暴露进程级 `molopconfig` 对象，用于控制进度显示、并行任务上限、日志和结构恢复选项。
单次调用的行为优先使用显式函数参数；只有需要影响后续调用的全局策略时才修改全局配置。

## 导入副作用与显式初始化

导入 `molop` 不会创建 `molop.log`、修改宿主进程的 Python 递归上限，也不会默认关闭 RDKit 或
Open Babel 的原生日志。`rdkit-dof` 也只在需要景深绘制或显式调用相关设置时加载。需要文件日志或
安静的原生日志时，显式创建配置并启用对应选项：

```python
from molop.config import MolOPConfig

config = MolOPConfig(
    log_to_file=True,
    log_file_path="run.log",
    suppress_rdkit_logs=True,
    suppress_openbabel_logs=True,
)
config.configure_native_logging()
```

`enable_file_logging()`、`disable_file_logging()` 负责文件 handler 的生命周期；
`set_max_recursion_depth()` 只在调用方明确请求时修改当前进程的递归上限。

主进程拓扑预热由 `prewarm_topologies` 控制，默认值为 `False`。默认情况下，依赖分子图的操作会在
使用结果的 spawn-like `loky` worker 中按需惰性重建，避免单独的主进程预热步骤。需要确定性主进程
缓存的工作流可以显式开启：

```python
from molop import molopconfig

molopconfig.prewarm_topologies = True
```

MolGR 重建失败默认采用严格模式。需要保留错误结构以便复核时，可以开启原始 fallback 保留：

```python
molopconfig.reconstruction_failure_policy = "return_suspicious"
```

保留的分子图会标记为 `topology_reconstruction_status == "suspicious_fallback"`，不能作为可信
化学结构直接使用。

## 并行默认值

`max_jobs` 默认是 `None`，表示自动检测当前进程可用的 CPU。MolOP 会在 joblib/loky 可检测的
调度器和容器限制、`os.cpu_count()`，以及操作系统支持时的当前进程 affinity 中取最严格的
有效值。MolOP 将 psutil 作为直接依赖，以启用 joblib 在 Windows 等受支持平台上的 affinity
fallback。因此 `n_jobs=-1` 会使用当前进程可用的 CPU，而不超过已检测到的分配额度。

Linux 支持 affinity 和 cgroup quota 检测。Windows 的进程 affinity 通过 psutil 获取，joblib 还会
应用 Windows 进程池的安全 worker 上限。macOS 没有可移植的进程级 CPU affinity API，因此自动
值使用当前进程可见 CPU 和 joblib 可检测的限制。若运行平台无法向进程暴露实际配额，所有平台
都可以通过 `max_jobs`、`MOLOP_MAX_JOBS` 或 `LOKY_MAX_CPU_COUNT` 显式设置上限。

需要设置进程级上限时传入正整数：

```python
from molop import molopconfig

molopconfig.max_jobs = 32
```

由环境创建默认配置时，可设置 `MOLOP_MAX_JOBS=32`。即使环境变量存在，显式构造
`MolOPConfig(max_jobs=None)` 仍会保持自动检测。正数 `n_jobs` 仍表示单次调用上限，并受
`effective_max_jobs` 限制；`n_jobs=1` 继续表示串行排查模式。CLI 可通过
`molop --max-jobs 32 parse ...` 为单次命令设置相同的进程级上限。

可能进入 MolGR 的任务使用独立的 `effective_molgr_max_jobs` 上限。自动值和显式值都会先受
`floor(2 * available_cpu_count / 3)` 限制（最少保留 1 个串行 worker），再受 `max_jobs` 限制。
该策略覆盖拓扑重建、依赖分子图的摘要、导出、格式转换以及 `filter_custom` 等用户回调；
文件解析以及不调用 MolGR 的其他操作仍使用完整的 `effective_max_jobs` 上限。

::: molop.config
