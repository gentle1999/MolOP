# Reader/Writer 插件开发

MolOP 的 IO 扩展分为两类：仓库内置 codec 和独立发行的第三方插件。新增格式前，先确认现有专用 codec 或 OpenBabel fallback 不能满足需求。

完整的 File/Frame 类层级、source span、解析诊断和 renderer 约束见 [Parser 契约](../contracts/parser.md)。本页只说明如何把一个 codec 接入 registry。

## 选择注册方式

| 场景 | 注册方式 |
| --- | --- |
| MolOP 仓库内置格式 | 在 `src/molop/io/logic/` 下的公开 `*FileParser.py` 或 `*File.py` 模块中定义 `register(registry)`；catalog 会自动扫描。 |
| 独立安装的第三方包 | 在包自身的 `pyproject.toml` 中声明 `molop.codecs` entry point。 |
| 宿主应用私有集成 | 在启动阶段取得目标 `Registry`，显式调用插件的 `register(registry)`。 |

内置格式不需要修改 MolOP 的 `pyproject.toml`，也不应声明 entry point。第三方包不在 builtin 扫描范围内，必须使用 entry point 或由宿主应用显式注册。

## 最小第三方包

推荐将解析、模型和注册分开：

```text
molop-myfmt/
  pyproject.toml
  src/
    molop_myfmt/
      __init__.py
      plugin.py
      reader.py
      writer.py
```

在插件包的 `pyproject.toml` 中声明依赖和 entry point：

```toml
[project]
name = "molop-myfmt"
version = "0.1.0"
dependencies = ["molop"]

[project.entry-points."molop.codecs"]
myfmt = "molop_myfmt.plugin:register"
```

示例中的无约束依赖只用于展示最小包结构。发布插件前，必须把它替换成实际测试过的 MolOP 版本范围。开发版文档可能领先于 PyPI；若插件依赖顶部横幅标记为未发布的接口，开发阶段应固定横幅中的 commit，等接口进入正式版本后再发布插件并设置依赖下限。详见[文档版本](../release/versioning.md)。

entry point 的目标必须是可调用的 `register(registry)`，也可以是带有同名方法的模块对象。加载发生在 registry 首次激活时，不应在模块导入阶段解析文件、加载大型模型或修改全局 registry。

```python
from __future__ import annotations

from typing import TYPE_CHECKING

from .reader import MyFmtReader

if TYPE_CHECKING:
    from molop.io.codec_registry import Registry


def register(registry: Registry) -> None:
    @registry.reader_factory(
        format_id="myfmt",
        extensions={".myfmt"},
        priority=100,
    )
    def make_reader() -> MyFmtReader:
        return MyFmtReader()
```

同一 entry point 可以注册多个 reader 和 writer。factory 应保持轻量；只有实际选择 codec 时才构造昂贵资源。

## Reader 契约

Reader 至少提供 `format_id`、`extensions`、`priority` 和 `read(path, **kwargs)`，并返回 `ParseResult`。`read()` 返回值必须是可加入 MolOP disk batch 的文件模型，而不是裸字典、frame 或外部工具对象。

```python
from __future__ import annotations

from pathlib import Path
from typing import Any

from molop.io.codec_types import ParseOptions, ParseResult, StructureLevel

from .parser import MyFmtFileParserDisk


class MyFmtReader:
    format_id = "myfmt"
    extensions = frozenset({".myfmt"})
    priority = 100

    def probe_file_format(self, path: str | Path) -> bool:
        return MyFmtFileParserDisk.probe_file_format(path)

    def read(self, path: str | Path, **kwargs: Any) -> ParseResult[object]:
        supplied = kwargs.get("parse_options")
        if supplied is not None and not isinstance(supplied, ParseOptions):
            raise TypeError("parse_options must be a ParseOptions instance")
        options = (supplied or ParseOptions()).resolved()

        parser = MyFmtFileParserDisk(
            forced_charge=options.total_charge,
            forced_multiplicity=options.total_multiplicity,
            only_extract_structure=options.only_extract_structure,
            only_last_frame=options.only_last_frame,
            capture_source_evidence=options.capture_source_evidence,
            source_encoding=options.source_encoding,
            parse_options=options,
        )
        value = parser.parse(
            str(path),
            total_charge=options.total_charge,
            total_multiplicity=options.total_multiplicity,
            release_file_content=options.release_file_content,
        )
        return ParseResult(
            value=value,
            level=StructureLevel.COORDS,
            detected_format=self.format_id,
        )
```

`ParseOptions` 是 batch、reader、file parser 和 frame parser 之间的不可变配置快照。插件必须优先使用传入的 `parse_options`，不能在每个 frame 中重新读取全局配置。只有 reader 被独立调用且没有收到该参数时，才自行构造并解析默认值。

`StructureLevel.COORDS` 表示最低保证为原子和坐标；只有 reader 保证生成合法分子图时才使用 `StructureLevel.GRAPH`。

### Probe 与异常

`probe_file_format(path)` 是推荐能力，但不是 `ReaderCodec` 的硬性要求。未实现时 registry 会把 reader 视为候选并直接进入正式解析，因此共享扩展名或 fallback reader 应实现快速、无副作用的 probe。

- 内容不属于当前格式时抛出 `FormatMismatchError`，batch 才会尝试下一个候选 reader。
- 文件属于当前格式但内容损坏或 parser 有缺陷时，抛出具体解析异常，不要伪装成格式不匹配。
- 可恢复的问题放入 `ParseResult.warnings`，不要只写日志。

```python
from molop.io.codec_types import ParseResult, ParseWarning, StructureLevel

return ParseResult(
    value=value,
    level=StructureLevel.COORDS,
    warnings=(
        ParseWarning(
            code="MYFMT.MISSING_OPTIONAL_BLOCK",
            message="Optional property block was not present.",
        ),
    ),
    detected_format="myfmt",
)
```

这些 warning 会保留在 `AutoParser(..., return_report=True)` 返回的逐文件结果中。默认 `AutoParser()` 仍只返回成功解析的 batch。

### 解析状态与源加载

Frame parser 必须使用不可变 `FrameParseContext` 读取当前 block、file/segment metadata 和 `ParseOptions`，不得把当前 block、游标或 metadata 保存到 parser 实例。组装 frame 时，context metadata 会覆盖格式 payload 中的同名字段。

基于 `BaseFileParserDisk` 的 parser 通常直接调用 `parse(path)`。如果上层已经读取源内容，可使用 `parse_bytes(..., file_path=...)` 或 `parse_decoded_source(..., file_path=...)`，在保留磁盘来源身份和精确 source offset 的同时避免再次打开文件。不要为了 probe 主动读取完整文件；probe/full-read 复用应由明确拥有预载源的调用层控制。

## Writer 契约

Writer 至少提供 `format_id`、`priority`、`required_level` 和 `write(value, **kwargs)`。通过 factory 注册时还要明确输出域和图策略：

```python
def register(registry: Registry) -> None:
    @registry.writer_factory(
        format_id="myfmt",
        required_level=StructureLevel.COORDS,
        domain="file",
        default_graph_policy="coords",
        priority=100,
    )
    def make_writer() -> MyFmtWriter:
        return MyFmtWriter()
```

- `domain="file"` 用于文件组合、多帧合并或拆分。
- `domain="frame"` 只用于格式确实具有独立单帧语义的情况。
- `default_graph_policy` 必须与格式要求一致；不要依赖隐式猜测。
- 重要的格式专用指令或原文区块应尽量保留，以支持 round-trip。

CLI 的动态 writer 选项由 `write()` 以及关联 file/frame 类 `_render()` 的具名参数推导。为参数提供准确类型、默认值和 docstring；`Literal` 会生成候选值，`bool` 会生成布尔开关。不要通过 `**kwargs` 隐藏需要暴露给用户的选项。

## 加载与优先级

Reader 和 writer 按 `priority` 从高到低选择；相同优先级按确定性的注册顺序排列。插件应只在确有覆盖意图时使用高于内置 codec 的优先级，并为扩展名冲突增加测试。

第三方 entry point 会按名称和值稳定排序后加载。单个插件导入失败、缺少可调用的 `register()` 或注册失败时，MolOP 会发出 `CodecPluginWarning` 并跳过该插件，不会阻止其他 codec 激活。插件测试应把这些 warning 当作安装或兼容性错误处理。

## 验证清单

至少覆盖以下行为：

1. entry point 安装后能够被默认 registry 发现，且仅注册一次。
2. 显式 `parser_detection="myfmt"` 和基于扩展名的自动选择都能解析 fixture。
3. probe 不匹配时能回退到下一个 reader，真实解析错误不会被吞掉。
4. `ParseOptions` 能到达 file/frame parser，reader warning 能出现在结构化报告中。
5. 返回的 disk file model 能加入 batch，并支持 summary、frame 选择和 `format_transform()`。
6. writer 覆盖字符串渲染、写盘、命名、file/frame domain 和失败行为。
7. 在最低支持 MolOP 版本和当前主线各运行一次契约测试。

仓库内置 codec 还需要同步更新中英文格式页和用户示例，并运行[质量门禁](../contributing/quality.md)。不要在文档中维护第二份格式清单；能力状态由 `tests/format_feature_coverage/support_matrix.py` 和格式页生成区段维护。
