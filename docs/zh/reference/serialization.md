# Source evidence 与序列化

把 MolOP 的文件级 metadata、frame 结果和 source evidence 导出给数据库或审计系统。普通分析使用
`AutoParser`、`to_summary_df` 或 frame 字段；需要跨系统保存解析事实时再使用本页接口。

## 导出 payload

```python
from molop import AutoParser

chem_file = AutoParser(
    "water_mp2.out",
    capture_source_evidence=True,
    release_file_content=True,
    n_jobs=1,
)[0]
file_payload = chem_file.to_unitless_dump_with_unit_keys(exclude_none=True)
frame_payloads = [
    frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    for frame in chem_file
]

print(file_payload["schema_version"], file_payload["source_format"])
print(len(frame_payloads), frame_payloads[0]["file_frame_index"])
```

??? example "真实输出"

    ```text
    molop-calculation-export-v1 orcaout
    1 0
    ```

文件 payload 和 frame payload 分开保存。文件 payload 负责来源与文件级 metadata；frame payload
负责坐标、拓扑、科学结果和 frame 级 source locator。

## 数值与单位

`to_unitless_dump_with_unit_keys()` 将 Pint quantity 转为 magnitude，并把规范单位写入 key，例如
`Energy.total_energy.hartree` 或 `coords (angstrom)`。数组默认序列化为 list；数据库使用数值
sidecar 时传入 `array_mode="ndarray"`。

| 对象 | 用途 | 稳定字段 |
| --- | --- | --- |
| 文件 payload | 文件 identity、来源、parser provenance | `schema_version`、`source_format`、`parser_provenance` |
| frame payload | 结构和计算结果 | `file_frame_index`、`source_span`、`parse_presence` |
| source evidence | 回溯原文位置和摘要 | `source_segments`、`source_block_sha256` |

## 索引与边界

- 用 `file_frame_index` 表示源文件中的稳定 frame 顺序。
- `frame_id` 属于当前保留的 frame 集合；筛选或重排后不应作为全局 identity。
- `capture_source_evidence=True` 才会请求额外的 source locator 和 parser provenance。
- MolOP 提供解析事实，不负责数据库 identity、payload 校验、数组编码、artifact 校验或 admission/QC。

详细的 source span 和 parser lifecycle 规则见[完整 Parser 契约](../developer/parser-contract.md)。

## 相关页面

- [API 契约](api_contracts.md)
- [按科学性质查找字段](model_fields.md)
- [格式支持概览](format_support.md)
