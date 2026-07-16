# 结构恢复

从只有元素和三维坐标的 frame 推断分子键、键级、形式电荷和自由基状态。

## 最短示例

```python
from molop import AutoParser

frame = AutoParser("water_mp2.out", n_jobs=1)[0][-1]
mol = frame.rdmol

if mol is None:
    raise RuntimeError("结构恢复失败")

print(mol.GetNumAtoms(), mol.GetNumBonds())
print(frame.smiles)
print(frame.topology_reconstruction_status)
```

共享样例输出：

```text
3 2
[H]O[H]
succeeded
```

访问 `frame.rdmol` 时，如果源格式没有提供分子图，MolOP 会按需调用 MolGR 从坐标恢复结构。

## 状态含义

| 状态 | 含义 |
| --- | --- |
| `provided` | 源格式已提供键、形式电荷和自由基信息 |
| `succeeded` | MolGR 从坐标得到正常候选 |
| `suspicious_fallback` | 得到可用候选，但属于需要复核的 fallback |
| `failed` | 未能得到 RDKit 分子，`rdmol` 为 `None` |
| `None` | 尚未触发结构访问，或当前对象没有状态 |

## 为什么结果需要复核

仅从几何恢复拓扑不是唯一问题。近距离接触、离子对、自由基、金属配合物和异常几何可能有
多个合理成键方案。`suspicious_fallback` 不等于无法使用，但不应无检查地进入高可信数据库。

## 指定电荷和多重度

源文件缺少可靠电荷信息时可在解析入口覆盖：

```python
batch = AutoParser(
    "radical.xyz",
    total_charge=0,
    total_multiplicity=2,
)
mol = batch[0][0].rdmol
```

覆盖值会影响结构恢复，不应只为“让算法成功”而猜测。

## 导出前检查

```python
frame = batch[0][-1]
mol = frame.rdmol

if mol is not None and frame.topology_reconstruction_status != "suspicious_fallback":
    print(frame.to_canonical_SMILES())
```

SDF、SMILES 和 CML 等图级 writer 依赖可用分子图；XYZ 和坐标模式的 Gaussian/ORCA input
可在不保留完整拓扑语义时仍导出坐标。

## 下一步

- [格式转换与导出](conversion.md)
- [SDF/MOL 格式](../reference/formats/sdf.md)
- [结构恢复 API](../reference/api/structure.md)
