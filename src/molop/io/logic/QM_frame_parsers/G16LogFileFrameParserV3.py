from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

from pydantic import PrivateAttr

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.QM_frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.QM_frame_models.G16V3Components import (
    G16V3BaseComponent,
    G16V3ComponentNode,
    G16V3ComponentTree,
    G16V3ComponentTreeBuilder,
    _has_meaningful_value,
)


class G16LogFileFrameParserV3Mixin:
    _last_component_tree: G16V3ComponentTree | None = PrivateAttr(default=None)

    @property
    def last_component_tree(self) -> G16V3ComponentTree | None:
        if self._last_component_tree is not None:
            self._ensure_synthetic_children(self._last_component_tree)
        return self._last_component_tree

    def build_component_tree(self, block: str) -> G16V3ComponentTree:
        return G16V3ComponentTreeBuilder.from_g16_text(
            block,
            only_extract_structure=cast(_HasParseMethod, self).only_extract_structure,
        )

    def _ensure_synthetic_children(self, tree: G16V3ComponentTree) -> G16V3ComponentTree:
        return G16V3ComponentTreeBuilder.expand_synthetic_children(tree)

    def _append_synthetic_child(
        self,
        parent_node: G16V3ComponentNode,
        component_cls: type[G16V3BaseComponent],
        payload: Mapping[str, Any],
        *,
        node_name: str | None = None,
    ) -> G16V3ComponentNode | None:
        normalized_payload = {
            key: value for key, value in payload.items() if _has_meaningful_value(value)
        }
        if not normalized_payload:
            return None
        parent_component = parent_node.component
        child_component = component_cls(
            span_start=parent_component.span_start if parent_component is not None else 0,
            span_end=parent_component.span_end if parent_component is not None else 0,
            raw_text="",
            payload=dict(normalized_payload),
        )
        child_node = G16V3ComponentNode(
            node_name=node_name or child_component.component_name,
            component=child_component,
            include_in_aggregation=False,
            include_in_render=False,
        )
        parent_node.children.append(child_node)
        return child_node

    def _expand_l716_children(self, tree: G16V3ComponentTree) -> None:
        for node in tree.root.children:
            component = node.component
            if component is None:
                continue
            for child_spec in component.build_synthetic_children():
                self._append_synthetic_child(
                    node,
                    child_spec.component_cls,
                    child_spec.payload,
                    node_name=child_spec.node_name,
                )

    def _parse_frame(self) -> Mapping[str, Any]:
        tree = self.build_component_tree(cast(_HasParseMethod, self)._block)
        for component in tree.iter_components():
            component.parse_payload(cast(_HasParseMethod, self)._block, tree)
        self._last_component_tree = tree
        contract_issues = tree.validate_contracts()
        if contract_issues:
            raise ValueError("Protocol tree contract violations: " + "; ".join(contract_issues))
        return {"component_tree_input": tree}


class G16LogFileFrameParserV3Memory(
    G16LogFileFrameParserV3Mixin, BaseFrameParser[G16LogFileFrameMemory]
):
    _file_frame_class_ = G16LogFileFrameMemory


class G16LogFileFrameParserV3Disk(
    G16LogFileFrameParserV3Mixin, BaseFrameParser[G16LogFileFrameDisk]
):
    _file_frame_class_ = G16LogFileFrameDisk
