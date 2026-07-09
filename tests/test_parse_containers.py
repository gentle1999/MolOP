from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext


class _SplitOncePattern:
    def __init__(self, marker: str) -> None:
        self.marker = marker

    def split_content(self, content: str) -> tuple[str, str]:
        head, separator, tail = content.partition(self.marker)
        if separator:
            return head, tail
        return "", content


def test_text_parse_context_advances_local_content_only() -> None:
    context = TextParseContext("header--body")

    focus = context.split(_SplitOncePattern("--"))

    assert focus == "header"
    assert context.content == "body"


def test_model_parse_result_returns_copy_and_merges_payload_fields() -> None:
    result = ModelParseResult({"software": "Gaussian", "empty": ""})

    result.set_missing_from({"software": "ORCA", "empty": "filled", "charge": 0})
    result.merge_payload_field(
        "energies",
        {"reference_energy": -1.0, "unused": None},
    )
    result.merge_payload_field(
        "energies",
        {"reference_energy": -2.0, "mp2_energy": -2.1},
    )

    data = result.model_data()
    data["software"] = "mutated"

    assert result.model_data()["software"] == "Gaussian"
    assert result.model_data()["empty"] == "filled"
    assert result.model_data()["charge"] == 0
    assert result.model_data()["energies"] == {
        "reference_energy": -1.0,
        "mp2_energy": -2.1,
    }


def test_model_parse_result_can_overwrite_payload_fields() -> None:
    result = ModelParseResult({"energies": {"reference_energy": -1.0}})

    result.merge_payload_field("energies", {"reference_energy": -2.0}, overwrite=True)

    assert result.model_data()["energies"] == {"reference_energy": -2.0}
