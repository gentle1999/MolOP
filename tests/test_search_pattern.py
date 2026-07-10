import pytest

from molop.io.base_models.SearchPattern import MolOPPattern


def test_molop_pattern_find_matches_preserves_named_groups() -> None:
    pattern = MolOPPattern(
        content_pattern=r"STATE\s+(?P<root>\d+):\s+(?P<ev>[-+]?\d+\.\d+)\s+eV",
        content_repeat=0,
    )

    matches = pattern.find_matches("STATE 1: 2.50 eV\nSTATE 2: 3.75 eV\n")

    assert pattern.group_names == ("root", "ev")
    assert pattern.has_named_group("ev") is True
    assert pattern.has_named_group("missing") is False
    assert [matched.group("root") for matched in matches] == ["1", "2"]
    assert matches[1].group("ev") == "3.75"
    assert matches[0].groupdict() == {"root": "1", "ev": "2.50"}
    assert pattern.get_named_group(matches[0], "ev") == "2.50"
    assert pattern.require_named_group(matches[1], "ev") == "3.75"
    assert pattern.named_group_dict(matches[0]) == {"root": "1", "ev": "2.50"}
    assert pattern.get_group(matches[0], "root") == "1"
    assert pattern.require_group(matches[0], "root") == "1"
    assert pattern.group_dict(matches[0]) == {"root": "1", "ev": "2.50"}
    assert pattern.find_named_group_values(
        "STATE 1: 2.50 eV\nSTATE 2: 3.75 eV\n",
        "ev",
    ) == [
        "2.50",
        "3.75",
    ]
    assert pattern.find_first_named_group_value(
        "STATE 1: 2.50 eV\nSTATE 2: 3.75 eV\n",
        "ev",
    ) == ("2.50")
    assert pattern.find_named_group_dicts("STATE 1: 2.50 eV\n")[0] == {
        "root": "1",
        "ev": "2.50",
    }
    assert pattern.find_group_values("STATE 1: 2.50 eV\nSTATE 2: 3.75 eV\n", "ev") == [
        "2.50",
        "3.75",
    ]
    assert pattern.find_first_group_value("STATE 1: 2.50 eV\nSTATE 2: 3.75 eV\n", "ev") == ("2.50")
    assert pattern.find_groupdicts("STATE 1: 2.50 eV\n")[0] == {"root": "1", "ev": "2.50"}
    assert pattern.find_group_dicts("STATE 1: 2.50 eV\n")[0] == {
        "root": "1",
        "ev": "2.50",
    }

    searched = pattern.search("prefix STATE 3: 4.10 eV")
    assert searched is not None
    assert searched.group("root") == "3"

    line_matched = pattern.match("STATE 4: 5.20 eV")
    assert line_matched is not None
    assert line_matched.group("ev") == "5.20"


def test_molop_pattern_named_group_helpers_use_regex_error_strategy() -> None:
    pattern = MolOPPattern(
        content_pattern=r"STATE[ \t]+(?P<root>\d+)(?:[ \t]+(?P<label>\w+))?",
        content_repeat=0,
    )

    assert pattern.find_named_group_values("STATE 1\nSTATE 2 excited\n", "label") == ["excited"]
    first_match = pattern.find_matches("STATE 1\n")[0]
    assert pattern.get_named_group(first_match, "label") is None
    with pytest.raises(ValueError, match="group did not match"):
        pattern.require_named_group(first_match, "label")
    with pytest.raises(IndexError, match="no such group"):
        pattern.find_named_group_values("STATE 1\n", "missing")
    with pytest.raises(IndexError, match="no such group"):
        pattern.find_first_named_group_value("unmatched text", "missing")


def test_molop_pattern_without_content_pattern_has_no_named_groups() -> None:
    pattern = MolOPPattern(start_pattern="BEGIN", end_pattern="END")

    assert pattern.group_names == ()
    assert pattern.find_matches("BEGIN\nvalue\nEND") == []
    with pytest.raises(IndexError, match="no such group"):
        pattern.find_named_group_values("BEGIN\nvalue\nEND", "value")


def test_molop_pattern_escapes_literal_text_for_dynamic_patterns() -> None:
    label = "Energy change (2.0)"
    pattern = MolOPPattern(
        content_pattern=rf"{MolOPPattern.escape_literal(label)}:\s+(?P<value>\d+)",
    )

    matched = pattern.search("Energy change (2.0): 5")

    assert matched is not None
    assert matched.group("value") == "5"


def test_molop_pattern_cursor_api_lives_on_base_pattern() -> None:
    pattern = MolOPPattern(
        start_pattern=r"^BEGIN$",
        end_pattern=r"^END$",
        content_pattern=r"value=(?P<value>\d+)",
        content_repeat=0,
    )
    content = "skip\nBEGIN\nvalue=1\nEND\nnoise\nBEGIN\nvalue=2\nEND\n"

    first_matches = pattern.find_content_matches(content)
    assert first_matches is not None
    assert [matched.group("value") for matched in first_matches] == ["1"]

    legacy_matches = pattern.get_matches("value=3\nvalue=4\n")
    assert legacy_matches is not None
    assert [matched.group("value") for matched in legacy_matches] == ["3", "4"]

    located_matches = pattern.match_content(content)
    assert located_matches is not None
    assert [matched.group("value") for matched in located_matches] == ["1"]

    first_block, cursor = pattern.split_content_from(content)
    assert "value=1" in first_block
    second_matches = pattern.find_content_matches(content, cursor + len("END"))
    assert second_matches is not None
    assert [matched.group("value") for matched in second_matches] == ["2"]


def test_molop_pattern_negative_repeat_preserves_regex_match_direction() -> None:
    pattern = MolOPPattern(
        content_pattern=r"(?r)^value=(?P<value>\d+)$",
        content_repeat=-2,
    )

    matches = pattern.find_matches("value=1\nvalue=2\nvalue=3\n")

    assert [matched.group("value") for matched in matches] == ["3", "2"]
