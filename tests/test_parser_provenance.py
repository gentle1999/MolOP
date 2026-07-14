from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import asdict, replace
from hashlib import sha256
from importlib.metadata import version as distribution_version

import pytest
from molgr.config import get_config as get_molgr_config
from molgr.config import set_config as set_molgr_config
from pydantic import ValidationError

from molop.config import molopconfig
from molop.io.base_models.ChemFile import BaseChemFile
from molop.io.base_models.ChemFileFrame import BaseChemFileFrame
from molop.io.base_models.source import ParserProvenance, canonical_json_sha256
from molop.io.logic.coords.parsers.XYZFileParser import XYZFileParserMemory


_XYZ_TEXT = "1\ncharge 0 multiplicity 1\nH 0.0 0.0 0.0\n"
_TWO_XYZ_TEXT = _XYZ_TEXT + "1\nsecond\nHe 1.0 0.0 0.0\n"


def test_file_schema_and_parser_provenance_stay_file_scoped() -> None:
    parsed = XYZFileParserMemory().parse(_XYZ_TEXT)

    assert BaseChemFile().schema_version == "molop-calculation-export-v1"
    assert parsed.schema_version == "molop-calculation-export-v1"
    assert parsed.model_dump()["schema_version"] == "molop-calculation-export-v1"
    assert parsed.parser_provenance is not None
    assert "schema_version" not in BaseChemFileFrame.model_fields
    assert "parser_provenance" not in BaseChemFileFrame.model_fields
    assert "parser_provenance" not in parsed[0].model_dump()


def test_parser_provenance_captures_complete_effective_config() -> None:
    parser = XYZFileParserMemory(
        forced_charge=-1,
        forced_multiplicity=2,
        only_extract_structure=True,
        only_last_frame=True,
        capture_source_evidence=True,
        source_encoding="UTF-8",
    )

    parsed = parser.parse(
        _TWO_XYZ_TEXT,
        total_charge=-2,
        total_multiplicity=3,
    )
    provenance = parsed.parser_provenance

    assert provenance is not None
    assert provenance.parser_id.endswith(".XYZFileParserMemory")
    assert provenance.parser_version == provenance.molop_version
    assert provenance.molop_version == distribution_version("molop")
    assert provenance.molgr_version == distribution_version("molgr")
    assert provenance.rdkit_version == distribution_version("rdkit")
    assert provenance.effective_config["parser"] == {
        "total_charge_override": -2,
        "total_multiplicity_override": 3,
        "only_extract_structure": True,
        "only_last_frame": True,
        "capture_source_evidence": True,
        "source_encoding": "utf-8",
    }
    assert provenance.effective_config["molop"] == {
        "force_unit_transform": molopconfig.force_unit_transform,
        "graph_reconstruction_backend": molopconfig.graph_reconstruction_backend,
        "make_dative_bonds": molopconfig.make_dative_bonds,
    }
    assert provenance.effective_config["molgr"] == asdict(get_molgr_config())


def test_effective_config_hash_is_strict_and_deterministic() -> None:
    first = XYZFileParserMemory().parse(_XYZ_TEXT).parser_provenance
    second = XYZFileParserMemory().parse(_XYZ_TEXT).parser_provenance

    assert first is not None
    assert second is not None
    assert first == second
    canonical_json = json.dumps(
        first.effective_config,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    assert first.effective_config_sha256 == sha256(canonical_json.encode("utf-8")).hexdigest()
    assert canonical_json_sha256({"z": 1, "a": {"y": 2, "x": 3}}) == (
        canonical_json_sha256({"a": {"x": 3, "y": 2}, "z": 1})
    )
    with pytest.raises(ValueError):
        canonical_json_sha256({"non_finite": float("nan")})
    with pytest.raises(ValidationError, match="does not match effective_config"):
        ParserProvenance(
            **(
                first.model_dump(exclude={"effective_config_sha256"})
                | {"effective_config_sha256": "0" * 64}
            )
        )


def test_parser_provenance_is_a_snapshot_of_global_configuration() -> None:
    original_force_unit_transform = molopconfig.force_unit_transform
    original_graph_backend = molopconfig.graph_reconstruction_backend
    original_make_dative_bonds = molopconfig.make_dative_bonds
    original_molgr_config = get_molgr_config()

    before = XYZFileParserMemory().parse(_XYZ_TEXT).parser_provenance
    assert before is not None
    before_config = deepcopy(before.effective_config)
    before_hash = before.effective_config_sha256

    try:
        molopconfig.force_unit_transform = not original_force_unit_transform
        molopconfig.graph_reconstruction_backend = (
            "python" if original_graph_backend == "cpp" else "cpp"
        )
        molopconfig.make_dative_bonds = not original_make_dative_bonds
        set_molgr_config(
            replace(
                original_molgr_config,
                resonance=replace(
                    original_molgr_config.resonance,
                    max_depth=original_molgr_config.resonance.max_depth + 1,
                ),
            )
        )

        after = XYZFileParserMemory().parse(_XYZ_TEXT).parser_provenance
        assert after is not None
        assert after.effective_config != before_config
        assert after.effective_config_sha256 != before_hash
        assert before.effective_config == before_config
        assert before.effective_config_sha256 == before_hash
    finally:
        molopconfig.force_unit_transform = original_force_unit_transform
        molopconfig.graph_reconstruction_backend = original_graph_backend
        molopconfig.make_dative_bonds = original_make_dative_bonds
        set_molgr_config(original_molgr_config)
