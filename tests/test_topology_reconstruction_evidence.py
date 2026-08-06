from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from molgr.config import CONFIG as MOLGR_CONFIG
from rdkit import Chem

from molop.config import molopconfig
from molop.io.base_models.Molecule import Molecule
from molop.io.base_models.source import canonical_json_sha256
from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserMemory
from molop.io.logic.orca.log.parsers.ORCALogFileParser import ORCALogFileParserMemory


FIXTURE_ROOT = Path(__file__).resolve().parent / "test_files"
CALCULATION_CASES = (
    (G16LogFileParserMemory, FIXTURE_ROOT / "g16log" / "H2O.log"),
    (
        ORCALogFileParserMemory,
        FIXTURE_ROOT / "orca" / "output_files" / "local" / "H2_sp_orca.out",
    ),
)
WATER_ATOMS = [8, 1, 1]
WATER_COORDS = np.asarray(
    [
        [0.0, 0.0, 0.0],
        [0.9572, 0.0, 0.0],
        [-0.2399872, 0.927297, 0.0],
    ]
)


@pytest.mark.parametrize(("parser_cls", "fixture"), CALCULATION_CASES)
def test_real_calc_frames_capture_trusted_topology_evidence(
    parser_cls: type[Any],
    fixture: Path,
) -> None:
    parsed = parser_cls(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]
    provenance = parsed.parser_provenance
    rdmol = frame.rdmol

    assert provenance is not None
    assert rdmol is not None
    assert [atom.GetAtomicNum() for atom in rdmol.GetAtoms()] == frame.atoms
    assert frame.source_to_topology_atom_permutation == list(range(len(frame.atoms)))
    assert frame.topology_reconstruction_status == "succeeded"
    assert (
        frame.topology_reconstruction_backend
        == (provenance.effective_config["molop"]["graph_reconstruction_backend"])
    )
    assert (
        frame.topology_make_dative_bonds
        == (provenance.effective_config["molop"]["make_dative_bonds"])
    )
    assert frame.topology_reconstruction_config_sha256 == canonical_json_sha256(
        {
            "backend": frame.topology_reconstruction_backend,
            "make_dative_bonds": frame.topology_make_dative_bonds,
            "molgr": provenance.effective_config["molgr"],
        }
    )
    assert frame.parse_presence["topology"] == "parsed"
    assert not any(
        diagnostic.code.startswith("MOL.PARSE.TOPOLOGY_") for diagnostic in frame.parse_diagnostics
    )

    payload = frame.to_unitless_dump_with_unit_keys(exclude_none=True)
    assert payload["source_to_topology_atom_permutation"] == list(range(len(frame.atoms)))
    assert payload["topology_reconstruction_backend"] == frame.topology_reconstruction_backend
    assert payload["topology_make_dative_bonds"] == frame.topology_make_dative_bonds
    assert payload["topology_reconstruction_config_sha256"] == (
        frame.topology_reconstruction_config_sha256
    )
    assert payload["topology_v3000_molblock"] == frame.topology_v3000_molblock
    assert "V3000" in payload["topology_v3000_molblock"]


@pytest.mark.parametrize("backend", ["cpp", "python"])
def test_molgr_backends_preserve_source_atom_order(
    backend: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", backend)
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)

    rdmol = molecule.rdmol

    assert rdmol is not None
    assert [atom.GetAtomicNum() for atom in rdmol.GetAtoms()] == WATER_ATOMS
    np.testing.assert_allclose(rdmol.GetConformer().GetPositions(), WATER_COORDS, atol=1e-5)
    assert molecule.source_to_topology_atom_permutation == [0, 1, 2]
    assert molecule.topology_reconstruction_backend == backend
    assert molecule.topology_reconstruction_status == "succeeded"


def test_molecule_reads_global_topology_config_at_reconstruction_time(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    captured: dict[str, Any] = {}

    def capture_config(
        xyz_block: str,
        _charge: int,
        _multiplicity: int,
        *,
        backend: str,
        make_dative_bonds: bool,
        config: Any,
    ) -> Chem.Mol:
        captured.update(
            backend=backend,
            make_dative_bonds=make_dative_bonds,
            config=config,
        )
        rdmol = Chem.MolFromXYZBlock(xyz_block)
        assert rdmol is not None
        return rdmol

    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    monkeypatch.setattr(molopconfig, "graph_reconstruction_backend", "python")
    monkeypatch.setattr(molopconfig, "make_dative_bonds", False)
    monkeypatch.setattr(
        MOLGR_CONFIG.resonance,
        "max_depth",
        MOLGR_CONFIG.resonance.max_depth + 1,
    )
    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", capture_config)

    assert molecule.rdmol is not None
    assert captured == {
        "backend": "python",
        "make_dative_bonds": False,
        "config": MOLGR_CONFIG,
    }
    assert molecule.topology_reconstruction_backend == "python"
    assert molecule.topology_make_dative_bonds is False


def test_topology_v3000_molblock_is_map_free_and_does_not_mutate_cached_rdmol() -> None:
    molecule = Molecule.from_coords(WATER_ATOMS, WATER_COORDS)
    rdmol = molecule.rdmol
    assert rdmol is not None
    for atom_index, atom in enumerate(rdmol.GetAtoms(), start=1):
        atom.SetAtomMapNum(atom_index)
    cached_maps = [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()]
    cached_conformer_count = rdmol.GetNumConformers()

    first = molecule.topology_v3000_molblock
    second = molecule.topology_v3000_molblock

    assert first is not None
    assert first == second
    assert "V3000" in first
    portable = Chem.MolFromMolBlock(
        first,
        sanitize=False,
        removeHs=False,
        strictParsing=True,
    )
    assert portable is not None
    assert [atom.GetAtomicNum() for atom in portable.GetAtoms()] == WATER_ATOMS
    assert [atom.GetAtomMapNum() for atom in portable.GetAtoms()] == [0, 0, 0]
    assert portable.GetNumBonds() == rdmol.GetNumBonds()
    assert rdmol.GetNumConformers() == cached_conformer_count
    assert [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()] == cached_maps


def test_calc_frame_topology_reconstruction_failure_is_structured(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")

    def fail_reconstruction(*_args: Any, **_kwargs: Any) -> Any:
        raise ValueError("synthetic reconstruction failure")

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", fail_reconstruction)
    fixture = FIXTURE_ROOT / "g16log" / "H2O.log"
    parsed = G16LogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]

    assert frame.rdmol is None
    assert frame.topology_v3000_molblock is None
    assert frame.topology_reconstruction_status == "failed"
    assert frame.topology_reconstruction_backend is not None
    assert frame.topology_reconstruction_config_sha256 is not None
    assert frame.source_to_topology_atom_permutation is None
    assert frame.parse_presence["topology"] == "parse_failed"
    assert frame.parse_completeness == "partial"
    assert parsed.parse_completeness == "partial"
    assert any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_RECONSTRUCTION_FAILED"
        for diagnostic in frame.parse_diagnostics
    )


def test_ambiguous_same_element_reordering_does_not_guess_permutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule_module = importlib.import_module("molop.io.base_models.Molecule")

    def reorder_equivalent_atoms(xyz_block: str, *_args: Any, **_kwargs: Any) -> Chem.Mol:
        rdmol = Chem.MolFromXYZBlock(xyz_block)
        assert rdmol is not None
        return Chem.RenumberAtoms(rdmol, [0, 2, 1])

    monkeypatch.setattr(molecule_module, "xyz_to_rdmol", reorder_equivalent_atoms)
    fixture = FIXTURE_ROOT / "g16log" / "H2O.log"
    parsed = G16LogFileParserMemory(
        capture_source_evidence=True,
        only_last_frame=True,
    ).parse(fixture.read_text(encoding="utf-8"))
    frame = parsed[-1]

    assert frame.rdmol is not None
    assert [atom.GetAtomicNum() for atom in frame.rdmol.GetAtoms()] == frame.atoms
    assert frame.source_to_topology_atom_permutation is None
    assert frame.topology_v3000_molblock is None
    assert frame.bonds == []
    assert frame.formal_charges == []
    assert frame.formal_num_radicals == []
    assert frame.parse_presence["topology"] == "parse_failed"
    assert frame.parse_completeness == "partial"
    assert parsed.parse_completeness == "partial"
    assert any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_ATOM_ORDER_MISMATCH"
        for diagnostic in frame.parse_diagnostics
    )
    assert not any(
        diagnostic.code == "MOL.PARSE.TOPOLOGY_RECONSTRUCTION_FAILED"
        for diagnostic in frame.parse_diagnostics
    )
