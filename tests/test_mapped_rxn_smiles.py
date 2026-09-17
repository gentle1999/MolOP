from __future__ import annotations

from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from molop.io.base_models.ChemFileFrame import BaseCalcFrame


def _molecule(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return molecule


def _atom_maps(smiles: str) -> list[int]:
    molecule = _molecule(smiles)
    return [atom.GetAtomMapNum() for atom in molecule.GetAtoms()]


def test_to_reaction_uses_one_based_source_atom_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    endpoints = (_molecule("C.C"), _molecule("CC"))
    monkeypatch.setattr(
        BaseCalcFrame,
        "possible_pre_post_ts",
        lambda self, **_kwargs: endpoints,
    )

    reaction = BaseCalcFrame().to_reaction()
    assert isinstance(reaction, rdChemReactions.ChemicalReaction)
    reactant_smiles, product_smiles = BaseCalcFrame().to_mapped_rxn_smiles().split(">>")

    assert _atom_maps(reactant_smiles) == [1, 2]
    assert _atom_maps(product_smiles) == [1, 2]
    assert [atom.GetAtomMapNum() for atom in endpoints[0].GetAtoms()] == [0, 0]
    assert [atom.GetAtomMapNum() for atom in endpoints[1].GetAtoms()] == [0, 0]


def test_to_mapped_rxn_smiles_exports_from_a_chemical_reaction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    endpoints = (_molecule("C.C"), _molecule("CC"))
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        BaseCalcFrame,
        "possible_pre_post_ts",
        lambda self, **_kwargs: endpoints,
    )

    def fake_reaction_to_smiles(
        reaction: rdChemReactions.ChemicalReaction,
        canonical: bool = True,
    ) -> str:
        captured["reaction"] = reaction
        captured["canonical"] = canonical
        assert reaction.GetNumReactantTemplates() == 2
        assert reaction.GetNumProductTemplates() == 1
        assert [
            atom.GetAtomMapNum()
            for template_index in range(reaction.GetNumReactantTemplates())
            for atom in reaction.GetReactantTemplate(template_index).GetAtoms()
        ] == [1, 2]
        return "exported-by-reaction"

    monkeypatch.setattr(rdChemReactions, "ReactionToSmiles", fake_reaction_to_smiles)

    assert BaseCalcFrame().to_mapped_rxn_smiles() == "exported-by-reaction"
    assert isinstance(captured["reaction"], rdChemReactions.ChemicalReaction)
    assert captured["canonical"] is False


def test_to_mapped_rxn_smiles_can_select_additional_endpoint_sampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = BaseCalcFrame()
    standard = (_molecule("C.C"), _molecule("CC"))
    additional = (_molecule("C.C"), _molecule("C=C"))
    calls: dict[str, Any] = {}

    def fake_possible(
        self: BaseCalcFrame,
        show_3D: bool = False,
        **_kwargs: Any,
    ) -> tuple[Chem.Mol, Chem.Mol]:
        assert self is frame
        calls["show_3D"] = show_3D
        return standard

    def fake_additional(
        self: BaseCalcFrame,
        pre: Chem.Mol,
        post: Chem.Mol,
        **_kwargs: Any,
    ) -> tuple[Chem.Mol, Chem.Mol]:
        assert self is frame
        assert (pre, post) == standard
        calls["additional"] = True
        return additional

    monkeypatch.setattr(BaseCalcFrame, "possible_pre_post_ts", fake_possible)
    monkeypatch.setattr(BaseCalcFrame, "additional_pre_post_ts", fake_additional)

    reaction_smiles = frame.to_mapped_rxn_smiles(additional=True)

    assert calls == {"show_3D": True, "additional": True}
    assert _atom_maps(reaction_smiles.split(">>")[0]) == [1, 2]
    assert _atom_maps(reaction_smiles.split(">>")[1]) == [1, 2]


def test_to_reaction_rejects_nonconserved_atom_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        BaseCalcFrame,
        "possible_pre_post_ts",
        lambda self, **_kwargs: (_molecule("CO"), _molecule("OC")),
    )

    with pytest.raises(ValueError, match="same atom order"):
        BaseCalcFrame().to_reaction()
