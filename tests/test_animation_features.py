from __future__ import annotations

from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest
from PIL import Image
from rdkit import Chem

from molop import AutoParser
from molop.config import molopconfig
from molop.io.base_models.ChemFile import BaseCalcFile, BaseCoordsFile
from molop.io.base_models.ChemFileFrame import BaseCalcFrame, BaseCoordsFrame
from molop.io.base_models.DataClasses import Energies, Status, Vibrations
from molop.unit import atom_ureg


TS_LOG = (
    Path(__file__).parent
    / "test_files"
    / "g16log"
    / "000000000000_000016928457_00_conf_01_ts.107c60f3cfcb.log"
)


@pytest.fixture(scope="module")
def ts_file() -> BaseCalcFile[BaseCalcFrame]:
    return cast(BaseCalcFile[BaseCalcFrame], AutoParser(TS_LOG, n_jobs=1)[0])


@pytest.fixture(scope="module")
def ts_frame(ts_file: BaseCalcFile[BaseCalcFrame]) -> BaseCalcFrame:
    frame = next((candidate for candidate in ts_file if candidate.is_TS), None)
    assert frame is not None
    return frame


def _water_frame(z_shift: float) -> BaseCoordsFrame:
    return BaseCoordsFrame(
        atoms=[8, 1, 1],
        coords=np.array(
            [
                [0.0, 0.0, z_shift],
                [0.0, 0.76, z_shift + 0.59],
                [0.0, -0.76, z_shift + 0.59],
            ]
        )
        * atom_ureg.angstrom,
        bonds=[(0, 1, 1, 0), (0, 2, 1, 0)],
        formal_charges=[0, 0, 0],
        formal_num_radicals=[0, 0, 0],
    )


def test_chem_file_draw_animation_renders_all_frames(tmp_path: Path) -> None:
    trajectory = BaseCoordsFile()
    trajectory.append(_water_frame(0.0))
    trajectory.append(BaseCoordsFrame())
    trajectory.append(_water_frame(0.2))

    gif_path = tmp_path / "trajectory.gif"
    gif_data = trajectory.draw_animation(
        file_path=gif_path,
        size=(180, 140),
        duration=80,
        return_image=False,
    )
    assert isinstance(gif_data, bytes)
    assert gif_data.startswith(b"GIF8")
    assert gif_path.read_bytes() == gif_data
    with Image.open(BytesIO(gif_data)) as image:
        assert image.n_frames == 2

    svg_text = trajectory.draw_animation(
        image_format="svg",
        size=(180, 140),
        return_image=False,
    )
    assert isinstance(svg_text, str)
    assert "<animate" in svg_text


def test_chem_file_animation_default_legends_follow_renderable_frames(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import rdkit_dof

    trajectory = BaseCoordsFile()
    trajectory.append(_water_frame(0.0))
    trajectory.append(BaseCoordsFrame())
    trajectory.append(_water_frame(0.2))
    captured: dict[str, object] = {}

    def fake_drawer(mols: object, **kwargs: object) -> bytes:
        captured["mols"] = mols
        captured["legends"] = kwargs["legends"]
        return b"GIF89a"

    monkeypatch.setattr(rdkit_dof, "MolsToDofGif", fake_drawer)

    assert trajectory.draw_animation(return_image=False) == b"GIF89a"
    rendered_mols = captured["mols"]
    assert isinstance(rendered_mols, list)
    assert len(rendered_mols) == 2
    assert captured["legends"] == ["Frame 0", "Frame 2"]


def test_chem_file_animation_filters_all_frame_aligned_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import rdkit_dof

    trajectory = BaseCoordsFile()
    trajectory.append(_water_frame(0.0))
    trajectory.append(BaseCoordsFrame())
    trajectory.append(_water_frame(0.2))
    captured: dict[str, object] = {}

    def fake_drawer(mols: object, **kwargs: object) -> bytes:
        captured["mols"] = mols
        captured.update(kwargs)
        return b"GIF89a"

    monkeypatch.setattr(rdkit_dof, "MolsToDofGif", fake_drawer)

    result = trajectory.draw_animation(
        duration=[60, 70, 80],
        legends=["first", "skipped", "last"],
        highlightAtomLists=[[0], [1], [2]],
        highlightBondLists=[[0], [1], []],
        return_image=False,
    )

    assert result == b"GIF89a"
    assert captured["duration"] == [60, 80]
    assert captured["legends"] == ["first", "last"]
    assert captured["highlightAtomLists"] == [[0], [2]]
    assert captured["highlightBondLists"] == [[0], []]


def test_ts_vibration_uses_the_unique_imaginary_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = BaseCalcFrame(
        atoms=[8, 1, 1],
        coords=np.array([[0.0, 0.0, 0.0], [0.0, 0.76, 0.59], [0.0, -0.76, 0.59]])
        * atom_ureg.angstrom,
        energies=Energies(electronic_energy=-1.0 * atom_ureg.hartree),
        status=Status(normal_terminated=True),
        vibrations=Vibrations(
            frequencies=np.array([100.0, -50.0, 200.0]) * atom_ureg.cm**-1,
            vibration_modes=[
                np.zeros((3, 3)) * atom_ureg.angstrom,
                np.ones((3, 3)) * atom_ureg.angstrom,
                np.zeros((3, 3)) * atom_ureg.angstrom,
            ],
        ),
    )
    captured: dict[str, int] = {}

    def fake_vibrate(
        self: BaseCalcFrame,
        vibration_id: int | None = None,
        **_kwargs: object,
    ) -> list[object]:
        captured["vibration_id"] = vibration_id if vibration_id is not None else -1
        return []

    monkeypatch.setattr(BaseCalcFrame, "vibrate", fake_vibrate)

    assert frame.is_TS is True
    assert frame.ts_vibration() == []
    assert captured["vibration_id"] == 1


def _marked_topology(smiles: str, marker: float) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for atom_index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(atom_index, (marker, float(atom_index), 0.0))
    molecule.AddConformer(conformer)
    return molecule


def test_possible_pre_post_ts_votes_on_each_displacement_side(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = BaseCalcFrame()
    topology_by_ratio = {
        -0.75: "CC",
        -1.0: "CC",
        -1.25: "C=C",
        -1.5: "CC",
        -1.75: "C=C",
        0.75: "C.C",
        1.0: "CC",
        1.25: "C.C",
        1.5: "CC",
        1.75: "C.C",
    }
    sampled_ratios: list[float] = []

    def fake_ts_vibration(
        self: BaseCalcFrame, ratio: float = 1.75, steps: int = 7
    ) -> list[SimpleNamespace]:
        assert self is frame
        assert steps == 1
        sampled_ratios.append(ratio)
        return [SimpleNamespace(rdmol=_marked_topology(topology_by_ratio[ratio], ratio))]

    monkeypatch.setattr(BaseCalcFrame, "ts_vibration", fake_ts_vibration)

    pre, post = frame.possible_pre_post_ts(show_3D=True, min_ratio=0.75, max_ratio=1.75, steps=5)

    assert sampled_ratios == [
        -0.75,
        -1.0,
        -1.25,
        -1.5,
        -1.75,
        0.75,
        1.0,
        1.25,
        1.5,
        1.75,
    ]
    assert Chem.MolToSmiles(pre) == "C.C"
    assert Chem.MolToSmiles(post) == "CC"
    assert pre.GetConformer().GetAtomPosition(0).x == pytest.approx(1.75)
    assert post.GetConformer().GetAtomPosition(0).x == pytest.approx(-1.5)


def test_additional_pre_post_ts_returns_endpoints_without_bond_changes() -> None:
    frame = BaseCalcFrame(
        atoms=[6, 6],
        coords=np.zeros((2, 3)) * atom_ureg.angstrom,
    )
    pre = Chem.MolFromSmiles("CC")
    post = Chem.MolFromSmiles("CC")
    assert pre is not None
    assert post is not None

    additional_pre, additional_post = frame.additional_pre_post_ts(pre, post)

    assert additional_pre is pre
    assert additional_post is post


def test_additional_pre_post_ts_fixes_bond_change_atoms_during_resampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import importlib

    frame = BaseCalcFrame(
        atoms=[6, 6, 6],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [0.0, 6.0, 0.0],
            ]
        )
        * atom_ureg.angstrom,
        vibrations=Vibrations(
            frequencies=np.array([-100.0, 100.0, 200.0]) * atom_ureg.cm**-1,
            vibration_modes=[
                np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                    ]
                )
                * atom_ureg.angstrom,
                np.zeros((3, 3)) * atom_ureg.angstrom,
                np.zeros((3, 3)) * atom_ureg.angstrom,
            ],
        ),
    )
    pre = _marked_topology("C.C.C", -1.8)
    post = _marked_topology("CC.C", 1.8)
    captured_coords: list[np.ndarray] = []
    frame_module = importlib.import_module("molop.io.base_models.ChemFileFrame")

    def fake_from_coords(*, coords: np.ndarray, **_kwargs: object) -> SimpleNamespace:
        captured_coords.append(np.array(coords, copy=True))
        return SimpleNamespace(rdmol=_marked_topology("C.C.C", float(coords[0, 0])))

    monkeypatch.setattr(frame_module, "check_crowding", lambda _rdmol: True)
    monkeypatch.setattr(frame_module.Molecule, "from_coords", staticmethod(fake_from_coords))
    monkeypatch.setattr(BaseCalcFrame, "_ts_vibration_id", lambda _self: 0)

    frame.additional_pre_post_ts(pre, post)

    assert len(captured_coords) == 18
    amplitudes = np.linspace(0.2, 1.8, num=9, endpoint=True)
    for direction, endpoint, side_coords in (
        (-1.0, pre, captured_coords[:9]),
        (1.0, post, captured_coords[9:]),
    ):
        endpoint_coords = endpoint.GetConformer().GetPositions()
        for coords, amplitude in zip(side_coords, amplitudes, strict=True):
            np.testing.assert_allclose(coords[:2], endpoint_coords[:2])
            np.testing.assert_allclose(coords[2], [direction * amplitude, 6.0, 0.0])


def test_possible_pre_post_ts_removes_conformers_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = BaseCalcFrame()

    def fake_ts_vibration(
        _self: BaseCalcFrame, ratio: float = 1.75, steps: int = 7
    ) -> list[SimpleNamespace]:
        assert steps == 1
        return [SimpleNamespace(rdmol=_marked_topology("CC", ratio))]

    monkeypatch.setattr(BaseCalcFrame, "ts_vibration", fake_ts_vibration)

    pre, post = frame.possible_pre_post_ts(steps=1)

    assert pre.GetNumConformers() == 0
    assert post.GetNumConformers() == 0


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"min_ratio": 0.0}, "min_ratio"),
        ({"min_ratio": float("nan")}, "min_ratio"),
        ({"max_ratio": 0.1}, "max_ratio"),
        ({"max_ratio": float("inf")}, "max_ratio"),
        ({"steps": 0}, "steps"),
    ],
)
def test_possible_pre_post_ts_rejects_invalid_sampling_ranges(
    kwargs: dict[str, float | int], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        BaseCalcFrame().possible_pre_post_ts(**kwargs)  # type: ignore[arg-type]


def _assert_sdf_endpoint(path: Path, atom_count: int) -> None:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecules = [molecule for molecule in supplier if molecule is not None]
    assert len(molecules) == 1
    assert molecules[0].GetNumAtoms() == atom_count
    assert molecules[0].GetNumConformers() == 1


def test_ts_animation_and_endpoint_export(ts_frame: BaseCalcFrame, tmp_path: Path) -> None:
    vibration_gif = ts_frame.draw_vibration_animation(
        vibration_id=0,
        steps=3,
        size=(180, 140),
        return_image=False,
    )
    ts_gif = ts_frame.draw_ts_vibration_animation(
        steps=3,
        size=(180, 140),
        return_image=False,
    )

    assert isinstance(vibration_gif, bytes)
    assert vibration_gif.startswith(b"GIF8")
    assert isinstance(ts_gif, bytes)
    assert ts_gif.startswith(b"GIF8")

    pre_path, post_path = ts_frame.save_pre_post_ts(tmp_path, prefix="candidate")
    assert pre_path.name == "candidate_pre.xyz"
    assert post_path.name == "candidate_post.xyz"
    for path in (pre_path, post_path):
        lines = path.read_text(encoding="utf-8").splitlines()
        assert int(lines[0]) == ts_frame.rdmol.GetNumAtoms()  # type: ignore[union-attr]
        assert len(lines) == int(lines[0]) + 2

    default_pre, default_post = ts_frame.save_pre_post_ts(tmp_path / "source-named")
    source_prefix = Path(ts_frame.filename).stem  # type: ignore[attr-defined]
    assert default_pre.name == f"{source_prefix}_frame_{ts_frame.frame_id:03d}_pre.xyz"
    assert default_post.name == f"{source_prefix}_frame_{ts_frame.frame_id:03d}_post.xyz"

    pre_sdf, post_sdf = ts_frame.save_pre_post_ts(tmp_path, prefix="candidate", format="sdf")
    assert pre_sdf.name == "candidate_pre.sdf"
    assert post_sdf.name == "candidate_post.sdf"
    for path in (pre_sdf, post_sdf):
        _assert_sdf_endpoint(path, ts_frame.rdmol.GetNumAtoms())  # type: ignore[union-attr]


def test_file_and_batch_endpoint_export(
    ts_file: BaseCalcFile[BaseCalcFrame], tmp_path: Path
) -> None:
    expected_frame_ids = {frame.frame_id for frame in ts_file if frame.is_TS}
    file_exports = ts_file.save_pre_post_ts(tmp_path / "single-file", format="sdf")
    assert set(file_exports) == expected_frame_ids
    for paths in file_exports.values():
        for path in paths:
            _assert_sdf_endpoint(path, ts_file[0].rdmol.GetNumAtoms())  # type: ignore[union-attr]

    first_copy = tmp_path / "first" / "source.ts.log"
    second_copy = tmp_path / "second" / "source.ts.log"
    first_copy.parent.mkdir()
    second_copy.parent.mkdir()
    first_copy.write_bytes(TS_LOG.read_bytes())
    second_copy.write_bytes(TS_LOG.read_bytes())
    batch = AutoParser([first_copy, second_copy], n_jobs=1)
    batch_exports = batch.save_pre_post_ts(tmp_path / "batch", format="sdf", n_jobs=1)

    assert set(batch_exports) == {str(first_copy), str(second_copy)}
    output_directories = {
        paths[0].parent for exports in batch_exports.values() for paths in exports.values()
    }
    assert len(output_directories) == 2
    for source_path, exports in batch_exports.items():
        assert set(exports) == expected_frame_ids
        source_directories = {paths[0].parent for paths in exports.values()}
        assert len(source_directories) == 1
        source_directory = source_directories.pop()
        assert source_directory.parent == tmp_path / "batch"
        assert source_directory.name.startswith(f"{Path(source_path).stem}-")
        for paths in exports.values():
            for path in paths:
                assert path.parent == source_directory
                _assert_sdf_endpoint(path, ts_file[0].rdmol.GetNumAtoms())  # type: ignore[union-attr]


def test_batch_endpoint_export_reconstructs_dynamic_candidates_in_task(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import importlib

    first_copy = tmp_path / "first" / "source.ts.log"
    second_copy = tmp_path / "second" / "source.ts.log"
    first_copy.parent.mkdir()
    second_copy.parent.mkdir()
    first_copy.write_bytes(TS_LOG.read_bytes())
    second_copy.write_bytes(TS_LOG.read_bytes())
    batch = AutoParser([first_copy, second_copy], n_jobs=1)

    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    original_iterator = molecule_module.iter_xyz_to_rdmol_batch
    calls: list[int] = []

    def counting_iterator(requests: Any, **kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append(len(request_list))
        return original_iterator(request_list, **kwargs)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", counting_iterator)

    exports = batch.save_pre_post_ts(tmp_path / "batch", format="sdf", n_jobs=1)

    assert set(exports) == {str(first_copy), str(second_copy)}
    # TS endpoint candidates are generated inside each file task.  They are
    # intentionally not collected and sent through the native batch helper.
    assert calls == []


def test_batch_summary_reconstructs_dynamic_ts_candidates_in_workers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(molopconfig, "prewarm_topologies", True)
    import importlib

    first_copy = tmp_path / "first" / "source.ts.log"
    second_copy = tmp_path / "second" / "source.ts.log"
    first_copy.parent.mkdir()
    second_copy.parent.mkdir()
    first_copy.write_bytes(TS_LOG.read_bytes())
    second_copy.write_bytes(TS_LOG.read_bytes())
    batch = AutoParser([first_copy, second_copy], n_jobs=1)

    molecule_module = importlib.import_module("molop.io.base_models.Molecule")
    original_iterator = molecule_module.iter_xyz_to_rdmol_batch
    calls: list[int] = []

    def counting_iterator(requests: Any, **kwargs: Any) -> Any:
        request_list = list(requests)
        calls.append(len(request_list))
        return original_iterator(request_list, **kwargs)

    monkeypatch.setattr(molecule_module, "iter_xyz_to_rdmol_batch", counting_iterator)

    summary = batch.to_summary_df(frame="all", n_jobs=2)

    assert len(summary) > 0
    # The parent only batches the selected source frames.  Temporary TS
    # candidates are generated and reconstructed in the fresh loky workers.
    assert calls == [68]
    ts_rows = summary[summary[("Status", "IsTS", "")].eq(True)]
    assert len(ts_rows) == 2
    assert ts_rows[("General", "PreCanonicalSMILES", "")].notna().all()
    assert ts_rows[("General", "PostCanonicalSMILES", "")].notna().all()


def test_batch_endpoint_export_skips_unsupported_files(tmp_path: Path) -> None:
    ts_copy = tmp_path / "candidate.log"
    xyz_path = tmp_path / "water.xyz"
    ts_copy.write_bytes(TS_LOG.read_bytes())
    xyz_path.write_text(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n",
        encoding="utf-8",
    )

    batch = AutoParser([ts_copy, xyz_path], n_jobs=1)
    exports = batch.save_pre_post_ts(tmp_path / "mixed", format="sdf", n_jobs=1)

    assert set(exports) == {str(ts_copy)}
    assert exports[str(ts_copy)]
