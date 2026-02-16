from pathlib import Path
from typing import cast

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import ORCAInpFileFrameDisk


def _ensure_str(rendered: str | list[str]) -> str:
    assert isinstance(rendered, str)
    assert rendered.strip()
    return rendered


def _write_orcainp(tmp_path: Path, content: str, name: str = "input.inp") -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_autoparser_orcainp_minimal_xyz_parses_atoms_and_coords(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! SP HF def2-SVP
* xyz 0 1
H 0.0 0.0 0.0
H 0.7 0.0 0.0
*
""",
    )
    batch = AutoParser(str(path))
    assert len(batch) > 0
    file_model = batch[0]
    assert len(file_model) > 0

    frame = file_model[-1]
    atoms = getattr(frame, "atoms", None)
    coords = getattr(frame, "coords", None)
    assert atoms is not None and len(atoms) >= 2
    assert coords is not None
    assert getattr(coords, "shape", None) is not None
    assert tuple(coords.shape) == (len(atoms), 3)


def test_autoparser_orcainp_raw_first_render_preserves_header_and_q_line(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! SP
* xyz 0 1
H 0.0 0.0 0.0
Q 0.50 0.000000 0.000000 1.000000
H 0.8 0.0 0.0
*
""",
    )
    batch = AutoParser(str(path))
    file_model = batch[0]

    rendered = _ensure_str(file_model.format_transform("orcainp"))
    assert "* xyz 0 1" in rendered
    assert "Q 0.50 0.000000 0.000000 1.000000" in rendered

    # Canonical ordering: atoms first, then Q lines
    lines = [line.strip() for line in rendered.splitlines()]
    atom_indices = [i for i, line in enumerate(lines) if line.startswith("H ")]
    q_indices = [i for i, line in enumerate(lines) if line.startswith("Q ")]

    assert atom_indices, "No atom lines found in rendered output"
    assert q_indices, "No Q lines found in rendered output"
    assert max(atom_indices) < min(q_indices), "Q lines must appear after all atom lines"


def test_autoparser_xyz_to_orcainp_structured_render_has_xyz_block() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    rendered = _ensure_str(file_model.format_transform("orcainp", graph_policy="prefer"))
    rendered_lines = [line.strip() for line in rendered.splitlines() if line.strip()]
    assert rendered_lines
    assert any(line.startswith("* xyz") for line in rendered_lines)
    assert rendered_lines[-1] == "*"


def test_autoparser_xyz_to_orcainp_render_overrides_include_keywords_and_resources() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    rendered = _ensure_str(
        file_model.format_transform(
            "orcainp",
            keywords="B3LYP def2-SVP",
            nprocs=4,
            maxcore=2000,
        )
    )
    assert "! B3LYP def2-SVP" in rendered
    assert "%pal" in rendered
    assert "nprocs 4" in rendered
    assert "%maxcore 2000" in rendered
    assert "* xyz" in rendered


def test_autoparser_orcainp_new_job_splits_frames_and_renders_delimiter(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! SP
* xyz 0 1
H 0.0 0.0 0.0
H 0.7 0.0 0.0
*
$new_job
! SP
* xyz 0 1
H 0.0 0.0 0.0
H 0.8 0.0 0.0
*
""",
    )
    batch = AutoParser(str(path))
    file_model = batch[0]
    assert len(file_model) == 2

    rendered = _ensure_str(file_model.format_transform("orcainp", frameID="all"))
    assert "\n$new_job\n" in rendered


def test_autoparser_orcainp_xyzfile_does_not_require_external_file(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! SP
* xyzfile 0 1 does_not_exist.xyz
""",
    )
    batch = AutoParser(str(path))
    assert len(batch) > 0
    file_model = batch[0]
    assert len(file_model) > 0

    frame = file_model[-1]
    atoms = getattr(frame, "atoms", None)
    assert atoms is not None
    assert len(atoms) == 0

    rendered = _ensure_str(file_model.format_transform("orcainp"))
    assert "* xyzfile 0 1 does_not_exist.xyz" in rendered


def test_autoparser_orcainp_metadata_population(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! B3LYP D3BJ def2-TZVP
%pal
  nprocs 4
end
%maxcore 2000
* xyz 0 1
H 0.0 0.0 0.0
H 0.7 0.0 0.0
*
""",
    )
    batch = AutoParser(str(path), parser_detection="orcainp")
    file_model = batch[0]
    frame = file_model[-1]
    typed_frame = cast(ORCAInpFileFrameDisk, frame)

    assert typed_frame.qm_software == "ORCA"
    assert typed_frame.method == "DFT"
    assert typed_frame.functional == "B3LYP"
    assert typed_frame.basis_set == "def2-TZVP"

    assert "%pal" in typed_frame.resources_raw
    assert "nprocs 4" in typed_frame.resources_raw
    assert "end" in typed_frame.resources_raw
    assert "%maxcore 2000" in typed_frame.resources_raw
