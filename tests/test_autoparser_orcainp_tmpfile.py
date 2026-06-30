from pathlib import Path
from typing import cast

import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.base_models.DataClasses import CoordinateContainer, CoordinateParameters
from molop.io.codec_exceptions import UnsupportedFormatError
from molop.io.codec_registry import get_supported_writer_formats
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import ORCAInpFileFrameDisk


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
    typed_frame = cast(ORCAInpFileFrameDisk, frame)
    assert typed_frame.keyword_lines[0].text == "SP HF def2-SVP"
    assert typed_frame.geometry is not None
    assert typed_frame.geometry.ctype == "xyz"
    assert isinstance(typed_frame.geometry, CoordinateContainer)
    assert len(typed_frame.geometry) == 2
    assert typed_frame.geometry.get_symbols() == ["H", "H"]


def test_autoparser_orcainp_point_charges_parse_into_frame_geometry(tmp_path: Path) -> None:
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
    frame = cast(ORCAInpFileFrameDisk, batch[0][-1])

    assert frame.geometry is not None
    assert len(frame.geometry.atoms) == 2
    assert frame.atoms == [1, 1]
    assert frame.geometry.point_charges == [
        {"charge": 0.5, "x": 0.0, "y": 0.0, "z": 1.0}
    ]


def test_orcainp_writer_is_not_registered() -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    assert "orcainp" not in get_supported_writer_formats()
    with pytest.raises(UnsupportedFormatError):
        file_model.format_transform("orcainp", graph_policy="prefer")


def test_autoparser_orcainp_new_job_splits_frames(tmp_path: Path) -> None:
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

    first_frame = cast(ORCAInpFileFrameDisk, file_model[0])
    second_frame = cast(ORCAInpFileFrameDisk, file_model[1])
    assert first_frame.frame_id == 0
    assert second_frame.frame_id == 1
    assert first_frame.geometry is not None
    assert second_frame.geometry is not None
    assert first_frame.geometry.atoms[1].x == 0.7
    assert second_frame.geometry.atoms[1].x == 0.8


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
    typed_frame = cast(ORCAInpFileFrameDisk, frame)
    assert typed_frame.geometry is not None
    assert typed_frame.geometry.ctype == "xyzfile"
    assert typed_frame.geometry.external_path == "does_not_exist.xyz"


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

    assert file_model.model_chemistry.functional == "B3LYP-D3BJ"
    assert file_model.task_requests[0].task_type == "sp"
    assert file_model.resource_request.num_cpu == 4

    assert typed_frame.qm_software == "ORCA"
    assert typed_frame.qm_software_version == "Any"
    assert typed_frame.method == "DFT"
    assert typed_frame.functional == "B3LYP-D3BJ"
    assert typed_frame.dispersion_correction == "D3BJ"
    assert typed_frame.basis_set == "def2-TZVP"
    assert typed_frame.keywords == "B3LYP D3BJ def2-TZVP"
    assert typed_frame.model_chemistry.method_family == "DFT"
    assert typed_frame.model_chemistry.functional == "B3LYP-D3BJ"
    assert typed_frame.model_chemistry.dispersion_correction == "D3BJ"
    assert typed_frame.model_chemistry.basis_set == "def2-TZVP"
    assert [basis.name for basis in typed_frame.model_chemistry.basis_sets] == ["def2-TZVP"]
    assert [task.task_type for task in typed_frame.task_requests] == ["sp"]
    assert typed_frame.resource_request.num_cpu == 4
    assert typed_frame.resource_request.memory is not None

    assert [block.name for block in typed_frame.blocks] == ["pal", "maxcore"]
    assert typed_frame.blocks[0].lines[0].text.strip() == "nprocs 4"
    assert typed_frame.blocks[1].lines[0].text.strip() == "2000"
    assert typed_frame.resources_raw == "%pal\n  nprocs 4\nend\n%maxcore 2000"
    assert typed_frame.request_num_cpu == 4
    assert typed_frame.request_memory is not None


def test_autoparser_orcainp_metadata_population_recognizes_d4(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! PBE0 RIJCOSX D4 def2-SVP
* xyz 0 1
H 0.0 0.0 0.0
H 0.7 0.0 0.0
*
""",
    )
    batch = AutoParser(str(path), parser_detection="orcainp")
    typed_frame = cast(ORCAInpFileFrameDisk, batch[0][-1])

    assert typed_frame.functional == "PBE0-D4"
    assert typed_frame.dispersion_correction == "D4"
    assert typed_frame.basis_set == "def2-SVP"


def test_autoparser_orcainp_mixed_basis_and_output_print_settings(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! OPBE opt def2-SV(P) def2-SVP/C TightSCF UNO

%output
Print[ P_UNO_OccNum ] = 1
end

* xyz 2 5
  Fe    -0.000094   -0.001538    0.015314  newgto "def2-TZVP" end
  N     -1.301044    0.998903    1.550485  newgto "def2-TZVP" end
  C      2.460016   -0.340681    1.821585
*
""",
        name="phen3.inp",
    )
    batch = AutoParser(str(path), parser_detection="orcainp")
    typed_frame = cast(ORCAInpFileFrameDisk, batch[0][-1])

    assert typed_frame.charge == 2
    assert typed_frame.multiplicity == 5
    assert typed_frame.method == "DFT"
    assert typed_frame.functional == "OPBE"
    assert typed_frame.basis_set == "def2-SV(P)"
    assert typed_frame.auxiliary_basis_set == "def2-SVP/C"
    assert typed_frame.has_mixed_basis is True
    assert typed_frame.model_chemistry.method_family == "DFT"
    assert typed_frame.model_chemistry.functional == "OPBE"
    assert typed_frame.model_chemistry.auxiliary_basis_set == "def2-SVP/C"
    assert typed_frame.model_chemistry.options["has_mixed_basis"] is True
    atom_basis_sets = [
        basis
        for basis in typed_frame.model_chemistry.basis_sets
        if basis.scope == "atom" and basis.role == "orbital"
    ]
    assert [basis.name for basis in atom_basis_sets] == ["def2-TZVP", "def2-TZVP"]
    assert [basis.atom_indices for basis in atom_basis_sets] == [[0], [1]]
    assert [task.task_type for task in typed_frame.task_requests] == ["opt"]
    assert typed_frame.geometry is not None
    assert typed_frame.atoms == [26, 7, 6]
    assert tuple(typed_frame.coords.shape) == (3, 3)
    assert typed_frame.geometry.atoms[0].basis_set == "def2-TZVP"
    assert typed_frame.geometry.atoms[1].basis_overrides[0].kind == "newgto"
    assert typed_frame.geometry.atoms[2].basis_overrides == []
    assert typed_frame.output_print_settings[0].target == "P_UNO_OccNum"
    assert typed_frame.output_print_settings[0].value == "1"
    assert typed_frame.blocks[0].name == "output"
    assert typed_frame.resources_raw == "%output\nPrint[ P_UNO_OccNum ] = 1\nend"


def test_autoparser_orcainp_mrci_newblock_is_not_truncated(tmp_path: Path) -> None:
    path = _write_orcainp(
        tmp_path,
        """! ano-pVDZ TightSCF

%casscf
 nel     7
 norb    6
 nroots  3
 mult    2
end

%mrci
 tsel       0
 tpre       0
 newblock 2 *
  nroots 3
  excitations none
  refs
   cas(7,6)
  end
 end
end

* int 1 2
 O     0   0   0   0.000000     0.000     0.000
 H     1   0   0   1.012277     0.000     0.000
 H     1   2   0   1.012177   109.288     0.000
end
""",
    )
    typed_frame = cast(
        ORCAInpFileFrameDisk, AutoParser(str(path), parser_detection="orcainp")[0][-1]
    )

    assert typed_frame.method == "MRCI"
    assert typed_frame.multi_reference_semantic.enabled is True
    assert len(typed_frame.multi_reference_semantic.new_blocks) == 1
    assert typed_frame.multi_reference_semantic.new_blocks[0].refs == "cas(7,6)"
    assert [task.task_type for task in typed_frame.task_requests] == ["multi_reference"]
    assert len(typed_frame.multireference_requests) == 1
    request = typed_frame.multireference_requests[0]
    assert request.enabled is True
    assert request.method == "MRCI"
    assert request.active_space is not None
    assert request.active_space.electrons == 7
    assert request.active_space.orbitals == 6
    assert len(request.state_blocks) == 1
    assert request.state_blocks[0].active_space is not None
    assert request.state_blocks[0].active_space.electrons == 7
    assert typed_frame.resources_raw.endswith("  end\n end\nend")


def test_autoparser_orcainp_paras_structures_scan_and_resolves_cartesian_variables(
    tmp_path: Path,
) -> None:
    path = _write_orcainp(
        tmp_path,
        """! ano-pVDZ VeryTightSCF NoPop Conv MRCI+Q

%paras  R = 0.85,1.1,7
        end

* xyz 0 1
F  0 0 0
H  0 0 {R}
*
""",
    )
    typed_frame = cast(
        ORCAInpFileFrameDisk, AutoParser(str(path), parser_detection="orcainp")[0][-1]
    )

    assert typed_frame.geometry is not None
    assert typed_frame.geometry.ctype == "xyz"
    assert isinstance(typed_frame.geometry.coordinate_parameters, CoordinateParameters)
    assert len(typed_frame.geometry.coordinate_parameters) == 1
    parameter = typed_frame.geometry.coordinate_parameters[0]
    assert parameter.name == "R"
    assert parameter.start == 0.85
    assert parameter.stop == 1.1
    assert parameter.steps == 7
    assert parameter.is_scan is True
    assert typed_frame.atoms == [9, 1]
    assert tuple(typed_frame.coords.shape) == (2, 3)
    assert typed_frame.coords.magnitude[1, 2] == pytest.approx(0.85)
