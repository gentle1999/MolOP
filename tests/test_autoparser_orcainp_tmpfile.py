from pathlib import Path
from typing import cast

import numpy as np
import pytest

from molop.io import AutoParser  # type: ignore[reportMissingImports]
from molop.io.base_models.DataClasses import CoordinateContainer, CoordinateParameters
from molop.io.base_models.FrameParser import FrameParseContext
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.codec_registry import get_supported_writer_formats
from molop.io.logic.coords.frame_models.XYZFileFrame import XYZFileFrameMemory
from molop.io.logic.orca.input.frame_models.ORCAInpFileFrame import ORCAInpFileFrameDisk
from molop.io.logic.orca.input.frame_parsers.ORCAInpFileFrameParser import (
    ORCAInpFileFrameParserMemory,
    parse_orca_input_frame_result,
)
from molop.io.logic.orca.input.parsers._orca_inp_metadata import parse_orca_input_metadata
from molop.io.logic.orca.input.parsers.ORCAInpFileParser import ORCAInpFileParserMemory
from molop.unit import atom_ureg


def _write_orcainp(tmp_path: Path, content: str, name: str = "input.inp") -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_orcainp_file_parser_artifact_metadata_is_model_ready() -> None:
    parser = ORCAInpFileParserMemory()
    metadata = parser._parse_artifact_metadata("! SP HF def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n")

    assert metadata == {
        "qm_software": "ORCA",
        "qm_software_version": "Any",
    }


def test_orcainp_frame_parser_returns_model_ready_semantics() -> None:
    parser = ORCAInpFileFrameParserMemory()
    block = """! B3LYP D3BJ def2-TZVP
%pal
  nprocs 4
end
%maxcore 2000
%output
Print[ P_UNO_OccNum ] = 1
end
* xyz 0 1
H 0.0 0.0 0.0
H 0.7 0.0 0.0
*
"""
    result = parse_orca_input_frame_result(block)
    assert isinstance(result, ModelParseResult)
    assert result.has_value("model_chemistry") is True
    metadata = parse_orca_input_metadata(block)
    assert metadata["functional"] == "B3LYP-D3BJ"
    assert metadata["request_num_cpu"] == 4
    assert metadata["request_memory"].magnitude == pytest.approx(2000.0)
    payload = parser._parse_frame(block, context=FrameParseContext(additional_data={}))

    assert payload["keywords"] == "B3LYP D3BJ def2-TZVP"
    assert payload["method"] == "DFT"
    assert payload["functional"] == "B3LYP-D3BJ"
    assert payload["basis_set"] == "def2-TZVP"
    assert payload["request_num_cpu"] == 4
    assert payload["request_memory"].magnitude == pytest.approx(2000.0)
    assert payload["model_chemistry"].functional == "B3LYP-D3BJ"
    assert payload["task_requests"][0].task_type == "sp"
    assert payload["output_print_settings"][0].target == "P_UNO_OccNum"


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
    assert frame.geometry.point_charges == [{"charge": 0.5, "x": 0.0, "y": 0.0, "z": 1.0}]


def test_orcainp_writer_builds_common_sp_input() -> None:
    atom_coords = [
        (6, -2.171223, 0.546866, -0.754713),
        (6, -2.039335, 1.221513, 0.473138),
        (1, -3.149621, 0.162675, -1.050236),
        (1, -1.515817, 0.807716, -1.589491),
        (1, -1.280986, 2.001154, 0.581962),
        (1, -2.920192, 1.352219, 1.105139),
        (6, -1.062025, -0.223728, 1.751743),
        (16, -2.049826, -1.592815, 1.218779),
        (6, -1.282622, -1.395756, -0.360778),
        (6, 0.045510, -0.925588, -0.164056),
        (6, 0.172494, -0.263103, 1.049059),
        (8, 1.232979, 0.531606, 1.362523),
        (6, 2.317203, 0.371785, 0.453516),
        (6, 1.822191, 0.302852, -0.986418),
        (8, 0.956597, -0.814109, -1.173951),
        (1, -1.146877, 0.125758, 2.781787),
        (1, -1.565218, -2.057787, -1.180760),
        (1, 2.984903, 1.233542, 0.596005),
        (1, 2.871089, -0.553312, 0.697110),
        (1, 1.288049, 1.236942, -1.240989),
        (1, 2.665526, 0.177568, -1.680268),
    ]
    frame = XYZFileFrameMemory(
        atoms=[atom for atom, *_ in atom_coords],
        coords=np.asarray([coords for _, *coords in atom_coords]) * atom_ureg.angstrom,
        charge=0,
        multiplicity=1,
    )

    rendered = frame.format_transform(
        "orcainp",
        keywords=("wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP"),
        nprocs=16,
        maxcore=4000,
        blocks={"scf": {"MaxIter": 300, "STABPerform": True}},
    )

    assert "orcainp" in get_supported_writer_formats()
    assert rendered.startswith(
        "! wB97M-V def2-TZVPP RIJCOSX def2/J TightSCF DefGrid3 NoAutoStart SP\n\n"
        "%pal\n  nprocs 16\nend\n\n"
        "%maxcore 4000\n\n"
        "%scf\n  MaxIter 300\n  STABPerform true\nend\n\n"
        "* xyz 0 1\n"
    )
    coordinate_lines = rendered.splitlines()[14:-1]
    assert len(coordinate_lines) == 21
    assert coordinate_lines[0] == "C     -2.1712230000     0.5468660000    -0.7547130000"
    assert coordinate_lines[7] == "S     -2.0498260000    -1.5928150000     1.2187790000"
    assert coordinate_lines[-1] == "H      2.6655260000     0.1775680000    -1.6802680000"
    assert rendered.endswith("\n*")


def test_orcainp_writer_uses_inp_output_extension(tmp_path: Path) -> None:
    fixture_path = Path(__file__).resolve().parent / "test_files" / "xyz" / "dsgdb9nsd_004015-7.xyz"
    batch = AutoParser(str(fixture_path))
    assert len(batch) > 0
    file_model = batch[0]

    rendered = file_model.format_transform(
        "orcainp",
        file_path=tmp_path / "calculation.placeholder",
        write_to_disk=True,
        keywords="HF def2-SVP SP",
        nprocs=2,
        maxcore=1000,
    )

    output_path = tmp_path / "calculation.inp"
    assert output_path.read_text(encoding="utf-8") == rendered
    assert not (tmp_path / "calculation.orcainp").exists()


def test_orcainp_writer_requires_explicit_keywords_for_cross_format_input() -> None:
    frame = XYZFileFrameMemory(
        atoms=[1],
        coords=np.zeros((1, 3)) * atom_ureg.angstrom,
    )

    with pytest.raises(ValueError, match="requires at least one keyword line"):
        frame.format_transform("orcainp")


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

    rendered = file_model.format_transform("orcainp", frame="all")
    assert rendered.count("$new_job") == 1
    assert rendered.count("! SP") == 2


def test_orcainp_writer_canonicalizes_percent_coords_without_duplication(
    tmp_path: Path,
) -> None:
    path = _write_orcainp(
        tmp_path,
        """! SP HF def2-SVP
%coords
  CTyp xyz
  Charge 0
  Mult 1
  coords
    H 0.0 0.0 0.0
    H 0.0 0.0 0.7
  end
end
""",
    )
    file_model = AutoParser(path, parser_detection="orcainp")[0]

    rendered = file_model.format_transform("orcainp")

    assert "%coords" not in rendered.lower()
    assert rendered.count("* xyz 0 1") == 1
    assert rendered.count("H ") == 2


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
