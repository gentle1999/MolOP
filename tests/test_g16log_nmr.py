from pathlib import Path

import numpy as np
import pytest

from molop.io.logic.gaussian.log.parsers.G16LogFileParser import G16LogFileParserDisk


FIXTURE = Path(__file__).resolve().parent / "test_files" / "g16log" / "test_nmr.log"
COUPLING_FIXTURE = (
    Path(__file__).resolve().parent / "test_files" / "g16log" / "test_nmr_coupling.log"
)


def test_g16log_nmr_shielding_tensors_are_attached_to_property_frame() -> None:
    parsed = G16LogFileParserDisk().parse(str(FIXTURE), release_file_content=True)

    assert len(parsed.frames) == 5
    assert all(frame.nmr is None for frame in parsed.frames[:-1])
    frame = parsed.frames[-1]
    assert [task.task_type for task in frame.task_requests] == ["property"]
    assert frame.nmr is not None
    assert frame.nmr.gauge == "GIAO"
    assert [item.atom_index for item in frame.nmr.shielding_tensors] == [0, 1, 2, 3, 4]
    assert [item.atom_symbol for item in frame.nmr.shielding_tensors] == ["C", "H", "H", "H", "H"]

    carbon = frame.nmr.shielding_tensors[0]
    assert carbon.isotropic is not None
    assert carbon.isotropic.m_as("ppm") == pytest.approx(191.1871)
    assert carbon.anisotropy is not None
    assert carbon.anisotropy.m_as("ppm") == pytest.approx(0.0022)
    assert carbon.principal_values is not None
    np.testing.assert_allclose(
        carbon.principal_values.m_as("ppm"),
        [191.1857, 191.1871, 191.1886],
    )
    np.testing.assert_allclose(
        carbon.shielding_tensor.m_as("ppm"),
        [
            [191.1880, 0.0000, 0.0014],
            [0.0000, 191.1871, -0.0000],
            [0.0010, -0.0000, 191.1864],
        ],
    )

    hydrogen = frame.nmr.shielding_tensors[2]
    assert hydrogen.isotropic is not None
    assert hydrogen.isotropic.m_as("ppm") == pytest.approx(31.2266)
    assert hydrogen.anisotropy is not None
    assert hydrogen.anisotropy.m_as("ppm") == pytest.approx(8.9013)


def test_g16log_nmr_shielding_is_embedded_in_rdkit_atom_properties() -> None:
    frame = G16LogFileParserDisk().parse(str(FIXTURE), release_file_content=True).frames[-1]

    mol = frame.qm_embedded_rdmol()
    assert mol is not None
    assert mol.GetProp("NMR_GAUGE_BY_GAUSSIAN") == "GIAO"

    carbon = mol.GetAtomWithIdx(0)
    isotropic_key = "NMR_SHIELDING_ISOTROPIC_PPM_BY_GAUSSIAN"
    tensor_xz_key = "NMR_SHIELDING_TENSOR_XZ_PPM_BY_GAUSSIAN"
    principal_key = "NMR_SHIELDING_PRINCIPAL_VALUE_1_PPM_BY_GAUSSIAN"
    convention_key = "NMR_SHIELDING_ANISOTROPY_CONVENTION_BY_GAUSSIAN"
    orientation_key = "NMR_SHIELDING_ORIENTATION_BY_GAUSSIAN"
    assert carbon.GetDoubleProp(isotropic_key) == pytest.approx(191.1871)
    assert carbon.GetDoubleProp(tensor_xz_key) == pytest.approx(0.0014)
    assert carbon.GetDoubleProp(principal_key) == pytest.approx(191.1857)
    assert carbon.GetProp(convention_key) == "Gaussian"
    assert carbon.GetProp(orientation_key) == "unknown"
    assert mol.HasProp(f"atom.dprop.{isotropic_key}")
    assert mol.HasProp(f"atom.prop.{convention_key}")

    sdf_block = frame.to_population_embedded_SDF_block()
    assert f"atom.dprop.{isotropic_key}" in sdf_block
    assert f"atom.prop.{convention_key}" in sdf_block

    unembedded = frame.qm_embedded_rdmol(embed_nmr=False)
    assert unembedded is not None
    assert not unembedded.GetAtomWithIdx(0).HasProp(isotropic_key)


def test_g16log_nmr_total_spin_spin_coupling_matrices_are_structured() -> None:
    parsed = G16LogFileParserDisk().parse(str(COUPLING_FIXTURE), release_file_content=True)
    frame = parsed.frames[-1]

    assert frame.nmr is not None
    assert frame.nmr.gauge == "GIAO"
    assert len(frame.nmr.shielding_tensors) == 9
    assert frame.nmr.coupling_atom_indices == list(range(9))
    assert frame.nmr.spin_spin_coupling_k is not None
    assert frame.nmr.spin_spin_coupling_j is not None
    assert frame.nmr.spin_spin_coupling_k.shape == (9, 9)
    assert frame.nmr.spin_spin_coupling_j.shape == (9, 9)
    assert set(frame.nmr.spin_spin_coupling_k_components) == {"FC", "SD", "PSO", "DSO"}
    assert set(frame.nmr.spin_spin_coupling_j_components) == {"FC", "SD", "PSO", "DSO"}
    assert frame.nmr.spin_spin_coupling_k[1, 0].m_as("Hz") == pytest.approx(20.634)
    assert frame.nmr.spin_spin_coupling_j[1, 0].m_as("Hz") == pytest.approx(40.7217)
    assert frame.nmr.spin_spin_coupling_j[8, 2].m_as("Hz") == pytest.approx(-63.4582)
    expected_k_10 = {"FC": 20.0849, "SD": 0.554805, "PSO": -0.0953026, "DSO": 0.0896238}
    expected_j_10 = {"FC": 39.6380, "SD": 1.09492, "PSO": -0.188082, "DSO": 0.176875}
    for component_name in ("FC", "SD", "PSO", "DSO"):
        k_component = frame.nmr.spin_spin_coupling_k_components[component_name]
        j_component = frame.nmr.spin_spin_coupling_j_components[component_name]
        assert k_component.shape == (9, 9)
        assert j_component.shape == (9, 9)
        assert k_component[1, 0].m_as("Hz") == pytest.approx(expected_k_10[component_name])
        assert j_component[1, 0].m_as("Hz") == pytest.approx(expected_j_10[component_name])
        np.testing.assert_allclose(k_component.m_as("Hz"), k_component.m_as("Hz").T)
        np.testing.assert_allclose(j_component.m_as("Hz"), j_component.m_as("Hz").T)
    np.testing.assert_allclose(
        frame.nmr.spin_spin_coupling_j.m_as("Hz"),
        frame.nmr.spin_spin_coupling_j.m_as("Hz").T,
    )
    np.testing.assert_allclose(
        sum(
            component.m_as("Hz") for component in frame.nmr.spin_spin_coupling_k_components.values()
        ),
        frame.nmr.spin_spin_coupling_k.m_as("Hz"),
        atol=1e-4,
    )
    np.testing.assert_allclose(
        sum(
            component.m_as("Hz") for component in frame.nmr.spin_spin_coupling_j_components.values()
        ),
        frame.nmr.spin_spin_coupling_j.m_as("Hz"),
        atol=1e-3,
    )
