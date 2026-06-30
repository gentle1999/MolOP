import numpy as np
import pytest
from pydantic import ValidationError

from molop.io.base_models.Bases import (
    PropertyBundle,
    PropertyPoint,
    PropertySeries,
    PropertyTable,
    SpectralBand,
    Spectrum,
    TensorProperty,
)
from molop.io.base_models.ChemFileFrame import BaseQMInputFrame
from molop.io.base_models.DataClasses import (
    ActiveSpace,
    BondOrders,
    ChargeSpinPopulations,
    Dispersions,
    ElectronicConfiguration,
    ElectronicState,
    ElectronicStates,
    Energies,
    ExcitedStateRequest,
    MolecularOrbitals,
    MultireferenceRequest,
    MultireferenceResult,
    MultireferenceStateBlock,
    Polarizability,
    QMBasisSet,
    QMModelChemistry,
    QMResourceRequest,
    QMTaskRequest,
    SinglePointProperties,
    ThermalInformations,
    TotalSpin,
    Vibration,
    Vibrations,
)
from molop.io.base_models.Molecule import Molecule
from molop.io.logic.qminput_frame_models.GJFFileFrame import GJFFileFrameMemory
from molop.io.logic.qminput_frame_models.ORCAInpFileFrame import ORCAInpFileFrameMemory
from molop.unit import atom_ureg


def test_molecule_total_electrons_accounts_for_charge() -> None:
    assert Molecule(atoms=[1], charge=0).total_electrons == 1
    assert Molecule(atoms=[1], charge=1).total_electrons == 0
    assert Molecule(atoms=[8, 1, 1], charge=-1).total_electrons == 11


def test_molecular_orbitals_support_fractional_occupancies() -> None:
    orbitals = MolecularOrbitals(
        alpha_energies=np.array([-0.5, -0.2, 0.1]) * atom_ureg.hartree,
        beta_energies=np.array([-0.45, -0.1, 0.2]) * atom_ureg.hartree,
        alpha_occupancies=[1.0, 0.5, 0.0],
        beta_occupancies=[1.0, 0.0, 0.0],
    )

    assert orbitals.HOMO_id == 1
    assert orbitals.beta_HOMO_id == 0
    assert orbitals.SOMO_ids == [1]
    assert orbitals.HOMO is not None
    assert orbitals.HOMO.alpha_occupancy == 0.5


def test_resource_request_normalizes_memory_unit() -> None:
    request = QMResourceRequest(num_cpu=8, memory=1 * atom_ureg.gigabyte)

    assert request.num_cpu == 8
    assert request.memory is not None
    assert request.memory.units == atom_ureg.megabyte
    assert request.memory.magnitude == pytest.approx(1000.0)


def test_model_chemistry_and_task_request_capture_common_input_semantics() -> None:
    model = QMModelChemistry(
        method_family="DFT",
        method="B3LYP",
        functional="B3LYP-GD3BJ",
        basis_set="def2-SVP",
        auxiliary_basis_set="def2/J",
        dispersion_correction="GD3BJ",
        basis_sets=[
            QMBasisSet(name="def2-SVP", role="orbital", scope="global"),
            QMBasisSet(name="def2/J", role="auxiliary", scope="global"),
            QMBasisSet(name="def2-TZVP", role="orbital", scope="atom", atom_indices=[0]),
        ],
        solvation_model="smd",
        solvent="water",
    )
    task = QMTaskRequest(
        task_type="opt",
        derivative_order=1,
        transition_state=True,
        scan=True,
        source_keywords=["OptTS"],
        source_blocks=["geom"],
    )

    assert model.functional == "B3LYP-GD3BJ"
    assert model.basis_sets[2].atom_indices == [0]
    assert task.task_type == "opt"
    assert task.transition_state is True
    assert task.scan is True


def test_excited_state_request_and_result_states_are_sequence_like() -> None:
    request = ExcitedStateRequest(
        enabled=True,
        family="TDDFT",
        nroots=5,
        root=2,
        triplets=True,
        properties=["nacme"],
    )
    states = ElectronicStates(
        states=[
            ElectronicState(state_index=0, label="S0", energy=-100 * atom_ureg.hartree),
            ElectronicState(
                state_index=1,
                root=1,
                label="T1",
                multiplicity=3,
                excitation_energy=1.5 * atom_ureg.eV,
                oscillator_strength=0.0,
                configurations=[
                    ElectronicConfiguration(
                        label="HOMO->LUMO",
                        coefficient=0.7,
                        weight=0.49,
                        orbital_indices=[4, 5],
                    )
                ],
            ),
        ]
    )

    assert request.enabled is True
    assert request.family == "TDDFT"
    assert len(states) == 2
    assert states.ground_state is not None
    assert states.ground_state.label == "S0"
    assert [state.label for state in states.excited_states] == ["T1"]
    assert states[1].configurations[0].orbital_indices == [4, 5]


def test_electronic_states_project_to_state_table_and_transition_spectrum() -> None:
    states = ElectronicStates(
        states=[
            ElectronicState(state_index=0, label="S0", energy=-100 * atom_ureg.hartree),
            ElectronicState(
                state_index=1,
                root=1,
                label="S1",
                excitation_energy=2.5 * atom_ureg.eV,
                oscillator_strength=0.12,
                transition_dipole=np.array([0.1, 0.2, 0.3]) * atom_ureg.debye,
            ),
        ]
    )

    table = states.to_state_table()
    assert isinstance(table, PropertyTable)
    assert table["excitation_energy"].units == atom_ureg.eV
    assert table.row(1)["label"] == "S1"

    transitions = states.to_transitions()
    assert len(transitions) == 1
    assert transitions[0].energy.units == atom_ureg.eV
    assert transitions[0].oscillator_strength == pytest.approx(0.12)

    spectrum = states.to_spectrum()
    assert isinstance(spectrum, Spectrum)
    assert len(spectrum.transitions) == 1

    bundle = states.to_property_bundle()
    assert "electronic_states" in bundle.tables
    assert "electronic_transitions" in bundle.spectra
    assert "electronic" in bundle.transitions


def test_multireference_request_and_result_capture_active_space_and_states() -> None:
    active_space = ActiveSpace(electrons=8, orbitals=6, raw="cas(8,6)")
    block = MultireferenceStateBlock(
        multiplicity=3,
        irrep="1",
        nroots=2,
        refs="cas(8,6)",
        active_space=active_space,
    )
    request = MultireferenceRequest(
        enabled=True,
        method="MRCI",
        reference_method="CASSCF",
        ci_type="MRCI",
        active_space=active_space,
        state_blocks=[block],
        thresholds={"tsel": 1e-7},
        corrections=["Davidson:DAV"],
    )
    result = MultireferenceResult(
        method="MRCI",
        reference_method="CASSCF",
        ci_type="MRCI",
        active_space=active_space,
        electronic_states=ElectronicStates(states=[ElectronicState(root=0, label="Triplet")]),
        corrections={"Davidson": -0.01 * atom_ureg.hartree},
    )

    assert request.active_space is not None
    assert request.active_space.electrons == 8
    assert request.state_blocks[0].multiplicity == 3
    assert request.thresholds["tsel"] == pytest.approx(1e-7)
    assert result.electronic_states is not None
    assert result.electronic_states[0].label == "Triplet"
    bundle = result.to_property_bundle()
    assert bundle.scalar_properties["Davidson_correction"].units == atom_ureg.hartree
    assert "electronic_states" in bundle.tables


def test_structured_qm_input_projection_is_authoritative() -> None:
    frame = BaseQMInputFrame(
        keywords="legacy keywords",
        method="legacy-method",
        functional="legacy-functional",
        basis_set="legacy-basis",
        resources_raw="%old",
        request_num_cpu=2,
        request_memory=256 * atom_ureg.megabyte,
        model_chemistry=QMModelChemistry(
            method_family="DFT",
            method="B3LYP",
            functional="B3LYP-D3BJ",
            basis_set="def2-SVP",
            raw_keywords="! B3LYP D3BJ def2-SVP",
        ),
        resource_request=QMResourceRequest(
            num_cpu=8,
            memory=1 * atom_ureg.gigabyte,
            raw="%pal nprocs 8 end",
        ),
    )

    frame.project_common_qm_fields()

    assert frame.keywords == "! B3LYP D3BJ def2-SVP"
    assert frame.method == "DFT"
    assert frame.functional == "B3LYP-D3BJ"
    assert frame.basis_set == "def2-SVP"
    assert frame.request_num_cpu == 8
    assert frame.request_memory is not None
    assert frame.request_memory.magnitude == pytest.approx(1000.0)
    assert frame.resources_raw == "%pal nprocs 8 end"


def test_legacy_qm_input_backfill_does_not_overwrite_structured_values() -> None:
    frame = BaseQMInputFrame(
        keywords="legacy keywords",
        method="legacy-method",
        functional="legacy-functional",
        basis_set="legacy-basis",
        resources_raw="%old",
        request_num_cpu=2,
        model_chemistry=QMModelChemistry(
            method_family="DFT",
            method="B3LYP",
            functional="B3LYP-D3BJ",
            basis_set="def2-SVP",
            raw_keywords="! B3LYP D3BJ def2-SVP",
        ),
        resource_request=QMResourceRequest(num_cpu=8, raw="%pal nprocs 8 end"),
    )

    frame.refresh_common_qm_containers()

    assert frame.model_chemistry.method == "B3LYP"
    assert frame.model_chemistry.functional == "B3LYP-D3BJ"
    assert frame.model_chemistry.basis_set == "def2-SVP"


def test_property_bundle_and_spectrum_store_generic_property_shapes() -> None:
    table = PropertyTable(
        columns={
            "orbital_index": np.array([0, 1]),
            "energy": np.array([-0.5, 0.1]) * atom_ureg.hartree,
        },
        row_labels=["0", "1"],
    )
    series = PropertySeries(
        axis_label="wavenumber",
        value_label="intensity",
        points=[
            PropertyPoint(label="p1", value=1.0 * atom_ureg.cm**-1),
            PropertyPoint(label="p2", value=2.0 * atom_ureg.cm**-1),
        ],
    )
    spectrum = Spectrum(
        label="IR",
        x_label="wavenumber",
        y_label="intensity",
        bands=[SpectralBand(label="band1", center=1000.0 * atom_ureg.cm**-1, intensity=0.8)],
    )
    bundle = PropertyBundle(
        scalar_properties={"energy": -100.0},
        tables={"orbitals": table},
        series={"ir": series},
        spectra={"ir": spectrum},
        metadata={"source": "orca"},
    )

    assert len(table) == 2
    assert table["energy"].units == atom_ureg.hartree
    assert table.row(0)["orbital_index"] == 0
    assert len(series) == 2
    assert series[0].label == "p1"
    assert "orbitals" in bundle.tables
    assert bundle.spectra["ir"].bands[0].label == "band1"
    dump = bundle.model_dump()
    assert dump["scalar_properties"]["energy"] == -100.0
    assert dump["metadata"]["source"] == "orca"


def test_molecular_orbitals_project_to_columnar_tables_without_eager_records() -> None:
    orbitals = MolecularOrbitals(
        electronic_state="S0",
        alpha_energies=np.array([-0.5, -0.2, 0.1]) * atom_ureg.hartree,
        beta_energies=np.array([-0.45, -0.1, 0.2]) * atom_ureg.hartree,
        alpha_occupancies=[1.0, 0.5, 0.0],
        beta_occupancies=[1.0, 0.0, 0.0],
        alpha_symmetries=["a", "b", "c"],
    )

    alpha_table = orbitals.to_orbital_table("alpha")
    assert isinstance(alpha_table, PropertyTable)
    assert alpha_table.metadata["source"] == "MolecularOrbitals"
    assert alpha_table.metadata["spin"] == "alpha"
    assert alpha_table["energy"].units == atom_ureg.hartree
    assert alpha_table.row(1)["occupancy"] == 0.5
    assert alpha_table.row(1)["symmetry"] == "b"

    bundle = orbitals.to_property_bundle()
    assert "alpha_orbitals" in bundle.tables
    assert "beta_orbitals" in bundle.tables
    assert bundle.scalar_properties["HOMO_id"] == 1
    assert bundle.scalar_properties["HOMO_energy"].units == atom_ureg.hartree


def test_charge_spin_populations_project_to_atomic_population_table() -> None:
    populations = ChargeSpinPopulations(
        mulliken_charges=[-0.2, 0.1, 0.1],
        mulliken_spins=[0.0, 0.5, -0.5],
        hirshfeld_charges=[-0.1, 0.05, 0.05],
    )

    table = populations.to_population_table(atom_symbols=["O", "H", "H"])
    assert isinstance(table, PropertyTable)
    assert table.column_names == [
        "atom_index",
        "atom_symbol",
        "mulliken_charges",
        "mulliken_spins",
        "hirshfeld_charges",
    ]
    assert table.row(0)["atom_symbol"] == "O"
    assert table.row(1)["mulliken_spins"] == 0.5

    bundle = populations.to_property_bundle(atom_symbols=["O", "H", "H"])
    assert "atomic_populations" in bundle.tables


def test_polarizability_projects_to_generic_bundle_without_losing_units() -> None:
    polarizability = Polarizability(
        isotropic_polarizability=10.0 * atom_ureg.bohr**3,
        anisotropic_polarizability=2.0 * atom_ureg.bohr**3,
        electric_dipole_moment=np.array([0.0, 1.0, 0.0]) * atom_ureg.debye,
        polarizability_tensor=np.eye(3) * atom_ureg.bohr**3,
    )

    bundle = polarizability.to_property_bundle()
    assert bundle.scalar_properties["isotropic_polarizability"].units == atom_ureg.bohr**3
    assert bundle.vector_properties["electric_dipole_moment"].units == atom_ureg.debye
    assert isinstance(bundle.tensor_properties["polarizability_tensor"], TensorProperty)
    assert bundle.tensor_properties["polarizability_tensor"].tensor.units == atom_ureg.bohr**3


def test_scalar_and_matrix_properties_project_to_generic_bundles() -> None:
    energies = Energies(
        electronic_energy=-100.0 * atom_ureg.hartree,
        reference_energy=-99.0 * atom_ureg.hartree,
    )
    energy_bundle = energies.to_property_bundle()
    assert energy_bundle.scalar_properties["total_energy"].units == atom_ureg.hartree
    assert energy_bundle.scalar_properties["electronic_energy"].magnitude == pytest.approx(-100.0)
    assert energy_bundle.scalar_properties["reference_energy"].magnitude == pytest.approx(-99.0)

    assert "scf_energy" not in energies.model_dump()
    with pytest.raises(ValidationError, match="scf_energy"):
        Energies.model_validate({"scf_energy": -98.0 * atom_ureg.hartree})

    thermal = ThermalInformations(
        ZPVE=1.0 * atom_ureg.kcal / atom_ureg.mol,
        S=2.0 * atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin,
    )
    thermal_bundle = thermal.to_property_bundle()
    assert thermal_bundle.scalar_properties["ZPVE"].units == atom_ureg.kcal / atom_ureg.mol
    assert thermal_bundle.scalar_properties["S"].units == (
        atom_ureg.calorie / atom_ureg.mol / atom_ureg.kelvin
    )

    spin = TotalSpin(spin_square=0.75, spin_quantum_number=0.5)
    spin_bundle = spin.to_property_bundle()
    assert spin_bundle.scalar_properties["spin_square"] == 0.75

    bond_orders = BondOrders(wiberg_bond_order=np.eye(2))
    bond_bundle = bond_orders.to_property_bundle()
    assert isinstance(bond_bundle.tensor_properties["wiberg_bond_order"], TensorProperty)
    assert np.allclose(bond_bundle.tensor_properties["wiberg_bond_order"].tensor, np.eye(2))

    dispersions = Dispersions(C6AA=1.0 * atom_ureg.bohr**6)
    dispersion_bundle = dispersions.to_property_bundle()
    assert dispersion_bundle.scalar_properties["C6AA"].units == atom_ureg.bohr**6


def test_single_point_properties_project_scalar_and_atomic_columns() -> None:
    properties = SinglePointProperties(
        vip=8.0 * atom_ureg.Unit("eV / particle"),
        fukui_positive=[0.1, 0.2],
        fukui_negative=[0.3, 0.4],
        fod=[0.5, 0.6],
    )

    table = properties.to_atomic_property_table(atom_symbols=["C", "O"])
    assert table.column_names == [
        "atom_index",
        "atom_symbol",
        "fukui_positive",
        "fukui_negative",
        "fod",
    ]
    assert table.row(0)["atom_symbol"] == "C"
    assert table.row(1)["fod"] == 0.6

    bundle = properties.to_property_bundle(atom_symbols=["C", "O"])
    assert bundle.scalar_properties["vip"].units == atom_ureg.Unit("eV / particle")
    assert "atomic_properties" in bundle.tables


def test_vibrations_keep_matrix_storage_with_lazy_views_and_spectrum_projection() -> None:
    vibrations = Vibrations(
        frequencies=np.array([-10.0, 100.0]) * atom_ureg.cm**-1,
        reduced_masses=np.array([1.0, 2.0]) * atom_ureg.amu,
        force_constants=np.array([0.5, 1.0]) * atom_ureg.Unit("mdyne/angstrom"),
        IR_intensities=np.array([3.0, 4.0]) * atom_ureg.Unit("km/mol"),
        vibration_modes=[
            np.array([[1.0, 0.0, 0.0]]) * atom_ureg.angstrom,
            np.array([[0.0, 1.0, 0.0]]) * atom_ureg.angstrom,
        ],
    )

    first_mode = vibrations[0]
    assert isinstance(first_mode, Vibration)
    assert first_mode.is_imaginary is True
    assert vibrations.frequencies.shape == (2,)
    assert vibrations.num_imaginary == 1
    assert isinstance(vibrations.imaginary_vibrations, Vibrations)
    assert len(vibrations.imaginary_vibrations) == 1

    spectrum = vibrations.to_ir_spectrum()
    assert isinstance(spectrum, Spectrum)
    assert len(spectrum.bands) == 2
    assert spectrum.bands[0].center is not None
    assert spectrum.bands[0].center.units == atom_ureg.cm**-1
    assert spectrum.bands[0].intensity is not None
    assert spectrum.bands[0].intensity.units == atom_ureg.Unit("km/mol")
    assert spectrum.bands[0].metadata["mode_index"] == 0


def test_legacy_qm_input_backfill_populates_empty_structured_values() -> None:
    frame = BaseQMInputFrame(
        keywords="! HF STO-3G",
        method="HF",
        basis_set="STO-3G",
        resources_raw="%pal nprocs 2 end",
        request_num_cpu=2,
    )

    frame.refresh_common_qm_containers()

    assert frame.model_chemistry.method == "HF"
    assert frame.model_chemistry.method_family == "HF"
    assert frame.model_chemistry.basis_set == "STO-3G"
    assert frame.model_chemistry.raw_keywords == "! HF STO-3G"
    assert frame.resource_request.num_cpu == 2
    assert frame.resource_request.raw == "%pal nprocs 2 end"


def test_gaussian_common_projection_ignores_stale_legacy_functional() -> None:
    frame = GJFFileFrameMemory.model_validate(
        {
            "route_section": {"route": "#p b3lyp/def2svp em=gd3bj opt"},
            "title_card": {"title": "projection"},
            "molecule_specifications": "0 1\nH 0.0 0.0 0.0",
            "functional": "stale-functional",
            "basis_set": "stale-basis",
            "method": "stale-method",
        }
    )

    assert frame.model_chemistry.functional == "b3lyp-GD3BJ"
    assert frame.model_chemistry.basis_set == "def2svp"
    assert frame.functional == "b3lyp-GD3BJ"
    assert frame.basis_set == "def2svp"
    assert frame.method == "DFT"


def test_gaussian_hf_projection_rejects_stale_legacy_functional() -> None:
    frame = GJFFileFrameMemory.model_validate(
        {
            "route_section": {"route": "#p hf/3-21g sp"},
            "title_card": {"title": "hf projection"},
            "molecule_specifications": "0 1\nH 0.0 0.0 0.0",
            "functional": "stale-functional",
            "basis_set": "stale-basis",
            "method": "stale-method",
        }
    )

    assert frame.model_chemistry.method == "hf"
    assert frame.model_chemistry.method_family == "HF"
    assert frame.model_chemistry.functional is None
    assert frame.method == "HF"
    assert frame.functional == ""
    assert frame.basis_set == "3-21g"


def test_orca_common_projection_keeps_structured_semantics_authoritative() -> None:
    frame = ORCAInpFileFrameMemory.model_validate(
        {
            "keyword_lines": [{"text": "B3LYP D3BJ def2-SVP"}],
            "blocks": [{"name": "pal", "raw_header": "%pal", "lines": [{"text": "nprocs 4"}]}],
            "geometry": {
                "ctype": "xyz",
                "charge": 0,
                "multiplicity": 1,
                "items": [
                    {
                        "symbol": "H",
                        "atomic_number": 1,
                        "x": 0.0,
                        "y": 0.0,
                        "z": 0.0,
                    }
                ],
            },
            "method": "stale-method",
            "functional": "stale-functional",
            "basis_set": "stale-basis",
            "request_num_cpu": 99,
        }
    )

    assert frame.model_chemistry.method_family == "DFT"
    assert frame.model_chemistry.functional == "B3LYP-D3BJ"
    assert frame.model_chemistry.basis_set == "def2-SVP"
    assert frame.model_chemistry.dispersion_correction == "D3BJ"
    assert frame.method == "DFT"
    assert frame.functional == "B3LYP-D3BJ"
    assert frame.basis_set == "def2-SVP"
    assert frame.request_num_cpu == 4
    assert frame.resource_request.num_cpu == 4
