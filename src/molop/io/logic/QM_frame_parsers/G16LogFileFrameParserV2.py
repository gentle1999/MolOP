from __future__ import annotations

from collections.abc import Mapping
from enum import Enum, auto
from typing import Any, cast

from molop.io.base_models.DataClasses import (
    ChargeSpinPopulations,
    Energies,
    MolecularOrbitals,
    Polarizability,
    ThermalInformations,
    TotalSpin,
    Vibrations,
)
from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.logic.QM_frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.QM_frame_models.G16V3Components import (
    _parse_running_time,
    _temperature_and_pressure_from_block,
)
from molop.io.logic.QM_frame_parsers._g16_v2_extractors import (
    ParseState,
    extract_berny_from_state,
    extract_electric_dipole_and_polarizability_from_state,
    extract_energies_and_total_spin_from_state,
    extract_forces_from_state,
    extract_hessian_from_state,
    extract_input_coords_from_state,
    extract_polarizability_from_state,
    extract_populations_from_state,
    extract_rotation_consts_from_state,
    extract_standard_coords_from_state,
    extract_tail_energies_from_state,
    extract_tail_hessian_from_state,
    extract_tail_metadata_from_state,
    extract_tail_polarizability_from_state,
    extract_tail_thermal_infos_from_state,
    extract_thermal_infos_from_state,
    extract_vibrations_from_state,
)
from molop.io.patterns.G16Patterns import g16_log_patterns
from molop.utils.functions import merge_models


class G16ParsePhase(Enum):
    """Explicit stages for the sequential Gaussian frame parser.

    Transition table
    ----------------
    HEADER -> ORIENTATION
    ORIENTATION -> STRUCTURE_ONLY_CHECK
    STRUCTURE_ONLY_CHECK -> DONE | ROTATION
    ROTATION -> SCF
    SCF -> ISOTROPIC_POLARIZABILITY
    ISOTROPIC_POLARIZABILITY -> POPULATION
    POPULATION -> FREQUENCY
    FREQUENCY -> THERMOCHEM
    THERMOCHEM -> FORCES
    FORCES -> HESSIAN
    HESSIAN -> OPTIMIZATION
    OPTIMIZATION -> ELECTRIC_RESPONSE
    ELECTRIC_RESPONSE -> ARCHIVE_TAIL
    ARCHIVE_TAIL -> DONE

    Notes
    -----
    - The order is intentionally biased toward the original fast V1 scan path.
    - Each phase consumes from the shared ``ParseState`` and advances the cursor.
    - V3 component trees are no longer part of the extraction path; they are rebuilt
      lazily from frame data for fakeG rendering and inspection.
    """

    HEADER = auto()
    ORIENTATION = auto()
    STRUCTURE_ONLY_CHECK = auto()
    ROTATION = auto()
    SCF = auto()
    ISOTROPIC_POLARIZABILITY = auto()
    POPULATION = auto()
    FREQUENCY = auto()
    THERMOCHEM = auto()
    FORCES = auto()
    HESSIAN = auto()
    OPTIMIZATION = auto()
    ELECTRIC_RESPONSE = auto()
    ARCHIVE_TAIL = auto()
    DONE = auto()


class G16LogFileFrameParserV2Mixin:
    """Explicit-state Gaussian frame parser.

    Design goals
    ------------
    - Preserve the sequential extraction performance characteristics of V1.
    - Replace implicit ``self._block`` mutation chains with an explicit ``ParseState``.
    - Keep phase ordering first-class so later maintenance can reason about transitions.
    - Leave V3 component trees to render/inspection paths instead of hot extraction.
    """

    def _run_header_phase(self, block: str, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse frame-global header data that does not depend on cursor state."""
        if charge_multiplicity := g16_log_patterns.CHARGE_MULTIPLICITY.match_content(block):
            infos["charge"] = int(charge_multiplicity[0][0])
            infos["multiplicity"] = int(charge_multiplicity[0][1])
        if running_time := _parse_running_time(block):
            infos["running_time"] = running_time
        return G16ParsePhase.ORIENTATION

    def _run_orientation_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Extract input and standard orientations, establishing core atom/coord payload."""
        atoms, coords = extract_input_coords_from_state(state)
        if atoms and coords is not None:
            infos["atoms"] = atoms
            infos["coords"] = coords

        atoms, standard_coords = extract_standard_coords_from_state(state)
        if atoms and standard_coords is not None:
            infos["atoms"] = atoms
            infos["standard_coords"] = standard_coords

        return G16ParsePhase.STRUCTURE_ONLY_CHECK

    def _run_structure_only_check(self) -> G16ParsePhase:
        """Stop early for structure-only mode before any expensive electronic parsing."""
        return (
            G16ParsePhase.DONE
            if cast(_HasParseMethod, self).only_extract_structure
            else G16ParsePhase.ROTATION
        )

    def _run_rotation_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse lightweight rotational metadata that may appear before SCF/population sections."""
        if (rotation_consts := extract_rotation_consts_from_state(state)) is not None:
            infos["rotation_constants"] = rotation_consts
        return G16ParsePhase.SCF

    def _run_scf_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse SCF/energy and spin-state information from the current cursor position."""
        energies_dict, total_spin_dict = extract_energies_and_total_spin_from_state(state)
        if energies_dict:
            infos["energies"] = Energies.model_validate(energies_dict)
        if total_spin_dict:
            infos["total_spin"] = TotalSpin.model_validate(total_spin_dict)
        return G16ParsePhase.ISOTROPIC_POLARIZABILITY

    def _run_isotropic_polarizability_phase(
        self, state: ParseState, infos: dict[str, Any]
    ) -> G16ParsePhase:
        """Parse early scalar polarizability data before the larger population section."""
        if (polarizability := extract_polarizability_from_state(state)) is not None:
            infos["polarizability"] = Polarizability.model_validate(polarizability)
        return G16ParsePhase.POPULATION

    def _run_population_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse population analysis, orbitals, and any coupled response data in that block."""
        if populations := extract_populations_from_state(state):
            if (molecular_orbitals := populations.get("molecular_orbitals")) is not None:
                populations["molecular_orbitals"] = MolecularOrbitals.model_validate(
                    molecular_orbitals
                )
            if (charge_spin_populations := populations.get("charge_spin_populations")) is not None:
                populations["charge_spin_populations"] = ChargeSpinPopulations.model_validate(
                    charge_spin_populations
                )
            if (polarizability := populations.get("polarizability")) is not None:
                populations["polarizability"] = Polarizability.model_validate(polarizability)
            infos.update(populations)
        return G16ParsePhase.FREQUENCY

    def _run_frequency_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse vibrational frequency blocks if present."""
        if vibrations := extract_vibrations_from_state(state):
            infos["vibrations"] = Vibrations.model_validate(vibrations)
        return G16ParsePhase.THERMOCHEM

    def _run_thermochem_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse thermochemistry summaries that often follow frequency sections."""
        cursor_before = state.cursor
        if (thermal_info := extract_thermal_infos_from_state(state)) is not None:
            infos["thermal_informations"] = ThermalInformations.model_validate(thermal_info)
        if temp_pressure := _temperature_and_pressure_from_block(state.content[cursor_before:]):
            infos.update(temp_pressure)
        return G16ParsePhase.FORCES

    def _run_forces_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse Cartesian forces when present."""
        if (forces := extract_forces_from_state(state)) is not None:
            infos["forces"] = forces
        return G16ParsePhase.HESSIAN

    def _run_hessian_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse Hessian / second-derivative blocks."""
        if (hessian := extract_hessian_from_state(state)) is not None:
            infos["hessian"] = hessian
        return G16ParsePhase.OPTIMIZATION

    def _run_optimization_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Parse Berny optimization state summaries."""
        if (berny := extract_berny_from_state(state)) is not None:
            infos["geometry_optimization_status"] = berny
        return G16ParsePhase.ELECTRIC_RESPONSE

    def _run_electric_response_phase(
        self, state: ParseState, infos: dict[str, Any]
    ) -> G16ParsePhase:
        """Merge late electric-dipole/polarizability sections into existing response data."""
        if (
            polarizability := extract_electric_dipole_and_polarizability_from_state(state)
        ) is not None:
            polarizability_model = Polarizability.model_validate(polarizability)
            if "polarizability" in infos:
                infos["polarizability"] = merge_models(
                    cast(Polarizability, infos["polarizability"]),
                    polarizability_model,
                    force_update=True,
                )
            else:
                infos["polarizability"] = polarizability_model
        return G16ParsePhase.ARCHIVE_TAIL

    def _run_archive_tail_phase(self, state: ParseState, infos: dict[str, Any]) -> G16ParsePhase:
        """Use archive-tail data as the final fallback/augmentation stage."""
        if tail := extract_tail_metadata_from_state(state):
            for key, value in tail.items():
                if key not in infos:
                    infos[key] = value

            if tail_energies := extract_tail_energies_from_state(state):
                tail_energies_model = Energies.model_validate(tail_energies)
                if "energies" in infos:
                    infos["energies"] = merge_models(
                        cast(Energies, infos["energies"]), tail_energies_model
                    )
                else:
                    infos["energies"] = tail_energies_model

            if tail_thermal_info := extract_tail_thermal_infos_from_state(state):
                tail_thermal_model = ThermalInformations.model_validate(tail_thermal_info)
                if "thermal_informations" in infos:
                    infos["thermal_informations"] = merge_models(
                        cast(ThermalInformations, infos["thermal_informations"]),
                        tail_thermal_model,
                    )
                else:
                    infos["thermal_informations"] = tail_thermal_model

            if tail_polarizability := extract_tail_polarizability_from_state(state):
                tail_polarizability_model = Polarizability.model_validate(tail_polarizability)
                if "polarizability" in infos:
                    infos["polarizability"] = merge_models(
                        cast(Polarizability, infos["polarizability"]),
                        tail_polarizability_model,
                    )
                else:
                    infos["polarizability"] = tail_polarizability_model

            if (
                tail_hessian := extract_tail_hessian_from_state(state)
            ) is not None and "hessian" not in infos:
                infos["hessian"] = tail_hessian

        return G16ParsePhase.DONE

    def _parse_frame(self) -> Mapping[str, Any]:
        """Execute the explicit phase machine until all extractors have run."""
        block = cast(_HasParseMethod, self)._block
        state = ParseState(block)
        infos: dict[str, Any] = {"qm_software": "Gaussian"}

        phase = G16ParsePhase.HEADER
        while phase is not G16ParsePhase.DONE:
            if phase is G16ParsePhase.HEADER:
                phase = self._run_header_phase(block, infos)
            elif phase is G16ParsePhase.ORIENTATION:
                phase = self._run_orientation_phase(state, infos)
            elif phase is G16ParsePhase.STRUCTURE_ONLY_CHECK:
                phase = self._run_structure_only_check()
            elif phase is G16ParsePhase.ROTATION:
                phase = self._run_rotation_phase(state, infos)
            elif phase is G16ParsePhase.SCF:
                phase = self._run_scf_phase(state, infos)
            elif phase is G16ParsePhase.ISOTROPIC_POLARIZABILITY:
                phase = self._run_isotropic_polarizability_phase(state, infos)
            elif phase is G16ParsePhase.POPULATION:
                phase = self._run_population_phase(state, infos)
            elif phase is G16ParsePhase.FREQUENCY:
                phase = self._run_frequency_phase(state, infos)
            elif phase is G16ParsePhase.THERMOCHEM:
                phase = self._run_thermochem_phase(state, infos)
            elif phase is G16ParsePhase.FORCES:
                phase = self._run_forces_phase(state, infos)
            elif phase is G16ParsePhase.HESSIAN:
                phase = self._run_hessian_phase(state, infos)
            elif phase is G16ParsePhase.OPTIMIZATION:
                phase = self._run_optimization_phase(state, infos)
            elif phase is G16ParsePhase.ELECTRIC_RESPONSE:
                phase = self._run_electric_response_phase(state, infos)
            elif phase is G16ParsePhase.ARCHIVE_TAIL:
                phase = self._run_archive_tail_phase(state, infos)

        return infos


class G16LogFileFrameParserV2Memory(
    G16LogFileFrameParserV2Mixin, BaseFrameParser[G16LogFileFrameMemory]
):
    _file_frame_class_ = G16LogFileFrameMemory


class G16LogFileFrameParserV2Disk(
    G16LogFileFrameParserV2Mixin, BaseFrameParser[G16LogFileFrameDisk]
):
    _file_frame_class_ = G16LogFileFrameDisk
