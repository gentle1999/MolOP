from __future__ import annotations

from collections.abc import Mapping
from enum import Enum, auto
from typing import Any, cast

from molop.io.base_models.FrameParser import BaseFrameParser, _HasParseMethod
from molop.io.base_models.ParseContainers import ModelParseResult, TextParseContext
from molop.io.logic.gaussian.log.frame_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.logic.gaussian.log.frame_parsers._g16_extractors import (
    ParseState,
    extract_archive_tail_payload_from_state,
    extract_berny_from_state,
    extract_electric_dipole_and_polarizability_from_state,
    extract_energies_and_total_spin_from_state,
    extract_forces_from_state,
    extract_hessian_from_state,
    extract_input_coords_from_state,
    extract_nmr_from_state,
    extract_polarizability_from_state,
    extract_populations_from_state,
    extract_rotation_consts_from_state,
    extract_standard_coords_from_state,
    extract_thermal_infos_from_state,
    extract_vibrations_from_state,
    merge_g16_energy_payloads,
)
from molop.io.logic.gaussian.log.frame_parsers._g16_shared import (
    _parse_running_time,
    _temperature_and_pressure_from_block,
)
from molop.io.logic.gaussian.log.parsers._g16_log_file_extractors import (
    extract_g16_termination_status,
)
from molop.io.logic.gaussian.log.parsers._g16_log_patterns import g16_log_patterns


class G16ParsePhase(Enum):
    """Explicit stages for the sequential Gaussian frame parser.

    Transition table
    ----------------
    HEADER -> ORIENTATION
    ORIENTATION -> STRUCTURE_ONLY_CHECK
    STRUCTURE_ONLY_CHECK -> DONE | ROTATION
    ROTATION -> SCF
    SCF -> ISOTROPIC_POLARIZABILITY
    ISOTROPIC_POLARIZABILITY -> NMR
    NMR -> POPULATION
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
    - The order is intentionally biased toward the historical fast scan path.
    - Each phase consumes from the shared ``ParseState`` and advances the cursor.
    - Component trees are rebuilt lazily from frame data for fakeG rendering and
      inspection instead of participating in extraction.
    """

    HEADER = auto()
    ORIENTATION = auto()
    STRUCTURE_ONLY_CHECK = auto()
    ROTATION = auto()
    SCF = auto()
    ISOTROPIC_POLARIZABILITY = auto()
    NMR = auto()
    POPULATION = auto()
    FREQUENCY = auto()
    THERMOCHEM = auto()
    FORCES = auto()
    HESSIAN = auto()
    OPTIMIZATION = auto()
    ELECTRIC_RESPONSE = auto()
    ARCHIVE_TAIL = auto()
    DONE = auto()


class G16LogFileFrameParserMixin:
    """Explicit-state Gaussian frame parser.

    Design goals
    ------------
    - Preserve the sequential extraction performance characteristics of the legacy parser.
    - Replace implicit ``self._block`` mutation chains with an explicit ``ParseState``.
    - Keep phase ordering first-class so later maintenance can reason about transitions.
    - Leave component trees to render/inspection paths instead of hot extraction.
    """

    def _run_header_phase(self, block: str, result: ModelParseResult) -> G16ParsePhase:
        """Parse frame-global header data that does not depend on cursor state."""
        if matched := g16_log_patterns.CHARGE_MULTIPLICITY.search(block):
            result.set("charge", int(matched.group("charge")))
            result.set("multiplicity", int(matched.group("multiplicity")))
        if running_time := _parse_running_time(block):
            result.set("running_time", running_time)
        if not cast(_HasParseMethod, self).only_extract_structure:
            status = extract_g16_termination_status(
                TextParseContext(block),
                include_termination=False,
            )
            if status is not None:
                result.set("status", status)
        return G16ParsePhase.ORIENTATION

    def _run_orientation_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Extract input and standard orientations, establishing core atom/coord payload."""
        atoms, coords = extract_input_coords_from_state(state)
        if atoms and coords is not None:
            result.set("atoms", atoms)
            result.set("coords", coords)

        atoms, standard_coords = extract_standard_coords_from_state(state)
        if atoms and standard_coords is not None:
            result.set("atoms", atoms)
            result.set("standard_coords", standard_coords)

        return G16ParsePhase.STRUCTURE_ONLY_CHECK

    def _run_structure_only_check(self) -> G16ParsePhase:
        """Stop early for structure-only mode before any expensive electronic parsing."""
        return (
            G16ParsePhase.DONE
            if cast(_HasParseMethod, self).only_extract_structure
            else G16ParsePhase.ROTATION
        )

    def _run_rotation_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse lightweight rotational metadata that may appear before SCF/population sections."""
        if (rotation_consts := extract_rotation_consts_from_state(state)) is not None:
            result.set("rotation_constants", rotation_consts)
        return G16ParsePhase.SCF

    def _run_scf_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse SCF/energy and spin-state information from the current cursor position."""
        energies_dict, total_spin_dict = extract_energies_and_total_spin_from_state(
            state,
            capture_source_evidence=cast(_HasParseMethod, self).capture_source_evidence,
        )
        if energies_dict:
            result.set("energies", energies_dict)
        if total_spin_dict:
            result.set("total_spin", total_spin_dict)
        return G16ParsePhase.ISOTROPIC_POLARIZABILITY

    def _run_isotropic_polarizability_phase(
        self, state: ParseState, result: ModelParseResult
    ) -> G16ParsePhase:
        """Parse early scalar polarizability data before the larger population section."""
        if (polarizability := extract_polarizability_from_state(state)) is not None:
            result.set("polarizability", polarizability)
        return G16ParsePhase.NMR

    @staticmethod
    def _run_nmr_phase(state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse per-nucleus magnetic shielding tensors before population analysis."""
        if (nmr := extract_nmr_from_state(state)) is not None:
            result.set("nmr", nmr)
        return G16ParsePhase.POPULATION

    def _run_population_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse population analysis, orbitals, and any coupled response data in that block."""
        if populations := extract_populations_from_state(state):
            result.update(populations)
        return G16ParsePhase.FREQUENCY

    def _run_frequency_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse vibrational frequency blocks if present."""
        if vibrations := extract_vibrations_from_state(state):
            result.set("vibrations", vibrations)
        return G16ParsePhase.THERMOCHEM

    def _run_thermochem_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse thermochemistry summaries that often follow frequency sections."""
        cursor_before = state.cursor
        if (thermal_info := extract_thermal_infos_from_state(state)) is not None:
            vibrations = result.fields.get("vibrations")
            frequencies = vibrations.get("frequencies") if isinstance(vibrations, Mapping) else None
            temperatures = thermal_info.get("vibrational_temperatures")
            if frequencies is not None and temperatures is not None:
                positive_mode_indices = [
                    mode_index
                    for mode_index, frequency in enumerate(frequencies)
                    if float(frequency.magnitude) > 0
                ]
                if len(positive_mode_indices) == len(temperatures):
                    thermal_info["vibrational_temperature_mode_indices"] = positive_mode_indices
                elif len(frequencies) == len(temperatures):
                    thermal_info["vibrational_temperature_mode_indices"] = list(
                        range(len(frequencies))
                    )
            result.set("thermal_informations", thermal_info)
        if temp_pressure := _temperature_and_pressure_from_block(state.content[cursor_before:]):
            result.update(temp_pressure)
        return G16ParsePhase.FORCES

    def _run_forces_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse Cartesian forces when present."""
        if (forces := extract_forces_from_state(state)) is not None:
            result.set("forces", forces)
            result.set("forces_axis_order", ("atom", "cartesian"))
            result.set("forces_atom_order", "source")
            result.set("forces_orientation", "unknown")
        return G16ParsePhase.HESSIAN

    def _run_hessian_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse Hessian / second-derivative blocks."""
        if (hessian := extract_hessian_from_state(state)) is not None:
            result.set("hessian", hessian)
            result.set("hessian_axis_order", ("atom_cartesian", "atom_cartesian"))
            result.set("hessian_atom_order", "source")
            result.set("hessian_orientation", "unknown")
        return G16ParsePhase.OPTIMIZATION

    def _run_optimization_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Parse Berny optimization state summaries."""
        if (
            berny := extract_berny_from_state(
                state,
                capture_source_evidence=cast(_HasParseMethod, self).capture_source_evidence,
            )
        ) is not None:
            result.set("geometry_optimization_status", berny)
        return G16ParsePhase.ELECTRIC_RESPONSE

    def _run_electric_response_phase(
        self, state: ParseState, result: ModelParseResult
    ) -> G16ParsePhase:
        """Merge late electric-dipole/polarizability sections into existing response data."""
        if (
            polarizability := extract_electric_dipole_and_polarizability_from_state(state)
        ) is not None:
            result.merge_payload_field("polarizability", polarizability, overwrite=True)
        return G16ParsePhase.ARCHIVE_TAIL

    def _run_archive_tail_phase(self, state: ParseState, result: ModelParseResult) -> G16ParsePhase:
        """Use archive-tail data as the final fallback/augmentation stage."""
        archive_payload = extract_archive_tail_payload_from_state(
            state,
            capture_source_evidence=cast(_HasParseMethod, self).capture_source_evidence,
        )
        if tail := archive_payload.get("metadata"):
            result.set_missing_from(tail)

        if tail_energies := archive_payload.get("energies"):
            live_energies = result.fields.get("energies")
            result.set(
                "energies",
                merge_g16_energy_payloads(
                    live_energies if isinstance(live_energies, Mapping) else None,
                    tail_energies,
                ),
            )

        if tail_thermal_info := archive_payload.get("thermal_informations"):
            result.merge_payload_field("thermal_informations", tail_thermal_info)

        if tail_polarizability := archive_payload.get("polarizability"):
            result.merge_payload_field("polarizability", tail_polarizability)

        if (tail_hessian := archive_payload.get("hessian")) is not None and not result.has_value(
            "hessian"
        ):
            result.set("hessian", tail_hessian)
            result.set("hessian_axis_order", ("atom_cartesian", "atom_cartesian"))
            result.set("hessian_atom_order", "source")
            result.set("hessian_orientation", "unknown")

        return G16ParsePhase.DONE

    def _parse_block_to_result(self, block: str) -> ModelParseResult:
        """Execute the explicit phase machine until all extractors have run."""
        state = ParseState(block)
        result = ModelParseResult({"qm_software": "Gaussian"})

        phase = G16ParsePhase.HEADER
        while phase is not G16ParsePhase.DONE:
            if phase is G16ParsePhase.HEADER:
                phase = self._run_header_phase(block, result)
            elif phase is G16ParsePhase.ORIENTATION:
                phase = self._run_orientation_phase(state, result)
            elif phase is G16ParsePhase.STRUCTURE_ONLY_CHECK:
                phase = self._run_structure_only_check()
            elif phase is G16ParsePhase.ROTATION:
                phase = self._run_rotation_phase(state, result)
            elif phase is G16ParsePhase.SCF:
                phase = self._run_scf_phase(state, result)
            elif phase is G16ParsePhase.ISOTROPIC_POLARIZABILITY:
                phase = self._run_isotropic_polarizability_phase(state, result)
            elif phase is G16ParsePhase.NMR:
                phase = self._run_nmr_phase(state, result)
            elif phase is G16ParsePhase.POPULATION:
                phase = self._run_population_phase(state, result)
            elif phase is G16ParsePhase.FREQUENCY:
                phase = self._run_frequency_phase(state, result)
            elif phase is G16ParsePhase.THERMOCHEM:
                phase = self._run_thermochem_phase(state, result)
            elif phase is G16ParsePhase.FORCES:
                phase = self._run_forces_phase(state, result)
            elif phase is G16ParsePhase.HESSIAN:
                phase = self._run_hessian_phase(state, result)
            elif phase is G16ParsePhase.OPTIMIZATION:
                phase = self._run_optimization_phase(state, result)
            elif phase is G16ParsePhase.ELECTRIC_RESPONSE:
                phase = self._run_electric_response_phase(state, result)
            elif phase is G16ParsePhase.ARCHIVE_TAIL:
                phase = self._run_archive_tail_phase(state, result)
            else:
                raise AssertionError(f"Unexpected G16 frame parse phase: {phase!r}")

        return result

    def _parse_frame(self) -> Mapping[str, Any]:
        """Return model-ready frame fields from the canonical state-machine result."""
        typed_self = cast(_HasParseMethod, self)
        result = self._parse_block_to_result(typed_self._block)
        if typed_self.capture_source_evidence:
            result.set("coordinate_source", "observed")
            result.set("coordinate_provenance", "Gaussian source geometry via frame.coords")
            if result.has_value("forces"):
                result.set("force_source_field", "forces")
                result.set("force_transformation", "Gaussian Cartesian forces in source atom order")
        return result.model_data()


class G16LogFileFrameParserMemory(
    G16LogFileFrameParserMixin, BaseFrameParser[G16LogFileFrameMemory]
):
    _file_frame_class_ = G16LogFileFrameMemory


class G16LogFileFrameParserDisk(G16LogFileFrameParserMixin, BaseFrameParser[G16LogFileFrameDisk]):
    _file_frame_class_ = G16LogFileFrameDisk
