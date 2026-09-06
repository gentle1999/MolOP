from __future__ import annotations

from collections.abc import Mapping
from enum import Enum, auto
from typing import Any, cast

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext, _HasParseMethod
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.orca.log.frame_models.ORCALogFileFrame import (
    ORCALogFileFrameDisk,
    ORCALogFileFrameMemory,
)
from molop.io.logic.orca.log.frame_parsers._orca_extractors import (
    extract_orca_atomic_masses,
    extract_orca_coords,
    extract_orca_electronic_states,
    extract_orca_energies,
    extract_orca_forces,
    extract_orca_geometry_optimization_status,
    extract_orca_polarizability,
    extract_orca_populations,
    extract_orca_solvent,
    extract_orca_vibrations,
)
from molop.io.logic.orca.log.parsers._orca_log_shared import (
    extract_orca_running_time,
    extract_orca_status,
)


class ORCALogParsePhase(Enum):
    """Explicit stages for ORCA output frame parsing."""

    STRUCTURE = auto()
    ENERGY = auto()
    STRUCTURE_ONLY_CHECK = auto()
    GRADIENT = auto()
    VIBRATION = auto()
    POPULATION = auto()
    POLARIZABILITY = auto()
    STATUS = auto()
    OPTIMIZATION = auto()
    SOLVENT = auto()
    ELECTRONIC_STATES = auto()
    DONE = auto()


class ORCALogFileFrameParserMixin:
    """Explicit-state ORCA output frame parser."""

    @staticmethod
    def _num_atoms_from_result(result: ModelParseResult) -> int | None:
        atoms = result.fields.get("atoms")
        if isinstance(atoms, list):
            return len(atoms)
        return None

    def _run_structure_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        capture_source_evidence = cast(_HasParseMethod, self).capture_source_evidence
        atoms, coords, coordinate_decimal_places = extract_orca_coords(
            text,
            capture_source_evidence=capture_source_evidence,
        )
        if atoms is not None and coords is not None:
            result.set("atoms", atoms)
            result.set("coords", coords)
            if (
                atomic_masses := extract_orca_atomic_masses(
                    text,
                    expected_atom_count=len(atoms),
                )
            ) is not None:
                result.set("atomic_masses", atomic_masses)
                result.set("atomic_masses_source", "orca_cartesian_au_mass")
            if capture_source_evidence:
                result.set("coordinate_source", "observed")
                result.set(
                    "coordinate_provenance",
                    "ORCA CARTESIAN COORDINATES (ANGSTROEM) in source atom order",
                )
                if coordinate_decimal_places is not None:
                    result.set("coordinate_decimal_places", coordinate_decimal_places)
        return ORCALogParsePhase.ENERGY

    def _run_energy_phase(
        self,
        text: str,
        result: ModelParseResult,
        context: FrameParseContext,
    ) -> ORCALogParsePhase:
        typed_self = cast(_HasParseMethod, self)
        energies = extract_orca_energies(
            text,
            capture_source_evidence=typed_self.capture_source_evidence,
            model_chemistry=context.additional_data.get("model_chemistry"),
        )
        if energies is not None:
            result.set("energies", energies)
        return ORCALogParsePhase.STRUCTURE_ONLY_CHECK

    def _run_structure_only_check(self) -> ORCALogParsePhase:
        return (
            ORCALogParsePhase.DONE
            if cast(_HasParseMethod, self).only_extract_structure
            else ORCALogParsePhase.GRADIENT
        )

    def _run_gradient_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        forces = extract_orca_forces(text, self._num_atoms_from_result(result))
        if forces is not None:
            result.set("forces", forces)
            result.set("forces_axis_order", ("atom", "cartesian"))
            result.set("forces_atom_order", "source")
            result.set("forces_orientation", "source")
            if cast(_HasParseMethod, self).capture_source_evidence:
                result.set("force_source_field", "gradient")
                result.set(
                    "force_transformation",
                    "Elementwise negation of the ORCA Cartesian gradient in source atom order",
                )
        return ORCALogParsePhase.VIBRATION

    def _run_vibration_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        vibrations = extract_orca_vibrations(text, self._num_atoms_from_result(result))
        if vibrations is not None:
            result.set("vibrations", vibrations)
        return ORCALogParsePhase.POPULATION

    def _run_population_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        populations = extract_orca_populations(text)
        if populations is not None:
            result.set("charge_spin_populations", populations)
        return ORCALogParsePhase.POLARIZABILITY

    def _run_polarizability_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        polarizability = extract_orca_polarizability(text)
        if polarizability is not None:
            result.set("polarizability", polarizability)
        return ORCALogParsePhase.STATUS

    def _run_status_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        status = extract_orca_status(text, include_termination=False)
        if status is not None:
            result.set("status", status)
        running_time = extract_orca_running_time(text)
        if running_time is not None:
            result.set("running_time", running_time)
        return ORCALogParsePhase.OPTIMIZATION

    def _run_optimization_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        opt_status = extract_orca_geometry_optimization_status(
            text,
            capture_source_evidence=cast(_HasParseMethod, self).capture_source_evidence,
        )
        if opt_status is not None:
            result.set("geometry_optimization_status", opt_status)
        return ORCALogParsePhase.SOLVENT

    def _run_solvent_phase(self, text: str, result: ModelParseResult) -> ORCALogParsePhase:
        solvent = extract_orca_solvent(text)
        if solvent is not None:
            result.set("solvent", solvent)
        return ORCALogParsePhase.ELECTRONIC_STATES

    def _run_electronic_states_phase(
        self, text: str, result: ModelParseResult
    ) -> ORCALogParsePhase:
        method = getattr(self, "_orca_method", None)
        electronic_states = extract_orca_electronic_states(text, method)
        if electronic_states is not None:
            result.set("electronic_states", electronic_states)
        return ORCALogParsePhase.DONE

    def _parse_block_to_result(
        self,
        text: str,
        context: FrameParseContext | None = None,
    ) -> ModelParseResult:
        context = context or FrameParseContext(additional_data={})
        result = ModelParseResult({"qm_software": "ORCA"})
        phase = ORCALogParsePhase.STRUCTURE
        while phase is not ORCALogParsePhase.DONE:
            if phase is ORCALogParsePhase.STRUCTURE:
                phase = self._run_structure_phase(text, result)
            elif phase is ORCALogParsePhase.ENERGY:
                phase = self._run_energy_phase(text, result, context)
            elif phase is ORCALogParsePhase.STRUCTURE_ONLY_CHECK:
                phase = self._run_structure_only_check()
            elif phase is ORCALogParsePhase.GRADIENT:
                phase = self._run_gradient_phase(text, result)
            elif phase is ORCALogParsePhase.VIBRATION:
                phase = self._run_vibration_phase(text, result)
            elif phase is ORCALogParsePhase.POPULATION:
                phase = self._run_population_phase(text, result)
            elif phase is ORCALogParsePhase.POLARIZABILITY:
                phase = self._run_polarizability_phase(text, result)
            elif phase is ORCALogParsePhase.STATUS:
                phase = self._run_status_phase(text, result)
            elif phase is ORCALogParsePhase.OPTIMIZATION:
                phase = self._run_optimization_phase(text, result)
            elif phase is ORCALogParsePhase.SOLVENT:
                phase = self._run_solvent_phase(text, result)
            elif phase is ORCALogParsePhase.ELECTRONIC_STATES:
                phase = self._run_electronic_states_phase(text, result)
            else:
                raise AssertionError(f"Unexpected ORCA log frame parse phase: {phase!r}")
        return result

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        return self._parse_block_to_result(block, context).model_data()


class ORCALogFileFrameParserMemory(
    ORCALogFileFrameParserMixin, BaseFrameParser[ORCALogFileFrameMemory]
):
    _file_frame_class_ = ORCALogFileFrameMemory


class ORCALogFileFrameParserDisk(
    ORCALogFileFrameParserMixin, BaseFrameParser[ORCALogFileFrameDisk]
):
    _file_frame_class_ = ORCALogFileFrameDisk
