from __future__ import annotations

from collections.abc import Mapping
from enum import Enum, auto
from typing import Any, cast

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext, _HasParseMethod
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.xtb.output.frame_models.XTBOutputFileFrame import (
    XTBOutputFileFrameDisk,
    XTBOutputFileFrameMemory,
)
from molop.io.logic.xtb.output.parsers._xtb_output_extractors import (
    extract_xtb_coords,
    extract_xtb_dipole,
    extract_xtb_final_energy,
    extract_xtb_optimization_status,
    extract_xtb_orbitals,
    extract_xtb_populations,
    extract_xtb_rotation_constants,
    extract_xtb_running_time,
    extract_xtb_single_point_properties,
    extract_xtb_status,
    extract_xtb_thermal_information,
    extract_xtb_vibrations,
)


class XTBOutputParsePhase(Enum):
    STRUCTURE = auto()
    STRUCTURE_ONLY_CHECK = auto()
    ENERGY = auto()
    ORBITALS = auto()
    POPULATIONS = auto()
    RESPONSE = auto()
    VIBRATIONS = auto()
    THERMOCHEMISTRY = auto()
    OPTIMIZATION = auto()
    STATUS = auto()
    DONE = auto()


class XTBOutputFileFrameParserMixin:
    def _run_structure_phase(
        self,
        text: str,
        result: ModelParseResult,
    ) -> XTBOutputParsePhase:
        atoms, coords, precision, source_label = extract_xtb_coords(text)
        if atoms is not None and coords is not None:
            result.set("atoms", atoms)
            result.set("coords", coords)
            if cast(_HasParseMethod, self).capture_source_evidence:
                result.set("coordinate_source", "observed")
                result.set("coordinate_provenance", f"xTB {source_label} in source atom order")
                if precision is not None:
                    result.set("coordinate_decimal_places", precision)
        return XTBOutputParsePhase.STRUCTURE_ONLY_CHECK

    def _run_structure_only_check(self) -> XTBOutputParsePhase:
        if cast(_HasParseMethod, self).only_extract_structure:
            return XTBOutputParsePhase.DONE
        return XTBOutputParsePhase.ENERGY

    def _run_energy_phase(
        self,
        text: str,
        result: ModelParseResult,
        context: FrameParseContext,
    ) -> XTBOutputParsePhase:
        typed_self = cast(_HasParseMethod, self)
        model_chemistry = context.additional_data.get("model_chemistry")
        method = getattr(model_chemistry, "method", None)
        energies = extract_xtb_final_energy(
            text,
            method=method,
            capture_source_evidence=typed_self.capture_source_evidence,
        )
        if energies is not None:
            result.set("energies", energies)
        return XTBOutputParsePhase.ORBITALS

    @staticmethod
    def _run_orbitals_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        if orbitals := extract_xtb_orbitals(text):
            result.set("molecular_orbitals", orbitals)
        return XTBOutputParsePhase.POPULATIONS

    @staticmethod
    def _run_populations_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        if populations := extract_xtb_populations(text):
            result.set("charge_spin_populations", populations)
        if properties := extract_xtb_single_point_properties(text):
            result.set("single_point_properties", properties)
        return XTBOutputParsePhase.RESPONSE

    @staticmethod
    def _run_response_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        if dipole := extract_xtb_dipole(text):
            result.set("polarizability", dipole)
        rotation_constants = extract_xtb_rotation_constants(text)
        if rotation_constants is not None:
            result.set("rotation_constants", rotation_constants)
        return XTBOutputParsePhase.VIBRATIONS

    @staticmethod
    def _run_vibration_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        if vibrations := extract_xtb_vibrations(text):
            result.set("vibrations", vibrations)
        return XTBOutputParsePhase.THERMOCHEMISTRY

    @staticmethod
    def _run_thermochemistry_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        if thermal := extract_xtb_thermal_information(text):
            result.set("thermal_informations", thermal)
        return XTBOutputParsePhase.OPTIMIZATION

    @staticmethod
    def _run_optimization_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        optimization, gradient_norm, gradient_threshold = extract_xtb_optimization_status(text)
        if optimization is not None:
            result.set("geometry_optimization_status", optimization)
        if gradient_norm is not None:
            result.set("gradient_norm", gradient_norm)
        if gradient_threshold is not None:
            result.set("gradient_norm_threshold", gradient_threshold)
        return XTBOutputParsePhase.STATUS

    @staticmethod
    def _run_status_phase(text: str, result: ModelParseResult) -> XTBOutputParsePhase:
        result.set("status", extract_xtb_status(text))
        if running_time := extract_xtb_running_time(text):
            result.set("running_time", running_time)
        return XTBOutputParsePhase.DONE

    def _parse_block_to_result(
        self,
        text: str,
        context: FrameParseContext | None = None,
    ) -> ModelParseResult:
        context = context or FrameParseContext(additional_data={})
        result = ModelParseResult({"qm_software": "xTB"})
        phase = XTBOutputParsePhase.STRUCTURE
        while phase is not XTBOutputParsePhase.DONE:
            if phase is XTBOutputParsePhase.STRUCTURE:
                phase = self._run_structure_phase(text, result)
            elif phase is XTBOutputParsePhase.STRUCTURE_ONLY_CHECK:
                phase = self._run_structure_only_check()
            elif phase is XTBOutputParsePhase.ENERGY:
                phase = self._run_energy_phase(text, result, context)
            elif phase is XTBOutputParsePhase.ORBITALS:
                phase = self._run_orbitals_phase(text, result)
            elif phase is XTBOutputParsePhase.POPULATIONS:
                phase = self._run_populations_phase(text, result)
            elif phase is XTBOutputParsePhase.RESPONSE:
                phase = self._run_response_phase(text, result)
            elif phase is XTBOutputParsePhase.VIBRATIONS:
                phase = self._run_vibration_phase(text, result)
            elif phase is XTBOutputParsePhase.THERMOCHEMISTRY:
                phase = self._run_thermochemistry_phase(text, result)
            elif phase is XTBOutputParsePhase.OPTIMIZATION:
                phase = self._run_optimization_phase(text, result)
            elif phase is XTBOutputParsePhase.STATUS:
                phase = self._run_status_phase(text, result)
            else:
                raise AssertionError(f"Unexpected xTB output parse phase: {phase!r}")
        return result

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        return self._parse_block_to_result(block, context).model_data()


class XTBOutputFileFrameParserMemory(
    XTBOutputFileFrameParserMixin,
    BaseFrameParser[XTBOutputFileFrameMemory],
):
    _file_frame_class_ = XTBOutputFileFrameMemory


class XTBOutputFileFrameParserDisk(
    XTBOutputFileFrameParserMixin,
    BaseFrameParser[XTBOutputFileFrameDisk],
):
    _file_frame_class_ = XTBOutputFileFrameDisk


__all__ = ["XTBOutputFileFrameParserDisk", "XTBOutputFileFrameParserMemory"]
