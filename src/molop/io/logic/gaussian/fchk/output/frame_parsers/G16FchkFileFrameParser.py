from __future__ import annotations

from collections.abc import Mapping
from enum import Enum, auto
from typing import Any, cast

from molop.io.base_models.FrameParser import BaseFrameParser, FrameParseContext, _HasParseMethod
from molop.io.base_models.ParseContainers import ModelParseResult
from molop.io.logic.gaussian.fchk.output.frame_models.G16FchkFileFrame import (
    G16FchkFileFrameDisk,
    G16FchkFileFrameMemory,
)
from molop.io.logic.gaussian.fchk.output.parsers._fchk_extractors import (
    extract_fchk_energies,
    extract_fchk_forces,
    extract_fchk_hessian,
    extract_fchk_nmr,
    extract_fchk_optimization_status,
    extract_fchk_orbitals,
    extract_fchk_polarizability,
    extract_fchk_populations,
    extract_fchk_status,
    extract_fchk_structure,
    extract_fchk_thermal_information,
    extract_fchk_total_spin,
    extract_fchk_vibrations,
    parse_fchk_frame_records,
)


class G16FchkParsePhase(Enum):
    STRUCTURE = auto()
    STRUCTURE_ONLY_CHECK = auto()
    ENERGY = auto()
    DERIVATIVES = auto()
    ORBITALS = auto()
    POPULATIONS = auto()
    RESPONSE = auto()
    NMR = auto()
    VIBRATIONS = auto()
    THERMOCHEMISTRY = auto()
    STATUS = auto()
    DONE = auto()


class G16FchkFileFrameParserMixin:
    @staticmethod
    def _num_atoms(result: ModelParseResult) -> int:
        atoms = result.fields.get("atoms")
        return len(atoms) if isinstance(atoms, list) else 0

    def _run_structure_phase(
        self, records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        atoms, coords = extract_fchk_structure(records)
        if atoms is not None and coords is not None:
            result.set("atoms", atoms)
            result.set("coords", coords)
            if cast(_HasParseMethod, self).capture_source_evidence:
                result.set("coordinate_source", "observed")
                result.set(
                    "coordinate_provenance",
                    "Gaussian fchk Current cartesian coordinates in source atom order",
                )
        return G16FchkParsePhase.STRUCTURE_ONLY_CHECK

    def _run_structure_only_check(self) -> G16FchkParsePhase:
        if cast(_HasParseMethod, self).only_extract_structure:
            return G16FchkParsePhase.DONE
        return G16FchkParsePhase.ENERGY

    def _run_energy_phase(
        self,
        records: Mapping[str, Any],
        result: ModelParseResult,
        context: FrameParseContext,
    ) -> G16FchkParsePhase:
        typed_self = cast(_HasParseMethod, self)
        model_chemistry = context.additional_data.get("model_chemistry")
        method = getattr(model_chemistry, "method", None)
        energies = extract_fchk_energies(
            records,
            method=method,
            capture_source_evidence=typed_self.capture_source_evidence,
        )
        if energies is not None:
            result.set("energies", energies)
        return G16FchkParsePhase.DERIVATIVES

    def _run_derivatives_phase(
        self, records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        num_atoms = self._num_atoms(result)
        forces = extract_fchk_forces(records, num_atoms)
        if forces is not None:
            result.set("forces", forces)
            result.set("forces_axis_order", ("atom", "cartesian"))
            result.set("forces_atom_order", "source")
            result.set("forces_orientation", "source")
            if cast(_HasParseMethod, self).capture_source_evidence:
                result.set("force_source_field", "Cartesian Gradient")
                result.set(
                    "force_transformation",
                    "Elementwise negation of the Gaussian fchk Cartesian Gradient",
                )
        hessian = extract_fchk_hessian(records, num_atoms)
        if hessian is not None:
            result.set("hessian", hessian)
            result.set("hessian_axis_order", ("atom_cartesian", "atom_cartesian"))
            result.set("hessian_atom_order", "source")
            result.set("hessian_orientation", "source")
        return G16FchkParsePhase.ORBITALS

    @staticmethod
    def _run_orbitals_phase(
        records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        if orbitals := extract_fchk_orbitals(records):
            result.set("molecular_orbitals", orbitals)
        if total_spin := extract_fchk_total_spin(records):
            result.set("total_spin", total_spin)
        return G16FchkParsePhase.POPULATIONS

    def _run_populations_phase(
        self, records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        if populations := extract_fchk_populations(records, self._num_atoms(result)):
            result.set("charge_spin_populations", populations)
        return G16FchkParsePhase.RESPONSE

    @staticmethod
    def _run_response_phase(
        records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        if polarizability := extract_fchk_polarizability(records):
            result.set("polarizability", polarizability)
        return G16FchkParsePhase.NMR

    def _run_nmr_phase(
        self,
        records: Mapping[str, Any],
        result: ModelParseResult,
        context: FrameParseContext,
    ) -> G16FchkParsePhase:
        atoms = result.fields.get("atoms")
        route = context.additional_data.get("keywords", "")
        if (
            isinstance(atoms, list)
            and (nmr := extract_fchk_nmr(records, atoms, route=str(route))) is not None
        ):
            result.set("nmr", nmr)
        return G16FchkParsePhase.VIBRATIONS

    def _run_vibrations_phase(
        self, records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        if vibrations := extract_fchk_vibrations(records, self._num_atoms(result)):
            result.set("vibrations", vibrations)
        return G16FchkParsePhase.THERMOCHEMISTRY

    @staticmethod
    def _run_thermochemistry_phase(
        records: Mapping[str, Any], result: ModelParseResult
    ) -> G16FchkParsePhase:
        if thermal := extract_fchk_thermal_information(records):
            result.set("thermal_informations", thermal)
        return G16FchkParsePhase.STATUS

    def _run_status_phase(
        self,
        records: Mapping[str, Any],
        result: ModelParseResult,
        context: FrameParseContext,
    ) -> G16FchkParsePhase:
        result.set("status", extract_fchk_status(records))
        task_requests = context.additional_data.get("task_requests", [])
        optimization = extract_fchk_optimization_status(records, list(task_requests))
        if optimization is not None:
            result.set("geometry_optimization_status", optimization)
        return G16FchkParsePhase.DONE

    def _parse_block_to_result(
        self,
        text: str,
        context: FrameParseContext | None = None,
    ) -> ModelParseResult:
        context = context or FrameParseContext(additional_data={})
        records = parse_fchk_frame_records(text)
        result = ModelParseResult({"qm_software": "Gaussian"})
        phase = G16FchkParsePhase.STRUCTURE
        while phase is not G16FchkParsePhase.DONE:
            if phase is G16FchkParsePhase.STRUCTURE:
                phase = self._run_structure_phase(records, result)
            elif phase is G16FchkParsePhase.STRUCTURE_ONLY_CHECK:
                phase = self._run_structure_only_check()
            elif phase is G16FchkParsePhase.ENERGY:
                phase = self._run_energy_phase(records, result, context)
            elif phase is G16FchkParsePhase.DERIVATIVES:
                phase = self._run_derivatives_phase(records, result)
            elif phase is G16FchkParsePhase.ORBITALS:
                phase = self._run_orbitals_phase(records, result)
            elif phase is G16FchkParsePhase.POPULATIONS:
                phase = self._run_populations_phase(records, result)
            elif phase is G16FchkParsePhase.RESPONSE:
                phase = self._run_response_phase(records, result)
            elif phase is G16FchkParsePhase.NMR:
                phase = self._run_nmr_phase(records, result, context)
            elif phase is G16FchkParsePhase.VIBRATIONS:
                phase = self._run_vibrations_phase(records, result)
            elif phase is G16FchkParsePhase.THERMOCHEMISTRY:
                phase = self._run_thermochemistry_phase(records, result)
            elif phase is G16FchkParsePhase.STATUS:
                phase = self._run_status_phase(records, result, context)
            else:
                raise AssertionError(f"Unexpected Gaussian fchk parse phase: {phase!r}")
        return result

    def _parse_frame(self, block: str, *, context: FrameParseContext) -> Mapping[str, Any]:
        return self._parse_block_to_result(block, context).model_data()


class G16FchkFileFrameParserMemory(
    G16FchkFileFrameParserMixin,
    BaseFrameParser[G16FchkFileFrameMemory],
):
    _file_frame_class_ = G16FchkFileFrameMemory


class G16FchkFileFrameParserDisk(
    G16FchkFileFrameParserMixin,
    BaseFrameParser[G16FchkFileFrameDisk],
):
    _file_frame_class_ = G16FchkFileFrameDisk


__all__ = ["G16FchkFileFrameParserDisk", "G16FchkFileFrameParserMemory"]
