"""
Author: TMJ
Date: 2025-08-01 16:13:58
LastEditors: TMJ
LastEditTime: 2025-11-20 13:58:47
Description: 请填写简介
"""

from typing import Any, Optional, Sequence

import numpy as np
from pint.facets.plain import PlainQuantity

from molop.config import moloplogger
from molop.io.base_models.DataClasses import ImplicitSolvation, Status
from molop.io.base_models.FileParser import BaseFileParserDisk, BaseFileParserMemory
from molop.io.patterns.G16Patterns import MolOPPattern, g16_log_patterns
from molop.io.QM_models.G16LogFile import G16LogFileDisk, G16LogFileMemory
from molop.io.QM_models.G16LogFileFrame import (
    G16LogFileFrameDisk,
    G16LogFileFrameMemory,
)
from molop.io.QM_parsers.G16LogFileFrameParser import (
    G16LogFileFrameParserDisk,
    G16LogFileFrameParserMemory,
)
from molop.unit import atom_ureg
from molop.utils.functions import find_rigid_transform

split_pattern = "Input orientation:"
split_pattern_2 = "Standard orientation:"


class G16LogFileParserMixin:
    _file_content: str
    _file_path: str
    only_extract_structure: bool

    def _parse_metadata(self, file_content: str) -> dict[str, Any]:
        metadata: dict[str, Any] = {"qm_software": "Gaussian"}
        self._file_content = file_content
        if version := self._parse_version():
            metadata["qm_software_version"] = version
        if options := self._parse_options():
            metadata["options"] = options
        if route := self._parse_keywords():
            metadata["keywords"] = route
        if title := self._parse_title():
            metadata["title_card"] = title
        if charge_multiplicity := self._parse_charge_multiplicity():
            metadata["charge"], metadata["multiplicity"] = charge_multiplicity
        transform_matrix = self._parse_standard_orientation_transformation_matrix()
        if transform_matrix is not None:
            metadata["standard_orientation_transformation_matrix"] = transform_matrix
        if self.only_extract_structure:
            return metadata
        if solvent := self._parse_solvent():
            metadata["solvent"] = solvent
        temperature, pressure = self._parse_temperature_and_pressure()
        if temperature:
            metadata["temperature"] = temperature
        if pressure:
            metadata["pressure"] = pressure
        if tail := self._parse_tail():
            metadata.update(tail)
        if running_time := self._parse_running_time():
            metadata["running_time"] = running_time
        if status := self._parse_termination_status():
            metadata["status"] = status
        return metadata

    def _split_file(self, file_content: str) -> Sequence[str]:
        if split_pattern in file_content:
            split_ = split_pattern
        elif split_pattern_2 in file_content:
            split_ = split_pattern_2
        fragments = file_content.split(split_)[1:]
        frame_contents = [f"{split_}\n{fragment}" for fragment in fragments]
        return frame_contents

    def _parse_version(self) -> str | None:
        focus_content, self._file_content = g16_log_patterns.VERSION.split_content(
            self._file_content
        )
        if matches := g16_log_patterns.VERSION.get_matches(focus_content):
            return matches[0][0]
        return None

    def _parse_options(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.OPTIONS.split_content(
            self._file_content
        )
        if matches := g16_log_patterns.OPTIONS.find_iter(focus_content):
            options = "\n".join(
                [f"{match.groups()[0]}={match.groups()[1]}" for match in matches]
            )
            self._file_content = continued_content
            return options
        return None

    def _parse_keywords(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.KEYWORDS.split_content(
            self._file_content
        )
        if len(keyword_lines := focus_content.splitlines()) >= 3:
            self._file_content = continued_content
            return "\n".join(keyword_lines[1:-1]).replace("\n ", "")
        return None

    def _parse_title(self) -> str | None:
        focus_content, continued_content = g16_log_patterns.TITLE.split_content(
            self._file_content
        )
        if len(title_lines := focus_content.splitlines()) >= 3:
            self._file_content = continued_content
            return "\n".join(title_lines[1:-1]).replace("\n ", "")
        return None

    def _parse_charge_multiplicity(self) -> tuple[int, int] | tuple[None, None]:
        if match := g16_log_patterns.CHARGE_MULTIPLICITY.match_content(
            self._file_content
        ):
            charge, multiplicity = int(match[0][0]), int(match[0][1])
            return charge, multiplicity
        return None, None

    def _parse_standard_orientation_transformation_matrix(self) -> Optional[np.ndarray]:
        try:
            coords = self._parse_coordinates(g16_log_patterns.INITIAL_INPUT_COORDS)
            if coords is None:
                return None
            standard_coords = self._parse_coordinates(g16_log_patterns.STANDARD_COORDS)
            if standard_coords is None:
                return None
            return find_rigid_transform(coords, standard_coords)
        except Exception as e:
            raise e

    def _parse_coordinates(self, pattern: MolOPPattern) -> Optional[np.ndarray]:
        focus_content, self._file_content = pattern.split_content(self._file_content)
        if matches := pattern.match_content(focus_content):
            return np.array(
                [list(map(float, match[1:])) for match in matches],
            )
        return None

    def _parse_solvent(self) -> Optional[ImplicitSolvation]:
        solvent_dict: dict[str, Any] = {}
        focus_content, self._file_content = (
            g16_log_patterns.SOLVENT_PARAMETERS.split_content(self._file_content)
        )
        if matches := g16_log_patterns.SOLVENT_MODEL.get_matches(focus_content):
            solvent_dict["solvent_model"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_ATOM_RADII.get_matches(focus_content):
            solvent_dict["atomic_radii"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_TYPE.get_matches(focus_content):
            solvent_dict["solvent"] = matches[0][0]
        if matches := g16_log_patterns.SOLVENT_EPS.get_matches(focus_content):
            solvent_dict["solvent_epsilon"] = float(matches[0][0])
        if matches := g16_log_patterns.SOLVENT_EPS_INF.get_matches(focus_content):
            solvent_dict["solvent_epsilon_infinite"] = float(matches[0][-1])
        if solvent_dict:
            return ImplicitSolvation.model_validate(solvent_dict)
        return None

    def _parse_temperature_and_pressure(
        self,
    ) -> tuple[PlainQuantity, PlainQuantity] | tuple[None, None]:
        index = self._file_content.find(" - Thermochemistry -")
        if index == -1:
            return None, None
        focus_content = self._file_content[index:]
        # self._file_content = self._file_content[index:]
        if match := g16_log_patterns.TEMPEREATURE_PRESSURE.match_content(focus_content):
            temperature, pressure = (
                float(match[0][0]) * atom_ureg.K,
                float(match[0][1]) * atom_ureg.atm,
            )
            return temperature, pressure
        return None, None

    def _parse_tail(self) -> dict[str, Any]:
        tail_dict: dict[str, Any] = {}
        start_index = self._file_content.find("1\\1\\GINC")
        if start_index == -1:
            start_index = self._file_content.find("1|1|GINC")
        if start_index == -1:
            return tail_dict
        end_index = self._file_content.find("\\@\n")
        if end_index == -1:
            end_index = self._file_content.find("@\n")
        if end_index == -1:
            return tail_dict
        focus_content = self._file_content[start_index:end_index]
        self._file_content = self._file_content[end_index:]
        focus_content = focus_content.replace("\n ", "")

        def parse_and_update(pattern: MolOPPattern, key: str):
            nonlocal focus_content
            try:
                sub_focus_content, sub_continued_content = pattern.split_content(
                    focus_content
                )
                moloplogger.debug(
                    f"{key} focus content: \n{sub_focus_content}\n{key} focus content end"
                    f"\ncorresponding file: {self._file_path}"
                )
                if matches := pattern.get_matches(sub_focus_content):
                    tail_dict[key] = matches[0][0]
                    focus_content = sub_continued_content
            except Exception as e:
                moloplogger.error(f"Error in parsing {key}: {e}")

        parse_and_update(g16_log_patterns.JOB_TYPE_IN_ARCHIVE_TAIL, "job_type")
        parse_and_update(g16_log_patterns.FUNCTIONAL_IN_ARCHIVE_TAIL, "functional")
        parse_and_update(g16_log_patterns.BASIS_SET_IN_ARCHIVE_TAIL, "basis")
        parse_and_update(g16_log_patterns.KEYWORDS_IN_ARCHIVE_TAIL, "keywords")
        parse_and_update(g16_log_patterns.TITLE_IN_ARCHIVE_TAIL, "title_card")

        sub_focus_content, sub_continued_content = (
            g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.split_content(
                focus_content
            )
        )
        if (
            matches
            := g16_log_patterns.CHARGE_SPIN_MULTIPLICITY_IN_ARCHIVE_TAIL.get_matches(
                sub_focus_content
            )
        ):
            charge_multiplicity = matches[0]
            tail_dict["charge"] = int(charge_multiplicity[0])
            tail_dict["multiplicity"] = int(charge_multiplicity[1])
            focus_content = sub_continued_content

        parse_and_update(
            g16_log_patterns.VERSION_IN_ARCHIVE_TAIL, "qm_software_version"
        )

        moloplogger.debug(f"parsed tail_dict: {tail_dict}")
        return tail_dict

    def _parse_running_time(self) -> Optional[PlainQuantity]:
        if matches := g16_log_patterns.JOB_TIME.match_content(self._file_content):
            total_seconds = 0
            for match in matches:
                days, hours, minutes, seconds = [float(i) for i in match[1:]]
                total_seconds += (days * 24 + hours) * 3600 + minutes * 60 + seconds
            return total_seconds * atom_ureg.second
        return None

    def _parse_termination_status(self) -> Optional[Status]:
        status_dict: dict[str, Any] = {}
        if matches := g16_log_patterns.TERMINATION_STATUS.match_content(
            self._file_content
        ):
            status = matches[-1][0]
            if status == "Normal":
                status_dict["normal_terminated"] = True
                status_dict["scf_converged"] = True
            else:
                status_dict["normal_terminated"] = False
        if status_dict:
            return Status.model_validate(status_dict)
        return None


class G16LogFileParserMemory(
    G16LogFileParserMixin,
    BaseFileParserMemory[
        G16LogFileMemory, G16LogFileFrameMemory, G16LogFileFrameParserMemory
    ],
):
    _chem_file = G16LogFileMemory
    _frame_parser = G16LogFileFrameParserMemory


class G16LogFileParserDisk(
    G16LogFileParserMixin,
    BaseFileParserDisk[G16LogFileDisk, G16LogFileFrameDisk, G16LogFileFrameParserDisk],
):
    _chem_file = G16LogFileDisk
    _frame_parser = G16LogFileFrameParserDisk
