"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2026-06-18 19:45:06
Description: 请填写简介
"""

import logging  # noqa: I001
import os
import sys
from contextlib import suppress
from typing import Any, Literal

from joblib import cpu_count as joblib_cpu_count

# RDKit must initialize before Open Babel; loading Open Babel first causes an
# ELF symbol collision between the bundled native libraries.
# isort: off
from pydantic import BaseModel, ConfigDict, Field, ValidationError
from rdkit import RDLogger
from rdkit_dof import dofconfig
from openbabel import pybel
# isort: on


RDLogger.DisableLog("rdApp.*")  # type: ignore
pybel.ob.obErrorLog.StopLogging()
moloplogger = logging.getLogger("molop")
moloplogger.propagate = False
file_handler = logging.FileHandler("molop.log")
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
sh_formatter = logging.Formatter("%(levelname)s - %(message)s")
stream_handler.setFormatter(sh_formatter)

moloplogger.setLevel(logging.INFO)

MAX_JOBS_ENV_VAR = "MOLOP_MAX_JOBS"


def available_cpu_count() -> int:
    """Return the logical CPUs currently available to this process."""

    counts = [joblib_cpu_count(), os.cpu_count()]
    if hasattr(os, "sched_getaffinity"):
        with suppress(OSError):
            counts.append(len(os.sched_getaffinity(0)))
    positive_counts = [int(count) for count in counts if count is not None and count > 0]
    return min(positive_counts, default=1)


class MolOPConfig(BaseModel):
    """
    Configuration class for MolOP operations.
    Used to manage settings related to molecule processing, fingerprint generation, and logging.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, validate_assignment=True)

    # --- General Settings ---
    show_progress_bar: bool = Field(default=True, description="Whether to display the progress bar")
    max_jobs: int | None = Field(
        default=None,
        ge=1,
        description="Maximum parallel jobs; None uses all CPUs available to the process",
    )

    # --- Advanced Settings ---
    graph_reconstruction_backend: Literal["cpp", "python"] = Field(
        default="cpp", description="Backend for graph reconstruction"
    )
    make_dative_bonds: bool = Field(default=True, description="Whether to make dative bonds")
    force_unit_transform: bool = Field(
        default=False, description="Whether to force unit conversion"
    )
    parallel_max_size: int = Field(
        default=8 * 1024**2,
        description="Maximum data size for parallel processing (bytes)",
    )
    max_recursion_depth: int = Field(default=3000, description="Maximum recursion depth for Python")

    # --- Log File Control ---
    log_to_file: bool = Field(default=False, description="Whether to write log messages to a file")

    # --- DOF Effect Drawer Control ---
    use_dof_effect_drawer: bool = Field(
        default=True, description="Whether to use DOF effect drawer"
    )

    def __init__(self, **data: Any):
        """
        Initializes the configuration object and configures the logger based on current settings.
        """
        if "max_jobs" not in data and (env_max_jobs := os.environ.get(MAX_JOBS_ENV_VAR)):
            data["max_jobs"] = env_max_jobs
        super().__init__(**data)
        sys.setrecursionlimit(self.max_recursion_depth)
        # Set log state based on initial configuration values
        if self.show_progress_bar:
            self.verbose()
        else:
            self.quiet()

        if self.log_to_file:
            self.enable_file_logging()
        else:
            self.disable_file_logging()

        if self.use_dof_effect_drawer:
            self.set_dof_effect_drawer(enable=True)
        else:
            self.set_dof_effect_drawer(enable=False)

    def quiet(self):
        """
        Disables the progress bar and console log output.
        This allows the program to run silently in the background.
        """
        self.show_progress_bar = False
        if stream_handler in moloplogger.handlers:
            moloplogger.removeHandler(stream_handler)

    def verbose(self):
        """
        Enables the progress bar and console log output.
        """
        self.show_progress_bar = True
        if stream_handler not in moloplogger.handlers:
            moloplogger.addHandler(stream_handler)

    def enable_file_logging(self):
        """Enable logging to a file."""
        self.log_to_file = True
        if file_handler not in moloplogger.handlers:
            moloplogger.addHandler(file_handler)
        logging.info(
            f"File logging enabled. Logs will be written to: {getattr(file_handler, 'baseFilename', 'N/A')}"
        )

    def disable_file_logging(self):
        """Disable logging to a file."""
        self.log_to_file = False
        if file_handler in moloplogger.handlers:
            moloplogger.removeHandler(file_handler)

    def set_log_level(self, level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]):
        """
        Sets the level for the molop logger.

        Args:
            level (str): The logging level, must be one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'.
        """
        try:
            moloplogger.setLevel(level)
            logging.info(f"Log level set to {level}")
        except ValueError as e:
            logging.error(f"Error setting log level: Invalid level '{level}'")
            raise ValueError(f"Invalid log level: {level}") from e

    def set_max_recursion_depth(self, depth: int):
        """
        Set the maximum recursion depth for Python.

        Args:
            depth (int): The maximum recursion depth.
        """
        sys.setrecursionlimit(depth)
        self.max_recursion_depth = depth
        logging.info(f"Maximum recursion depth set to {depth}")

    @property
    def effective_max_jobs(self) -> int:
        """Current process-aware worker limit after applying ``max_jobs``."""

        available_jobs = available_cpu_count()
        return available_jobs if self.max_jobs is None else min(available_jobs, self.max_jobs)

    def set_n_jobs(self, n_jobs: int) -> int:
        """Resolve automatic or explicit parallelism within the effective worker limit."""

        return self.effective_max_jobs if n_jobs <= 0 else min(n_jobs, self.effective_max_jobs)

    def set_dof_effect_drawer(self, enable: bool):
        """
        Set whether to use the DOF effect drawer.

        Args:
            enable (bool): Whether to use the DOF effect drawer.
        """
        dofconfig.enable_ipython_integration(enable)


# --- Global Configuration Instance ---
# Create a globally available configuration instance
try:
    molopconfig = MolOPConfig()
    molopconfig.set_log_level("INFO")
except ValidationError as e:
    logging.error(f"Configuration validation failed: {e}")
    molopconfig = MolOPConfig(max_jobs=None)
    molopconfig.quiet()
    molopconfig.disable_file_logging()
    logging.critical(
        "Default safe configuration has been used. Please check and fix your custom configuration."
    )
