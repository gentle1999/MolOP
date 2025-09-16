"""
Author: TMJ
Date: 2025-01-15 23:01:22
LastEditors: TMJ
LastEditTime: 2025-09-15 22:48:01
Description: 请填写简介
"""

import logging
import multiprocessing
import sys
from typing import Literal

from openbabel import pybel
from pydantic import BaseModel, ConfigDict, Field, ValidationError
from rdkit import RDLogger
from rdkit.Chem.rdFingerprintGenerator import (
    FingerprintGenerator64,
    GetAtomPairGenerator,
    GetMorganGenerator,
    GetRDKitFPGenerator,
    GetTopologicalTorsionGenerator,
)

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


class MolOPConfig(BaseModel):
    """
    Configuration class for MolOP operations.
    Used to manage settings related to molecule processing, fingerprint generation, and logging.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # --- General Settings ---
    show_progress_bar: bool = Field(
        default=True, description="Whether to display the progress bar"
    )
    max_jobs: int = Field(default=16, description="Maximum number of parallel jobs")

    # --- Fingerprint Generator Settings ---
    morgan_fpgen: FingerprintGenerator64 = Field(
        default=GetMorganGenerator(radius=3, fpSize=1024, includeChirality=True),
        description="Morgan fingerprint generator",
    )
    atompair_fpgen: FingerprintGenerator64 = Field(
        default=GetAtomPairGenerator(fpSize=1024, includeChirality=True),
        description="AtomPair fingerprint generator",
    )
    rdkit_fpgen: FingerprintGenerator64 = Field(
        default=GetRDKitFPGenerator(fpSize=1024),
        description="RDKit fingerprint generator",
    )
    topological_torsion_fpgen: FingerprintGenerator64 = Field(
        default=GetTopologicalTorsionGenerator(fpSize=1024),
        description="TopologicalTorsion fingerprint generator",
    )

    # --- Advanced Settings ---
    strict_structure_recovery: bool = Field(
        default=False, description="Whether to perform strict structure recovery"
    )
    max_structure_recovery_time: int = Field(
        default=10, description="Maximum structure recovery time (seconds)"
    )
    allow_spin_change: bool = Field(
        default=False, description="Whether to allow spin changes"
    )
    force_unit_transform: bool = Field(
        default=False, description="Whether to force unit conversion"
    )
    parallel_max_size: int = Field(
        default=8 * 1024**2,
        description="Maximum data size for parallel processing (bytes)",
    )
    max_recursion_depth: int = Field(
        default=3000, description="Maximum recursion depth for Python"
    )

    # --- New: Log File Control ---
    log_to_file: bool = Field(
        default=False, description="Whether to write log messages to a file"
    )

    def __init__(self, **data):
        """
        Initializes the configuration object and configures the logger based on current settings.
        """
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

    def set_log_level(
        self, level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    ):
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

    def set_n_jobs(self, n_jobs: int):
        return (
            min(n_jobs, multiprocessing.cpu_count())
            if n_jobs > 0
            else min(multiprocessing.cpu_count(), molopconfig.max_jobs)
        )


# --- Global Configuration Instance ---
# Create a globally available configuration instance
try:
    molopconfig = MolOPConfig()
    molopconfig.set_log_level("INFO")
except ValidationError as e:
    logging.error(f"Configuration validation failed: {e}")
    molopconfig = MolOPConfig()
    molopconfig.quiet()
    molopconfig.disable_file_logging()
    logging.critical(
        "Default safe configuration has been used. Please check and fix your custom configuration."
    )
