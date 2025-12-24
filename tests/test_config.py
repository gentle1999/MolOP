import logging

import pytest

from molop.config import MolOPConfig, molopconfig


def test_config_defaults():
    config = MolOPConfig()
    assert config.show_progress_bar is True
    assert config.max_jobs == 16
    assert config.max_structure_recovery_time == 10.0
    assert config.log_to_file is False


def test_logging_control():
    config = MolOPConfig()

    # Test quiet/verbose
    config.quiet()
    assert config.show_progress_bar is False
    # Verifying handler removal is tricky without mocking, but we can check state

    config.verbose()
    assert config.show_progress_bar is True


def test_file_logging():
    config = MolOPConfig()
    config.enable_file_logging()
    assert config.log_to_file is True

    config.disable_file_logging()
    assert config.log_to_file is False


def test_log_level():
    config = MolOPConfig()
    config.set_log_level("DEBUG")
    assert logging.getLogger("molop").level == logging.DEBUG

    with pytest.raises(ValueError):
        config.set_log_level("INVALID")  # type: ignore


def test_recursion_depth():
    config = MolOPConfig()
    import sys

    original_limit = sys.getrecursionlimit()

    try:
        config.set_max_recursion_depth(5000)
        assert sys.getrecursionlimit() == 5000
    finally:
        sys.setrecursionlimit(original_limit)


def test_global_config_instance():
    assert isinstance(molopconfig, MolOPConfig)
