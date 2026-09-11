import importlib.metadata
from types import SimpleNamespace

import pytest
from pydantic import ValidationError

import molop.config as config_module


def _count_handler(logger, handler):
    return sum(1 for item in logger.handlers if item is handler)


def _remove_all(logger, handler):
    while handler in logger.handlers:
        logger.removeHandler(handler)


def _add_once(logger, handler):
    if handler not in logger.handlers:
        logger.addHandler(handler)


def test_config_init_false_branches(monkeypatch):
    logger = config_module.moloplogger
    stream_before = config_module.stream_handler in logger.handlers
    file_before = config_module.file_handler in logger.handlers
    dof_calls = []
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, enable: dof_calls.append(enable),
    )

    config = config_module.MolOPConfig(
        show_progress_bar=False,
        log_to_file=True,
        use_dof_effect_drawer=False,
    )

    assert config.show_progress_bar is False
    assert config.log_to_file is True
    assert config_module.stream_handler not in logger.handlers
    assert config_module.file_handler in logger.handlers
    assert dof_calls == [False]

    _remove_all(logger, config_module.stream_handler)
    _remove_all(logger, config_module.file_handler)
    if stream_before:
        _add_once(logger, config_module.stream_handler)
    if file_before:
        _add_once(logger, config_module.file_handler)


def test_native_logging_policy_defaults_to_suppressed(monkeypatch):
    calls = []

    class FakeRDLogger:
        @staticmethod
        def DisableLog(pattern):
            calls.append(("rdkit", "disable", pattern))

        @staticmethod
        def EnableLog(pattern):
            calls.append(("rdkit", "enable", pattern))

    class FakeErrorLog:
        def StopLogging(self):
            calls.append(("openbabel", "stop"))

        def StartLogging(self):
            calls.append(("openbabel", "start"))

    monkeypatch.setattr(config_module, "RDLogger", FakeRDLogger)
    monkeypatch.setattr(
        config_module,
        "pybel",
        SimpleNamespace(ob=SimpleNamespace(obErrorLog=FakeErrorLog())),
    )

    config = config_module.MolOPConfig()

    assert config.suppress_rdkit_logs is True
    assert config.suppress_openbabel_logs is True
    assert calls == [("rdkit", "disable", "rdApp.*"), ("openbabel", "stop")]


def test_native_logging_policy_can_be_reenabled(monkeypatch):
    calls = []

    class FakeRDLogger:
        @staticmethod
        def DisableLog(pattern):
            calls.append(("rdkit", "disable", pattern))

        @staticmethod
        def EnableLog(pattern):
            calls.append(("rdkit", "enable", pattern))

    class FakeErrorLog:
        def StopLogging(self):
            calls.append(("openbabel", "stop"))

        def StartLogging(self):
            calls.append(("openbabel", "start"))

    monkeypatch.setattr(config_module, "RDLogger", FakeRDLogger)
    monkeypatch.setattr(
        config_module,
        "pybel",
        SimpleNamespace(ob=SimpleNamespace(obErrorLog=FakeErrorLog())),
    )

    config_module.MolOPConfig(
        suppress_rdkit_logs=False,
        suppress_openbabel_logs=False,
    )

    assert calls == [("rdkit", "enable", "rdApp.*"), ("openbabel", "start")]


def test_set_n_jobs_boundaries(monkeypatch):
    logger = config_module.moloplogger
    stream_before = config_module.stream_handler in logger.handlers
    file_before = config_module.file_handler in logger.handlers
    monkeypatch.setattr(config_module, "available_cpu_count", lambda: 8)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig(max_jobs=6)

    assert config.set_n_jobs(3) == 3
    assert config.set_n_jobs(99) == 6
    assert config.set_n_jobs(0) == 6
    assert config.set_n_jobs(-4) == 6

    _remove_all(logger, config_module.stream_handler)
    _remove_all(logger, config_module.file_handler)
    if stream_before:
        _add_once(logger, config_module.stream_handler)
    if file_before:
        _add_once(logger, config_module.file_handler)


def test_auto_max_jobs_uses_process_available_cpu_count(monkeypatch):
    monkeypatch.setattr(config_module, "available_cpu_count", lambda: 24)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig(max_jobs=None)

    assert config.effective_max_jobs == 24
    assert config.set_n_jobs(-1) == 24
    assert config.set_n_jobs(100) == 24


def test_reconstruction_failure_policy_is_exposed_and_applied(monkeypatch):
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )
    monkeypatch.setattr(
        config_module.MOLGR_CONFIG.interface,
        "reconstruction_failure_policy",
        "raise",
    )
    config = config_module.MolOPConfig(reconstruction_failure_policy="return_suspicious")

    assert config.reconstruction_failure_policy == "return_suspicious"
    assert config.apply_molgr_reconstruction_policy() == "return_suspicious"
    assert config_module.MOLGR_CONFIG.interface.reconstruction_failure_policy == (
        "return_suspicious"
    )


@pytest.mark.parametrize(
    ("available", "expected"),
    [(1, 1), (2, 1), (3, 2), (4, 2), (8, 5), (24, 16)],
)
def test_molgr_parallelism_uses_two_thirds_of_available_cpu(monkeypatch, available, expected):
    monkeypatch.setattr(config_module, "available_cpu_count", lambda: available)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig(max_jobs=None)

    assert config.effective_molgr_max_jobs == expected
    assert config.set_molgr_n_jobs(-1) == expected
    assert config.set_molgr_n_jobs(99) == expected


def test_molgr_parallelism_also_respects_max_jobs(monkeypatch):
    monkeypatch.setattr(config_module, "available_cpu_count", lambda: 12)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig(max_jobs=4)

    assert config.effective_molgr_max_jobs == 4
    assert config.set_molgr_n_jobs(99) == 4


def test_available_cpu_count_uses_joblib_process_limit(monkeypatch):
    monkeypatch.setattr(config_module, "joblib_cpu_count", lambda: 8)

    assert config_module.available_cpu_count() == 8


def test_available_cpu_count_never_falls_below_one(monkeypatch):
    monkeypatch.setattr(config_module, "joblib_cpu_count", lambda: 0)

    assert config_module.available_cpu_count() == 1


def test_available_cpu_count_matches_joblib_runtime_limit():
    assert config_module.available_cpu_count() == config_module.joblib_cpu_count()


def test_runtime_metadata_declares_psutil_for_affinity_fallback():
    requirements = importlib.metadata.requires("molop") or []

    assert any(requirement.lower().startswith("psutil") for requirement in requirements)


def test_max_jobs_environment_default_and_explicit_override(monkeypatch):
    monkeypatch.setenv(config_module.MAX_JOBS_ENV_VAR, "4")
    monkeypatch.setattr(config_module, "available_cpu_count", lambda: 12)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    environment_config = config_module.MolOPConfig()
    automatic_config = config_module.MolOPConfig(max_jobs=None)
    explicit_config = config_module.MolOPConfig(max_jobs=7)

    assert environment_config.max_jobs == 4
    assert environment_config.effective_max_jobs == 4
    assert automatic_config.max_jobs is None
    assert automatic_config.effective_max_jobs == 12
    assert explicit_config.effective_max_jobs == 7


def test_max_jobs_rejects_non_positive_values(monkeypatch):
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    with pytest.raises(ValidationError):
        config_module.MolOPConfig(max_jobs=0)

    config = config_module.MolOPConfig(max_jobs=2)
    with pytest.raises(ValidationError):
        config.max_jobs = 0


def test_quiet_verbose_handler_toggling_without_duplicates(monkeypatch):
    logger = config_module.moloplogger
    stream_before = config_module.stream_handler in logger.handlers
    file_before = config_module.file_handler in logger.handlers
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig()
    _remove_all(logger, config_module.stream_handler)

    config.verbose()
    config.verbose()
    assert _count_handler(logger, config_module.stream_handler) == 1

    config.quiet()
    config.quiet()
    assert _count_handler(logger, config_module.stream_handler) == 0

    _remove_all(logger, config_module.stream_handler)
    _remove_all(logger, config_module.file_handler)
    if stream_before:
        _add_once(logger, config_module.stream_handler)
    if file_before:
        _add_once(logger, config_module.file_handler)


def test_file_logging_toggling_without_duplicates(monkeypatch):
    logger = config_module.moloplogger
    stream_before = config_module.stream_handler in logger.handlers
    file_before = config_module.file_handler in logger.handlers
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig()
    _remove_all(logger, config_module.file_handler)

    config.enable_file_logging()
    config.enable_file_logging()
    assert _count_handler(logger, config_module.file_handler) == 1

    config.disable_file_logging()
    config.disable_file_logging()
    assert _count_handler(logger, config_module.file_handler) == 0

    _remove_all(logger, config_module.stream_handler)
    _remove_all(logger, config_module.file_handler)
    if stream_before:
        _add_once(logger, config_module.stream_handler)
    if file_before:
        _add_once(logger, config_module.file_handler)
