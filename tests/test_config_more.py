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


def test_set_n_jobs_boundaries(monkeypatch):
    logger = config_module.moloplogger
    stream_before = config_module.stream_handler in logger.handlers
    file_before = config_module.file_handler in logger.handlers
    monkeypatch.setattr(config_module.multiprocessing, "cpu_count", lambda: 8)
    monkeypatch.setattr(
        type(config_module.dofconfig),
        "enable_ipython_integration",
        lambda _self, _enable: None,
    )

    config = config_module.MolOPConfig(max_jobs=6)

    assert config.set_n_jobs(3) == 3
    assert config.set_n_jobs(99) == 8
    assert config.set_n_jobs(0) == 6
    assert config.set_n_jobs(-4) == 6

    _remove_all(logger, config_module.stream_handler)
    _remove_all(logger, config_module.file_handler)
    if stream_before:
        _add_once(logger, config_module.stream_handler)
    if file_before:
        _add_once(logger, config_module.file_handler)


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
