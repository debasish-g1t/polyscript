"""Tests for ``polyscript.utils.logger``."""

from __future__ import annotations

import io
import os
import tempfile

import pytest
from polyscript.utils.logger import (
    CRITICAL,
    DEBUG,
    ERROR,
    INFO,
    SILENT,
    VERBOSITY_TO_LEVEL,
    WARNING,
    PolyLogger,
    get_logger,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def buf() -> io.StringIO:
    return io.StringIO()


@pytest.fixture
def log() -> PolyLogger:
    return PolyLogger("test", level=DEBUG)


# ---------------------------------------------------------------------------
# Verbosity integer → level mapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "verbosity, expected",
    [(0, SILENT), (1, WARNING), (2, INFO), (3, DEBUG)],
)
def test_verbosity_mapping(verbosity, expected):
    log = PolyLogger("v", level=verbosity)
    assert log.level == expected


def test_verbosity_out_of_range_passthrough():
    log = PolyLogger("raw", level=99)
    assert log.level == 99
    assert log.level > CRITICAL  # behaves like SILENT


# ---------------------------------------------------------------------------
# Level gating
# ---------------------------------------------------------------------------


def test_silent_blocks_everything(buf):
    log = PolyLogger("s", SILENT)
    log.add_stream(buf, level=DEBUG)
    log.debug("d")
    log.info("i")
    log.warning("w")
    log.error("e")
    log.critical("c")
    assert buf.getvalue() == ""


def test_warning_blocks_debug_info(buf):
    log = PolyLogger("w", WARNING)
    log.add_stream(buf, level=DEBUG)
    log.debug("d")
    log.info("i")
    log.warning("w")
    log.error("e")
    out = buf.getvalue()
    assert "d" not in out
    assert "i" not in out
    assert "w" in out
    assert "e" in out


def test_info_blocks_debug(buf):
    log = PolyLogger("i", INFO)
    log.add_stream(buf, level=DEBUG)
    log.debug("d")
    log.info("i")
    log.warning("w")
    out = buf.getvalue()
    assert "d" not in out
    assert "i" in out
    assert "w" in out


def test_debug_passes_everything(buf):
    log = PolyLogger("d", DEBUG)
    log.add_stream(buf, level=DEBUG)
    log.debug("d")
    log.info("i")
    out = buf.getvalue()
    assert "d" in out
    assert "i" in out


# ---------------------------------------------------------------------------
# String level names
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, expected",
    [
        ("SILENT", SILENT),
        ("DEBUG", DEBUG),
        ("INFO", INFO),
        ("WARNING", WARNING),
        ("ERROR", ERROR),
        ("CRITICAL", CRITICAL),
    ],
)
def test_string_level_names(name, expected):
    log = PolyLogger(name, level=name)
    assert log.level == expected


def test_unknown_string_falls_back_to_debug():
    log = PolyLogger("u", level="NONEXISTENT")
    assert log.level == DEBUG


# ---------------------------------------------------------------------------
# Dynamic set_level / set_verbosity
# ---------------------------------------------------------------------------


def test_set_level_changes_threshold(buf):
    log = PolyLogger("dyn", level=INFO)
    log.add_stream(buf, level=DEBUG)
    log.debug("d")
    assert "d" not in buf.getvalue()

    log.set_level(DEBUG)
    buf2 = io.StringIO()
    log.add_stream(buf2, level=DEBUG)
    log.debug("d2")
    assert "d2" in buf2.getvalue()


def test_set_verbosity(buf):
    log = PolyLogger("v", level=DEBUG)
    log.add_stream(buf, level=DEBUG)

    log.set_verbosity(0)  # SILENT
    log.warning("w0")
    assert buf.getvalue() == ""

    log.set_verbosity(1)  # WARNING
    log.info("i1")
    log.warning("w1")
    assert "i1" not in buf.getvalue()
    assert "w1" in buf.getvalue()

    log.set_verbosity(3)  # DEBUG
    log.debug("d3")
    assert "d3" in buf.getvalue()


# ---------------------------------------------------------------------------
# Multiple handlers — different "versions" of output
# ---------------------------------------------------------------------------


def test_handlers_independent_levels():
    buf_verbose = io.StringIO()
    buf_concise = io.StringIO()

    log = PolyLogger("multi", level=DEBUG)
    log.add_stream(buf_verbose, level=DEBUG, fmt="%(levelname)s|%(message)s")
    log.add_stream(buf_concise, level=WARNING, fmt="%(levelname)s|%(message)s")

    log.debug("detail")
    log.warning("alert")

    assert "detail" in buf_verbose.getvalue()
    assert "alert" in buf_verbose.getvalue()

    assert "detail" not in buf_concise.getvalue()
    assert "alert" in buf_concise.getvalue()


def test_handlers_independent_formats():
    buf_a = io.StringIO()
    buf_b = io.StringIO()

    log = PolyLogger("fmt", level=DEBUG)
    log.add_stream(buf_a, fmt="A:%(message)s")
    log.add_stream(buf_b, fmt="B:%(message)s")

    log.info("hello")
    assert "A:hello" in buf_a.getvalue()
    assert "B:hello" in buf_b.getvalue()


# ---------------------------------------------------------------------------
# File handler
# ---------------------------------------------------------------------------


def test_file_handler():
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log") as tf:
        path = tf.name

    try:
        log = PolyLogger("file", level=DEBUG)
        log.add_file(path, level=DEBUG, fmt="%(message)s")
        log.info("file-logged")
        log.close()

        with open(path) as fh:
            content = fh.read()
        assert "file-logged" in content
    finally:
        os.unlink(path)


# ---------------------------------------------------------------------------
# get_logger convenience
# ---------------------------------------------------------------------------


def test_get_logger_defaults():
    log = get_logger("quick", level=WARNING)
    assert log.level == WARNING
    assert len(log._logger.handlers) == 1  # console only


def test_get_logger_with_file():
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log") as tf:
        path = tf.name

    try:
        log = get_logger("qf", level=DEBUG, log_file=path)
        assert len(log._logger.handlers) == 2  # console + file
        log.info("quick-file")
        log.close()

        with open(path) as fh:
            assert "quick-file" in fh.read()
    finally:
        os.unlink(path)


def test_get_logger_no_console():
    log = get_logger("nc", console=False)
    assert len(log._logger.handlers) == 0


# ---------------------------------------------------------------------------
# Level-check properties
# ---------------------------------------------------------------------------


def test_is_properties():
    log = PolyLogger("props", level=INFO)
    assert not log.is_debug
    assert log.is_info
    assert log.is_warning
    assert log.is_error

    log.set_level(DEBUG)
    assert log.is_debug


# ---------------------------------------------------------------------------
# __repr__
# ---------------------------------------------------------------------------


def test_repr(log):
    r = repr(log)
    assert "PolyLogger" in r
    assert "test" in r


# ---------------------------------------------------------------------------
# add_handler with custom handler
# ---------------------------------------------------------------------------


def test_add_handler_custom():
    import logging

    buf = io.StringIO()
    handler = logging.StreamHandler(buf)
    handler.setFormatter(logging.Formatter("%(message)s"))

    log = PolyLogger("custom")
    log.add_handler(handler, level=INFO)
    log.info("custom-msg")
    assert "custom-msg" in buf.getvalue()


# ---------------------------------------------------------------------------
# std_logger compatibility
# ---------------------------------------------------------------------------


def test_std_logger():
    log = PolyLogger("compat")
    assert isinstance(log.std_logger, __import__("logging").Logger)


# ---------------------------------------------------------------------------
# Chainable handler methods
# ---------------------------------------------------------------------------


def test_handler_chaining():
    buf = io.StringIO()
    log = (
        PolyLogger("chain", level=DEBUG)
        .add_stream(buf, level=DEBUG, fmt="%(message)s")
        .add_console(level=SILENT)  # shouldn't interfere
    )
    log.info("chained")
    assert "chained" in buf.getvalue()


# ---------------------------------------------------------------------------
# close  /  warn alias
# ---------------------------------------------------------------------------


def test_close_removes_handlers():
    log = PolyLogger("closing")
    log.add_console()
    assert len(log._logger.handlers) == 1
    log.close()
    assert len(log._logger.handlers) == 0


def test_warn_alias(buf):
    log = PolyLogger("warn", WARNING)
    log.add_stream(buf, level=DEBUG, fmt="%(message)s")
    log.warn("via-warn")
    assert "via-warn" in buf.getvalue()


def test_flush_does_not_remove_handlers():
    log = PolyLogger("flush-test")
    log.add_console()
    count_before = len(log._logger.handlers)
    log.flush()
    assert len(log._logger.handlers) == count_before


# ===========================================================================
# BaseExecutor
# ===========================================================================

from polyscript.utils.executor import BaseExecutor


class _FakeRunner(BaseExecutor):
    """Minimal subclass for testing."""

    _log_level = DEBUG

    def run(self, msg: str) -> None:
        self.logger.info(msg)


class _FakeSilentRunner(BaseExecutor):
    _log_level = SILENT

    def run(self, msg: str) -> None:
        self.logger.info(msg)


class TestBaseExecutor:
    def test_auto_logger_created(self):
        runner = _FakeRunner()
        assert isinstance(runner.logger, PolyLogger)
        assert runner.logger.level == DEBUG

    def test_explicit_logger_passed(self):
        custom = PolyLogger("custom", level=SILENT)
        runner = _FakeRunner(logger=custom)
        assert runner.logger is custom
        assert runner.logger.level == SILENT

    def test_level_propagates_to_instance(self, buf):
        runner = _FakeRunner()
        runner.logger.add_stream(buf, level=DEBUG, fmt="%(message)s")
        runner.run("hello")
        assert "hello" in buf.getvalue()

    def test_silent_subclass(self, buf):
        runner = _FakeSilentRunner()
        runner.logger.add_stream(buf, level=DEBUG)
        runner.run("should not appear")
        assert buf.getvalue() == ""

    def test_set_log_level_classmethod(self):
        original = _FakeRunner._log_level
        try:
            _FakeRunner.set_log_level(WARNING)
            runner = _FakeRunner()
            assert runner.logger.level == WARNING
        finally:
            _FakeRunner._log_level = original

    def test_set_log_level_does_not_affect_other_classes(self):
        orig_a = _FakeRunner._log_level
        orig_b = _FakeSilentRunner._log_level
        try:
            _FakeRunner.set_log_level(WARNING)
            runner_b = _FakeSilentRunner()
            assert runner_b.logger.level == SILENT  # unchanged
        finally:
            _FakeRunner._log_level = orig_a
            _FakeSilentRunner._log_level = orig_b
