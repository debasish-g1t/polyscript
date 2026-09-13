"""Tests for BaseExecutor."""

import pytest
from polyscript.utils.executor import BaseExecutor
from polyscript.utils.logger import INFO, WARNING, PolyLogger


# ---------------------------------------------------------------------------
# Test subclasses
# ---------------------------------------------------------------------------


class _DefaultExecutor(BaseExecutor):
    """Subclass with default _log_level (INFO)."""


class _WarnExecutor(BaseExecutor):
    """Subclass that overrides _log_level to WARNING."""
    _log_level = WARNING


class _MultiInherit(BaseExecutor):
    """Tests that self.logger is not clobbered when a parent sets it."""

    def __init__(self, **kwargs):
        self.logger = PolyLogger(name="pre_set", level=WARNING)
        super().__init__(**kwargs)


class _KwargsConsumer:
    """Cooperative parent that consumes extra kwargs."""

    def __init__(self, *, extra=None, **kwargs):
        self._caught = {"extra": extra}
        super().__init__(**kwargs)


class _KwargsPassthrough(_KwargsConsumer, BaseExecutor):
    """Tests that **kwargs flow through the cooperative MRO."""


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestBaseExecutor:
    """Tests for the BaseExecutor mixin."""

    # ── auto-created logger ────────────────────────────────────────────

    def test_auto_creates_logger_with_class_name(self):
        """When no logger is passed, one is created with the class name."""
        ex = _DefaultExecutor()
        assert ex.logger._logger.name == "_defaultexecutor"
        assert ex.logger.level == INFO
        # Console handler should be attached
        assert len(ex.logger._logger.handlers) == 1

    def test_respects_subclass_log_level(self):
        """Subclass _log_level override is honoured."""
        ex = _WarnExecutor()
        assert ex.logger.level == WARNING

    def test_explicit_logger_bypasses_default(self):
        """A passed-in logger is used as-is."""
        custom = PolyLogger(name="custom", level=WARNING)
        ex = _DefaultExecutor(logger=custom)
        assert ex.logger is custom           # same object, not a copy
        assert ex.logger.level == WARNING    # subclass _log_level ignored
        assert ex.logger._logger.name == "custom"

    # ── hasattr guard ──────────────────────────────────────────────────

    def test_does_not_clobber_preset_logger(self):
        """If self.logger is already set, BaseExecutor leaves it alone."""
        ex = _MultiInherit()
        assert ex.logger._logger.name == "pre_set"
        assert ex.logger.level == WARNING

    # ── kwargs passthrough ─────────────────────────────────────────────

    def test_kwargs_forwarded_to_super_init(self):
        """Extra **kwargs reach super().__init__()."""
        ex = _KwargsPassthrough(extra="value")
        assert ex._caught == {"extra": "value"}

    # ── class-level level control ──────────────────────────────────────

    def test_set_log_level_classmethod(self):
        """set_log_level changes _log_level for future instances."""

        class _Tmp(BaseExecutor):
            _log_level = INFO

        assert _Tmp._log_level == INFO
        _Tmp.set_log_level(WARNING)
        assert _Tmp._log_level == WARNING

        ex = _Tmp()
        assert ex.logger.level == WARNING
