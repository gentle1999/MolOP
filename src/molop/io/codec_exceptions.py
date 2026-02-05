"""Lightweight codec exception definitions."""


class UnsupportedFormatError(ValueError):
    """Raised when no codec supports the requested format."""


class ParseError(ValueError):
    """Raised when a reader codec cannot parse input data."""


class ConversionError(ValueError):
    """Raised when converting between structure levels fails."""


class MissingOptionalDependencyError(ImportError):
    """Raised when a codec requires an unavailable optional dependency."""
