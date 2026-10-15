"""Public API for backend.validation."""

from .checks import (
    section_scalar,
    validate_patched_text,
)

__all__ = ['section_scalar', 'validate_patched_text']
