"""Public read-only CIF parsing and fingerprint API."""

from .service import (
    parse_structure,
    read_block,
    read_text,
    sha256,
)

__all__ = ["sha256", "read_text", "parse_structure", "read_block"]
