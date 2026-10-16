"""Comparison result contract shared by algorithms and workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class ComparisonResult:
    local_path: str
    reference_path: str
    byte_exact: bool
    text_exact: bool
    structure_equal: bool
    reference_hash_valid: bool | None
    local_sha256: str
    reference_sha256: str
    details: dict[str, Any]

    @property
    def passed(self) -> bool:
        """Semantic equality is the default; formatting is reported separately."""
        return self.structure_equal and self.reference_hash_valid is not False
