"""Comparison result contract shared by algorithms and workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal


@dataclass(frozen=True, slots=True)
class VerificationDecision:
    passed: bool
    status: Literal["pass", "fail", "inconclusive"]
    semantic_match: bool
    reference_verified: bool
    strict: bool
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class VerificationPolicy:
    """Strict acceptance requires semantics, trusted reference hash, and bytes.

    A normal semantic comparison does not certify reference provenance. Missing
    trust evidence is distinct from a measured mismatch; neither passes strict
    acceptance. This policy does not change crystal matching or its tolerances.
    """

    strict: bool = False

    def evaluate(self, result: ComparisonResult) -> VerificationDecision:
        issues = []
        if not result.structure_equal:
            issues.append("CIF scientific semantics differ")
        if result.reference_hash_valid is False:
            issues.append("reference SHA-256 does not match expected hash")
        if self.strict and not result.byte_exact:
            issues.append("strict verification requires byte-identical files")
        status: Literal["pass", "fail", "inconclusive"] = "fail" if issues else "pass"
        if self.strict and result.reference_hash_valid is None:
            issues.append("strict verification requires a trusted reference SHA-256")
            if status == "pass":
                status = "inconclusive"
        return VerificationDecision(
            passed=status == "pass", status=status, semantic_match=result.structure_equal,
            reference_verified=result.reference_hash_valid is True, strict=self.strict,
            issues=tuple(issues),
        )


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
        return VerificationPolicy().evaluate(self).passed
