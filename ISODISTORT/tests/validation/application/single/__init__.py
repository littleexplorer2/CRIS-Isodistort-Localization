"""Public single-CIF comparison workflow."""

from .service import (
    ComparisonResult,
    compare_cif,
)
from ...backend.results import VerificationDecision, VerificationPolicy

__all__ = ["compare_cif", "ComparisonResult", "VerificationDecision", "VerificationPolicy"]
