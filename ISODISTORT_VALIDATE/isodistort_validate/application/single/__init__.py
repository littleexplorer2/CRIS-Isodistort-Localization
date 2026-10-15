"""Public single-CIF comparison workflow."""

from .service import (
    ComparisonResult,
    compare_cif,
)

__all__ = ["compare_cif", "ComparisonResult"]
