"""Serialize a restricted set of official ISO commands; no scientific solver."""

import re

from .runtime import InputError


def build_iso_query(parent, irrep, opd, wyckoff=None):
    """Prepare explicit user choices for the official ISOTROPY command parser.

    No k parameters or atomic coordinates are inferred. The caller must verify
    that the selected IR/OPD/Wyckoff labels apply to the selected parent.
    """
    if isinstance(parent, bool) or not isinstance(parent, int) or not 1 <= parent <= 230:
        raise InputError("Parent space-group number must be an integer from 1 to 230.")
    if not isinstance(irrep, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9+\-]{0,31}", irrep):
        raise InputError("IR label must use letters, digits, + or - with no spaces.")
    if not isinstance(opd, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9]{0,15}", opd):
        raise InputError("OPD label must use letters and digits with no spaces.")
    if wyckoff is not None and (
        not isinstance(wyckoff, str) or not re.fullmatch(r"[a-z]", wyckoff)
    ):
        raise InputError("Optional Wyckoff label must be one lowercase letter.")
    lines = [
        "PAGE NOBREAK", f"VALUE PARENT {parent}", f"VALUE IRREP {irrep}",
        f"VALUE DIRECTION {opd}", "SHOW PARENT", "SHOW IRREP", "SHOW DIRECTION",
        "SHOW SUBGROUP", "SHOW BASIS", "SHOW ORIGIN", "SHOW SIZE", "SHOW INDEX",
        "DISPLAY ISOTROPY",
    ]
    if wyckoff is not None:
        lines.extend([
            f"VALUE WYCKOFF {wyckoff}", "SHOW MICROSCOPIC VECTOR", "DISPLAY DISTORTION",
        ])
    lines.append("QUIT")
    return ("\n".join(lines) + "\n").encode("ascii")
