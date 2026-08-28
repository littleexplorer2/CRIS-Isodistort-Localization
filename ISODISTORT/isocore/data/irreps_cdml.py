"""CDML / Kovalev irrep tags from ISOTROPY ``data_irreps.txt`` (read-only).

Maps Miller-Love irrep labels (e.g. GM1+) to Kovalev tags (e.g. k14t1)
using the iso suite database under ``isobyu/``. This is crystallographic
metadata for the parent space group + k point, not a per-subgroup answer key.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

_ISO_DATA = Path(__file__).resolve().parents[2] / "isobyu" / "data_irreps.txt"


def _section_tokens(text: str, start: str, end: str) -> list[str]:
    lines = text.splitlines()
    s = e = None
    for i, line in enumerate(lines):
        if line.strip() == start:
            s = i + 1
        elif line.strip() == end and s is not None:
            e = i
            break
    if s is None or e is None:
        return []
    blob = "\n".join(lines[s:e])
    return [tok.strip() for tok in re.findall(r'"([^"]*)"', blob)]


@lru_cache(maxsize=1)
def _ml_kov_tables() -> tuple[tuple[str, ...], tuple[str, ...]]:
    if not _ISO_DATA.is_file():
        return (), ()
    text = _ISO_DATA.read_text(encoding="latin-1", errors="replace")
    ml = _section_tokens(text, "irrep_label", "irrep_label_bc")
    kov = _section_tokens(text, "irrep_label_kov", "irrep_label_kov_id")
    n = min(len(ml), len(kov))
    return tuple(ml[:n]), tuple(kov[:n])


def lookup_irrep_kovalev(
    irrep_label: str,
    *,
    k_kovalev: str | None = None,
    companion_irreps: list[str] | None = None,
) -> str | None:
    """Return Kovalev irrep tag (e.g. ``k14t1``) for a Miller-Love label.

    When several database rows share the same ML label, prefer the row whose
    Kovalev prefix matches ``k_kovalev`` (from the k-point table) and whose
    local neighborhood best covers ``companion_irreps`` (other IRs at that k).
    """
    label = (irrep_label or "").strip()
    if not label:
        return None
    ml, kov = _ml_kov_tables()
    if not ml:
        return None
    k_prefix = (k_kovalev or "").strip().lower()
    if k_prefix and not k_prefix.startswith("k"):
        k_prefix = f"k{k_prefix}"
    companions = {c.strip() for c in (companion_irreps or []) if c and c.strip()}
    companions.discard(label)

    candidates: list[int] = []
    for i, name in enumerate(ml):
        if name != label:
            continue
        tag = kov[i]
        if k_prefix and not tag.lower().startswith(k_prefix + "t"):
            continue
        candidates.append(i)
    if not candidates and k_prefix:
        # Fallback: any row with this ML label.
        candidates = [i for i, name in enumerate(ml) if name == label]
    if not candidates:
        return None
    if len(candidates) == 1 or not companions:
        return kov[candidates[0]]

    best_i = candidates[0]
    best_score = -1
    for i in candidates:
        lo, hi = max(0, i - 24), min(len(ml), i + 25)
        nearby = set(ml[lo:hi])
        score = len(companions & nearby)
        tag = kov[i]
        if k_prefix and tag.lower().startswith(k_prefix + "t"):
            score += 5
        if score > best_score:
            best_score = score
            best_i = i
    return kov[best_i]
