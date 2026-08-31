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

# Parametric / line irreps (LD, SM, …) are often absent from data_irreps.txt.
# When lookup fails, synthesize ``{k}t{τ}`` using CDML τ remaps where known.
# Keyed by (parent_sg, k_label) → {Miller-Love index → Kovalev τ index}.
_LINE_IRREP_TAU: dict[tuple[int, str], dict[int, int]] = {
    # I4/mmm (#139) Λ / LD: LD1→τ1, LD2→τ3, LD3→τ2, LD4→τ4, LD5→τ5
    (139, "LD"): {1: 1, 2: 3, 3: 2, 4: 4, 5: 5},
}


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


def _section_ints(text: str, start: str, end: str) -> list[int]:
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
    return [int(x) for x in re.findall(r"-?\d+", blob)]


@lru_cache(maxsize=1)
def _ml_kov_tables() -> tuple[tuple[str, ...], tuple[str, ...], tuple[int, ...]]:
    if not _ISO_DATA.is_file():
        return (), (), ()
    text = _ISO_DATA.read_text(encoding="latin-1", errors="replace")
    ml = _section_tokens(text, "irrep_label", "irrep_label_bc")
    kov = _section_tokens(text, "irrep_label_kov", "irrep_label_kov_id")
    sgs = _section_ints(text, "irrep_space_group", "irrep_space_group_pointer")
    n = min(len(ml), len(kov), len(sgs) or len(ml))
    if not sgs:
        sgs = [0] * n
    return tuple(ml[:n]), tuple(kov[:n]), tuple(sgs[:n])


def _normalize_k_prefix(k_kovalev: str | None) -> str:
    k_prefix = (k_kovalev or "").strip().lower()
    if k_prefix and not k_prefix.startswith("k"):
        k_prefix = f"k{k_prefix}"
    return k_prefix


def _synthesize_line_tag(
    label: str,
    *,
    k_kovalev: str | None,
    parent_sg: int | None,
    k_point_label: str | None,
) -> str | None:
    """Build ``k10t1``-style tags for line irreps missing from data_irreps."""
    k_prefix = _normalize_k_prefix(k_kovalev)
    if not k_prefix:
        return None
    m = re.fullmatch(r"([A-Za-z]+)(\d+)([+-]?)", label)
    if not m:
        return None
    letter, num_s, parity = m.group(1), m.group(2), m.group(3)
    k_lab = (k_point_label or "").strip()
    if k_lab and letter.upper() != k_lab.upper():
        return None
    n = int(num_s)
    remap = _LINE_IRREP_TAU.get((int(parent_sg or 0), letter.upper()), {})
    tau = remap.get(n, n)
    return f"{k_prefix}t{tau}{parity}"


def lookup_irrep_kovalev(
    irrep_label: str,
    *,
    k_kovalev: str | None = None,
    parent_sg: int | None = None,
    k_point_label: str | None = None,
    companion_irreps: list[str] | None = None,
) -> str | None:
    """Return Kovalev irrep tag (e.g. ``k14t1``) for a Miller-Love label.

    When several database rows share the same ML label, prefer the row for
    ``parent_sg`` whose Kovalev prefix matches ``k_kovalev`` (from the k-point
    table). ``companion_irreps`` is kept for API compatibility but ignored —
    global companion windows previously biased the wrong SG block.
    """
    del companion_irreps  # deprecated; parent_sg disambiguates correctly
    label = (irrep_label or "").strip()
    if not label:
        return None
    ml, kov, sgs = _ml_kov_tables()
    if not ml:
        return _synthesize_line_tag(
            label,
            k_kovalev=k_kovalev,
            parent_sg=parent_sg,
            k_point_label=k_point_label,
        )
    k_prefix = _normalize_k_prefix(k_kovalev)
    parent = int(parent_sg) if parent_sg else 0

    def _collect(*, require_sg: bool) -> list[int]:
        out: list[int] = []
        for i, name in enumerate(ml):
            if name != label:
                continue
            if require_sg and parent and sgs[i] != parent:
                continue
            tag = kov[i]
            if k_prefix and not tag.lower().startswith(k_prefix + "t"):
                continue
            out.append(i)
        return out

    candidates = _collect(require_sg=True) if parent else []
    if not candidates:
        candidates = _collect(require_sg=False)
    if not candidates and not k_prefix:
        # No k-point hint: any row with this ML label (optionally SG-scoped).
        candidates = [
            i
            for i, name in enumerate(ml)
            if name == label and (not parent or sgs[i] == parent)
        ]
        if not candidates and parent:
            candidates = [i for i, name in enumerate(ml) if name == label]
    if candidates:
        return kov[candidates[0]]
    # Prefer synthesize over a wrong-k DB row (parametric Y/LD often absent).
    return _synthesize_line_tag(
        label,
        k_kovalev=k_kovalev,
        parent_sg=parent_sg,
        k_point_label=k_point_label,
    )
