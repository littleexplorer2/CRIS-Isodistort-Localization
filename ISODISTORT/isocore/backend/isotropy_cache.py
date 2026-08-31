"""List / delete iso-generated isotropy-subgroup database files (``i*.iso``).

When Method 2 generates a missing parametric-k database, ``iso`` writes files
like ``i0665800.iso`` into the WSL staging directory (``~/.id/tmp``) after
``cd`` there. These are reusable caches, not the read-only ``isobyu/data_*.txt``
suite databases.
"""
from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .base_wrapper import BaseWrapper

_ISO_NAME_RE = re.compile(r"^i\d+\.iso$", re.I)
_HEADER_RE = re.compile(
    r"Space\s+Group\s+(\d+)\s*,\s*Irrep\s+(\S+)(?:\s*,\s*kparam\s+(\S+))?",
    re.I,
)


@dataclass(frozen=True)
class IsotropyCacheEntry:
    """One generated ``i*.iso`` cache file."""

    name: str
    size: int
    mtime_iso: str
    parent_sg: int | None = None
    irrep: str | None = None
    kparam: str | None = None
    header: str = ""

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "size": self.size,
            "mtime": self.mtime_iso,
            "parent_sg": self.parent_sg,
            "irrep": self.irrep,
            "kparam": self.kparam,
            "header": self.header,
        }


def _stage_tmp(wrapper: "BaseWrapper") -> str:
    if not getattr(wrapper, "_stage_dir", ""):
        # Ensure WSL staging exists (lazy init for callers that only manage cache).
        if hasattr(wrapper, "_init_wsl_environment") and wrapper._mode == "wsl":
            wrapper._init_wsl_environment()
    stage = getattr(wrapper, "_stage_dir", "") or ""
    if not stage:
        raise RuntimeError("WSL/native staging directory is not available")
    return stage


def list_isotropy_cache(wrapper: "BaseWrapper") -> list[IsotropyCacheEntry]:
    """Return sorted cache entries (newest first)."""
    stage = _stage_tmp(wrapper)
    # name|size|mtime_epoch|first_line
    cmd = (
        f"cd {shlex.quote(stage)} && "
        r"for f in i*.iso; do "
        r"[ -f \"$f\" ] || continue; "
        r"sz=$(wc -c < \"$f\" | tr -d ' '); "
        r"mt=$(stat -c %Y \"$f\" 2>/dev/null || stat -f %m \"$f\"); "
        r"hd=$(head -n 1 \"$f\" | tr '\t' ' '); "
        r"printf '%s\t%s\t%s\t%s\n' \"$f\" \"$sz\" \"$mt\" \"$hd\"; "
        r"done"
    )
    result = wrapper._wsl(cmd)
    if result.returncode != 0:
        return []
    entries: list[IsotropyCacheEntry] = []
    for line in result.stdout.splitlines():
        parts = line.split("\t", 3)
        if len(parts) < 3:
            continue
        name, size_s, mt_s = parts[0].strip(), parts[1].strip(), parts[2].strip()
        header = parts[3].strip() if len(parts) > 3 else ""
        if not _ISO_NAME_RE.match(name):
            continue
        try:
            size = int(size_s)
            mt = int(float(mt_s))
        except ValueError:
            continue
        m = _HEADER_RE.search(header)
        entries.append(
            IsotropyCacheEntry(
                name=name,
                size=size,
                mtime_iso=datetime.fromtimestamp(mt, tz=timezone.utc)
                .astimezone()
                .strftime("%Y-%m-%d %H:%M:%S"),
                parent_sg=int(m.group(1)) if m else None,
                irrep=m.group(2) if m else None,
                kparam=m.group(3) if m else None,
                header=header,
            )
        )
    entries.sort(key=lambda e: e.mtime_iso, reverse=True)
    return entries


def delete_isotropy_cache(
    wrapper: "BaseWrapper", names: list[str]
) -> dict[str, list[str]]:
    """Delete selected ``i*.iso`` files. Returns ``deleted`` / ``skipped`` lists."""
    stage = _stage_tmp(wrapper)
    deleted: list[str] = []
    skipped: list[str] = []
    for raw in names:
        name = _safe_name(raw)
        if not _ISO_NAME_RE.match(name):
            skipped.append(raw)
            continue
        dst = f"{stage}/{name}"
        result = wrapper._wsl(
            f"if [ -f {shlex.quote(dst)} ]; then rm -f {shlex.quote(dst)} && echo DEL; "
            f"else echo MISS; fi"
        )
        out = (result.stdout or "").strip()
        if result.returncode == 0 and out.startswith("DEL"):
            deleted.append(name)
        else:
            skipped.append(name)
    return {"deleted": deleted, "skipped": skipped}


def _safe_name(raw: str) -> str:
    """Basename-only sanitize for cache file names."""
    text = (raw or "").replace("\\", "/").strip()
    return text.rsplit("/", 1)[-1]
