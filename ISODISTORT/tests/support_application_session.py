"""Tests for support_application_session; migrated without changing scientific expectations."""
from __future__ import annotations

from tests.support_paths import TESTS_ROOT
import pytest
import re
import shutil
import subprocess


# Original coverage: test_robustness.py.

_robustness_CIFS_DIR = TESTS_ROOT / "cifs_30"

# Original coverage: test_wsl.py.

def _wsl_wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(  # noqa: PLW1510
            ["wsl.exe", "--status"],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False

_wsl_required = pytest.mark.skipif(
    not _wsl_wsl_available(), reason="WSL 不可用，跳过真实二进制 / 金标准测试"
)

def _mode_core_token(label: str) -> str:
    """Identity of a complete-mode label ignoring parent-site display names.

    Official IsoVIZ compact labels drop the k prefix for Gamma *and* for
    I4/mmm M=(1,1,1) (written ``M1+[Al:d:dsp]…``).  Local pretty labels keep
    ``[1,1,1]``.  Both must hash to the same core token.
    """
    text = str(label)

    def _looks_like_k(token: str) -> bool:
        return "," in token or "/" in token

    m = re.search(
        r"\[([^]]+)\]([A-Za-z0-9+-]+)(?:\([^)]*\))?"
        r"\[[^:]+:([a-z]):dsp\]([A-Za-z0-9_]+(?:\([^)]*\))?)",
        text,
    )
    if m and _looks_like_k(m.group(1)):
        return f"{m.group(1)}|{m.group(2)}|{m.group(3)}|{m.group(4)}"
    m = re.search(
        r"([A-Za-z0-9+-]+)\[[^:]+:([a-z]):dsp\]([A-Za-z0-9_]+(?:\([^)]*\))?)",
        text,
    )
    if not m:
        return text.strip()
    irrep, letter, site = m.group(1), m.group(2), m.group(3)
    if irrep.startswith("GM"):
        k = "0,0,0"
    elif irrep.startswith("M"):
        k = "1,1,1"
    else:
        k = "0,0,0"
    return f"{k}|{irrep}|{letter}|{site}"
