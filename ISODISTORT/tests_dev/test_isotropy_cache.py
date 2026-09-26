"""Tests for iso-generated isotropy-subgroup cache list/delete."""
from __future__ import annotations

import re

import pytest

from isocore.backend.iso_wrapper import IsoWrapper
from isocore.backend.isotropy_cache import (
    _safe_name,
    delete_isotropy_cache,
    list_isotropy_cache,
)


def test_safe_name_strips_path():
    assert _safe_name(r"a/b/i0123400.iso") == "i0123400.iso"
    assert _safe_name(r"C:\tmp\i0123400.iso") == "i0123400.iso"


def test_list_isotropy_cache_smoke():
    try:
        wrapper = IsoWrapper()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"WSL/iso unavailable: {exc}")
    entries = list_isotropy_cache(wrapper)
    assert isinstance(entries, list)
    for e in entries:
        assert re.match(r"^i\d+\.iso$", e.name, re.I)
        assert e.size >= 0
        d = e.to_dict()
        assert d["name"] == e.name


def test_delete_rejects_non_iso_names():
    try:
        wrapper = IsoWrapper()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"WSL/iso unavailable: {exc}")
    result = delete_isotropy_cache(wrapper, ["../etc/passwd", "data_isotropy.txt"])
    assert result["deleted"] == []
    assert len(result["skipped"]) == 2
