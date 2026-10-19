"""Tests for backend.messages; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.messages import MESSAGES
from backend.messages import t


# Original coverage: test_basics.py.

# --- English message formatting ---

def test_messages_are_flat_english():
    assert isinstance(MESSAGES, dict)
    assert "zh" not in MESSAGES
    assert "en" not in MESSAGES
    assert t("load.done", sg=139, sym="I4/mmm", n=10).startswith("[Loaded]")
    assert "Space group" in t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "139" in t("load.done", sg=139, sym="I4/mmm", n=10)

def test_unknown_key_returns_key():
    assert t("no.such.key") == "no.such.key"

def test_no_language_switch_keys():
    assert "ui.menu.language" not in MESSAGES
    assert "ui.lang.current" not in MESSAGES
    assert "dist.gen" not in MESSAGES
    assert "dist.domainsBtn" not in MESSAGES

def test_no_terminal_ui_message_keys():
    for key in (
        "ui.banner.title",
        "ui.menu.search",
        "prefs.terminalBlock",
        "m2.nmodRemoved",
        "dist.computeModesAsk",
    ):
        assert key not in MESSAGES

def test_distortion_and_method2_help_keys():
    assert "m2.genDbHelp" in MESSAGES
    assert "isotropy-subgroup" in MESSAGES["m2.genDbHelp"].lower() or (
        "isotropy" in MESSAGES["m2.genDbHelp"].lower()
        and "subgroup" in MESSAGES["m2.genDbHelp"].lower()
    )
    assert "cache" in MESSAGES["m2.genDbHelp"].lower()
    assert "\n" not in MESSAGES["m2.genDbHelp"]
    assert "m2.genDbManage" in MESSAGES
    assert "m2.genDbDelete" in MESSAGES
    assert "dist.method4" in MESSAGES
    assert "dist.tableLabel" in MESSAGES
    assert "not limited by the filter" not in MESSAGES["dist.tableNote"]
    assert "dist.zipNotM4" in MESSAGES
    assert "dist.zipWait" in MESSAGES
    assert "dist.zipDone" in MESSAGES
    assert "st.elapsed" in MESSAGES
    assert "st.busyHint" in MESSAGES
    assert "independent incommensurate" in MESSAGES["l.nmod"].lower()
    assert "ok.nmod" in MESSAGES
    assert "ss.title" not in MESSAGES
    assert "ui.menu.superspace" not in MESSAGES

def test_numbered_web_placeholders():
    text = t("ok.nsup", 3)
    assert "3" in text
