"""English UI strings (isocore.i18n)."""

from isocore.i18n import MESSAGES, t


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


def test_distortion_and_method2_help_keys():
    assert "m2.genDbHelp" in MESSAGES
    assert "Generate isotropy subgroups" in MESSAGES["m2.genDbHelp"]
    assert "\n" not in MESSAGES["m2.genDbHelp"]
    assert "dist.method4" in MESSAGES
    assert "dist.tableLabel" in MESSAGES
    assert "not limited by the filter" not in MESSAGES["dist.tableNote"]
    assert "dist.zipNotM4" in MESSAGES
    assert "ui.dist.table" in MESSAGES


def test_numbered_web_placeholders():
    text = t("ok.nsup", 3)
    assert "3" in text
