"""
测试：国际化（isocore.i18n）—— zh / en 两种模式
"""
import pytest

from isocore.i18n import (
    LANGUAGES,
    MESSAGES,
    TERMS_EN2ZH,
    get_language,
    set_language,
    t,
    term_en2zh,
    term_zh2en,
    translate_term,
)


def test_default_language_is_en():
    assert get_language() == "en"


def test_supported_languages():
    assert LANGUAGES == ("zh", "en")


def test_set_and_get_language():
    set_language("en")
    assert get_language() == "en"
    set_language("zh")
    assert get_language() == "zh"


def test_invalid_language_raises():
    with pytest.raises(ValueError):
        set_language("fr")
    with pytest.raises(ValueError):
        set_language("mixed")  # mixed 模式已移除


def test_message_translation_and_format():
    set_language("zh")
    zh = t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "空间群" in zh and "139" in zh

    set_language("en")
    en = t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "Space group" in en and "139" in en

    set_language("en")


def test_zh_messages_are_chinese_only():
    """中文文案应为纯中文（修复中文模式下仍显示英文的问题）。"""
    set_language("zh")
    # 网页端文案不应出现“中文 / English”双语文案
    for key in ("hStatus", "hCif", "hTypes", "btn.load", "ok.loaded", "st.wait"):
        text = t(key)
        assert "/" not in text, f"zh[{key}] 含双语文案: {text}"
        assert "English" not in text, f"zh[{key}] 含英文: {text}"
    # 终端主菜单标题
    assert t("ui.search_page") == "搜索页"
    assert t("ui.menu.exit") == "0. 退出"
    assert t("ui.msg.state") == "当前会话状态"
    set_language("en")


def test_unknown_key_returns_key():
    assert t("no.such.key") == "no.such.key"


def test_terms_dictionary():
    # 中文术语以《晶体学名词》为准
    assert TERMS_EN2ZH["space group"] == "空间群"
    assert TERMS_EN2ZH["Wyckoff position"] == "Wyckoff 位置"
    assert TERMS_EN2ZH["irreducible representation (IR)"] == "不可约表示 (IR)"
    assert term_en2zh("crystal system") == "晶系"
    assert term_zh2en("序参量") == "order parameter"


def test_translate_term_by_language():
    assert translate_term("space group", "zh") == "空间群"
    assert translate_term("space group", "en") == "space group"


def test_message_catalogs_parallel():
    """中英文案目录键集合必须一致。"""
    assert set(MESSAGES["zh"].keys()) == set(MESSAGES["en"].keys())
