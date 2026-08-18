"""
测试：国际化（isocore.i18n）—— zh / en / mixed 三种模式
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
    assert LANGUAGES == ("zh", "en", "mixed")


def test_set_and_get_language():
    set_language("en")
    assert get_language() == "en"
    set_language("zh")
    assert get_language() == "zh"
    set_language("mixed")
    assert get_language() == "mixed"


def test_invalid_language_raises():
    with pytest.raises(ValueError):
        set_language("fr")


def test_message_translation_and_format():
    set_language("zh")
    zh = t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "空间群" in zh and "139" in zh

    set_language("en")
    en = t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "Space group" in en and "139" in en

    set_language("en")


def test_mixed_mode_keeps_chinese_connectors_but_english_terms():
    """mixed 模式：核心专有名词英文，其余中文衔接。"""
    set_language("mixed")
    msg = t("subgroups.found", n=5)
    # 专有名词转英文（术语表为单数形式）
    assert "isotropy subgroup" in msg
    # 中文衔接语保留
    assert "共找到" in msg
    # 不包含纯中文版本的全名
    assert "各向同性子群" not in msg


def test_mixed_mode_load_message():
    set_language("mixed")
    msg = t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "space group" in msg
    assert "个原子" in msg
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
