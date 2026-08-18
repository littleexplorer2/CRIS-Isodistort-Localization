"""
isocore.i18n - 中英双语 + 中英混杂国际化支持

三种输出模式：
- "zh"    全中文
- "en"    全英文
- "mixed" 中英混杂：核心科学专有名词用英文（如 space group / Wyckoff position /
          order parameter / supercell / isotropy subgroup），其余用户输入提示
          与输出提示用中文衔接。

三种使用方式均可随时切换语言：
1. 网页端：页面右上角 中/EN/中+EN 切换按钮（前端自己渲染，术语表经 /api/i18n 下发）
2. 终端（main_terminal.py）：主菜单“切换语言”选项，或配置文件 runtime.language
3. Python API：``IsoDistort(language="en")`` / ``iso.set_language("mixed")``

实现：
- 语言状态为进程级全局（线程安全）；网页端每次请求按 ?lang= 设置。
- ``t(key, **kwargs)``：查界面文案目录（messages.py）
- ``term_en2zh()`` / ``term_zh2en()`` / ``translate_term()``：科学术语对照（terms.py）
- mixed 模式：取中文文案后，按术语表把已收录的核心专有名词替换为英文
  （按中文长度降序匹配，避免短词先替换导致错误）
"""

from __future__ import annotations

import threading

from .messages import MESSAGES
from .terms import (
    TERMS_EN2ZH,
    TERMS_ZH2EN,
    term_en2zh,
    term_zh2en,
    translate_term,
)

__all__ = [
    "LANGUAGES",
    "MESSAGES",
    "TERMS_EN2ZH",
    "TERMS_ZH2EN",
    "get_language",
    "set_language",
    "t",
    "t_lang",
    "term_en2zh",
    "term_zh2en",
    "translate_term",
]

# 支持的语言：zh 中文 / en English / mixed 中英混杂
LANGUAGES = ("zh", "en", "mixed")

# 进程级全局语言（线程安全；网页端按请求设置）
_language = "en"
_language_lock = threading.Lock()

# mixed 模式专用：中文术语按长度降序（长词先替换，避免短词误替换）
_MIXED_TERMS = sorted(TERMS_EN2ZH.values(), key=len, reverse=True)


def _to_mixed(text: str) -> str:
    """把中文文案中的核心专有名词替换为英文（mixed 模式）。"""
    for zh_term in _MIXED_TERMS:
        if zh_term in text:
            text = text.replace(zh_term, TERMS_ZH2EN.get(zh_term, zh_term))
    return text


def set_language(language: str) -> None:
    """设置全局语言（"zh" / "en" / "mixed"）。"""
    global _language  # noqa: PLW0603 - 模块级全局语言状态的标准写法（线程锁保护）
    if language not in LANGUAGES:
        raise ValueError(f"不支持的语言: {language}，可选 {LANGUAGES}")
    with _language_lock:
        _language = language


def get_language() -> str:
    """获取当前全局语言。"""
    with _language_lock:
        return _language


def _resolve_template(key: str, lang: str) -> str | None:
    """取文案模板，支持跨语言回退。"""
    template = MESSAGES.get(lang, {}).get(key)
    if template is None:
        fallback = "zh" if lang in ("en", "mixed") else "en"
        template = MESSAGES.get(fallback, {}).get(key)
    return template


def _render(template: str, language: str, **kwargs) -> str:
    try:
        text = template.format(**kwargs)
    except (KeyError, IndexError, ValueError):
        text = template
    if language == "mixed":
        text = _to_mixed(text)
    return text


def t(key: str, **kwargs) -> str:
    """
    按当前语言取界面文案并格式化。

    Args:
        key: 文案键（见 messages.py）
        **kwargs: str.format 占位符

    Returns:
        str: 当前语言的文案；缺失时回退到另一语言，再缺失返回键名本身
    """
    lang = get_language()
    template = _resolve_template(key, lang)
    if template is None:
        return key
    return _render(template, lang, **kwargs)


def t_lang(key: str, language: str, **kwargs) -> str:
    """按指定语言取文案（不修改全局语言；网页端/多线程场景使用）。"""
    template = _resolve_template(key, language)
    if template is None:
        return key
    return _render(template, language, **kwargs)
