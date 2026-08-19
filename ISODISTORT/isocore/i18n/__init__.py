"""
isocore.i18n - 中英双语国际化支持

两种输出模式：
- "zh"    全中文
- "en"    全英文

三种使用方式均可随时切换语言：
1. 网页端：页面右上角语言下拉菜单（zh/en），选中即切换
2. 终端（main_terminal.py）：主菜单“切换语言”选项，或配置文件 runtime.language
3. Python API：``IsoDistort(language="en")`` / ``iso.set_language("zh")``

实现：
- 语言状态为进程级全局（线程安全）；网页端每次请求按 ?lang= 设置。
- ``t(key, **kwargs)``：查界面文案目录（messages.py）
- ``term_en2zh()`` / ``term_zh2en()`` / ``translate_term()``：科学术语对照（terms.py）
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

# 支持的语言：zh 中文 / en English
LANGUAGES = ("zh", "en")

# 进程级全局语言（线程安全；网页端按请求设置）
_language = "en"
_language_lock = threading.Lock()


def set_language(language: str) -> None:
    """设置全局语言（"zh" 中文 / "en" English）。"""
    global _language  # noqa: PLW0603 - 模块级全局语言状态的标准写法（线程锁保护）
    if language not in LANGUAGES:
        raise ValueError(f"不支持的语言: {language}，可选 {LANGUAGES}")
    with _language_lock:
        _language = language


def get_language() -> str:
    """获取当前全局语言。"""
    with _language_lock:
        return _language


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
    template = MESSAGES.get(lang, {}).get(key)
    if template is None:
        fallback = "zh" if lang == "en" else "en"
        template = MESSAGES.get(fallback, {}).get(key)
    if template is None:
        return key
    try:
        return template.format(**kwargs)
    except (KeyError, IndexError, ValueError):
        return template


def t_lang(key: str, language: str, **kwargs) -> str:
    """按指定语言取文案（不修改全局语言；网页端/多线程场景使用）。"""
    template = MESSAGES.get(language, {}).get(key)
    if template is None:
        fallback = "zh" if language == "en" else "en"
        template = MESSAGES.get(fallback, {}).get(key)
    if template is None:
        return key
    try:
        return template.format(**kwargs)
    except (KeyError, IndexError, ValueError):
        return template
