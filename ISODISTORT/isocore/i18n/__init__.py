"""English UI strings for the terminal, web page, and API console output.

``t(key, *args, **kwargs)`` looks up ``MESSAGES`` and formats placeholders.
"""

from __future__ import annotations

from .messages import MESSAGES

__all__ = [
    "MESSAGES",
    "t",
]


def t(key: str, *args, **kwargs) -> str:
    """Return the English string for ``key``, formatted with ``args`` / ``kwargs``.

    Unknown keys are returned unchanged. Format errors fall back to the raw template.
    """
    template = MESSAGES.get(key, key)
    try:
        if args or kwargs:
            return template.format(*args, **kwargs)
        return template
    except (KeyError, IndexError, ValueError):
        return template
