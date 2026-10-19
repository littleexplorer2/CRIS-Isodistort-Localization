"""The shared UTF-8 JSON contract for workflow DTOs and HTTP adapters."""
from __future__ import annotations

import json


def encode_json(value) -> bytes:
    """Encode valid JSON without silently replacing invalid Unicode or numbers."""
    try:
        return json.dumps(
            value, ensure_ascii=False, default=str, allow_nan=False,
        ).encode("utf-8")
    except UnicodeEncodeError as exc:
        raise ValueError("JSON text must contain only Unicode scalar values") from exc


def response_data(value):
    """Detach a DTO into its actual JSON values before committing its action.

    Resolve the existing string fallback here, so an adapter cannot encounter
    a second, stateful conversion after the transaction has already committed.
    """
    return json.loads(encode_json(value))
