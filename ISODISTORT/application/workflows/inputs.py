"""Validate transport inputs before constructing workers or changing state."""
from __future__ import annotations

import copy
import math
import re
import uuid

from backend.models.phase_path import DISTORTION_TYPES, validate_modulation_count

from .wire import encode_json


def nonnegative_integer(value, name: str, *, text: bool = False) -> int:
    if text:
        if not isinstance(value, str) or re.fullmatch(r"[0-9]+", value) is None:
            raise ValueError(f"{name} must be a non-negative integer")
        return int(value)
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a non-negative JSON integer")
    return value


def action_input(action: str, data: dict | None) -> dict:
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise ValueError("request body must be a JSON object")
    # JSON escapes can decode to isolated surrogates; reject them before any
    # worker construction or upload write, including in otherwise unused keys.
    encode_json(data)
    values = copy.deepcopy(data)
    for key in ("generate", "maximal_subgroup_only"):
        if key in values and type(values[key]) is not bool:
            raise ValueError(f"{key} must be a JSON boolean")
    if action == "set_types":
        types = values.get("types", ["strain"])
        if not isinstance(types, list) or any(
            not isinstance(value, str) or value not in DISTORTION_TYPES for value in types
        ):
            raise ValueError("types must be a list of supported distortion types")
        values["types"] = types or ["strain"]
        scope = values.get("scope")
        if scope is not None:
            if not isinstance(scope, dict):
                raise ValueError("scope must be an object mapping distortion types to species")
            normalized = {}
            for name, selected in scope.items():
                if name not in DISTORTION_TYPES:
                    raise ValueError(f"unsupported scope distortion type: {name}")
                if isinstance(selected, str) and selected in ("*", "all"):
                    normalized[name] = ["*"]
                elif isinstance(selected, list) and all(
                    isinstance(species, str) and species for species in selected
                ):
                    normalized[name] = selected
                else:
                    raise ValueError(f"scope.{name} must be a list of species or '*'")
            values["scope"] = normalized
    if action == "method2":
        values["subgroup_idx"] = nonnegative_integer(values.get("subgroup_idx"), "subgroup_idx")
        nmod = values.get("nmod", values.get("number_of_independent_modulations", 0))
        values["nmod"] = validate_modulation_count(nonnegative_integer(nmod, "nmod"))
        if "number_of_independent_modulations" in values:
            alias = nonnegative_integer(values["number_of_independent_modulations"], "nmod")
            if alias != values["nmod"]:
                raise ValueError("nmod and number_of_independent_modulations disagree")
    if action in ("load_cif", "method4"):
        if not isinstance(values.get("content"), str) or not values["content"].strip():
            raise ValueError("CIF content must be a nonempty string")
        if "filename" in values and not isinstance(values["filename"], str):
            raise ValueError("filename must be a string")
    if action == "method4":
        distance = values.get("robust_distance_threshold", 0.25)
        if type(distance) not in (float, int) or not math.isfinite(distance) or distance <= 0:
            raise ValueError("robust_distance_threshold must be a positive finite number")
    return values


def query_value(values: dict[str, list[str]], name: str, default: str | None = None):
    if name not in values:
        return default
    entries = values[name]
    if not isinstance(entries, list) or len(entries) != 1 or not isinstance(entries[0], str):
        raise ValueError(f"{name} must be provided exactly once")
    return entries[0]


def export_job_id(value: str | None) -> str:
    if value is None:
        return str(uuid.uuid4())
    try:
        parsed = uuid.UUID(value)
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError("job_id must be a UUID") from exc
    if str(parsed) != value.lower():
        raise ValueError("job_id must use the canonical UUID format")
    return str(parsed)
