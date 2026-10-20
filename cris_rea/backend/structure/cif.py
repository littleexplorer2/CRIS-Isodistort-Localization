"""Bounded CIF 1.1 data-file reader for the explicit ordinary-3D contract."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import re

from .model import (
    AtomicSite, CrystalStructure, Lattice, METRIC_RELATIVE_TOLERANCE,
    StructureError, expand_structure,
)
from .symmetry import parse_symmetry_operation

MAX_CIF_BYTES = 1024 * 1024
_NUMBER = re.compile(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))(?:[eE]([+-]?\d+))?(?:\((\d+)\))?\Z")


@dataclass(frozen=True)
class CifNumber:
    value: Fraction
    standard_uncertainty: Fraction | None
    raw: str


def parse_cif_number(value: str) -> CifNumber:
    """CIF numeric value with optional trailing standard uncertainty, exactly."""
    match = _NUMBER.fullmatch(value) if isinstance(value, str) and len(value) <= 128 else None
    if match is None:
        raise StructureError("invalid_number", "Expected a finite CIF number; missing values and fractions are not numeric CIF tokens.",
                             value=value)
    mantissa, exponent_text, uncertainty = match.groups()
    exponent = int(exponent_text or 0)
    if abs(exponent) > 300 or len(mantissa) > 64 or (uncertainty and len(uncertainty) > 64):
        raise StructureError("number_limit", "Numeric magnitude or precision exceeds this reader's bounds.")
    scale = Fraction(10) ** exponent
    result = Fraction(mantissa) * scale
    decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    su = Fraction(int(uncertainty), 10 ** decimals) * scale if uncertainty else None
    return CifNumber(result, su, value)


@dataclass(frozen=True)
class _Token:
    text: str
    quoted: bool = False


def _tokens(text: str) -> list[_Token]:
    # Preserve raw input separately for provenance; normalize only line endings here.
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    tokens = []
    i = 0
    while i < len(text):
        if text[i] in " \t\n":
            i += 1
            continue
        if text[i] == "#":
            end = text.find("\n", i)
            i = len(text) if end < 0 else end + 1
            continue
        if text[i] == ";" and (i == 0 or text[i - 1] == "\n"):
            end = text.find("\n;", i + 1)
            if end < 0:
                raise StructureError("cif_syntax", "Unterminated semicolon text field.")
            closing = end + 2
            if closing < len(text) and text[closing] not in " \t\n":
                raise StructureError("cif_syntax", "Closing semicolon must be followed by whitespace.")
            tokens.append(_Token(text[i + 1:end], True))
            i = closing
            continue
        if text[i] in "'\"":
            quote = text[i]
            start = i + 1
            i = start
            while i < len(text) and text[i] != "\n":
                if text[i] == quote and (i + 1 == len(text) or text[i + 1] in " \t\n"):
                    tokens.append(_Token(text[start:i], True))
                    i += 1
                    break
                i += 1
            else:
                raise StructureError("cif_syntax", "Unterminated or multiline quoted CIF token.")
            continue
        start = i
        while i < len(text) and text[i] not in " \t\n":
            i += 1
        value = text[start:i]
        if value.startswith(("[", "]", "$")):
            raise StructureError("unsupported_cif_syntax", "Reserved or CIF2 value syntax is unsupported.")
        tokens.append(_Token(value))
    return tokens


def _control(token: _Token) -> bool:
    if token.quoted:
        return False
    value = token.text.lower()
    return (value.startswith(("_", "data_", "save_"))
            or value in ("loop_", "stop_", "global_"))


def _read_block(text: str) -> tuple[str, dict[str, _Token], list[dict[str, tuple[_Token, ...]]]]:
    tokens = _tokens(text)
    if not tokens or tokens[0].quoted or not tokens[0].text.lower().startswith("data_") or len(tokens[0].text) == 5:
        raise StructureError("cif_syntax", "One named data_ block is required.")
    name = tokens[0].text[5:]
    scalars: dict[str, _Token] = {}
    loops: list[dict[str, tuple[_Token, ...]]] = []
    seen: set[str] = set()
    i = 1

    def claim(tag: str) -> None:
        if tag in seen:
            raise StructureError("duplicate_cif_tag", "A CIF tag occurs more than once.", tag=tag)
        seen.add(tag)

    while i < len(tokens):
        token = tokens[i]
        key = token.text.lower()
        if not token.quoted and key.startswith("data_"):
            raise StructureError("multiple_data_blocks", "Select one data block explicitly before parsing.")
        if not token.quoted and key == "loop_":
            i += 1
            tags = []
            while i < len(tokens) and not tokens[i].quoted and tokens[i].text.startswith("_"):
                tag = tokens[i].text.lower()
                claim(tag)
                tags.append(tag)
                i += 1
            if not tags:
                raise StructureError("cif_syntax", "A loop must contain column tags.")
            values = []
            while i < len(tokens) and not _control(tokens[i]):
                values.append(tokens[i])
                i += 1
            if not values or len(values) % len(tags):
                raise StructureError("cif_loop_cardinality", "A loop must contain a non-empty whole number of rows.")
            loops.append({tag: tuple(values[col::len(tags)]) for col, tag in enumerate(tags)})
        elif not token.quoted and key.startswith("_"):
            claim(key)
            if i + 1 >= len(tokens) or _control(tokens[i + 1]):
                raise StructureError("cif_syntax", "Missing scalar value.", tag=key)
            scalars[key] = tokens[i + 1]
            i += 2
        else:
            raise StructureError("unsupported_cif_syntax", "Unexpected value, save frame or reserved control word.",
                                 token=token.text)
    return name, scalars, loops


def _numeric(token: _Token, tag: str) -> CifNumber:
    if token.quoted:
        raise StructureError("quoted_numeric_value", "CIF numeric fields must not be quoted strings.", tag=tag)
    return parse_cif_number(token.text)


def _alias(mapping: dict, names: tuple[str, ...], *, required: bool = False):
    present = [name for name in names if name in mapping]
    if len(present) > 1:
        raise StructureError("ambiguous_cif_alias", "Provide only one spelling of each semantic field.", tags=present)
    if not present:
        if required:
            raise StructureError("missing_cif_field", "Required CIF field is missing.", alternatives=names)
        return None
    return mapping[present[0]]


def parse_cif(
    content: str | bytes, *, metric_relative_tolerance: float = METRIC_RELATIVE_TOLERANCE,
) -> CrystalStructure:
    """Read one independent input; no files, network, database or program calls."""
    if isinstance(content, str):
        try:
            raw = content.encode("ascii")
        except UnicodeEncodeError as exc:
            raise StructureError("unsupported_cif_encoding", "This CIF 1.1 reader requires ASCII.") from exc
    elif isinstance(content, bytes):
        raw = content
    else:
        raise StructureError("invalid_input_type", "Pass CIF text or bytes, not a path or file object.")
    if not raw or len(raw) > MAX_CIF_BYTES:
        raise StructureError("cif_size_limit", "CIF input must be nonempty and at most 1 MiB.")
    if any(value not in (9, 10, 13) and not 32 <= value <= 126 for value in raw):
        raise StructureError("unsupported_cif_encoding", "CIF 1.1 permits only printable ASCII, tab and newlines.")
    text = raw.decode("ascii")
    if re.match(r"\s*#\\#CIF_2", text, re.IGNORECASE):
        raise StructureError("unsupported_cif_version", "CIF2 input is not supported.")
    name, scalars, loops = _read_block(text)
    all_fields = {**scalars}
    for loop in loops:
        all_fields.update(loop)
    for tag, values in all_fields.items():
        if any(marker in tag for marker in (
                "_magn", "_moment", "_ssg", "_superspace", "_modulation", "_mod_",
                "_wave_vector", "_fourier", "_atom_site_displace_", "_atom_site_occ_",
                "_atom_site_rot_")):
            raise StructureError("magnetic_or_modulated_unsupported", "Magnetic and modulated structures are unsupported.", tag=tag)
        if tag in ("_atom_site_disorder_assembly", "_atom_site_disorder_group"):
            entries = values if isinstance(values, tuple) else (values,)
            if any(item.text not in (".", "?") for item in entries):
                raise StructureError("disorder_unsupported", "Explicit disorder models are unsupported.", tag=tag)
    lattice_tags = tuple("_cell_" + suffix for suffix in (
        "length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"))
    numbers = [_numeric(_alias(scalars, (tag,), required=True), tag) for tag in lattice_tags]
    lattice = Lattice(*(number.value for number in numbers),
                      standard_uncertainties=tuple(number.standard_uncertainty for number in numbers))
    operation_aliases = ("_space_group_symop_operation_xyz", "_space_group_symop.operation_xyz",
                         "_symmetry_equiv_pos_as_xyz")
    operation_columns = [loop[tag] for loop in loops for tag in operation_aliases if tag in loop]
    if not operation_columns:
        raise StructureError("missing_symmetry_operations",
                             "An explicit xyz operation loop is required. A space-group number alone is insufficient.")
    if len(operation_columns) != 1:
        raise StructureError("ambiguous_symmetry_operations", "More than one symmetry-operation column was supplied.")
    operations = tuple(parse_symmetry_operation(token.text) for token in operation_columns[0])
    atom_loops = [loop for loop in loops if "_atom_site_label" in loop]
    if len(atom_loops) != 1:
        raise StructureError("missing_or_ambiguous_atom_loop", "Exactly one atom loop with labels is required.")
    atoms = atom_loops[0]
    for tag in ("_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"):
        _alias(atoms, (tag,), required=True)
    sites = []
    notices = []
    for row, label_token in enumerate(atoms["_atom_site_label"]):
        if label_token.text in (".", "?") or atoms["_atom_site_type_symbol"][row].text in (".", "?"):
            raise StructureError("invalid_atom_identity", "Atom identity cannot be unknown.")
        coords = [_numeric(atoms["_atom_site_fract_" + axis][row], "_atom_site_fract_" + axis) for axis in "xyz"]
        occupancy = (_numeric(atoms["_atom_site_occupancy"][row], "_atom_site_occupancy")
                     if "_atom_site_occupancy" in atoms else CifNumber(Fraction(1), None, "1"))
        wyckoff_tokens = _alias(atoms, ("_atom_site_wyckoff_symbol",))
        wyckoff = wyckoff_tokens[row].text if wyckoff_tokens else None
        wyckoff = None if wyckoff in (".", "?") else wyckoff
        multiplicity_tokens = _alias(atoms, ("_atom_site_symmetry_multiplicity", "_atom_site_site_symmetry_multiplicity"))
        multiplicity_source = next((tag for tag in ("_atom_site_symmetry_multiplicity",
                                                   "_atom_site_site_symmetry_multiplicity") if tag in atoms), None)
        multiplicity = None
        if multiplicity_tokens and multiplicity_tokens[row].text not in (".", "?"):
            number = _numeric(multiplicity_tokens[row], "_atom_site_symmetry_multiplicity").value
            if number.denominator != 1:
                raise StructureError("invalid_multiplicity", "Multiplicity must be an integer.")
            multiplicity = int(number)
        sites.append(AtomicSite(label_token.text, atoms["_atom_site_type_symbol"][row].text,
                                tuple(number.value for number in coords), occupancy.value, wyckoff, multiplicity,
                                tuple(number.standard_uncertainty for number in coords), occupancy.standard_uncertainty,
                                multiplicity_source))
        if occupancy.value < 1:
            notices.append(f"partial_single_component_occupancy:{label_token.text}")
        if multiplicity_source == "_atom_site_symmetry_multiplicity":
            notices.append(f"legacy_multiplicity_tag:{label_token.text}")
    number_token = _alias(scalars, ("_space_group_it_number", "_symmetry_int_tables_number", "_space_group.it_number"))
    group_number = None
    if number_token and number_token.text not in (".", "?"):
        number = _numeric(number_token, "_space_group_IT_number").value
        if number.denominator != 1 or not 1 <= number <= 230:
            raise StructureError("invalid_space_group_number", "An ordinary 3D group number must be an integer from 1 to 230.")
        group_number = int(number)
    name_token = _alias(scalars, ("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m", "_space_group.name_h-m_alt"))
    group_name = name_token.text if name_token and name_token.text not in (".", "?") else None
    result = CrystalStructure(name, lattice, tuple(sites), operations, group_number, group_name,
                              hashlib.sha256(raw).hexdigest(), tuple(notices))
    # This also checks mixed/overlapping orbits and any declared multiplicity.
    expand_structure(result, metric_relative_tolerance=metric_relative_tolerance)
    return result
