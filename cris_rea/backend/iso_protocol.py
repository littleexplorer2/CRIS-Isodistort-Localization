"""Exact protocol boundary for official ISO 9.6.1 space-group queries.

Scientific contract: matrices act on column vectors in parent conventional
fractional coordinates. Basis matrices have the printed basis vectors as
columns. Seitz translations retain their full rational values, including
integers; reducing them modulo the parent cell would lose supercell phases.
Printed OPD vectors are metadata only because ISO rounds some coefficients.
See docs/ISO_COMMAND_CONTRACT.md and evidence/phase2/iso-contract/.
"""
from dataclasses import dataclass
from fractions import Fraction
import re

Vector = tuple[Fraction, Fraction, Fraction]
Matrix = tuple[Vector, Vector, Vector]


class IsoProtocolError(ValueError):
    """The official output or requested query does not meet this contract."""


@dataclass(frozen=True)
class SeitzOperation:
    rotation: Matrix
    translation: Vector


@dataclass(frozen=True)
class ParentGroup:
    number: int
    symbol: str
    primitive_basis: Matrix
    operations: tuple[SeitzOperation, ...]


@dataclass(frozen=True)
class IsotropyResult:
    parent_number: int
    parent_symbol: str
    irrep: str
    opd: str
    opd_vector: str
    subgroup_number: int
    subgroup_symbol: str
    index: int
    size: int
    basis: Matrix
    origin: Vector
    operations: tuple[SeitzOperation, ...]


@dataclass(frozen=True)
class IrrepMetadata:
    parent_number: int
    parent_symbol: str
    irrep: str
    wavevector_expression: str
    parameter_symbols: tuple[str, ...]
    dimension: int
    # Each row's columns correspond to canonical parameters (a, b, g).
    wavevector_coefficients: Matrix
    wavevector_offset: Vector

    def k_vector(self, k_parameters=()):
        """Evaluate the nominal star representative, not the selected OPD arm."""
        if not isinstance(k_parameters, (tuple, list)) \
                or len(k_parameters) != len(self.parameter_symbols):
            raise IsoProtocolError(
                f"This irrep requires {len(self.parameter_symbols)} k parameters "
                f"in order {self.parameter_symbols}."
            )
        provided = dict(zip(self.parameter_symbols, map(_fraction, k_parameters)))
        return tuple(self.wavevector_offset[i] + sum(
            self.wavevector_coefficients[i][j] * provided.get(symbol, Fraction(0))
            for j, symbol in enumerate(("a", "b", "g"))
        ) for i in range(3))


# SCREEN only controls presentation; parsers also accept wrapped operations.
# 1000 was verified against the retained binary; avoid untested extreme widths.
_PREFIX = (
    "PAGE NOBREAK\nSCREEN 1000\n"
    "LABEL VECTOR CONVENTIONAL\nLABEL ELEMENT INTERNATIONAL\n"
)
_RATIONAL = re.compile(r"[+-]?\d+(?:/[1-9]\d*)?\Z")
_LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9+\-]{0,31}\Z")
_OPD = re.compile(r"[A-Za-z0-9]{1,16}\Z")


def _group(value):
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 230:
        raise IsoProtocolError("Space-group number must be an integer from 1 to 230.")
    return value


def _label(value, *, direction=False):
    pattern = _OPD if direction else _LABEL
    if not isinstance(value, str) or not pattern.fullmatch(value):
        raise IsoProtocolError("Use one Miller-Love irrep or one named OPD; command text is not accepted.")
    return value.upper()


def _fraction(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction)):
        raise IsoProtocolError("Exact rationals must be integers, Fraction objects, or p/q strings.")
    if isinstance(value, str) and not _RATIONAL.fullmatch(value):
        raise IsoProtocolError("Exact rational text must be an integer or p/q, without decimals.")
    try:
        return Fraction(value)
    except (ZeroDivisionError, ValueError) as exc:
        raise IsoProtocolError("Invalid rational value.") from exc


def _vector(values):
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise IsoProtocolError("A vector must have three exact rational components.")
    return tuple(_fraction(value) for value in values)


def _matrix(values):
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise IsoProtocolError("A matrix must contain three rows.")
    result = tuple(_vector(row) for row in values)
    if _det(result) == 0:
        raise IsoProtocolError("Basis matrix must be invertible.")
    return result


def _det(m):
    return sum(m[0][j] * (
        m[1][(j + 1) % 3] * m[2][(j + 2) % 3]
        - m[1][(j + 2) % 3] * m[2][(j + 1) % 3]
    ) for j in range(3))


def _columns_to_matrix(columns):
    return tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))


def _format_vector(vector):
    return ",".join(str(x) for x in vector)


def build_parent_query(parent):
    return (
        _PREFIX + f"VALUE PARENT {_group(parent)}\n"
        "SHOW PARENT\nSHOW BASIS\nSHOW ELEMENTS\nDISPLAY PARENT\nQUIT\n"
    ).encode("ascii")


def build_irrep_query(parent, irrep):
    """Ask ISO for symbolic k coordinates and dimension before setting k values."""
    return (
        _PREFIX + f"VALUE PARENT {_group(parent)}\n"
        f"VALUE IRREP {_label(irrep)}\n"
        "SHOW PARENT\nSHOW IRREP\nSHOW KPOINT\nSHOW DIMENSION\n"
        "DISPLAY IRREP\nQUIT\n"
    ).encode("ascii")


def build_isotropy_query(parent, irrep, opd, k_parameters=(), allow_generate=True):
    if not isinstance(k_parameters, (tuple, list)) or len(k_parameters) > 3:
        raise IsoProtocolError("Provide zero to three exact k parameters.")
    values = tuple(_fraction(x) for x in k_parameters)
    if not isinstance(allow_generate, bool):
        raise IsoProtocolError("allow_generate must be a boolean.")
    kline = ""
    if values:
        kline = f"VALUE KVALUE {len(values)}," + ",".join(map(str, values)) + "\n"
    # An empty input line confirms creation of a missing parameter-k database.
    # If no confirmation is needed it is a harmless empty ISO command.
    confirmation = "\n" if allow_generate else ""
    ending = "QUIT\n" if allow_generate else "QUIT\nQUIT\n"
    return (
        _PREFIX + f"VALUE PARENT {_group(parent)}\n"
        f"VALUE IRREP {_label(irrep)}\n" + kline
        + f"VALUE DIRECTION {_label(opd, direction=True)}\n"
        "SHOW PARENT\nSHOW IRREP\nSHOW DIRECTION VECTOR\nSHOW SUBGROUP\n"
        "SHOW INDEX\nSHOW SIZE\nSHOW BASIS\nSHOW ORIGIN\nSHOW ELEMENTS\n"
        "DISPLAY ISOTROPY\n" + confirmation + ending
    ).encode("ascii")


def build_direction_query(parent, subgroup, basis, origin):
    matrix = _matrix(basis)
    columns = _columns_to_matrix(matrix)
    point = _vector(origin)
    return (
        _PREFIX + f"VALUE PARENT {_group(parent)}\n"
        f"VALUE SUBGROUP {_group(subgroup)}\n"
        "VALUE BASIS " + " ".join(_format_vector(v) for v in columns) + "\n"
        "VALUE ORIGIN " + _format_vector(point) + "\n"
        "SHOW IRREP\nSHOW KPOINT\nSHOW DIRECTION VECTOR\nDISPLAY DIRECTION\nQUIT\n"
    ).encode("ascii")


def _clean(text):
    if not isinstance(text, str) or len(text) > 20 * 1024 * 1024:
        raise IsoProtocolError("ISO output must be bounded text.")
    if "Current setting is International (new ed.) with conventional basis vectors." not in text:
        raise IsoProtocolError("ISO did not confirm the supported initial conventional setting.")
    if re.search(
        r"(?im)^\s*\**\s*(?:Syntax error|Error(?:\s|:)|This program has bombed|"
        r"Fortran runtime error|Quit display)\b", text
    ):
        raise IsoProtocolError("ISO reported a failed or aborted query.")
    lines = []
    for line in text.splitlines():
        # A star after non-prompt content is ISO's documented truncation marker.
        if line.rstrip().endswith("*") and line.strip("* \t"):
            raise IsoProtocolError("ISO truncated a display line.")
        line = line.lstrip("*").strip()
        if not line:
            continue
        if line.startswith(("Parent ", "Irrep (", "Frequency", "Elements")):
            continue
        lines.append(line)
    return lines


def _body(text, row_pattern):
    lines = _clean(text)
    starts = [i for i, line in enumerate(lines) if re.match(row_pattern, line)]
    if len(starts) != 1:
        raise IsoProtocolError("Expected exactly one ISO data row; received zero or multiple rows.")
    # Command builders issue one display and then quit. Banner and generated-table
    # messages precede the data row. Continuations may contain only the final field.
    return " ".join(lines[starts[0]:])


def _triples(tail):
    items = re.findall(r"\(([^()]*)\)", tail)
    remainder = re.sub(r"\([^()]*\)", "", tail)
    if remainder.strip(" ,\t\r\n"):
        raise IsoProtocolError("Unexpected or incomplete text in ISO coordinate fields.")
    if any(len(item.split(",")) != 3 for item in items):
        raise IsoProtocolError("ISO coordinate fields must be triples.")
    return items


def _numeric_triple(text):
    return _vector([item.strip() for item in text.split(",")])


def parse_seitz_operation(text):
    """Parse ISO International notation without eval or modulo reduction."""
    parts = text.strip().strip("()").split(",")
    if len(parts) != 3:
        raise IsoProtocolError("A Seitz operator needs three coordinate expressions.")
    rotation, translation = [], []
    for part in parts:
        expression = re.sub(r"\s+", "", part).lower()
        terms = re.findall(r"[+-]?[^+-]+", expression)
        if not expression or "".join(terms) != expression:
            raise IsoProtocolError("Malformed Seitz coordinate expression.")
        coefficients = [Fraction(0)] * 3
        offset = Fraction(0)
        for term in terms:
            sign = -1 if term.startswith("-") else 1
            unsigned = term.lstrip("+-")
            variable = re.fullmatch(r"(\d+(?:/[1-9]\d*)?)?([xyz])", unsigned)
            if variable:
                coefficient = _fraction(variable.group(1) or 1)
                coefficients["xyz".index(variable.group(2))] += sign * coefficient
            else:
                offset += sign * _fraction(unsigned)
        rotation.append(tuple(coefficients))
        translation.append(offset)
    matrix = tuple(rotation)
    if abs(_det(matrix)) != 1:
        raise IsoProtocolError("Seitz rotation determinant must be +1 or -1.")
    return SeitzOperation(matrix, tuple(translation))


def _operations(items):
    operations = tuple(parse_seitz_operation(item) for item in items)
    if not operations or len(set(operations)) != len(operations):
        raise IsoProtocolError("ISO must provide nonempty distinct Seitz representatives.")
    return operations


def parse_parent_output(text, expected_parent=None):
    body = _body(text, r"^\d+\s+\S+\s+\(")
    match = re.match(r"^(\d+)\s+(\S+)\s+(.+)\Z", body)
    if not match:
        raise IsoProtocolError("Malformed parent-group row.")
    number = _group(int(match[1]))
    if expected_parent is not None and number != _group(expected_parent):
        raise IsoProtocolError("ISO returned a different parent group.")
    fields = _triples(match[3])
    if len(fields) < 4:
        raise IsoProtocolError("Missing primitive basis or parent operations.")
    primitive = _matrix(_columns_to_matrix([_numeric_triple(t) for t in fields[:3]]))
    return ParentGroup(number, match[2], primitive, _operations(fields[3:]))


def _wavevector_expression(text):
    parts = text[1:-1].split(",")
    if len(parts) != 3:
        raise IsoProtocolError("The nominal wavevector must have three components.")
    rows, offsets = [], []
    for part in parts:
        expression = re.sub(r"\s+", "", part).lower()
        terms = re.findall(r"[+-]?[^+-]+", expression)
        if not expression or "".join(terms) != expression:
            raise IsoProtocolError("Malformed wavevector expression.")
        row, offset = [Fraction(0)] * 3, Fraction(0)
        for term in terms:
            sign = -1 if term.startswith("-") else 1
            unsigned = term.lstrip("+-")
            parameter = re.fullmatch(r"(\d+(?:/[1-9]\d*)?)?\*?([abg])", unsigned)
            divided = re.fullmatch(r"([abg])/([1-9]\d*)", unsigned)
            if parameter:
                symbol = parameter.group(2)
                coefficient = _fraction(parameter.group(1) or 1)
            elif divided:
                symbol = divided.group(1)
                coefficient = Fraction(1, int(divided.group(2)))
            else:
                offset += sign * _fraction(unsigned)
                continue
            row[("a", "b", "g").index(symbol)] += sign * coefficient
        rows.append(tuple(row))
        offsets.append(offset)
    # Canonical alpha/beta/gamma order, independent of where symbols appear.
    symbols = tuple(symbol for j, symbol in enumerate(("a", "b", "g"))
                    if any(row[j] for row in rows))
    return tuple(rows), tuple(offsets), symbols


def parse_irrep_output(text, expected_parent=None, expected_irrep=None):
    body = _body(text, r"^\d+\s+\S+\s+\S+\s+\(")
    match = re.match(r"^(\d+)\s+(\S+)\s+(\S+)\s+(\([^()]*\))\s+(\d+)\Z", body)
    if not match:
        raise IsoProtocolError("Malformed irrep metadata row.")
    parent, irrep = _group(int(match[1])), _label(match[3])
    if expected_parent is not None and parent != _group(expected_parent):
        raise IsoProtocolError("ISO returned a different parent group.")
    if expected_irrep is not None and irrep != _label(expected_irrep):
        raise IsoProtocolError("ISO returned a different irrep.")
    dimension = int(match[5])
    if dimension < 1:
        raise IsoProtocolError("Irrep dimension must be positive.")
    coefficients, offset, symbols = _wavevector_expression(match[4])
    return IrrepMetadata(parent, match[2], irrep, match[4], symbols, dimension,
                         coefficients, offset)


def parse_isotropy_output(text, expected_parent=None, expected_irrep=None, expected_opd=None):
    body = _body(text, r"^\d+\s+\S+\s+\S+\s+\d+\s+\S+\s+\d+\s+\d+\s+\S+")
    match = re.match(
        r"^(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)"
        r"\s+(\S+)\s+(\([^()]*\))\s+(.+)\Z", body
    )
    if not match:
        raise IsoProtocolError("Malformed isotropy-subgroup row.")
    parent, subgroup = _group(int(match[1])), _group(int(match[4]))
    irrep, opd = _label(match[3]), _label(match[8], direction=True)
    if expected_parent is not None and parent != _group(expected_parent):
        raise IsoProtocolError("ISO returned a different parent group.")
    if expected_irrep is not None and irrep != _label(expected_irrep):
        raise IsoProtocolError("ISO returned a different irrep.")
    if expected_opd is not None and opd != _label(expected_opd, direction=True):
        raise IsoProtocolError("ISO returned a different OPD.")
    index, size = int(match[6]), int(match[7])
    if index < 1 or size < 1:
        raise IsoProtocolError("Index and primitive-cell size ratio must be positive.")
    fields = _triples(match[10])
    if len(fields) < 5:
        raise IsoProtocolError("Missing basis, origin, or subgroup operations.")
    basis = _matrix(_columns_to_matrix([_numeric_triple(t) for t in fields[:3]]))
    origin = _numeric_triple(fields[3])
    return IsotropyResult(parent, match[2], irrep, opd, match[9], subgroup, match[5],
                          index, size, basis, origin, _operations(fields[4:]))
