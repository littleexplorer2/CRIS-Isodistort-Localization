"""Exact fixed-space tests for inverse-Landau subgroup feasibility.

For a finite group ``G``, a real representation ``D`` and a candidate
subgroup ``H``, the Reynolds projector

``P_H = (1 / |H|) * sum(D(h), h in H)``

projects the order-parameter space onto ``Fix(H)``.  A non-zero order
parameter whose isotropy subgroup is *exactly* ``H`` exists if and only if

1. ``Fix(H)`` is non-zero, and
2. no element of ``G \\ H`` acts as the identity on all of ``Fix(H)``.

The second condition says that the pointwise stabilizer of ``Fix(H)`` is
``H``.  It is sufficient because a finite union of proper linear subspaces
cannot cover a non-zero real vector space: a generic vector in ``Fix(H)``
then avoids every larger stabilizer.  This is the representation-theoretic
fixed-space criterion used in inverse Landau analysis; it establishes
symmetry feasibility, not energetic or thermodynamic stability.

All group, projector and row-reduction operations in this module use
:class:`fractions.Fraction`.  Matrices act on column order parameters,
``q' = D(g) q``, and group multiplication ``g*h`` means "apply ``h``, then
``g``".  The module deliberately has no dependency on Method-3 result rows
or ISO output and is not connected to a product search path yet.

The homogeneous-strain constructor uses Cartesian rotations and the tensor
law ``epsilon' = R epsilon R.T`` with component order
``(xx, yy, zz, yz, xz, xy)`` (tensor shear, not engineering shear).  Hence
pure translations act trivially and strain alone cannot lower translation
symmetry.  The site-displacement constructor is limited to a finite
commensurate site orbit in a common Cartesian frame; Bloch phases and
incommensurate/superspace representations must instead be supplied as an
explicit, independently derived matrix representation.

Scientific basis: the group-averaging/Reynolds projector and isotropy
subgroup criterion are standard finite-group representation theory; their
use for crystallographic order parameters follows the isotropy-subgroup and
coupled-order-parameter framework of Hatch & Stokes, Phys. Rev. B 65,
014113 (2002), and Stokes & Campbell, Acta Cryst. A73, 4-13 (2017).
"""

from __future__ import annotations

from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import TypeAlias

from ..utils.lattice import as_fraction

GroupElement: TypeAlias = Hashable
RationalVectorN: TypeAlias = tuple[Fraction, ...]
RationalMatrixN: TypeAlias = tuple[tuple[Fraction, ...], ...]
MatrixInput: TypeAlias = Sequence[Sequence[str | int | float | Fraction]]


def _rational_matrix(
    values: MatrixInput,
    *,
    rows: int | None = None,
    columns: int | None = None,
    description: str = "matrix",
) -> RationalMatrixN:
    matrix = tuple(tuple(as_fraction(value) for value in row) for row in values)
    if not matrix:
        raise ValueError(f"{description} must have at least one row")
    width = len(matrix[0])
    if width == 0 or any(len(row) != width for row in matrix):
        raise ValueError(f"{description} must be a non-empty rectangular matrix")
    if rows is not None and len(matrix) != rows:
        raise ValueError(f"{description} must have {rows} rows, got {len(matrix)}")
    if columns is not None and width != columns:
        raise ValueError(f"{description} must have {columns} columns, got {width}")
    return matrix


def _zero_matrix(rows: int, columns: int) -> RationalMatrixN:
    return tuple(tuple(Fraction(0) for _ in range(columns)) for _ in range(rows))


def _identity_matrix(dimension: int) -> RationalMatrixN:
    return tuple(
        tuple(Fraction(row == column) for column in range(dimension))
        for row in range(dimension)
    )


def _transpose(matrix: RationalMatrixN) -> RationalMatrixN:
    return tuple(
        tuple(matrix[row][column] for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    )


def _matrix_multiply(
    left: RationalMatrixN, right: RationalMatrixN
) -> RationalMatrixN:
    if len(left[0]) != len(right):
        raise ValueError("matrix dimensions are incompatible for multiplication")
    return tuple(
        tuple(
            sum(
                left[row][inner] * right[inner][column]
                for inner in range(len(right))
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def _matrix_vector(
    matrix: RationalMatrixN, vector: RationalVectorN
) -> RationalVectorN:
    if len(matrix[0]) != len(vector):
        raise ValueError("matrix and vector dimensions are incompatible")
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(len(vector)))
        for row in range(len(matrix))
    )


def _block_diagonal(matrices: Sequence[RationalMatrixN]) -> RationalMatrixN:
    dimension = sum(len(matrix) for matrix in matrices)
    result = [list(row) for row in _zero_matrix(dimension, dimension)]
    offset = 0
    for matrix in matrices:
        size = len(matrix)
        for row in range(size):
            for column in range(size):
                result[offset + row][offset + column] = matrix[row][column]
        offset += size
    return tuple(tuple(row) for row in result)


def _column_space_basis(matrix: RationalMatrixN) -> tuple[RationalVectorN, ...]:
    """Return independent original columns using exact RREF pivot columns."""
    work = [list(row) for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        source = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if source is None:
            continue
        work[pivot_row], work[source] = work[source], work[pivot_row]
        pivot = work[pivot_row][column]
        work[pivot_row] = [value / pivot for value in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[pivot_row][index]
                for index in range(column_count)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return tuple(
        tuple(matrix[row][column] for row in range(row_count))
        for column in pivot_columns
    )


@dataclass(frozen=True)
class FiniteGroup:
    """A finite group stored as an exact multiplication table.

    The table contains indices into ``elements``.  Use :meth:`from_operation`
    rather than constructing it manually in normal callers; all group axioms
    are then checked before any scientific calculation is performed.
    """

    elements: tuple[GroupElement, ...]
    multiplication_table: tuple[tuple[int, ...], ...]
    identity_index: int
    name: str = ""

    def __post_init__(self) -> None:
        count = len(self.elements)
        if count == 0:
            raise ValueError("a finite group must contain at least one element")
        try:
            unique_count = len(set(self.elements))
        except TypeError as exc:
            raise TypeError("finite-group elements must be hashable") from exc
        if unique_count != count:
            raise ValueError("finite-group elements must be unique")
        if not 0 <= self.identity_index < count:
            raise ValueError("identity_index is outside the element table")
        if len(self.multiplication_table) != count or any(
            len(row) != count for row in self.multiplication_table
        ):
            raise ValueError("multiplication table must be square and match elements")
        if any(
            not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < count
            for row in self.multiplication_table
            for value in row
        ):
            raise ValueError("multiplication table entries must be valid integer indices")

        identity = self.identity_index
        for index in range(count):
            if (
                self.multiplication_table[identity][index] != index
                or self.multiplication_table[index][identity] != index
            ):
                raise ValueError("declared identity does not act as a two-sided identity")
            if not any(
                self.multiplication_table[index][candidate] == identity
                and self.multiplication_table[candidate][index] == identity
                for candidate in range(count)
            ):
                raise ValueError(f"group element {self.elements[index]!r} has no inverse")

        for left in range(count):
            for middle in range(count):
                left_middle = self.multiplication_table[left][middle]
                for right in range(count):
                    lhs = self.multiplication_table[left_middle][right]
                    middle_right = self.multiplication_table[middle][right]
                    rhs = self.multiplication_table[left][middle_right]
                    if lhs != rhs:
                        raise ValueError("multiplication table is not associative")

    @classmethod
    def from_operation(
        cls,
        elements: Iterable[GroupElement],
        identity: GroupElement,
        multiply: Callable[[GroupElement, GroupElement], GroupElement],
        *,
        name: str = "",
    ) -> FiniteGroup:
        """Build and validate a finite group from labelled elements."""
        labels = tuple(elements)
        try:
            indices = {element: index for index, element in enumerate(labels)}
        except TypeError as exc:
            raise TypeError("finite-group elements must be hashable") from exc
        if len(indices) != len(labels):
            raise ValueError("finite-group elements must be unique")
        if identity not in indices:
            raise ValueError("declared identity is not in the element list")
        table: list[tuple[int, ...]] = []
        for left in labels:
            row: list[int] = []
            for right in labels:
                product = multiply(left, right)
                if product not in indices:
                    raise ValueError(
                        f"group is not closed: {left!r} * {right!r} = {product!r}"
                    )
                row.append(indices[product])
            table.append(tuple(row))
        return cls(labels, tuple(table), indices[identity], name=name)

    @property
    def identity(self) -> GroupElement:
        return self.elements[self.identity_index]

    def index(self, element: GroupElement) -> int:
        try:
            return self.elements.index(element)
        except ValueError as exc:
            raise KeyError(f"unknown group element {element!r}") from exc

    def multiply(
        self, left: GroupElement, right: GroupElement
    ) -> GroupElement:
        product = self.multiplication_table[self.index(left)][self.index(right)]
        return self.elements[product]

    def inverse(self, element: GroupElement) -> GroupElement:
        index = self.index(element)
        for candidate in range(len(self.elements)):
            if (
                self.multiplication_table[index][candidate] == self.identity_index
                and self.multiplication_table[candidate][index] == self.identity_index
            ):
                return self.elements[candidate]
        raise RuntimeError("validated finite group unexpectedly lacks an inverse")

    def subgroup(self, elements: Iterable[GroupElement]) -> tuple[GroupElement, ...]:
        """Validate and return a subgroup in canonical parent-group order."""
        supplied = tuple(elements)
        try:
            supplied_set = set(supplied)
        except TypeError as exc:
            raise TypeError("subgroup elements must be hashable") from exc
        if not supplied:
            raise ValueError("a subgroup cannot be empty")
        if len(supplied_set) != len(supplied):
            raise ValueError("subgroup elements must not be repeated")
        unknown = supplied_set.difference(self.elements)
        if unknown:
            raise ValueError(f"subgroup contains unknown elements: {unknown!r}")
        if self.identity not in supplied_set:
            raise ValueError("subgroup must contain the parent-group identity")
        for left in supplied_set:
            if self.inverse(left) not in supplied_set:
                raise ValueError("candidate subgroup is not closed under inverses")
            for right in supplied_set:
                if self.multiply(left, right) not in supplied_set:
                    raise ValueError("candidate subgroup is not closed under multiplication")
        return tuple(element for element in self.elements if element in supplied_set)


@dataclass(frozen=True)
class RationalRepresentation:
    """An exact finite-dimensional rational representation of a finite group."""

    group: FiniteGroup
    matrices: tuple[RationalMatrixN, ...]
    name: str = ""
    verification_method: str = "full_matrix_homomorphism"

    def __post_init__(self) -> None:
        self._validate_shape_and_identity()
        self._validate_homomorphism()

    def _validate_shape_and_identity(self) -> None:
        count = len(self.group.elements)
        if len(self.matrices) != count:
            raise ValueError("one representation matrix is required per group element")
        if not self.matrices:
            raise ValueError("a representation must contain matrices")
        dimension = len(self.matrices[0])
        if dimension == 0:
            raise ValueError("representation dimension must be positive")
        for matrix in self.matrices:
            if len(matrix) != dimension or any(len(row) != dimension for row in matrix):
                raise ValueError("all representation matrices must be square and equal-sized")
        identity = _identity_matrix(dimension)
        if self.matrices[self.group.identity_index] != identity:
            raise ValueError("the group identity must map to the identity matrix")

    def _validate_homomorphism(self) -> None:
        count = len(self.group.elements)
        for left in range(count):
            for right in range(count):
                product = self.group.multiplication_table[left][right]
                if (
                    _matrix_multiply(self.matrices[left], self.matrices[right])
                    != self.matrices[product]
                ):
                    raise ValueError(
                        "matrices do not form a representation: "
                        f"D({self.group.elements[left]!r}) D({self.group.elements[right]!r}) "
                        f"!= D({self.group.elements[product]!r})"
                    )

    @classmethod
    def _from_exact_construction(
        cls,
        group: FiniteGroup,
        matrices: tuple[RationalMatrixN, ...],
        *,
        name: str,
        proof: str,
    ) -> RationalRepresentation:
        """Construct from an exact functor/intertwiner proof.

        This private path never accepts arbitrary caller matrices.  Its users
        first verify the source representations and discrete actions exactly;
        tensor products, direct sums and invariant quotients then preserve the
        homomorphism algebraically.  Shape and identity are still checked on
        the materialised matrices.
        """
        if not proof:
            raise ValueError("an exact-construction proof label is required")
        instance = object.__new__(cls)
        object.__setattr__(instance, "group", group)
        object.__setattr__(instance, "matrices", matrices)
        object.__setattr__(instance, "name", name)
        object.__setattr__(instance, "verification_method", proof)
        instance._validate_shape_and_identity()
        return instance

    @classmethod
    def from_matrices(
        cls,
        group: FiniteGroup,
        matrices: Mapping[GroupElement, MatrixInput] | Sequence[MatrixInput],
        *,
        name: str = "",
    ) -> RationalRepresentation:
        """Parse exact matrices and verify the representation homomorphism."""
        if isinstance(matrices, Mapping):
            missing = [element for element in group.elements if element not in matrices]
            extras = [element for element in matrices if element not in group.elements]
            if missing or extras:
                raise ValueError(
                    f"representation matrix keys mismatch group; missing={missing!r}, "
                    f"extra={extras!r}"
                )
            ordered = [matrices[element] for element in group.elements]
        else:
            ordered = list(matrices)
            if len(ordered) != len(group.elements):
                raise ValueError("one representation matrix is required per group element")
        parsed = tuple(
            _rational_matrix(matrix, description="representation matrix")
            for matrix in ordered
        )
        return cls(group, parsed, name=name)

    @property
    def dimension(self) -> int:
        return len(self.matrices[0])

    def matrix(self, element: GroupElement) -> RationalMatrixN:
        return self.matrices[self.group.index(element)]


@dataclass(frozen=True, init=False)
class ExactCharacterRepresentation:
    """Exact character of a representation with a construction proof.

    Character data is sufficient for fixed-space dimensions and the
    generated-subgroup stabilizer test, while avoiding dense matrices for
    large finite-supercell permutation representations.  Instances can only
    be created by the verified constructors in this module; arbitrary caller
    character tables are deliberately not accepted.
    """

    group: FiniteGroup
    characters: tuple[Fraction, ...]
    dimension: int
    name: str
    verification_method: str

    @classmethod
    def _from_exact_construction(
        cls,
        group: FiniteGroup,
        characters: Sequence[Fraction],
        *,
        dimension: int,
        name: str,
        proof: str,
    ) -> ExactCharacterRepresentation:
        parsed = tuple(as_fraction(value) for value in characters)
        if len(parsed) != len(group.elements):
            raise ValueError("one character value is required per group element")
        if dimension < 1:
            raise ValueError("character representation dimension must be positive")
        if parsed[group.identity_index] != dimension:
            raise ValueError("identity character must equal representation dimension")
        if not proof:
            raise ValueError("an exact character-construction proof is required")
        instance = object.__new__(cls)
        object.__setattr__(instance, "group", group)
        object.__setattr__(instance, "characters", parsed)
        object.__setattr__(instance, "dimension", dimension)
        object.__setattr__(instance, "name", name)
        object.__setattr__(instance, "verification_method", proof)
        return instance


def character_representation(
    representation: RationalRepresentation,
    *,
    name: str | None = None,
) -> ExactCharacterRepresentation:
    """Compress a verified matrix representation to its exact character."""

    characters = tuple(
        sum(matrix[index][index] for index in range(representation.dimension))
        for matrix in representation.matrices
    )
    return ExactCharacterRepresentation._from_exact_construction(
        representation.group,
        characters,
        dimension=representation.dimension,
        name=representation.name if name is None else name,
        proof=f"exact_trace_of:{representation.verification_method}",
    )


def direct_sum_character_representations(
    representations: Sequence[ExactCharacterRepresentation],
    *,
    name: str = "coupled character",
) -> ExactCharacterRepresentation:
    """Return the exact character of a direct sum."""

    components = tuple(representations)
    if not components:
        raise ValueError("at least one character representation is required")
    group = components[0].group
    if any(component.group != group for component in components[1:]):
        raise ValueError("all character representations must use the same group")
    return ExactCharacterRepresentation._from_exact_construction(
        group,
        tuple(
            sum(component.characters[index] for component in components)
            for index in range(len(group.elements))
        ),
        dimension=sum(component.dimension for component in components),
        name=name,
        proof="exact_direct_sum_of_verified_characters",
    )


@dataclass(frozen=True)
class FixedSpaceFeasibility:
    """Exact result of the inverse-Landau fixed-space criterion."""

    subgroup: tuple[GroupElement, ...]
    projector: RationalMatrixN
    fixed_basis: tuple[RationalVectorN, ...]
    pointwise_stabilizer: tuple[GroupElement, ...]
    is_reachable: bool
    obstruction: str | None

    @property
    def fixed_dimension(self) -> int:
        return len(self.fixed_basis)


@dataclass(frozen=True)
class FixedSpaceFeasibilitySummary:
    """Basis-free exact result of the same inverse-Landau criterion.

    This summary is intended for large finite-supercell representations where
    materialising and row-reducing a dense Reynolds projector is unnecessary.
    The dimension is obtained from the exact character average

    ``dim Fix(H) = |H|^-1 sum_h trace(D(h))``.

    An element ``g`` fixes all of ``Fix(H)`` exactly when
    ``dim Fix(<H,g>) == dim Fix(H)``.  The latter follows because
    ``Fix(<H,g>) = Fix(H) intersect Fix(g)``.  Thus this route evaluates the
    identical pointwise-stabilizer criterion without a numerical rank or a
    hidden floating tolerance.
    """

    subgroup: tuple[GroupElement, ...]
    fixed_dimension: int
    pointwise_stabilizer: tuple[GroupElement, ...]
    is_reachable: bool
    obstruction: str | None
    verification_method: str = (
        "exact_character_average_and_generated_subgroup_dimension"
    )


class FixedSpaceFeasibilityAnalyzer:
    """Cache exact character/subgroup calculations for one representation."""

    def __init__(
        self,
        representation: RationalRepresentation | ExactCharacterRepresentation,
    ) -> None:
        self.representation = representation
        self._element_indices = {
            element: index
            for index, element in enumerate(representation.group.elements)
        }
        if isinstance(representation, ExactCharacterRepresentation):
            self._characters = representation.characters
        else:
            self._characters = tuple(
                sum(matrix[index][index] for index in range(representation.dimension))
                for matrix in representation.matrices
            )
        self._dimension_cache: dict[frozenset[int], int] = {}
        self._closure_cache: dict[tuple[frozenset[int], int], frozenset[int]] = {}

    def _fixed_dimension(self, subgroup_indices: frozenset[int]) -> int:
        cached = self._dimension_cache.get(subgroup_indices)
        if cached is not None:
            return cached
        average = sum(
            (self._characters[index] for index in subgroup_indices),
            start=Fraction(0),
        ) / len(subgroup_indices)
        if average.denominator != 1:
            raise RuntimeError(
                "exact subgroup character average is unexpectedly non-integral"
            )
        dimension = int(average)
        if not 0 <= dimension <= self.representation.dimension:
            raise RuntimeError("exact fixed-space dimension is outside valid bounds")
        self._dimension_cache[subgroup_indices] = dimension
        return dimension

    def _generated_subgroup(
        self,
        subgroup_indices: frozenset[int],
        generator_index: int,
    ) -> frozenset[int]:
        cache_key = (subgroup_indices, generator_index)
        cached = self._closure_cache.get(cache_key)
        if cached is not None:
            return cached
        closure = set(subgroup_indices)
        queue: list[int] = []
        if generator_index not in closure:
            closure.add(generator_index)
            queue.append(generator_index)
        table = self.representation.group.multiplication_table
        while queue:
            current = queue.pop()
            for existing in tuple(closure):
                for product in (table[current][existing], table[existing][current]):
                    if product not in closure:
                        closure.add(product)
                        queue.append(product)
        result = frozenset(closure)
        self._closure_cache[cache_key] = result
        return result

    def analyze(
        self,
        subgroup: Iterable[GroupElement],
    ) -> FixedSpaceFeasibilitySummary:
        """Apply the exact criterion using characters rather than dense RREF."""

        group = self.representation.group
        canonical = group.subgroup(subgroup)
        subgroup_indices = frozenset(
            self._element_indices[element] for element in canonical
        )
        fixed_dimension = self._fixed_dimension(subgroup_indices)
        stabilizer_indices = tuple(
            index
            for index in range(len(group.elements))
            if self._fixed_dimension(
                self._generated_subgroup(subgroup_indices, index)
            )
            == fixed_dimension
        )
        stabilizer = tuple(group.elements[index] for index in stabilizer_indices)
        if fixed_dimension == 0:
            obstruction = "Fix(H) is zero-dimensional; no non-zero order parameter exists"
            reachable = False
        else:
            canonical_set = set(canonical)
            extras = tuple(
                element for element in stabilizer if element not in canonical_set
            )
            if extras:
                obstruction = (
                    "elements outside H act identically on all of Fix(H): "
                    f"{extras!r}"
                )
                reachable = False
            else:
                obstruction = None
                reachable = True
        return FixedSpaceFeasibilitySummary(
            subgroup=canonical,
            fixed_dimension=fixed_dimension,
            pointwise_stabilizer=stabilizer,
            is_reachable=reachable,
            obstruction=obstruction,
        )


def fixed_space_projector(
    representation: RationalRepresentation,
    subgroup: Iterable[GroupElement],
) -> RationalMatrixN:
    """Return ``P_H = |H|^-1 sum_h D(h)`` using exact arithmetic."""
    canonical = representation.group.subgroup(subgroup)
    dimension = representation.dimension
    total = [list(row) for row in _zero_matrix(dimension, dimension)]
    for element in canonical:
        matrix = representation.matrix(element)
        for row in range(dimension):
            for column in range(dimension):
                total[row][column] += matrix[row][column]
    denominator = Fraction(len(canonical))
    projector = tuple(
        tuple(value / denominator for value in row) for row in total
    )
    if _matrix_multiply(projector, projector) != projector:
        raise RuntimeError("validated subgroup average is unexpectedly not idempotent")
    return projector


def fixed_space_basis(
    representation: RationalRepresentation,
    subgroup: Iterable[GroupElement],
) -> tuple[RationalVectorN, ...]:
    """Return an exact column basis for ``Fix(H)``."""
    return _column_space_basis(fixed_space_projector(representation, subgroup))


def pointwise_stabilizer(
    representation: RationalRepresentation,
    basis: Iterable[Sequence[str | int | float | Fraction]],
) -> tuple[GroupElement, ...]:
    """Return elements acting as identity on every vector in ``span(basis)``."""
    vectors = tuple(
        tuple(as_fraction(value) for value in vector) for vector in basis
    )
    if any(len(vector) != representation.dimension for vector in vectors):
        raise ValueError("fixed-space basis vector has the wrong dimension")
    return tuple(
        element
        for element in representation.group.elements
        if all(
            _matrix_vector(representation.matrix(element), vector) == vector
            for vector in vectors
        )
    )


def analyze_fixed_space(
    representation: RationalRepresentation,
    subgroup: Iterable[GroupElement],
) -> FixedSpaceFeasibility:
    """Apply the exact fixed-space/pointwise-stabilizer feasibility test."""
    canonical = representation.group.subgroup(subgroup)
    projector = fixed_space_projector(representation, canonical)
    basis = _column_space_basis(projector)
    stabilizer = pointwise_stabilizer(representation, basis)
    if not basis:
        obstruction = "Fix(H) is zero-dimensional; no non-zero order parameter exists"
        reachable = False
    else:
        extras = tuple(element for element in stabilizer if element not in canonical)
        if extras:
            obstruction = (
                "elements outside H act identically on all of Fix(H): "
                f"{extras!r}"
            )
            reachable = False
        else:
            obstruction = None
            reachable = True
    return FixedSpaceFeasibility(
        subgroup=canonical,
        projector=projector,
        fixed_basis=basis,
        pointwise_stabilizer=stabilizer,
        is_reachable=reachable,
        obstruction=obstruction,
    )


def direct_sum_representations(
    representations: Sequence[RationalRepresentation],
    *,
    name: str = "coupled",
) -> RationalRepresentation:
    """Return the direct sum used for coupled order parameters.

    For a subgroup ``H``, the fixed space of the direct sum is the direct sum
    of component fixed spaces.  Its pointwise stabilizer is therefore the
    intersection of their pointwise stabilizers, which is the required
    coupled-order-parameter criterion.
    """
    components = tuple(representations)
    if not components:
        raise ValueError("at least one representation is required")
    group = components[0].group
    if any(component.group != group for component in components[1:]):
        raise ValueError("all coupled representations must use the same finite group")
    matrices = tuple(
        _block_diagonal(
            tuple(component.matrices[index] for component in components)
        )
        for index in range(len(group.elements))
    )
    return RationalRepresentation._from_exact_construction(
        group,
        matrices,
        name=name,
        proof="exact_direct_sum_of_verified_representations",
    )


def _cartesian_rotations(
    group: FiniteGroup,
    rotations: Mapping[GroupElement, MatrixInput] | Sequence[MatrixInput],
) -> RationalRepresentation:
    representation = RationalRepresentation.from_matrices(
        group, rotations, name="Cartesian vector"
    )
    if representation.dimension != 3:
        raise ValueError("Cartesian rotations must be 3x3")
    identity = _identity_matrix(3)
    for element, rotation in zip(group.elements, representation.matrices, strict=True):
        if _matrix_multiply(rotation, _transpose(rotation)) != identity:
            raise ValueError(
                f"Cartesian rotation for {element!r} is not exactly orthogonal"
            )
    return representation


def polar_vector_representation(
    group: FiniteGroup,
    rotations: Mapping[GroupElement, MatrixInput] | Sequence[MatrixInput],
    *,
    name: str = "polar vector",
) -> RationalRepresentation:
    """Return the exact polar-vector representation ``D(g)=R_g``."""
    validated = _cartesian_rotations(group, rotations)
    return RationalRepresentation(group, validated.matrices, name=name)


def symmetric_square_representation(
    vector_representation: RationalRepresentation,
    *,
    name: str = "symmetric square",
) -> RationalRepresentation:
    """Return ``Sym^2(V)`` for an exact three-dimensional vector action.

    In an orthonormal Cartesian basis this is the physical symmetric-strain
    representation.  In crystallographic fractional axes it is an exactly
    conjugate tensor representation because the parent metric is invariant;
    fixed-space dimensions and stabilizers are therefore unchanged without
    rationalising generally irrational Cartesian rotation entries.
    """
    if vector_representation.dimension != 3:
        raise ValueError("a symmetric second-rank tensor requires a 3D vector action")
    components = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
    basis_tensors: list[RationalMatrixN] = []
    for first, second in components:
        tensor = [list(row) for row in _zero_matrix(3, 3)]
        tensor[first][second] = Fraction(1)
        tensor[second][first] = Fraction(1)
        basis_tensors.append(tuple(tuple(row) for row in tensor))

    tensor_matrices: list[RationalMatrixN] = []
    for action in vector_representation.matrices:
        columns: list[RationalVectorN] = []
        for tensor in basis_tensors:
            transformed = _matrix_multiply(
                _matrix_multiply(action, tensor), _transpose(action)
            )
            if transformed != _transpose(transformed):
                raise RuntimeError("transformed symmetric tensor lost symmetry")
            columns.append(
                tuple(transformed[row][column] for row, column in components)
            )
        tensor_matrices.append(
            tuple(
                tuple(columns[column][row] for column in range(6))
                for row in range(6)
            )
        )
    return RationalRepresentation._from_exact_construction(
        vector_representation.group,
        tuple(tensor_matrices),
        name=name,
        proof="exact_symmetric_square_of_verified_vector_representation",
    )


def homogeneous_strain_representation(
    group: FiniteGroup,
    rotations: Mapping[GroupElement, MatrixInput] | Sequence[MatrixInput],
    *,
    name: str = "homogeneous strain",
) -> RationalRepresentation:
    """Build the six-dimensional symmetric Cartesian strain representation."""
    validated = _cartesian_rotations(group, rotations)
    return symmetric_square_representation(validated, name=name)


def _site_permutations(
    group: FiniteGroup,
    permutations: Mapping[GroupElement, Sequence[int]] | Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    if isinstance(permutations, Mapping):
        missing = [element for element in group.elements if element not in permutations]
        extras = [element for element in permutations if element not in group.elements]
        if missing or extras:
            raise ValueError(
                f"site-permutation keys mismatch group; missing={missing!r}, "
                f"extra={extras!r}"
            )
        ordered = [permutations[element] for element in group.elements]
    else:
        ordered = list(permutations)
        if len(ordered) != len(group.elements):
            raise ValueError("one site permutation is required per group element")
    if any(
        not isinstance(value, int) or isinstance(value, bool)
        for permutation in ordered
        for value in permutation
    ):
        raise ValueError("site indices must be integers")
    parsed = tuple(tuple(permutation) for permutation in ordered)
    if not parsed or not parsed[0]:
        raise ValueError("a site orbit must contain at least one site")
    site_count = len(parsed[0])
    expected = set(range(site_count))
    if any(len(permutation) != site_count or set(permutation) != expected for permutation in parsed):
        raise ValueError("each site action must be a permutation of 0..N-1")
    if parsed[group.identity_index] != tuple(range(site_count)):
        raise ValueError("the group identity must fix every site index")
    for left in range(len(group.elements)):
        for right in range(len(group.elements)):
            product = group.multiplication_table[left][right]
            composed = tuple(parsed[left][parsed[right][site]] for site in range(site_count))
            if composed != parsed[product]:
                raise ValueError("site permutations do not respect group multiplication")
    return parsed


def site_vector_character_representation(
    group: FiniteGroup,
    permutations: Mapping[GroupElement, Sequence[int]] | Sequence[Sequence[int]],
    vector_representation: RationalRepresentation,
    *,
    remove_uniform_vectors: bool = False,
    name: str = "site-vector character",
) -> ExactCharacterRepresentation:
    """Build an exact site-vector character without dense site matrices.

    For the tensor product of a site permutation ``P_g`` and vector action
    ``V_g``, ``chi(g) = fixed_sites(g) * trace(V_g)``.  Uniform site vectors
    form an invariant copy of ``V``; the internal-displacement quotient thus
    has character ``(fixed_sites(g) - 1) * trace(V_g)``.  Site actions are
    first checked against the complete exact group multiplication table.
    """

    if vector_representation.group != group:
        raise ValueError("site and vector actions must use the same finite group")
    parsed = _site_permutations(group, permutations)
    site_count = len(parsed[0])
    if remove_uniform_vectors and site_count < 2:
        raise ValueError("at least two sites are required for internal displacements")
    vector_characters = tuple(
        sum(
            matrix[index][index]
            for index in range(vector_representation.dimension)
        )
        for matrix in vector_representation.matrices
    )
    offset = int(remove_uniform_vectors)
    characters = tuple(
        (sum(source == target for source, target in enumerate(permutation)) - offset)
        * vector_character
        for permutation, vector_character in zip(
            parsed, vector_characters, strict=True
        )
    )
    return ExactCharacterRepresentation._from_exact_construction(
        group,
        characters,
        dimension=(site_count - offset) * vector_representation.dimension,
        name=name,
        proof=(
            "exact_verified_site_permutation_tensor_vector_character"
            + ("_mod_uniform_vector_subrepresentation" if offset else "")
        ),
    )


def site_scalar_representation(
    group: FiniteGroup,
    permutations: Mapping[GroupElement, Sequence[int]] | Sequence[Sequence[int]],
    *,
    name: str = "site scalar",
) -> RationalRepresentation:
    """Build the scalar site-permutation representation for a finite orbit."""
    parsed = _site_permutations(group, permutations)
    site_count = len(parsed[0])
    matrices: list[RationalMatrixN] = []
    for permutation in parsed:
        matrix = [list(row) for row in _zero_matrix(site_count, site_count)]
        for source, target in enumerate(permutation):
            matrix[target][source] = Fraction(1)
        matrices.append(tuple(tuple(row) for row in matrix))
    return RationalRepresentation._from_exact_construction(
        group,
        tuple(matrices),
        name=name,
        proof="exact_verified_site_permutation",
    )


def site_vector_representation(
    group: FiniteGroup,
    permutations: Mapping[GroupElement, Sequence[int]] | Sequence[Sequence[int]],
    vector_representation: RationalRepresentation,
    *,
    name: str = "site vector",
) -> RationalRepresentation:
    """Tensor a finite site permutation with an exact vector action."""
    if vector_representation.group != group:
        raise ValueError("site and vector actions must use the same finite group")
    parsed = _site_permutations(group, permutations)
    site_count = len(parsed[0])
    vector_dimension = vector_representation.dimension
    dimension = vector_dimension * site_count
    matrices: list[RationalMatrixN] = []
    for permutation, vector_action in zip(
        parsed, vector_representation.matrices, strict=True
    ):
        matrix = [list(row) for row in _zero_matrix(dimension, dimension)]
        for source, target in enumerate(permutation):
            for output_component in range(vector_dimension):
                for input_component in range(vector_dimension):
                    matrix[
                        vector_dimension * target + output_component
                    ][vector_dimension * source + input_component] = vector_action[
                        output_component
                    ][input_component]
        matrices.append(tuple(tuple(row) for row in matrix))
    return RationalRepresentation._from_exact_construction(
        group,
        tuple(matrices),
        name=name,
        proof="exact_tensor_product_of_verified_site_and_vector_actions",
    )


def remove_uniform_site_vectors(
    representation: RationalRepresentation,
    *,
    site_count: int,
    vector_dimension: int = 3,
    name: str = "internal site vectors",
) -> RationalRepresentation:
    """Quotient a site-vector action by uniform rigid translations.

    The quotient coordinates are ``u_site - u_0`` for every site after the
    first.  The intertwining identity ``C D(g) = D_rel(g) C`` is verified
    exactly, so this function cannot silently remove a non-invariant subspace.
    """
    if site_count < 2:
        raise ValueError("at least two sites are required for internal displacements")
    if vector_dimension < 1 or representation.dimension != site_count * vector_dimension:
        raise ValueError("site/vector dimensions do not match the representation")
    relative_dimension = (site_count - 1) * vector_dimension
    full_dimension = site_count * vector_dimension

    projection = [
        [Fraction(0) for _ in range(full_dimension)]
        for _ in range(relative_dimension)
    ]
    section = [
        [Fraction(0) for _ in range(relative_dimension)]
        for _ in range(full_dimension)
    ]
    for site in range(1, site_count):
        for component in range(vector_dimension):
            relative = (site - 1) * vector_dimension + component
            projection[relative][site * vector_dimension + component] = Fraction(1)
            projection[relative][component] = Fraction(-1)
            section[site * vector_dimension + component][relative] = Fraction(1)
    projection_matrix = tuple(tuple(row) for row in projection)
    section_matrix = tuple(tuple(row) for row in section)

    relative_matrices: list[RationalMatrixN] = []
    for matrix in representation.matrices:
        relative = _matrix_multiply(
            _matrix_multiply(projection_matrix, matrix), section_matrix
        )
        if _matrix_multiply(projection_matrix, matrix) != _matrix_multiply(
            relative, projection_matrix
        ):
            raise ValueError("uniform site-vector subspace is not invariant")
        relative_matrices.append(relative)
    return RationalRepresentation._from_exact_construction(
        representation.group,
        tuple(relative_matrices),
        name=name,
        proof="exact_invariant_quotient_intertwiner",
    )


def site_displacement_representation(
    group: FiniteGroup,
    permutations: Mapping[GroupElement, Sequence[int]] | Sequence[Sequence[int]],
    rotations: Mapping[GroupElement, MatrixInput] | Sequence[MatrixInput],
    *,
    name: str = "commensurate site displacement",
) -> RationalRepresentation:
    """Build ``site permutation tensor polar vector`` in a Cartesian frame.

    ``permutation[source] = target`` and each site occupies three consecutive
    Cartesian components.  This intentionally excludes k-dependent phase
    factors; callers handling Bloch/superspace modes must derive and supply
    those matrices explicitly via :meth:`RationalRepresentation.from_matrices`.
    """
    vector = _cartesian_rotations(group, rotations)
    return site_vector_representation(
        group,
        permutations,
        vector,
        name=name,
    )


__all__ = [
    "ExactCharacterRepresentation",
    "FiniteGroup",
    "FixedSpaceFeasibility",
    "FixedSpaceFeasibilityAnalyzer",
    "FixedSpaceFeasibilitySummary",
    "GroupElement",
    "MatrixInput",
    "RationalMatrixN",
    "RationalRepresentation",
    "RationalVectorN",
    "analyze_fixed_space",
    "character_representation",
    "direct_sum_character_representations",
    "direct_sum_representations",
    "fixed_space_basis",
    "fixed_space_projector",
    "homogeneous_strain_representation",
    "pointwise_stabilizer",
    "polar_vector_representation",
    "remove_uniform_site_vectors",
    "site_displacement_representation",
    "site_scalar_representation",
    "site_vector_character_representation",
    "site_vector_representation",
    "symmetric_square_representation",
]
