"""Child-fixed displacive modes and (3+d) commensurate lock-ins.

Local ``iso DISPLAY BUSH`` only works at special k.  For a parametric k
(LD, DT, Y, …) the isotropy subgroup is still a 3D space group when the
wavevector is rational; its complete displacive modes are the identity
representation of that child group.  This module:

1. Builds the child supercell from the subgroup basis.
2. Collects parent k-points that fold to Gamma of that supercell
   (nmod=0) or harmonics of the independent modulation(s) (nmod>=1).
3. Lifts ``smodes`` cartesian patterns onto the supercell and keeps
   those that are invariant under the child space group.
4. Splits by parent Wyckoff and labels
   ``Parent[k]IR(opd)[Site:letter:dsp]siteIR(a)``.

The same fixed-space calculation also supplies secondary modes on a parent
Wyckoff orbit for which ``DISPLAY BUSH`` reports ``There is no root mode``.
In that special-k case BUSH remains the authoritative basis on every orbit
where it has a root mode; the Reynolds-projected ``smodes`` space is used only
on entirely rootless orbits.  This is the direct-sum decomposition by parent
Wyckoff orbit, not a compound- or irrep-specific answer table.

No per-compound answers are stored; LD1/C1 and a rootless Al2 orbit are both
instances of the same representation calculation.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from itertools import product

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup

from ..backend.iso_mode_models import (
    ModeIdentity,
    ModeSubspaceValidation,
    SignedColumnAlignment,
    validate_cartesian_subspaces,
    validate_microscopic_mode_source,
    validate_signed_cartesian_columns,
)
from ..backend.iso_wrapper import (
    BushMode,
    DistortionMode,
    MicroscopicDomainExtensionEvidence,
    SubgroupInfo,
    _component_label,
    microscopic_bush_column_digest,
    microscopic_bush_row_digest,
    microscopic_bush_source_float_column_digest,
    microscopic_extension_digest,
    microscopic_source_basis_digest,
)
from ..backend.smodes_wrapper import SmodesAtom, SmodesModeBlock, SmodesWrapper
from ..data.kpoints_official import (
    KPOINT_OFFICIAL,
    official_kparams_to_iso,
    official_special_k_coords,
)
from ..structure.coordinate_transform import build_supercell
from ..utils.config_loader import get_config
from ..utils.exceptions import IsodistortError
from ..utils.lattice import (
    centering_primitive_matrix,
    determinant,
    inverse,
    is_integral,
    multiply,
    rational_matrix,
)
from ..utils.opd_format import k_star_fraction_vectors
from ..utils.schoenflies import hm_symbol
from .search_methods import (
    _k_compatible_with_supercell,
    _param_value_candidates,
)


@dataclass
class ParametricModeResult:
    """Complete displacive modes of a parametric-k isotropy subgroup."""

    modes: list[DistortionMode]
    supercell_displacements: dict[str, np.ndarray]
    labels: dict[str, str]
    nmod: int
    note: str = ""


def _wyckoff_orbit_id(site: dict, ordinal: int = 0) -> str:
    """Return the stable physical-orbit identity, with a legacy fallback.

    Current ``symmetry_info`` producers provide a content-addressed
    ``orbit_id``.  The fallback keeps older API fixtures usable while still
    distinguishing repeated occurrences of the same Wyckoff letter.
    """
    identity = str(site.get("orbit_id") or "").strip()
    if identity:
        return identity
    letter = str(site.get("wyckoff_letter") or site.get("letter") or "")
    species = str(site.get("species") or site.get("element") or "")
    indices = ",".join(
        str(int(value)) for value in (site.get("equivalent_indices") or ())
    )
    representative = str(site.get("representative_index", ordinal))
    return f"legacy:{letter}:{species}:{representative}:{indices}"


def _orbit_key_token(orbit_id: str) -> str:
    """Filesystem/export-safe token for an internal amplitude key."""
    token = re.sub(r"[^A-Za-z0-9_.+-]+", "-", str(orbit_id)).strip("-")
    return token or "orbit"


def _instantiate_bush_modes_by_orbit(
    modes: Sequence[DistortionMode],
    wyckoff_sites: Sequence[dict],
    requested_letters: Sequence[str],
    requested_orbit_ids: Sequence[str] | None = None,
) -> list[DistortionMode]:
    """Give each independent occurrence of a Wyckoff position its own copy.

    ISO ``DISPLAY BUSH`` is parameterized by Wyckoff *letter*.  When the
    parent contains several independent sites on that position, the returned
    symbolic basis applies independently to every physical orbit.  Treating
    the letter as the orbit identity collapses those copies and can also mix
    incompatible free-coordinate values during child-cell mapping.
    """
    requested = {str(value) for value in requested_letters if str(value)}
    requested_orbit_sequence = [
        str(value) for value in (requested_orbit_ids or ()) if str(value)
    ]
    requested_orbits = set(requested_orbit_sequence)
    if len(requested_orbits) != len(requested_orbit_sequence):
        raise ValueError("requested physical orbit identities contain duplicates")
    all_requested_entries = [
        (ordinal, site)
        for ordinal, site in enumerate(wyckoff_sites)
        if str(site.get("wyckoff_letter") or site.get("letter") or "")
        in requested
    ]
    all_orbit_ids = [
        _wyckoff_orbit_id(site, ordinal) for ordinal, site in all_requested_entries
    ]
    if len(set(all_orbit_ids)) != len(all_orbit_ids):
        raise ValueError("physical Wyckoff orbit identities are not unique")
    if requested_orbits and not requested_orbits.issubset(set(all_orbit_ids)):
        missing = sorted(requested_orbits - set(all_orbit_ids))
        raise ValueError(
            "requested physical Wyckoff orbit metadata is missing: " + ", ".join(missing)
        )
    site_entries = [
        (ordinal, site)
        for ordinal, site in all_requested_entries
        if (
            not requested_orbits
            or _wyckoff_orbit_id(site, ordinal) in requested_orbits
        )
    ]
    if not wyckoff_sites:
        unresolved: list[DistortionMode] = []
        for mode in modes:
            identity = mode.mode_identity
            if identity is not None:
                identity = replace(
                    identity,
                    site_irrep=None,
                    copy_index=None,
                    component_index=None,
                    component_label=None,
                    source="unresolved",
                    status="unresolved",
                    reason="physical_orbit_metadata_unavailable",
                )
            unresolved.append(replace(
                mode,
                site_irrep="",
                mode_identity=identity,
            ))
        return unresolved
    if not site_entries:
        if modes:
            raise ValueError("no requested physical Wyckoff orbit can receive the modes")
        return []
    all_letter_counts: dict[str, int] = {}
    for site in wyckoff_sites:
        letter = str(site.get("wyckoff_letter") or site.get("letter") or "")
        all_letter_counts[letter] = all_letter_counts.get(letter, 0) + 1
    by_letter: dict[str, list[tuple[int, dict]]] = {}
    for ordinal, site in site_entries:
        letter = str(site.get("wyckoff_letter") or site.get("letter") or "")
        by_letter.setdefault(letter, []).append((ordinal, site))

    expanded: list[DistortionMode] = []
    for mode in modes:
        letter = str(mode.wyckoff_site or "")
        matches = by_letter.get(letter) or []
        if not matches:
            raise ValueError(
                f"mode Wyckoff letter {letter!r} has no physical orbit metadata"
            )
        repeated = all_letter_counts.get(letter, 0) > 1
        for ordinal, site in matches:
            orbit_id = _wyckoff_orbit_id(site, ordinal)
            key = str(mode.amplitude_key or mode.irrep_label)
            if repeated:
                marker = f"__{letter}__"
                replacement = f"__{letter}@{_orbit_key_token(orbit_id)}__"
                key = (
                    key.replace(marker, replacement, 1)
                    if marker in key
                    else f"{key}__{replacement.strip('_')}"
                )
            identity = getattr(mode, "mode_identity", None)
            if identity is not None:
                identity = replace(identity, orbit_id=orbit_id)
            expanded.append(replace(
                mode,
                amplitude_key=key,
                wyckoff_orbit_id=orbit_id,
                mode_identity=identity,
            ))
    return expanded


_AffineCoordinate = tuple[Fraction, Fraction, Fraction, Fraction]
_SymbolicPointKey = tuple[_AffineCoordinate, _AffineCoordinate, _AffineCoordinate]
_IntegerTranslation = tuple[int, int, int]
_SymbolicPointTransport = tuple[_SymbolicPointKey, _IntegerTranslation]


@dataclass(frozen=True)
class _SymbolicSubspaceProof:
    """The exact symbolic-row domain used for one whole-space comparison."""

    validation: ModeSubspaceValidation
    comparison_domain_kind: str = ""
    reference_points: tuple[_SymbolicPointKey, ...] = ()
    comparison_points: tuple[_SymbolicPointKey, ...] = ()
    reference_full: tuple[np.ndarray, ...] = ()
    reference_comparison: tuple[np.ndarray, ...] = ()
    candidate_comparison: tuple[np.ndarray, ...] = ()
    transports: tuple[_SymbolicPointTransport, ...] = ()
    ordered_reference_basis_digest: str = ""
    ordered_candidate_binding_digest: str = ""


def _ordered_symbolic_basis_digest(
    points: Sequence[_SymbolicPointKey],
    arrays: Sequence[np.ndarray],
) -> str:
    """Bind an ordered symbolic basis to its normalized point domain."""

    return microscopic_extension_digest({
        "points": tuple(points),
        "arrays": tuple(
            np.asarray(array, dtype=float).tolist() for array in arrays
        ),
    })


def _ordered_candidate_binding_digest(
    modes: Sequence[DistortionMode],
) -> str:
    """Bind proof columns to the exact ordered ISO source identities and rows."""

    payload: list[dict[str, object]] = []
    for mode in modes:
        identity = mode.mode_identity
        provenance = mode.microscopic_provenance
        if identity is None or provenance is None:
            raise ValueError("proof_candidate_identity_or_provenance_missing")
        if (
            identity.component_index != provenance.source_column_index
            or identity.global_irrep != provenance.query_irrep_label
            or identity.global_irrep != provenance.source_global_irrep
            or identity.wyckoff_letter != provenance.source_wyckoff_letter
            or identity.site_irrep != provenance.source_site_irrep
        ):
            raise ValueError("proof_candidate_identity_provenance_mismatch")
        if microscopic_bush_row_digest(mode.bush_modes) != provenance.exact_row_digest:
            raise ValueError("proof_candidate_source_row_digest_mismatch")
        if (
            microscopic_bush_source_float_column_digest(mode.bush_modes)
            != provenance.source_float_column_digest
        ):
            raise ValueError("proof_candidate_source_vector_digest_mismatch")
        payload.append({
            "source_order_key": provenance.source_order_key,
            "exact_source_token": provenance.exact_source_token,
            "identity_token": identity.stable_token,
            "exact_row_digest": provenance.exact_row_digest,
            "exact_vector_digest": provenance.exact_vector_digest,
            "source_float_column_digest": provenance.source_float_column_digest,
        })
    return microscopic_extension_digest(payload)


def _parse_affine_coordinate(expression: str) -> _AffineCoordinate:
    """Parse one exact affine ISO coordinate in ``x,y,z``.

    Only rational affine expressions are accepted.  Unsupported syntax fails
    closed rather than comparing two points through a numerical substitution
    that could accidentally identify distinct Wyckoff parameters.
    """

    compact = re.sub(r"\s+", "", str(expression))
    if not compact:
        raise ValueError("empty symbolic coordinate")
    coefficients = {"x": Fraction(0), "y": Fraction(0), "z": Fraction(0)}
    constant = Fraction(0)
    term_re = re.compile(
        r"(?P<sign>[+-]?)"
        r"(?:"
        r"(?:(?P<coefficient>\d+(?:/\d+)?)?\*?(?P<variable>[xyz])"
        r"(?:/(?P<divisor>\d+))?)"
        r"|(?P<constant>\d+(?:/\d+)?)"
        r")"
    )
    offset = 0
    while offset < len(compact):
        match = term_re.match(compact, offset)
        if match is None:
            raise ValueError(f"unsupported affine coordinate {expression!r}")
        sign = -1 if match.group("sign") == "-" else 1
        variable = match.group("variable")
        if variable is None:
            constant += sign * Fraction(match.group("constant"))
        else:
            coefficient = Fraction(match.group("coefficient") or 1)
            if match.group("divisor"):
                coefficient /= int(match.group("divisor"))
            coefficients[variable] += sign * coefficient
        offset = match.end()
    return (
        coefficients["x"],
        coefficients["y"],
        coefficients["z"],
        constant,
    )


def _symbolic_point_key(bush: BushMode) -> _SymbolicPointKey:
    raw = tuple(str(value) for value in bush.point_raw)
    if len(raw) != 3 or not all(value.strip() for value in raw):
        raise ValueError("mode row lacks a complete symbolic representative point")
    return tuple(_parse_affine_coordinate(value) for value in raw)  # type: ignore[return-value]


def _point_mod_integer_key(point: _SymbolicPointKey) -> _SymbolicPointKey:
    return tuple(
        (x, y, z, constant % 1)
        for x, y, z, constant in point
    )  # type: ignore[return-value]


def _integer_point_translation(
    reference: _SymbolicPointKey,
    candidate: _SymbolicPointKey,
) -> _IntegerTranslation:
    """Return the exact conventional-cell translation candidate-reference."""

    translation: list[int] = []
    for reference_coord, candidate_coord in zip(reference, candidate, strict=True):
        if reference_coord[:3] != candidate_coord[:3]:
            raise ValueError("symbolic points have different affine parameters")
        delta = candidate_coord[3] - reference_coord[3]
        if delta.denominator != 1:
            raise ValueError("symbolic points do not differ by an integer translation")
        translation.append(int(delta))
    return tuple(translation)  # type: ignore[return-value]


def _exact_mode_k_identity(
    mode: DistortionMode,
) -> tuple[int, tuple[Fraction, Fraction, Fraction]]:
    """Return one mode's exact parent SG and printed representative k vector."""

    identity = mode.mode_identity
    parent_sg = int(identity.parent_sg) if identity is not None else 0
    raw = tuple(identity.k_coordinates) if identity is not None else ()
    if not raw:
        raw = tuple(
            token.strip()
            for token in str(mode.k_coords_label or "").split(",")
            if token.strip()
        )
    if parent_sg <= 0 or len(raw) != 3:
        raise ValueError("mode lacks one exact parent-SG/k identity")
    try:
        k_coordinates = tuple(Fraction(str(token).strip()) for token in raw)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("mode k coordinates are not exact rational numbers") from exc
    return parent_sg, k_coordinates  # type: ignore[return-value]


def _common_exact_k_identity(
    candidate_modes: Sequence[DistortionMode],
) -> tuple[int, tuple[Fraction, Fraction, Fraction]]:
    """Return the microscopic columns' common exact k identity.

    A BUSH table can contain secondary irreps while its legacy mode metadata
    still carries the selected primary k.  The ISO microscopic block is the
    per-global-irrep identity source used by the verified-label gate, so only
    those candidate columns may authorize a nonzero symbolic transporter.
    """

    identities = [_exact_mode_k_identity(mode) for mode in candidate_modes]
    if not identities or any(identity != identities[0] for identity in identities[1:]):
        raise ValueError("mode columns do not share one exact parent-SG/k identity")
    return identities[0]


def _is_parent_gamma(
    parent_sg: int,
    k_coordinates: tuple[Fraction, Fraction, Fraction],
) -> bool:
    """Whether k is an exact reciprocal vector of the centered parent lattice."""

    if any(value.denominator != 1 for value in k_coordinates):
        return False
    symbol = (hm_symbol(int(parent_sg)) or "P").replace(" ", "")
    centering = next((char for char in symbol if char.isalpha()), "P").upper()
    half = Fraction(1, 2)
    third = Fraction(1, 3)
    extra_translations = {
        "I": ((half, half, half),),
        "A": ((0, half, half),),
        "B": ((half, 0, half),),
        "C": ((half, half, 0),),
        "F": ((0, half, half), (half, 0, half), (half, half, 0)),
        "R": ((2 * third, third, third), (third, 2 * third, 2 * third)),
    }
    return all(
        sum(
            (
                wavevector * translation
                for wavevector, translation in zip(
                    k_coordinates, centering_translation, strict=True,
                )
            ),
            Fraction(0),
        ).denominator == 1
        for centering_translation in extra_translations.get(centering, ())
    )


def _proven_real_transport_phase(
    translation: _IntegerTranslation,
    *,
    parent_sg: int,
    k_coordinates: tuple[Fraction, Fraction, Fraction],
) -> int:
    """Return the common real Bloch phase of every arm of the exact k star.

    ``ModeIdentity.k_coordinates`` does not bind individual columns to star
    arms.  Such a binding is unnecessary when the requested parent-lattice
    translation acts as the same real scalar on the *complete* star: every
    possible arm assignment then gives the same ``+1`` or ``-1``.  Different
    arm phases or a non-real phase remain unresolved and fail closed.
    """

    if not any(translation):
        return 1
    if _is_parent_gamma(parent_sg, k_coordinates):
        return 1
    phases: set[int] = set()
    for arm in k_star_fraction_vectors(k_coordinates, parent_sg):
        exponent = sum(
            component * shift
            for component, shift in zip(arm, translation, strict=True)
        ) % 1
        if exponent == 0:
            phases.add(1)
        elif exponent == Fraction(1, 2):
            phases.add(-1)
        else:
            raise ValueError(
                "non-Gamma integer translation lacks component-to-star-arm "
                "phase evidence"
            )
    if len(phases) != 1:
        raise ValueError(
            "non-Gamma integer translation lacks component-to-star-arm phase evidence"
        )
    return next(iter(phases))


def _mode_point_rows(mode: DistortionMode) -> dict[_SymbolicPointKey, np.ndarray]:
    rows: dict[_SymbolicPointKey, np.ndarray] = {}
    for bush in mode.bush_modes:
        point = _symbolic_point_key(bush)
        if point in rows:
            raise ValueError("mode contains a duplicate symbolic representative point")
        if len(bush.displacements) != 1:
            raise ValueError("mode columns must be split before subspace validation")
        vector = np.asarray(bush.displacements[0], dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("mode displacement must be one finite 3-vector")
        rows[point] = vector
    if not rows:
        raise ValueError("mode contains no symbolic representative points")
    return rows


def _common_point_rows(
    modes: Sequence[DistortionMode],
) -> tuple[list[_SymbolicPointKey], list[dict[_SymbolicPointKey, np.ndarray]]]:
    row_maps = [_mode_point_rows(mode) for mode in modes]
    if not row_maps:
        return [], []
    points = set(row_maps[0])
    if any(set(rows) != points for rows in row_maps[1:]):
        raise ValueError("mode columns do not share one symbolic point domain")
    return sorted(points), row_maps


def _symbolic_mode_arrays(
    modes: Sequence[DistortionMode],
    point_order: Sequence[_SymbolicPointKey],
) -> list[np.ndarray]:
    """Materialise split ISO columns on one proven symbolic point domain."""

    arrays: list[np.ndarray] = []
    for mode in modes:
        rows = _mode_point_rows(mode)
        if any(point not in rows for point in point_order):
            raise ValueError("mode point is absent from the requested symbolic domain")
        arrays.append(np.asarray([rows[point] for point in point_order], dtype=float))
    return arrays


def _match_candidate_points_to_reference(
    reference_points: Sequence[_SymbolicPointKey],
    reference_rows: Sequence[dict[_SymbolicPointKey, np.ndarray]],
    candidate_points: Sequence[_SymbolicPointKey],
    *,
    parent_sg: int,
    k_coordinates: tuple[Fraction, Fraction, Fraction],
) -> list[_SymbolicPointTransport]:
    """Match a sparse ISO representative domain to BUSH without guessing phase.

    Exact affine expressions win.  A nonzero integer transporter is retained,
    rather than erased modulo one, and must have a proven scalar real Bloch
    phase.  A non-Gamma transport is also proven when every exact k-star arm
    has the same real phase, so no component-to-arm assignment is needed.
    Multiple possible transporters are accepted only when every reference
    column gives the same phase-transported row.
    """

    by_mod: dict[_SymbolicPointKey, list[_SymbolicPointKey]] = {}
    for point in reference_points:
        by_mod.setdefault(_point_mod_integer_key(point), []).append(point)

    edges: dict[int, tuple[_SymbolicPointKey, ...]] = {}
    transports: dict[tuple[int, _SymbolicPointKey], _IntegerTranslation] = {}
    for candidate_index, point in enumerate(candidate_points):
        # A literal common representative is the strongest possible edge and
        # is locked to itself.  Nonliteral representatives may use every
        # modulo-integer BUSH image only after all such images are proven to
        # carry the same transported row on every reference column.
        choices = [point] if point in reference_points else sorted(
            by_mod.get(_point_mod_integer_key(point), [])
        )
        if not choices:
            raise ValueError("candidate point is absent from the BUSH point domain")
        transported_rows: list[tuple[_SymbolicPointKey, int]] = []
        for choice in choices:
            translation = _integer_point_translation(choice, point)
            try:
                phase = _proven_real_transport_phase(
                    translation,
                    parent_sg=parent_sg,
                    k_coordinates=k_coordinates,
                )
            except ValueError:
                # This image supplies no scalar proof for an unknown arm
                # assignment.  Another image can still be a valid edge when
                # its translation has one common real phase on the full star.
                continue
            transported_rows.append((choice, phase))
            transports[(candidate_index, choice)] = translation
        if not transported_rows:
            raise ValueError(
                "non-Gamma integer translation lacks component-to-star-arm "
                "phase evidence"
            )
        first_choice, first_phase = transported_rows[0]
        for other, other_phase in transported_rows[1:]:
            if any(
                not np.array_equal(
                    first_phase * rows[first_choice],
                    other_phase * rows[other],
                )
                for rows in reference_rows
            ):
                raise ValueError(
                    "integer-translation point match has an unresolved phase"
                )
        edges[candidate_index] = tuple(choice for choice, _phase in transported_rows)

    # Find a deterministic complete injection rather than greedily choosing
    # the first image in each modulo class.  Exact edges have width one and are
    # considered first; constrained nonexact rows precede flexible ones.
    reference_owner: dict[_SymbolicPointKey, int] = {}
    assignment: dict[int, _SymbolicPointKey] = {}

    def _augment(candidate_index: int, seen: set[_SymbolicPointKey]) -> bool:
        for choice in edges[candidate_index]:
            if choice in seen:
                continue
            seen.add(choice)
            incumbent = reference_owner.get(choice)
            if incumbent is None or _augment(incumbent, seen):
                reference_owner[choice] = candidate_index
                assignment[candidate_index] = choice
                return True
        return False

    candidate_order = sorted(
        range(len(candidate_points)),
        key=lambda index: (len(edges[index]), candidate_points[index]),
    )
    for candidate_index in candidate_order:
        if not _augment(candidate_index, set()):
            raise ValueError(
                "candidate points map non-injectively onto the BUSH domain"
            )
    return [
        (
            assignment[index],
            transports[(index, assignment[index])],
        )
        for index in range(len(candidate_points))
    ]


def _transported_symbolic_mode_arrays(
    modes: Sequence[DistortionMode],
    transports: Sequence[_SymbolicPointTransport],
    *,
    parent_sg: int,
    k_coordinates: tuple[Fraction, Fraction, Fraction],
) -> list[np.ndarray]:
    """Materialise reference columns on the candidate point representatives."""

    arrays: list[np.ndarray] = []
    for mode in modes:
        rows = _mode_point_rows(mode)
        transported_rows: list[np.ndarray] = []
        for point, translation in transports:
            if point not in rows:
                raise ValueError("mode point is absent from the transported domain")
            phase = _proven_real_transport_phase(
                translation,
                parent_sg=parent_sg,
                k_coordinates=k_coordinates,
            )
            transported_rows.append(phase * rows[point])
        arrays.append(np.asarray(transported_rows, dtype=float))
    return arrays


def _symbolic_subspace_proof(
    reference_modes: Sequence[DistortionMode],
    candidate_modes: Sequence[DistortionMode],
    lattice: np.ndarray,
) -> _SymbolicSubspaceProof:
    """Build the exact-row proof used for comparison and domain extension.

    ``DISPLAY BUSH`` commonly expands every child-cell image while
    ``DISPLAY DISTORTION`` prints a smaller representative set.  Comparison is
    therefore performed on the candidate's point domain, after proving that
    restricting the BUSH columns to that domain preserves their full rank.
    """

    reference_points, reference_rows = _common_point_rows(reference_modes)
    candidate_points, _candidate_rows = _common_point_rows(candidate_modes)
    if not reference_points or not candidate_points:
        return _SymbolicSubspaceProof(
            ModeSubspaceValidation(False, 0, 0, 0.0, "empty_point_domain")
        )
    parent_sg, k_coordinates = _common_exact_k_identity(candidate_modes)
    reference_full = _symbolic_mode_arrays(reference_modes, reference_points)
    candidate = _symbolic_mode_arrays(candidate_modes, candidate_points)
    reference_self = validate_cartesian_subspaces(
        reference_full, reference_full, lattice
    )
    full_rank = reference_self.reference_rank
    if full_rank != len(reference_modes):
        return _SymbolicSubspaceProof(
            ModeSubspaceValidation(
                False,
                full_rank,
                0,
                0.0,
                "reference_columns_not_independent",
            )
        )
    candidate_self = validate_cartesian_subspaces(candidate, candidate, lattice)
    if candidate_self.reference_rank != len(candidate_modes):
        return _SymbolicSubspaceProof(
            ModeSubspaceValidation(
                False,
                full_rank,
                candidate_self.reference_rank,
                0.0,
                "candidate_columns_not_independent",
            )
        )
    if len(reference_modes) != len(candidate_modes):
        return _SymbolicSubspaceProof(
            ModeSubspaceValidation(
                False,
                full_rank,
                candidate_self.reference_rank,
                0.0,
                "rank_mismatch",
            )
        )
    ordered_reference_basis_digest = _ordered_symbolic_basis_digest(
        reference_points, reference_full,
    )
    ordered_candidate_binding_digest = _ordered_candidate_binding_digest(
        candidate_modes,
    )

    # Prefer an exact common representative domain.  If restriction to those
    # literally identical symbolic rows is injective on both spaces, it fixes
    # every column coefficient without assigning any star arm or Bloch phase.
    # This is the strongest available proof for non-Gamma blocks whose two ISO
    # displays choose partly different integer-translated representatives.
    common_points = sorted(set(reference_points).intersection(candidate_points))
    if common_points:
        reference_common = _symbolic_mode_arrays(reference_modes, common_points)
        candidate_common = _symbolic_mode_arrays(candidate_modes, common_points)
        reference_common_rank = validate_cartesian_subspaces(
            reference_common, reference_common, lattice
        ).reference_rank
        candidate_common_rank = validate_cartesian_subspaces(
            candidate_common, candidate_common, lattice
        ).reference_rank
        if (
            reference_common_rank == full_rank
            and candidate_common_rank == candidate_self.reference_rank
        ):
            validation = validate_cartesian_subspaces(
                reference_common, candidate_common, lattice,
            )
            return _SymbolicSubspaceProof(
                validation=validation,
                comparison_domain_kind="exact_common_points",
                reference_points=tuple(reference_points),
                comparison_points=tuple(common_points),
                reference_full=tuple(reference_full),
                reference_comparison=tuple(reference_common),
                candidate_comparison=tuple(candidate_common),
                transports=tuple((point, (0, 0, 0)) for point in common_points),
                ordered_reference_basis_digest=ordered_reference_basis_digest,
                ordered_candidate_binding_digest=ordered_candidate_binding_digest,
            )

    # A rank-losing common domain cannot determine the basis change.  Only in
    # that case consider integer-translated rows, retaining the existing
    # fail-closed phase proof for non-Gamma star arms.
    matched_reference_points = _match_candidate_points_to_reference(
        reference_points,
        reference_rows,
        candidate_points,
        parent_sg=parent_sg,
        k_coordinates=k_coordinates,
    )
    reference_restricted = _transported_symbolic_mode_arrays(
        reference_modes,
        matched_reference_points,
        parent_sg=parent_sg,
        k_coordinates=k_coordinates,
    )
    restricted_rank = validate_cartesian_subspaces(
        reference_restricted, reference_restricted, lattice
    ).reference_rank
    if restricted_rank != full_rank:
        candidate_rank = validate_cartesian_subspaces(
            candidate, candidate, lattice
        ).reference_rank
        return _SymbolicSubspaceProof(
            ModeSubspaceValidation(
                False,
                full_rank,
                candidate_rank,
                0.0,
                "candidate_point_domain_loses_reference_rank",
            )
        )
    validation = validate_cartesian_subspaces(reference_restricted, candidate, lattice)
    return _SymbolicSubspaceProof(
        validation=validation,
        comparison_domain_kind="phase_proven_transported_points",
        reference_points=tuple(reference_points),
        comparison_points=tuple(candidate_points),
        reference_full=tuple(reference_full),
        reference_comparison=tuple(reference_restricted),
        candidate_comparison=tuple(candidate),
        transports=tuple(matched_reference_points),
        ordered_reference_basis_digest=ordered_reference_basis_digest,
        ordered_candidate_binding_digest=ordered_candidate_binding_digest,
    )


def validate_symbolic_mode_subspaces(
    reference_modes: Sequence[DistortionMode],
    candidate_modes: Sequence[DistortionMode],
    lattice: np.ndarray,
) -> ModeSubspaceValidation:
    """Compare BUSH and microscopic columns as complete Cartesian spaces."""

    return _symbolic_subspace_proof(
        reference_modes, candidate_modes, lattice,
    ).validation


def _cartesian_column_matrix(
    arrays: Sequence[np.ndarray],
    lattice: np.ndarray,
) -> np.ndarray:
    """Stack fractional displacement arrays as Cartesian mode columns."""

    lattice_array = np.asarray(lattice, dtype=float)
    if lattice_array.shape != (3, 3) or not np.all(np.isfinite(lattice_array)):
        raise ValueError("domain extension requires one finite 3x3 lattice")
    vectors: list[np.ndarray] = []
    row_count: int | None = None
    for item in arrays:
        array = np.asarray(item, dtype=float)
        if array.ndim != 2 or array.shape[1] != 3 or not np.all(np.isfinite(array)):
            raise ValueError("domain-extension modes must be finite (n_points, 3) arrays")
        if row_count is None:
            row_count = int(array.shape[0])
        elif int(array.shape[0]) != row_count:
            raise ValueError("domain-extension modes do not share one point domain")
        vectors.append((array @ lattice_array).reshape(-1))
    if not vectors:
        return np.zeros((0, 0), dtype=float)
    return np.column_stack(vectors)


def _matrix_rank_condition(
    matrix: np.ndarray,
    expected_rank: int,
    *,
    rtol: float,
    atol: float,
    subject: str,
) -> tuple[float, float]:
    """Return condition/threshold after a scale-aware full-rank proof."""

    singular = np.linalg.svd(np.asarray(matrix, dtype=float), compute_uv=False)
    largest = float(singular[0]) if singular.size else 0.0
    threshold = max(float(atol), float(rtol) * largest)
    machine_threshold = (
        np.finfo(float).eps * max(np.asarray(matrix).shape, default=1) * largest
    )
    numerical_rank = int(np.count_nonzero(singular > machine_threshold))
    if numerical_rank != expected_rank or singular.size < expected_rank:
        raise ValueError(f"{subject}_rank_loss")
    smallest = float(singular[expected_rank - 1])
    condition = largest / smallest
    maximum_condition = 1.0 / max(float(rtol), np.finfo(float).eps)
    if (
        not np.isfinite(condition)
        or condition > maximum_condition
        or smallest <= threshold
    ):
        raise ValueError(f"{subject}_ill_conditioned")
    return condition, threshold


def _extend_microscopic_modes_to_reference_domain(
    reference_modes: Sequence[DistortionMode],
    candidate_modes: Sequence[DistortionMode],
    lattice: np.ndarray,
    proof: _SymbolicSubspaceProof,
    *,
    rtol: float = 1e-8,
    atol_angstrom: float = 1e-10,
) -> list[DistortionMode]:
    """Uniquely extend sparse canonical ISO columns over the full BUSH domain.

    The direct ISO provenance remains attached to the sparse source column.
    A separate immutable evidence object records the only accepted operation:
    solve ``R_common @ T = C_common`` on the already proven comparison domain,
    then apply that same nonsingular ``T`` to ``R_full``.  No column ordering,
    missing row is guessed here.  A transported non-Gamma row is accepted only
    when its exact parent k star proves one common real scalar phase.
    """

    if not proof.validation.matched:
        raise ValueError("subspace_not_matched")
    if any(mode.microscopic_domain_extension is not None for mode in candidate_modes):
        raise ValueError("source_column_already_extended")
    if (
        not proof.ordered_candidate_binding_digest
        or _ordered_candidate_binding_digest(candidate_modes)
        != proof.ordered_candidate_binding_digest
    ):
        raise ValueError("proof_candidate_binding_mismatch")
    if not proof.ordered_reference_basis_digest:
        raise ValueError("proof_reference_basis_mismatch")
    if (
        _ordered_symbolic_basis_digest(
            proof.reference_points, proof.reference_full,
        )
        != proof.ordered_reference_basis_digest
    ):
        raise ValueError("proof_reference_basis_mismatch")
    current_reference_full = _symbolic_mode_arrays(
        reference_modes, proof.reference_points,
    )
    if (
        _ordered_symbolic_basis_digest(
            proof.reference_points, current_reference_full,
        )
        != proof.ordered_reference_basis_digest
    ):
        raise ValueError("proof_reference_basis_mismatch")
    candidate_points, _candidate_rows = _common_point_rows(candidate_modes)
    if (
        len(candidate_points) == len(proof.reference_points)
        and set(candidate_points) == set(proof.reference_points)
    ):
        return list(candidate_modes)

    mode_count = len(candidate_modes)
    if mode_count <= 0 or len(reference_modes) != mode_count:
        raise ValueError("mode_count_mismatch")
    if (
        len(proof.reference_comparison) != mode_count
        or len(proof.candidate_comparison) != mode_count
        or len(proof.reference_full) != mode_count
    ):
        raise ValueError("comparison_basis_width_mismatch")

    reference_common = _cartesian_column_matrix(
        proof.reference_comparison, lattice,
    )
    canonical_common = _cartesian_column_matrix(
        proof.candidate_comparison, lattice,
    )
    reference_condition, _reference_threshold = _matrix_rank_condition(
        reference_common,
        mode_count,
        rtol=rtol,
        atol=atol_angstrom,
        subject="reference_common",
    )
    canonical_condition, _canonical_threshold = _matrix_rank_condition(
        canonical_common,
        mode_count,
        rtol=rtol,
        atol=atol_angstrom,
        subject="canonical_common",
    )
    change_of_basis, _residuals, solved_rank, _singular = np.linalg.lstsq(
        reference_common,
        canonical_common,
        rcond=None,
    )
    if solved_rank != mode_count or change_of_basis.shape != (mode_count, mode_count):
        raise ValueError("change_of_basis_not_unique")
    change_condition, _change_threshold = _matrix_rank_condition(
        change_of_basis,
        mode_count,
        rtol=rtol,
        atol=np.finfo(float).eps,
        subject="change_of_basis",
    )
    reconstructed_common = reference_common @ change_of_basis
    common_residual = float(np.linalg.norm(
        reconstructed_common - canonical_common,
    ))
    common_tolerance = float(
        atol_angstrom
        + rtol * max(
            float(np.linalg.norm(reconstructed_common)),
            float(np.linalg.norm(canonical_common)),
        )
    )
    if not np.isfinite(common_residual) or common_residual > common_tolerance:
        raise ValueError("change_of_basis_residual_exceeds_tolerance")

    # ``T`` was solved in the Cartesian metric so the residual has physical
    # length units.  Right-side basis coefficients are unchanged by the common
    # invertible fractional-to-Cartesian map, hence the same ``T`` extends the
    # full fractional columns without pairing reference and canonical columns.
    reference_full_fractional = np.column_stack([
        np.asarray(array, dtype=float).reshape(-1)
        for array in proof.reference_full
    ])
    extended_matrix = reference_full_fractional @ change_of_basis
    if not np.all(np.isfinite(extended_matrix)):
        raise ValueError("extended_basis_not_finite")
    point_count = len(proof.reference_points)
    extended_arrays = tuple(
        extended_matrix[:, column].reshape(point_count, 3)
        for column in range(mode_count)
    )
    extended_self = validate_cartesian_subspaces(
        extended_arrays, extended_arrays, lattice,
        rtol=rtol, atol=atol_angstrom,
    )
    if extended_self.reference_rank != mode_count:
        raise ValueError("extended_basis_rank_loss")
    full_equivalence = validate_cartesian_subspaces(
        proof.reference_full, extended_arrays, lattice,
        rtol=rtol, atol=atol_angstrom,
    )
    if not full_equivalence.matched:
        raise ValueError("extended_basis_subspace_mismatch")

    templates = {
        _symbolic_point_key(bush): bush
        for bush in reference_modes[0].bush_modes
    }
    if set(templates) != set(proof.reference_points):
        raise ValueError("reference_point_template_mismatch")
    preliminary: list[DistortionMode] = []
    for mode, array in zip(candidate_modes, extended_arrays, strict=True):
        bushes = [
            replace(
                templates[point],
                irrep_label=mode.irrep_label,
                opd_symbol=mode.opd_symbol,
                wyckoff_letter=mode.wyckoff_site,
                displacements=[[float(value) for value in array[row]]],
            )
            for row, point in enumerate(proof.reference_points)
        ]
        preliminary.append(replace(
            mode,
            basis_vectors=[list(bush.displacements[0]) for bush in bushes],
            bush_modes=bushes,
            microscopic_domain_extension=None,
        ))

    parent_sg, k_coordinates = _common_exact_k_identity(candidate_modes)
    if proof.comparison_domain_kind == "exact_common_points":
        extended_comparison = _symbolic_mode_arrays(
            preliminary, proof.comparison_points,
        )
    elif proof.comparison_domain_kind == "phase_proven_transported_points":
        extended_comparison = _transported_symbolic_mode_arrays(
            preliminary,
            proof.transports,
            parent_sg=parent_sg,
            k_coordinates=k_coordinates,
        )
    else:
        raise ValueError("comparison_domain_unresolved")
    restricted_extended = _cartesian_column_matrix(extended_comparison, lattice)
    restriction_residual = float(np.linalg.norm(
        restricted_extended - canonical_common,
    ))
    restriction_tolerance = float(
        atol_angstrom
        + rtol * max(
            float(np.linalg.norm(restricted_extended)),
            float(np.linalg.norm(canonical_common)),
        )
    )
    if (
        not np.isfinite(restriction_residual)
        or restriction_residual > restriction_tolerance
    ):
        raise ValueError("extended_basis_restriction_mismatch")

    point_payload = tuple(proof.reference_points)
    comparison_point_payload = tuple(proof.comparison_points)
    transport_payload = tuple(
        {
            "candidate_point": candidate_point,
            "reference_point": reference_point,
            "integer_translation": translation,
            "proven_real_phase": _proven_real_transport_phase(
                translation,
                parent_sg=parent_sg,
                k_coordinates=k_coordinates,
            ),
        }
        for candidate_point, (reference_point, translation) in zip(
            proof.comparison_points, proof.transports, strict=True,
        )
    )
    reference_basis_digest = _ordered_symbolic_basis_digest(
        point_payload, proof.reference_full,
    )
    source_basis_digest = microscopic_source_basis_digest(candidate_modes)
    change_digest = microscopic_extension_digest(change_of_basis.tolist())
    comparison_point_digest = microscopic_extension_digest(comparison_point_payload)
    full_point_digest = microscopic_extension_digest(point_payload)
    transport_digest = microscopic_extension_digest(transport_payload)
    extended_column_digests = tuple(
        microscopic_bush_column_digest(mode.bush_modes) for mode in preliminary
    )
    extended_basis_digest = microscopic_extension_digest({
        "full_point_digest": full_point_digest,
        "ordered_column_digests": extended_column_digests,
    })
    maximum_residual = max(common_residual, restriction_residual)
    residual_tolerance = max(common_tolerance, restriction_tolerance)
    maximum_normalized_residual = max(
        common_residual / common_tolerance,
        restriction_residual / restriction_tolerance,
    )
    maximum_condition = 1.0 / max(float(rtol), np.finfo(float).eps)
    result: list[DistortionMode] = []
    for column, mode in enumerate(preliminary):
        provenance = mode.microscopic_provenance
        if provenance is None:
            raise ValueError("canonical_source_provenance_missing")
        evidence = MicroscopicDomainExtensionEvidence(
            schema_version="bush-domain-extension-v1",
            relation_status="verified_bush_domain_extension",
            comparison_domain_kind=proof.comparison_domain_kind,
            comparison_point_count=len(proof.comparison_points),
            full_point_count=len(proof.reference_points),
            mode_count=mode_count,
            reference_rank=mode_count,
            canonical_output_column_index=column,
            comparison_point_digest=comparison_point_digest,
            full_point_digest=full_point_digest,
            transport_digest=transport_digest,
            reference_basis_digest=reference_basis_digest,
            canonical_source_basis_digest=source_basis_digest,
            canonical_source_column_token=provenance.exact_source_token,
            canonical_source_order_key=provenance.source_order_key,
            change_of_basis_digest=change_digest,
            extended_basis_digest=extended_basis_digest,
            extended_column_digest=extended_column_digests[column],
            coefficient_column=tuple(
                float(value) for value in change_of_basis[:, column]
            ),
            reference_condition_number=reference_condition,
            canonical_condition_number=canonical_condition,
            change_of_basis_condition_number=change_condition,
            rank_relative_tolerance=float(rtol),
            cartesian_absolute_tolerance_angstrom=float(atol_angstrom),
            maximum_condition_number=maximum_condition,
            maximum_cartesian_residual_angstrom=maximum_residual,
            residual_tolerance_angstrom=residual_tolerance,
            maximum_normalized_residual=maximum_normalized_residual,
            residual_norm="cartesian_frobenius_angstrom",
            matrix_convention="R_common@T=C_common;C_full=R_full@T",
        )
        result.append(replace(mode, microscopic_domain_extension=evidence))
    return result


def validate_symbolic_mode_column_alignment(
    current_modes: Sequence[DistortionMode],
    canonical_modes: Sequence[DistortionMode],
    lattice: np.ndarray,
) -> SignedColumnAlignment:
    """Prove that each BUSH column is one scaled canonical ISO column.

    Point rows are first transported onto the canonical microscopic symbolic
    domain.  The subsequent Cartesian test accepts a column permutation and a
    positive or negative scalar, but rejects every non-diagonal basis change.
    """

    current_points, current_rows = _common_point_rows(current_modes)
    canonical_points, _canonical_rows = _common_point_rows(canonical_modes)
    if not current_points or not canonical_points:
        return SignedColumnAlignment(False, 0, 0, reason="empty_point_domain")
    parent_sg, k_coordinates = _common_exact_k_identity(canonical_modes)
    transports = _match_candidate_points_to_reference(
        current_points,
        current_rows,
        canonical_points,
        parent_sg=parent_sg,
        k_coordinates=k_coordinates,
    )
    current = _transported_symbolic_mode_arrays(
        current_modes,
        transports,
        parent_sg=parent_sg,
        k_coordinates=k_coordinates,
    )
    canonical = _symbolic_mode_arrays(canonical_modes, canonical_points)
    return validate_signed_cartesian_columns(canonical, current, lattice)


def _ordered_microscopic_source_columns(
    modes: Sequence[DistortionMode],
) -> tuple[list[DistortionMode], str | None]:
    """Validate and restore exact ISO query/block/column source order."""

    if not modes:
        return [], "microscopic_source_columns_missing"
    provenances = [mode.microscopic_provenance for mode in modes]
    if any(provenance is None for provenance in provenances):
        return [], "microscopic_column_provenance_missing"
    resolved = [provenance for provenance in provenances if provenance is not None]
    if any(provenance.query_order is None for provenance in resolved):
        return [], "microscopic_query_assignment_unresolved"
    if len({provenance.query_digest for provenance in resolved}) != 1:
        return [], "microscopic_query_provenance_mismatch"
    tokens = [provenance.stable_token for provenance in resolved]
    if len(tokens) != len(set(tokens)):
        return [], "microscopic_source_column_not_unique"
    positions = [
        (
            provenance.query_digest,
            provenance.source_block_order,
            provenance.source_column_index,
        )
        for provenance in resolved
    ]
    if len(positions) != len(set(positions)):
        return [], "microscopic_source_position_not_unique"

    by_block: dict[tuple[str, int], list] = {}
    for provenance in resolved:
        by_block.setdefault(
            (provenance.query_digest, provenance.source_block_order), []
        ).append(provenance)
    for block in by_block.values():
        counts = {provenance.source_block_column_count for provenance in block}
        if len(counts) != 1:
            return [], "microscopic_source_block_width_mismatch"
        count = next(iter(counts))
        if {provenance.source_column_index for provenance in block} != set(range(count)):
            return [], "microscopic_source_block_column_missing"

    for mode, provenance in zip(modes, resolved, strict=True):
        identity = mode.mode_identity
        if identity is None:
            return [], "microscopic_identity_provenance_mismatch"
        if identity.status == "unresolved":
            return [], identity.reason or "microscopic_identity_unresolved"
        try:
            exact_k = tuple(Fraction(token) for token in identity.k_coordinates)
        except (ValueError, ZeroDivisionError):
            return [], "microscopic_k_identity_not_exact"
        if len(exact_k) != 3:
            return [], identity.reason or "microscopic_k_identity_unresolved"
        if (
            identity.component_index != provenance.source_column_index
            or identity.global_irrep != provenance.query_irrep_label
            or identity.global_irrep != provenance.source_global_irrep
            or identity.wyckoff_letter != provenance.source_wyckoff_letter
            or identity.site_irrep != provenance.source_site_irrep
        ):
            return [], "microscopic_identity_provenance_mismatch"
    return sorted(
        modes,
        key=lambda mode: mode.microscopic_provenance.source_order_key,  # type: ignore[union-attr]
    ), None


def _unresolved_modes(
    modes: Sequence[DistortionMode],
    *,
    parent_sg: int,
    reason: str,
) -> list[DistortionMode]:
    unresolved: list[DistortionMode] = []
    for mode in modes:
        k_tokens = tuple(
            token for token in str(mode.k_coords_label or "").split(",") if token != ""
        )
        identity = ModeIdentity.unresolved(
            parent_sg=parent_sg,
            global_irrep=mode.irrep_label,
            k_coordinates=k_tokens,
            wyckoff_letter=mode.wyckoff_site,
            orbit_id=mode.wyckoff_orbit_id,
            reason=reason,
        )
        unresolved.append(replace(
            mode,
            site_irrep="",
            mode_identity=identity,
            microscopic_domain_extension=None,
        ))
    return unresolved


def _select_verified_microscopic_modes(
    parent: Structure,
    bush_modes: Sequence[DistortionMode],
    microscopic_modes: Sequence[DistortionMode],
    *,
    parent_sg: int,
) -> tuple[list[DistortionMode], dict[tuple[str, str], ModeSubspaceValidation]]:
    """Install direction-resolved ISO columns after whole-space validation.

    ``DISPLAY BUSH`` and direction-resolved ``DISPLAY DISTORTION`` are two
    bases of the same child-fixed space.  A multidimensional fixed space may
    use a non-diagonal change of basis, so signed one-column alignment is not
    a scientific requirement here.  Rank and principal-angle equality of the
    complete symbolic Cartesian spaces is required before the ordered exact
    ISO columns, their amplitude keys, and their provenance are installed.
    """

    bush_groups: dict[tuple[str, str], list[DistortionMode]] = {}
    microscopic_groups: dict[tuple[str, str], list[DistortionMode]] = {}
    for mode in bush_modes:
        bush_groups.setdefault((mode.irrep_label, mode.wyckoff_site), []).append(mode)
    for mode in microscopic_modes:
        microscopic_groups.setdefault((mode.irrep_label, mode.wyckoff_site), []).append(mode)

    selected: list[DistortionMode] = []
    diagnostics: dict[tuple[str, str], ModeSubspaceValidation] = {}
    for key, reference in bush_groups.items():
        candidate = microscopic_groups.get(key)
        if not candidate:
            selected.extend(_unresolved_modes(
                reference,
                parent_sg=parent_sg,
                reason="iso_microscopic_block_missing",
            ))
            continue
        candidate, provenance_reason = _ordered_microscopic_source_columns(candidate)
        if provenance_reason is not None:
            diagnostics[key] = ModeSubspaceValidation(
                False, 0, 0, 0.0, provenance_reason
            )
            selected.extend(_unresolved_modes(
                reference,
                parent_sg=parent_sg,
                reason=provenance_reason,
            ))
            continue
        lattice = np.asarray(parent.lattice.matrix, dtype=float)
        proof: _SymbolicSubspaceProof | None = None
        try:
            proof = _symbolic_subspace_proof(reference, candidate, lattice)
            validation = proof.validation
        except (ValueError, np.linalg.LinAlgError):
            validation = ModeSubspaceValidation(
                False, 0, 0, 0.0, "subspace_validation_error"
            )
        diagnostics[key] = validation
        if validation.matched:
            identities = [mode.mode_identity for mode in candidate]
            identity_metadata_complete = all(
                identity is not None
                and identity.source == "iso_microscopic"
                and identity.site_irrep
                and identity.component_label
                and identity.global_irrep == mode.irrep_label
                and identity.wyckoff_letter == mode.wyckoff_site
                for mode, identity in zip(candidate, identities, strict=True)
            )
            identity_tokens = [
                identity.stable_token
                for identity in identities
                if identity is not None
            ]
            if (
                not identity_metadata_complete
                or len(identity_tokens) != len(set(identity_tokens))
            ):
                reason = (
                    "microscopic_identity_metadata_missing"
                    if not identity_metadata_complete
                    else "microscopic_identity_not_unique"
                )
                diagnostics[key] = replace(validation, matched=False, reason=reason)
                selected.extend(_unresolved_modes(
                    reference,
                    parent_sg=parent_sg,
                    reason=reason,
                ))
                continue
            if proof is None:
                raise RuntimeError("matched microscopic subspace lacks its proof")
            try:
                candidate = _extend_microscopic_modes_to_reference_domain(
                    reference,
                    candidate,
                    lattice,
                    proof,
                )
            except (ValueError, np.linalg.LinAlgError) as exc:
                reason = f"microscopic_domain_extension_failed:{exc}"
                diagnostics[key] = replace(
                    validation, matched=False, reason=reason,
                )
                selected.extend(_unresolved_modes(
                    reference,
                    parent_sg=parent_sg,
                    reason=reason,
                ))
                continue
            selected.extend(
                replace(
                    mode,
                    mode_identity=replace(
                        identity,
                        status="verified",
                        reason=None,
                    ),
                )
                for mode, identity in zip(candidate, identities, strict=True)
                if identity is not None
            )
        else:
            selected.extend(_unresolved_modes(
                reference,
                parent_sg=parent_sg,
                reason=validation.reason or "iso_bush_subspace_mismatch",
            ))
    return selected, diagnostics


def _ordered_complete_microscopic_space(
    modes: Sequence[DistortionMode],
) -> list[DistortionMode]:
    """Validate every ISO query block and return one deterministic union.

    A child-fixed orbit can contain columns from several global irreps, hence
    several independent ``DISPLAY DISTORTION`` queries.  Completeness is
    checked inside each immutable query provenance before those query spaces
    are combined; grouping the whole orbit by one digest would incorrectly
    reject a legitimate direct sum.
    """

    by_query: dict[tuple[str, int], list[DistortionMode]] = {}
    for mode in modes:
        provenance = mode.microscopic_provenance
        if provenance is None:
            raise ValueError("rootless microscopic column provenance is missing")
        if provenance.query_order is None:
            raise ValueError("rootless microscopic query order is unresolved")
        by_query.setdefault(
            (provenance.query_digest, int(provenance.query_order)), []
        ).append(mode)
    ordered_queries: list[tuple[int, str, list[DistortionMode]]] = []
    for (digest, query_order), query_modes in by_query.items():
        query_irreps = {
            mode.microscopic_provenance.query_irrep_label  # type: ignore[union-attr]
            for mode in query_modes
        }
        if len(query_irreps) != 1 or None in query_irreps:
            raise ValueError(
                "one rootless microscopic query order names multiple irreps"
            )
        ordered, reason = _ordered_microscopic_source_columns(query_modes)
        if reason is not None:
            raise ValueError(
                f"rootless microscopic query {digest[:12]} is incomplete: {reason}"
            )
        ordered_queries.append((query_order, digest, ordered))
    query_positions = [item[0] for item in ordered_queries]
    if len(query_positions) != len(set(query_positions)):
        raise ValueError("rootless microscopic queries have duplicate source order")
    return [
        mode
        for _query_order, _digest, query_modes in sorted(ordered_queries)
        for mode in query_modes
    ]


def _replace_rootless_supplement_basis(
    supplement: ParametricModeResult,
    microscopic_modes: Sequence[DistortionMode],
    microscopic_displacements: Mapping[str, np.ndarray],
    child_lattice: np.ndarray,
) -> ParametricModeResult:
    """Replace a complete SMODES fixed space by canonical ISO columns.

    The comparison is intentionally per physical parent orbit and on the
    *whole* Cartesian space.  A non-diagonal basis change is valid; missing
    columns, rank loss, a cross-orbit column, or incomplete provenance rejects
    the replacement.  On success the returned amplitude keys and arrays are
    both the canonical ISO basis, keeping user amplitudes in one coordinate
    system.
    """

    reference_modes = list(supplement.modes)
    candidate_modes = list(microscopic_modes)
    if not reference_modes or not candidate_modes:
        raise ValueError("rootless basis replacement requires two nonempty spaces")
    reference_orbits = {
        str(mode.wyckoff_orbit_id or "") for mode in reference_modes
    }
    candidate_orbits = {
        str(mode.wyckoff_orbit_id or "") for mode in candidate_modes
    }
    if "" in reference_orbits or "" in candidate_orbits:
        raise ValueError("rootless basis replacement requires physical orbit IDs")
    if candidate_orbits != reference_orbits:
        raise ValueError(
            "rootless microscopic columns cross or omit physical orbits: "
            f"reference={sorted(reference_orbits)}, candidate={sorted(candidate_orbits)}"
        )

    reference_keys = [str(mode.amplitude_key or mode.irrep_label) for mode in reference_modes]
    candidate_keys = [str(mode.amplitude_key or mode.irrep_label) for mode in candidate_modes]
    if (
        any(not key for key in reference_keys + candidate_keys)
        or len(reference_keys) != len(set(reference_keys))
        or len(candidate_keys) != len(set(candidate_keys))
    ):
        raise ValueError("rootless mode amplitude keys must be nonempty and unique")
    if set(reference_keys) != set(supplement.supercell_displacements):
        raise ValueError("rootless SMODES arrays do not cover the complete reference basis")
    if set(candidate_keys) != set(microscopic_displacements):
        raise ValueError("rootless ISO arrays do not cover the complete canonical basis")

    orbit_order: list[str] = []
    for mode in reference_modes:
        orbit_id = str(mode.wyckoff_orbit_id)
        if orbit_id not in orbit_order:
            orbit_order.append(orbit_id)
    verified_modes: list[DistortionMode] = []
    verified_arrays: dict[str, np.ndarray] = {}
    for orbit_id in orbit_order:
        reference = [
            mode for mode in reference_modes
            if str(mode.wyckoff_orbit_id) == orbit_id
        ]
        candidate = _ordered_complete_microscopic_space([
            mode for mode in candidate_modes
            if str(mode.wyckoff_orbit_id) == orbit_id
        ])
        reference_arrays = [
            np.asarray(
                supplement.supercell_displacements[
                    str(mode.amplitude_key or mode.irrep_label)
                ],
                dtype=float,
            )
            for mode in reference
        ]
        candidate_arrays = [
            np.asarray(
                microscopic_displacements[
                    str(mode.amplitude_key or mode.irrep_label)
                ],
                dtype=float,
            )
            for mode in candidate
        ]
        shapes = {array.shape for array in reference_arrays + candidate_arrays}
        if len(shapes) != 1 or next(iter(shapes), ()) != (
            reference_arrays[0].shape[0],
            3,
        ):
            raise ValueError("rootless fixed-space columns use incompatible atom frames")
        validation = validate_cartesian_subspaces(
            reference_arrays,
            candidate_arrays,
            np.asarray(child_lattice, dtype=float),
        )
        if (
            not validation.matched
            or len(reference) != len(candidate)
            or validation.reference_rank != len(reference)
            or validation.candidate_rank != len(candidate)
        ):
            raise ValueError(
                "rootless ISO space is not the complete child-fixed space for "
                f"orbit {orbit_id!r}: {validation.reason or 'rank_or_dimension_mismatch'}"
            )
        orbit_tokens: set[str] = set()
        for mode, array in zip(candidate, candidate_arrays, strict=True):
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            if identity is None or provenance is None:
                raise ValueError("rootless microscopic identity is incomplete")
            if str(identity.orbit_id) != orbit_id:
                raise ValueError("rootless microscopic identity crosses physical orbits")
            verified_identity = replace(identity, status="verified", reason=None)
            validate_microscopic_mode_source(verified_identity, provenance)
            token = verified_identity.stable_token
            if token in orbit_tokens:
                raise ValueError("rootless microscopic mode identity is not unique")
            orbit_tokens.add(token)
            verified_mode = replace(mode, mode_identity=verified_identity)
            key = str(verified_mode.amplitude_key or verified_mode.irrep_label)
            verified_modes.append(verified_mode)
            verified_arrays[key] = array

    if len(verified_arrays) != len(verified_modes):
        raise ValueError("rootless canonical amplitude keys collide across physical orbits")
    return ParametricModeResult(
        modes=verified_modes,
        supercell_displacements=verified_arrays,
        labels={},
        nmod=supplement.nmod,
        note=(
            supplement.note + "; canonical ISO microscopic basis installed"
        ).strip("; "),
    )


def _canonicalize_rootless_supplement(
    parent: Structure,
    symmetry_info: dict,
    subgroup: SubgroupInfo,
    supplement: ParametricModeResult,
    raw_microscopic_modes: Sequence[DistortionMode],
    rootless_letters: Sequence[str],
    rootless_orbit_ids: Sequence[str],
) -> ParametricModeResult:
    """Map ISO rootless columns into the SMODES child frame and prove equality."""

    from .affine_embeddings import (  # noqa: PLC0415
        embedding_from_identity,
        parent_affine_group,
    )
    from .distortion_mapper import DistortionMapper  # noqa: PLC0415

    sites = list(symmetry_info.get("wyckoff_sites") or [])
    rootless_set = {str(value) for value in rootless_letters}
    relevant_raw = [
        mode for mode in raw_microscopic_modes
        if str(mode.wyckoff_site or "") in rootless_set
    ]
    candidates = _instantiate_bush_modes_by_orbit(
        relevant_raw,
        sites,
        rootless_letters,
        rootless_orbit_ids,
    )
    parent_group = parent_affine_group(parent)
    embedding = embedding_from_identity(subgroup, parent_group)
    basis = subgroup.basis_vectors or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        candidates,
        basis,
        cartesian_tolerance=get_config().symmetry_cartesian_tolerance_angstrom,
        subgroup_operations=embedding.operations,
        subgroup_translation_lattice=[
            [float(value) for value in row] for row in embedding.lattice
        ],
        subgroup_context=subgroup,
    )
    child = build_supercell(parent, basis)
    return _replace_rootless_supplement_basis(
        supplement,
        candidates,
        mapped,
        np.asarray(child.lattice.matrix, dtype=float),
    )


def compute_special_modes_with_rootless_supplement(
    parent: Structure,
    symmetry_info: dict,
    subgroup: SubgroupInfo,
    wyckoff_letters: Sequence[str],
    iso,
    smodes: SmodesWrapper | None,
    *,
    kpoints: Sequence | None = None,
    wyckoff_orbit_ids: Sequence[str] | None = None,
) -> ParametricModeResult:
    """Return complete special-k modes with audited microscopic identities.

    ``DISPLAY BUSH`` is conditioned on the selected primary order parameter.
    Consequently it can legitimately say ``There is no root mode`` for a
    parent Wyckoff orbit even though that orbit carries secondary modes fixed
    by the selected child group.  The displacement representation is a direct
    sum over parent Wyckoff orbits.  ISO microscopic columns are adopted only
    when their complete Cartesian column space matches the corresponding BUSH
    space by rank and principal angles.  BUSH columns remain numerically usable
    but explicitly unresolved when that evidence is missing or inconsistent.
    A wholly rootless orbit is completed from
    ``Fix(H) = image(|H|^-1 sum_h D(h))`` and likewise remains unresolved.

    The returned ``supercell_displacements`` contains only supplemented
    vectors.  Callers must map the retained BUSH modes with the exact subgroup
    Seitz operations and merge the two disjoint key sets.
    """
    requested = list(dict.fromkeys(
        str(value).strip() for value in wyckoff_letters if str(value).strip()
    ))
    requested_orbit_sequence = [
        str(value) for value in (wyckoff_orbit_ids or ()) if str(value)
    ]
    requested_orbits = set(requested_orbit_sequence)
    if len(requested_orbits) != len(requested_orbit_sequence):
        raise ValueError("requested physical orbit identities contain duplicates")
    parent_sg = int(
        symmetry_info.get("space_group_number") or subgroup.parent_sg or 0
    )
    raw_bush_modes = iso.calc_distortion_modes(
        parent_sg,
        subgroup,
        wyckoff_letters=requested,
    )
    microscopic_diagnostic = ""
    try:
        raw_microscopic_modes = iso.calc_microscopic_distortion_modes(
            parent_sg,
            subgroup,
            requested,
        )
    except (
        AttributeError,
        IsodistortError,
        ValueError,
        IndexError,
        np.linalg.LinAlgError,
    ) as exc:
        raw_microscopic_modes = []
        microscopic_diagnostic = (
            "ISO microscopic query failed: "
            f"{type(exc).__name__}: {exc}"
        )
    if not raw_microscopic_modes and not microscopic_diagnostic:
        microscopic_diagnostic = "ISO microscopic query returned no columns"
    if raw_microscopic_modes:
        selected_raw, _microscopic_diagnostics = _select_verified_microscopic_modes(
            parent,
            raw_bush_modes,
            raw_microscopic_modes,
            parent_sg=parent_sg,
        )
    else:
        selected_raw = _unresolved_modes(
            raw_bush_modes,
            parent_sg=parent_sg,
            reason="iso_microscopic_output_unavailable",
        )
    sites = list(symmetry_info.get("wyckoff_sites") or [])
    bush_modes = _instantiate_bush_modes_by_orbit(
        selected_raw, sites, requested, requested_orbit_sequence,
    )
    if sites:
        selected_sites = [
            (ordinal, site)
            for ordinal, site in enumerate(sites)
            if str(site.get("wyckoff_letter") or site.get("letter") or "")
            in requested
            and (
                not requested_orbits
                or _wyckoff_orbit_id(site, ordinal) in requested_orbits
            )
        ]
        rooted_orbits = {
            str(mode.wyckoff_orbit_id or "")
            for mode in bush_modes
            if str(mode.wyckoff_orbit_id or "")
        }
        rootless_sites = [
            (ordinal, site) for ordinal, site in selected_sites
            if _wyckoff_orbit_id(site, ordinal) not in rooted_orbits
        ]
        rootless = list(dict.fromkeys(
            str(site.get("wyckoff_letter") or site.get("letter") or "")
            for _ordinal, site in rootless_sites
        ))
        rootless_orbits = [
            _wyckoff_orbit_id(site, ordinal)
            for ordinal, site in rootless_sites
        ]
    else:
        rooted = {
            str(mode.wyckoff_site or "").strip()
            for mode in bush_modes
            if str(mode.wyckoff_site or "").strip()
        }
        rootless = [letter for letter in requested if letter not in rooted]
        rootless_orbits = []
    if not rootless or smodes is None:
        return ParametricModeResult(
            modes=bush_modes,
            supercell_displacements={},
            labels={},
            nmod=0,
            note=(
                "special k: every requested orbit has a BUSH root"
                if not rootless
                else "special k: rootless orbit supplement unavailable without smodes"
            ),
        )

    supplement = compute_parametric_modes(
        parent,
        symmetry_info,
        subgroup,
        rootless,
        smodes,
        kpoints=kpoints,
        nmod=0,
        wyckoff_orbit_ids=rootless_orbits or None,
    )
    if raw_microscopic_modes and sites and rootless_orbits:
        try:
            supplement = _canonicalize_rootless_supplement(
                parent,
                symmetry_info,
                subgroup,
                supplement,
                raw_microscopic_modes,
                rootless,
                rootless_orbits,
            )
        except (ValueError, IsodistortError, np.linalg.LinAlgError) as exc:
            # The Reynolds/SMODES basis remains a valid numerical fixed space,
            # but without a proven whole-space ISO replacement its individual
            # scientific identities stay explicitly unresolved and export
            # continues to fail closed.
            supplement = replace(
                supplement,
                note=(
                    supplement.note
                    + "; canonical ISO microscopic replacement rejected: "
                    + str(exc)
                ).strip("; "),
            )
    note_parts = [
        "special k: verified microscopic roots plus child-fixed supplement on "
        f"rootless orbit(s) {', '.join(rootless)}"
    ]
    if supplement.note:
        note_parts.append(supplement.note)
    if microscopic_diagnostic:
        note_parts.append(microscopic_diagnostic)
    return ParametricModeResult(
        modes=[*bush_modes, *supplement.modes],
        supercell_displacements=dict(supplement.supercell_displacements),
        labels=dict(supplement.labels),
        nmod=0,
        note="; ".join(note_parts),
    )


@dataclass
class _KCandidate:
    label: str
    coords: np.ndarray
    params: list[str] | None
    coord_tokens: tuple[str, ...]
    # ``coords`` is the website's representative arm (used in labels).
    # ``active_coords`` and ``arm_rotation`` identify a symmetry-equivalent
    # arm that is actually periodic in the requested child cell.  They differ
    # for cases such as P4/mmm DT=(0,b,0) in a 3a x 2b supercell, where the
    # compatible arm is along a*.
    active_coords: np.ndarray | None = None
    arm_rotation: np.ndarray | None = None


def format_k_token(value: float) -> str:
    """Render a k-component as a short fraction (``1/6``) when possible."""
    frac = Fraction(float(value)).limit_denominator(48)
    if abs(float(frac) - float(value)) > 1e-8:
        return f"{float(value):g}"
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def parse_k_token(token: str) -> float:
    return float(Fraction(str(token).strip()))


def opd_free_components(opd_dir_raw: str, opd_symbol: str = "") -> list[str]:
    """Letters that are free in the OPD, e.g. ``(a,b)`` → ``['a','b']``."""
    raw = (opd_dir_raw or "").strip()
    if raw.startswith("(") and ")" in raw:
        body = raw[1:raw.index(")")]
    else:
        body = raw
    letters: list[str] = []
    for part in re.split(r"[,;]", body):
        tok = part.strip()
        if re.fullmatch(r"[a-z]", tok):
            letters.append(tok)
    if letters:
        return letters
    if (opd_symbol or "").upper().startswith("C"):
        return ["a", "b"]
    return ["a"]


def is_harmonic_of(kvec: np.ndarray, fundamentals: Sequence[np.ndarray],
                   *, max_order: int = 12, atol: float = 1e-5) -> bool:
    """True if ``k`` is Gamma or ``n q + G`` for nonzero integer n.

    Gamma (the 3D average structure) is kept for ``nmod>=1``.  Other
    folding k such as M are *not* harmonics of an LD line along c*.
    """
    k = np.asarray(kvec, dtype=float)
    if np.allclose(k, 0.0, atol=atol):
        return True
    for fundamental in fundamentals:
        qvec = np.asarray(fundamental, dtype=float)
        qn = float(np.linalg.norm(qvec))
        if qn < atol:
            continue
        for n in range(-max_order, max_order + 1):
            if n == 0:
                continue
            residual = k - n * qvec
            residual -= np.round(residual)
            if np.allclose(residual, 0.0, atol=atol):
                return True
    return False


def _numeric_k(parent_sg: int, label: str, coords: Sequence[str],
               params: Sequence[str] | None) -> np.ndarray | None:
    tokens = official_special_k_coords(
        parent_sg, label, list(coords), list(params or []),
    )
    if not tokens:
        tokens = [str(c) for c in coords]
    try:
        return np.array([parse_k_token(t) for t in tokens], dtype=float)
    except (ValueError, ZeroDivisionError):
        return None


def _eval_affine_k_component(
    component: str,
    param_names: Sequence[str],
    param_values: Sequence[Fraction],
) -> Fraction | None:
    """Evaluate an affine CDML component such as ``-a+1`` or ``2g``.

    The k-point tables only use rational constants and integer multiples of
    named parameters.  A small parser keeps this deterministic and avoids
    evaluating arbitrary text.
    """
    text = str(component).replace(" ", "")
    if not text:
        return None
    values = dict(zip(param_names, param_values, strict=False))
    total = Fraction(0)
    terms = re.findall(r"[+-]?[^+-]+", text)
    if not terms or "".join(terms) != text:
        return None
    for term in terms:
        sign = -1 if term.startswith("-") else 1
        body = term[1:] if term[:1] in "+-" else term
        try:
            total += sign * Fraction(body)
            continue
        except (ValueError, ZeroDivisionError):
            pass
        match = re.fullmatch(r"(\d*)?([A-Za-z][A-Za-z0-9_]*)", body)
        if match is None or match.group(2) not in values:
            return None
        coeff = int(match.group(1) or "1")
        total += sign * coeff * values[match.group(2)]
    return total


def _parametric_kpoint_template(
    parent_sg: int,
    kp,
) -> tuple[list[str], list[str]] | None:
    override = KPOINT_OFFICIAL.get(int(parent_sg), {}).get(
        str(getattr(kp, "label", "") or "").strip()
    )
    if override is not None:
        _kovalev, coords, params = override
        if params:
            return [str(c) for c in coords], [str(p) for p in params]
        return None
    params = [str(p) for p in (getattr(kp, "parameters", None) or [])]
    coords = [str(c) for c in (getattr(kp, "coordinates", None) or [])]
    if not params:
        return None
    return coords, params


def folding_k_candidates(
    parent_sg: int,
    kpoints: Sequence,
    basis: np.ndarray,
    primary: _KCandidate | None,
    nmod: int,
    *,
    child_centering: str = "P",
) -> list[_KCandidate]:
    """Parent k-points that are Gamma of the child supercell (plus nmod filter)."""
    basis_exact = rational_matrix(basis)
    parent_symbol = (hm_symbol(int(parent_sg or 1)) or "P").replace(" ", "")
    parent_centering = next(
        (char for char in parent_symbol if char.isalpha()), "P"
    ).upper()
    child_primitive = multiply(
        centering_primitive_matrix(child_centering), basis_exact,
    )
    parent_primitive = centering_primitive_matrix(parent_centering)
    relative_primitive = multiply(child_primitive, inverse(parent_primitive))
    if not is_integral(relative_primitive):
        raise ValueError(
            "Subgroup primitive translations must form a sublattice of "
            "the parent primitive translation lattice"
        )
    det = abs(int(determinant(relative_primitive)))
    if det < 1:
        raise ValueError("Subgroup primitive translation lattice is singular")
    # Compatibility phases are evaluated in the parent conventional frame;
    # the finite quotient order/Smith exponent is evaluated relative to the
    # parent primitive lattice.  Keeping these two roles separate is essential
    # for centred parents and fractional conventional basis rows.
    m = np.asarray(
        [[float(value) for value in row] for row in child_primitive],
        dtype=float,
    )
    found: list[_KCandidate] = []
    seen: set[tuple[str, ...]] = set()

    try:
        raw_star_rotations = [
            np.asarray(op.rotation_matrix, dtype=float)
            for op in SpaceGroup.from_int_number(int(parent_sg)).symmetry_ops
        ]
    except (ValueError, KeyError, TypeError):
        raw_star_rotations = [np.eye(3)]
    star_rotations: list[np.ndarray] = []
    rotation_keys: set[tuple[int, ...]] = set()
    for rotation in raw_star_rotations:
        key = tuple(np.rint(rotation).astype(int).ravel())
        if key not in rotation_keys:
            rotation_keys.add(key)
            star_rotations.append(rotation)

    centering_vectors = _parent_centering_vectors(parent_sg)

    def _canonical_arm(kvec: np.ndarray) -> np.ndarray:
        """Reduce by the *parent* reciprocal lattice, including centering.

        Component-wise modulo one is only valid for a primitive conventional
        cell.  In an I lattice, for example, (1,1,1) is not Gamma because the
        integer shift has odd parity.
        """
        k = np.asarray(kvec, dtype=float)
        # Enumerate every parent-reciprocal translation near ``-k``.  A
        # component-wise ``rint`` is not a canonical reduction on a Brillouin-
        # zone boundary: at half-integers several equally short integer shifts
        # exist, and a centering condition can make more than one of them
        # admissible.  In an I lattice, for example,
        # (1/2,1/2,-1/2) and (-1/2,1/2,1/2) differ by the allowed reciprocal
        # vector (1,0,-1).  Missing that tie double-counts one P-star class and
        # can prematurely consume the finite-quotient stopping count.
        shift_ranges = [
            range(int(np.floor(-value)) - 1, int(np.ceil(-value)) + 2)
            for value in k
        ]
        candidates = [
            k + np.asarray(shift, dtype=float)
            for shift in product(*shift_ranges)
            if all(
                np.isclose(
                    float(np.dot(shift, centering)),
                    round(float(np.dot(shift, centering))),
                    atol=1e-8,
                )
                for centering in centering_vectors
            )
        ]
        if not candidates:
            return k
        return min(
            candidates,
            key=lambda value: (
                round(float(np.dot(value, value)), 12),
                tuple(round(abs(float(x)), 12) for x in value),
                tuple(round(float(x), 12) for x in value),
            ),
        )

    def _star_arms(kvec: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
        """Unique ``k``-star arms with a parent operation mapping each arm.

        ``SymmOp`` rotations act on direct fractional columns as ``R x``;
        reciprocal columns therefore transform as ``R^-T k``.  Opposite arms
        form one real sine/cosine space and are represented only once.
        """
        arms: list[tuple[np.ndarray, np.ndarray]] = []
        arm_keys: set[tuple[str, ...]] = set()
        for rot in star_rotations:
            try:
                arm = np.linalg.inv(rot).T @ np.asarray(kvec, dtype=float)
            except np.linalg.LinAlgError:
                continue
            arm = _canonical_arm(arm)
            opposite = _canonical_arm(-arm)
            key_a = tuple(format_k_token(x) for x in arm)
            key_b = tuple(format_k_token(x) for x in opposite)
            key = min(key_a, key_b)
            if key in arm_keys:
                continue
            arm_keys.add(key)
            arms.append((arm, rot))
        return arms or [(np.asarray(kvec, dtype=float), np.eye(3))]

    # Freeze the finite reciprocal quotient represented by this child cell.
    # This gives an exact stopping condition for the CDML stratum walk below,
    # so a 6x6 cell does not needlessly traverse every GP parameter triple
    # after all 20 real +/- k classes have already been classified.
    quotient_keys: set[tuple[str, ...]] = set()
    # Include the positive boundary representative.  It is redundant for a
    # primitive parent, but not for a centered conventional cell: an integer
    # conventional reciprocal shift can be outside the parent reciprocal
    # lattice (I: odd h+k+l).  For I4/mmm 1x1x6 this extra coset is the M
    # class; omitting it made the early-stop count consume LD 5/6.
    for numerators in product(range(det + 1), repeat=3):
        trial = np.asarray(numerators, dtype=float) / float(det)
        if not _k_compatible_with_supercell(trial, m):
            continue
        arm = _canonical_arm(trial)
        opposite = _canonical_arm(-arm)
        key_a = tuple(format_k_token(x) for x in arm)
        key_b = tuple(format_k_token(x) for x in opposite)
        quotient_keys.add(min(key_a, key_b))
    target_real_k_count = len(quotient_keys)

    def _add(cand: _KCandidate) -> None:
        if target_real_k_count and len(seen) >= target_real_k_count:
            return
        for arm, rot in _star_arms(cand.coords):
            if not _k_compatible_with_supercell(arm, m):
                continue
            arm_key = tuple(format_k_token(x) for x in arm)
            opposite_key = tuple(format_k_token(x) for x in _canonical_arm(-arm))
            # A physical +/- arm belongs to one canonical CDML stratum.  The
            # primary is inserted first, followed by special points, lines,
            # planes and GP, so boundary aliases never trigger duplicate
            # smodes calculations or duplicate mode spaces.
            key = min(arm_key, opposite_key)
            if key in seen:
                continue
            seen.add(key)
            found.append(replace(
                cand,
                active_coords=np.asarray(arm, dtype=float),
                arm_rotation=np.asarray(rot, dtype=float),
            ))
            if target_real_k_count and len(seen) >= target_real_k_count:
                break

    param_vals = _param_value_candidates(relative_primitive)
    extra = {Fraction(n, det) for n in range(1, det) if det % n == 0 or n <= det}
    param_vals = sorted(
        {Fraction(0)} | set(param_vals) | extra,
        key=lambda f: (f.denominator, f.numerator),
    )

    # The selected primary star owns its physical grid points even when a
    # lower-dimensional line template reaches the same coordinates.
    if primary is not None:
        _add(primary)

    def _kp_rank(kp) -> tuple[int, int, str]:
        if bool(getattr(kp, "is_special", False)):
            return (0, 0, str(getattr(kp, "label", "")))
        tmpl = _parametric_kpoint_template(parent_sg, kp)
        return (1, len(tmpl[1]) if tmpl else 99, str(getattr(kp, "label", "")))

    for kp in sorted(kpoints, key=_kp_rank):
        if target_real_k_count and len(seen) >= target_real_k_count:
            break
        label = str(getattr(kp, "label", "") or "").strip()
        if not label:
            continue
        coords = list(getattr(kp, "coordinates", None) or [])
        if getattr(kp, "is_special", False):
            kvec = _numeric_k(parent_sg, label, coords, None)
            if kvec is None:
                continue
            _add(_KCandidate(label, kvec, None, tuple(format_k_token(x) for x in kvec)))
            continue
        tmpl = _parametric_kpoint_template(parent_sg, kp)
        if tmpl is None:
            continue
        tmpl_coords, pnames = tmpl
        for pvals in product(param_vals, repeat=len(pnames)):
            if target_real_k_count and len(seen) >= target_real_k_count:
                break
            values = [
                _eval_affine_k_component(comp, pnames, pvals)
                for comp in tmpl_coords
            ]
            if any(value is None for value in values):
                continue
            kvec = np.asarray([float(value) for value in values], dtype=float)
            params = [
                f"{value.numerator}/{value.denominator}"
                if value.denominator != 1 else str(value.numerator)
                for value in pvals
            ]
            _add(_KCandidate(
                label, kvec, params,
                tuple(format_k_token(x) for x in kvec),
            ))

    if nmod <= 0 or primary is None:
        return found
    fundamentals = [primary.coords]
    return [
        cand for cand in found
        if is_harmonic_of(cand.coords, fundamentals)
    ]


def _child_ops(
    space_group_number: int,
    origin: Sequence[float] | None = None,
    structure: Structure | None = None,
) -> list[SymmOp]:
    """Return target-group operations in the actual child-cell origin.

    ``SpaceGroup.from_int_number`` uses its own standard origin/setting, while
    iso reports an embedding origin relative to the parent.  Depending on the
    Hall/origin choice already selected by pymatgen, applying that shift again
    can move most operations off the parent sites (notably I4/mmm M4- ->
    P4_2/nmc).  Try both conjugation conventions plus the standard origin and,
    when the child structure is available, select the first complete action.
    This is a geometric setting resolution, not a compound/IR answer table.
    """
    ops = list(_standard_child_ops(int(space_group_number)))
    origin_v = np.asarray(origin if origin is not None else [0.0, 0.0, 0.0],
                          dtype=float).reshape(3)
    if np.allclose(origin_v, 0.0) or structure is None:
        if np.allclose(origin_v, 0.0):
            return ops
        shifted: list[SymmOp] = []
        for op in ops:
            rot = np.asarray(op.rotation_matrix, dtype=float)
            trans = np.asarray(op.translation_vector, dtype=float)
            new_t = trans + origin_v - rot @ origin_v
            shifted.append(SymmOp.from_rotation_and_translation(rot, new_t))
        return shifted

    # If r_parent = r_child B + origin, conjugating a standard child
    # operation back onto the unshifted supercell gives t' = t + o - R o.
    # ``+origin`` is therefore the physical convention.  Keep the opposite
    # sign and zero as setting fallbacks, but do not let an accidentally
    # site-compatible opposite embedding win merely because it was tried
    # first (that projected out valid root modes for several nonzero origins).
    trial_origins = (origin_v, -origin_v, np.zeros(3, dtype=float))
    best_ops = ops
    best_score = -1
    probe = np.zeros((len(structure), 3), dtype=float)
    for trial_origin in trial_origins:
        shifted = []
        for op in ops:
            rot = np.asarray(op.rotation_matrix, dtype=float)
            trans = np.asarray(op.translation_vector, dtype=float)
            new_t = trans + trial_origin - rot @ trial_origin
            shifted.append(SymmOp.from_rotation_and_translation(rot, new_t))
        score = sum(
            _transform_displacement(probe, structure, op) is not None
            for op in shifted
        )
        if score > best_score:
            best_ops, best_score = shifted, score
        if score == len(ops):
            return shifted
    return best_ops


@lru_cache(maxsize=230)
def _standard_child_ops(space_group_number: int) -> tuple[SymmOp, ...]:
    """Operations in the same ITA setting used by Distortion exports.

    A space-group number alone can have several Hall settings.  ISODISTORT's
    fixed preferences select origin choice 2 when it exists, whereas
    ``pymatgen.SpaceGroup.from_int_number`` may choose another Hall entry.
    Build the choice-2 operation set from spglib's database so mode projection
    and CIF/ISOVIZ export use the same setting.  This is space-group metadata,
    not a per-IR answer table.
    """
    try:
        from spglib import (  # noqa: PLC0415 - optional runtime dependency
            get_spacegroup_type,
            get_symmetry_from_database,
        )

        first_hall = 0
        choice2_hall = 0
        for hall_number in range(1, 531):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                info = get_spacegroup_type(hall_number)
                info_number = int(info.number)
                info_choice = str(getattr(info, "choice", "") or "")
            if info_number != int(space_group_number):
                continue
            if not first_hall:
                first_hall = hall_number
            if info_choice == "2":
                choice2_hall = hall_number
                break
        hall = choice2_hall or first_hall
        if hall:
            data = get_symmetry_from_database(hall)
            operations: list[SymmOp] = []
            seen: set[tuple[float, ...]] = set()
            for rotation, translation in zip(
                data["rotations"], data["translations"], strict=False,
            ):
                rot = np.asarray(rotation, dtype=float)
                trans = np.mod(np.asarray(translation, dtype=float), 1.0)
                trans[np.isclose(trans, 1.0, atol=1e-10)] = 0.0
                key = tuple(np.round(np.concatenate((rot.ravel(), trans)), 10))
                if key in seen:
                    continue
                seen.add(key)
                operations.append(
                    SymmOp.from_rotation_and_translation(rot, trans)
                )
            if operations:
                return tuple(operations)
    except (ImportError, KeyError, TypeError, ValueError):
        pass
    try:
        return tuple(
            SpaceGroup.from_int_number(int(space_group_number)).symmetry_ops
        )
    except (ValueError, KeyError, TypeError):
        return (
            SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3)),
        )


def _child_origin_candidates(
    origin_parent: Sequence[float],
    basis: np.ndarray,
    parent_translations: np.ndarray,
) -> list[np.ndarray]:
    """Lift a parent-cell origin through the child translation quotient.

    ``iso`` can print an embedding origin outside the parent reference cell
    (for example ``9/4`` along an axis of a sixfold cell).  ``SubgroupInfo``
    intentionally stores coordinates modulo the *parent* cell, which loses
    that integer lift even though it is not an integer translation of the
    child cell.  The atoms already enumerated by :func:`build_supercell`
    provide a complete set of parent-translation cosets, so try those lifts
    rather than guessing a bounded integer search.
    """
    try:
        inv_basis = np.linalg.inv(np.asarray(basis, dtype=float))
    except np.linalg.LinAlgError:
        return [np.zeros(3, dtype=float)]
    origin = np.asarray(origin_parent, dtype=float).reshape(3)
    translations = np.asarray(parent_translations, dtype=float).reshape((-1, 3))
    if translations.size == 0:
        translations = np.zeros((1, 3), dtype=float)
    candidates: list[np.ndarray] = []
    seen: set[tuple[float, ...]] = set()
    for translation in np.vstack((np.zeros((1, 3)), translations)):
        child = (origin + translation) @ inv_basis
        child -= np.floor(child)
        child[np.isclose(child, 1.0, atol=1e-10)] = 0.0
        key = tuple(np.round(child, 10))
        if key in seen:
            continue
        seen.add(key)
        candidates.append(child)
    return candidates


def _element_symbol(species: str) -> str:
    m = re.match(r"([A-Z][a-z]?)", str(species or ""))
    return m.group(1) if m else str(species or "")


def _nearest_index(
    frac: np.ndarray,
    coords: np.ndarray,
    species: str,
    species_list: Sequence[str],
    lattice_matrix: np.ndarray | None = None,
) -> int:
    """Index of the nearest periodic same-species site.

    Fractional-coordinate Euclidean distance is not a physical metric for a
    non-cubic cell.  When a lattice is supplied, rank candidates by the
    minimum-image Cartesian distance ``||(x-y) A||`` in Å, where the rows
    of ``A`` are the direct-lattice vectors used by pymatgen.
    """
    best_i, best = -1, 1e9
    lattice = None if lattice_matrix is None else Lattice(lattice_matrix)
    for i, other in enumerate(coords):
        if species_list[i] != species:
            continue
        if lattice is None:
            d = frac - other
            d -= np.round(d)
            dist = float(np.linalg.norm(d))
        else:
            # Component-wise fractional wrapping is not the nearest-image
            # problem in a skew cell.  Pymatgen solves the direct-lattice
            # metric problem and returns the physical distance in Å.
            dist, _image = lattice.get_distance_and_image(frac, other)
            dist = float(dist)
        if dist < best:
            best, best_i = dist, i
    return best_i


def mode_is_invariant(
    disp: np.ndarray,
    structure: Structure,
    ops: Sequence[SymmOp],
    *,
    atol: float = 2e-3,
) -> bool:
    """True if ``disp`` transforms as the identity of ``ops``.

    ``atol`` is the dimensionless tolerance for normalized displacement
    components.  Atomic site matching is a separate physical question and
    uses the configured Cartesian symmetry tolerance in Å.
    """
    actions = _displacement_actions(structure, ops)
    if actions is None:
        return False
    for targets, rotation in actions:
        transformed = _apply_displacement_action(disp, targets, rotation)
        work = np.asarray(disp, dtype=float)
        nrm = float(np.linalg.norm(work))
        if nrm < 1e-12:
            return False
        if float(np.max(np.abs(transformed / nrm - work / nrm))) > atol:
            return False
    return True


def _site_bijection(
    candidates: Sequence[Sequence[tuple[int, float]]],
) -> np.ndarray | None:
    """Return a deterministic perfect source-to-target matching, if one exists."""
    n_sites = len(candidates)
    target_owner = np.full(n_sites, -1, dtype=int)
    source_target = np.full(n_sites, -1, dtype=int)

    def augment(source: int, seen_targets: set[int]) -> bool:
        for target, _distance in candidates[source]:
            if target in seen_targets:
                continue
            seen_targets.add(target)
            previous = int(target_owner[target])
            if previous < 0 or augment(previous, seen_targets):
                target_owner[target] = source
                source_target[source] = target
                return True
        return False

    source_order = sorted(
        range(n_sites), key=lambda source: (len(candidates[source]), source)
    )
    if any(not augment(source, set()) for source in source_order):
        return None
    if np.any(source_target < 0) or len(set(source_target.tolist())) != n_sites:
        raise RuntimeError("validated site matching did not produce a bijection")
    return source_target


def _periodic_species_bijection(
    source_fractional: np.ndarray,
    target_fractional: np.ndarray,
    source_species: Sequence[str],
    target_species: Sequence[str],
    lattice: Lattice,
    tolerance_angstrom: float,
) -> np.ndarray | None:
    """Match periodic sites by species with the true direct-lattice metric."""
    sources = np.asarray(source_fractional, dtype=float)
    targets = np.asarray(target_fractional, dtype=float)
    if (
        sources.ndim != 2
        or targets.ndim != 2
        or sources.shape[1:] != (3,)
        or targets.shape[1:] != (3,)
        or len(sources) != len(source_species)
        or len(targets) != len(target_species)
        or len(sources) != len(targets)
    ):
        raise ValueError("periodic site matching requires equal N x 3 site arrays")
    tolerance = float(tolerance_angstrom)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("site-matching tolerance must be finite and positive in Å")

    distances = np.asarray(
        lattice.get_all_distances(sources, targets), dtype=float
    )
    candidates: list[list[tuple[int, float]]] = []
    for source in range(len(sources)):
        options = sorted(
            (
                (target, float(distances[source, target]))
                for target in range(len(targets))
                if target_species[target] == source_species[source]
                and distances[source, target] <= tolerance
            ),
            key=lambda item: (item[1], item[0]),
        )
        if not options:
            return None
        candidates.append(options)
    return _site_bijection(candidates)


def _displacement_actions(
    structure: Structure,
    ops: Sequence[SymmOp],
    *,
    atol: float | None = None,
) -> list[tuple[np.ndarray, np.ndarray]] | None:
    """Precompute each operation's site permutation and vector rotation.

    Site matching depends only on the structure and operation, not on the mode
    vector.  Complete-mode searches can project hundreds of seeds on a large
    child cell, so repeating the O(N²) nearest-site search for every seed was
    the dominant cost.  This plan is mathematically identical to applying each
    ``SymmOp`` directly and is reused by all Reynolds projections.  ``atol``
    is a Cartesian site-matching tolerance in Å; by default it comes from
    ``defaults.symmetry_cartesian_tolerance_angstrom``.  This avoids treating
    the same fractional residual as equivalent in physically different cells.
    """
    coords = np.asarray(structure.frac_coords, dtype=float)
    site_tolerance = (
        get_config().symmetry_cartesian_tolerance_angstrom
        if atol is None else float(atol)
    )
    if not np.isfinite(site_tolerance) or site_tolerance <= 0.0:
        raise ValueError("site-matching tolerance must be finite and positive")
    species = [site.species_string for site in structure]
    actions: list[tuple[np.ndarray, np.ndarray]] = []
    for op in ops:
        rotation = np.asarray(op.rotation_matrix, dtype=float)
        translation = np.asarray(op.translation_vector, dtype=float)
        new_pos = (coords @ rotation.T + translation) % 1.0
        # A local nearest-neighbour choice is not a proof of a site
        # permutation: two sources can prefer the same target even though a
        # complete species-preserving assignment exists.  Solve the finite
        # bipartite matching problem exactly (for this tolerance graph), with
        # constrained sources first and deterministic distance/index order.
        targets = _periodic_species_bijection(
            new_pos,
            coords,
            species,
            species,
            structure.lattice,
            site_tolerance,
        )
        if targets is None:
            return None
        actions.append((targets, rotation))
    return actions


def _apply_displacement_action(
    disp: np.ndarray,
    targets: np.ndarray,
    rotation: np.ndarray,
) -> np.ndarray:
    """Apply a precomputed child-group action to one displacement field."""
    work = np.asarray(disp, dtype=float)
    mapped = np.zeros_like(work)
    mapped[targets] = work @ rotation.T
    return mapped


def _transform_displacement(
    disp: np.ndarray,
    structure: Structure,
    op: SymmOp,
    *,
    atol: float | None = None,
) -> np.ndarray | None:
    """Apply one child operation; ``atol`` is the site tolerance in Å."""
    actions = _displacement_actions(structure, [op], atol=atol)
    if actions is None:
        return None
    targets, rotation = actions[0]
    return _apply_displacement_action(disp, targets, rotation)


def _project_invariant(
    disp: np.ndarray,
    structure: Structure,
    ops: Sequence[SymmOp],
    *,
    atol: float | None = None,
    actions: Sequence[tuple[np.ndarray, np.ndarray]] | None = None,
) -> np.ndarray | None:
    """Reynolds-project a vector field onto the child identity irrep."""
    projected = np.zeros_like(np.asarray(disp, dtype=float))
    if not ops:
        return np.asarray(disp, dtype=float)
    if actions is None:
        actions = _displacement_actions(structure, ops, atol=atol)
    if actions is None:
        return None
    for targets, rotation in actions:
        projected += _apply_displacement_action(disp, targets, rotation)
    projected /= float(len(ops))
    if float(np.linalg.norm(projected)) < 1e-10:
        return None
    return projected


def _map_smodes_block_to_supercell(
    block: SmodesModeBlock,
    supercell: Structure,
) -> np.ndarray | None:
    """Tile smodes cartesian modes onto the child supercell.

    ``smodes`` writes each IR in *its* commensurate superlattice (often
    smaller than the isotropy-subgroup cell).  Atoms of the child cell are
    matched to the smodes pattern by reducing through that superlattice so
    a 1×1×2 Z-mode repeats correctly inside a 1×1×6 lock-in cell.
    """
    dense = block.dense_displacements()
    if dense.size == 0:
        return None
    A = np.asarray(block.superlattice_cart, dtype=float)
    inv_a: np.ndarray | None = None
    if A.shape == (3, 3) and abs(float(np.linalg.det(A))) > 1e-8:
        try:
            inv_a = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            inv_a = None
    sc_cart = np.asarray(supercell.cart_coords, dtype=float)
    inv_child = np.linalg.inv(np.asarray(supercell.lattice.matrix, dtype=float))
    sc_elem = [_element_symbol(site.species_string) for site in supercell]
    out = np.zeros((len(supercell), 3), dtype=float)
    used = 0
    for j, cart in enumerate(sc_cart):
        best = 1e9
        best_disp = None
        for atom, cart_disp in zip(block.atoms, dense, strict=False):
            if _element_symbol(atom.species) != sc_elem[j]:
                continue
            if inv_a is not None:
                frac = (cart - atom.cart) @ inv_a
                frac -= np.round(frac)
                dist = float(np.linalg.norm(frac @ A))
            else:
                dist = float(np.linalg.norm(cart - atom.cart))
            if dist < best:
                best = dist
                best_disp = cart_disp
        if best_disp is None or best > 1.35:
            continue
        out[j] = best_disp @ inv_child
        used += 1
    if used == 0 or float(np.max(np.abs(out))) < 1e-12:
        return None
    return out


def _parent_centering_vectors(space_group_number: int) -> list[np.ndarray]:
    """Conventional-cell centering translations of the parent lattice."""
    sym = (hm_symbol(int(space_group_number or 1)) or "P").replace(" ", "")
    letter = next((ch for ch in sym if ch.isalpha()), "P").upper()
    extras = {
        "I": [(0.5, 0.5, 0.5)],
        "A": [(0.0, 0.5, 0.5)],
        "B": [(0.5, 0.0, 0.5)],
        "C": [(0.5, 0.5, 0.0)],
        "F": [(0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0)],
        "R": [(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)],
    }
    out = [np.zeros(3, dtype=float)]
    for vec in extras.get(letter, ()):
        out.append(np.asarray(vec, dtype=float))
    return out


def _quadrature_supercell_shift(
    kvec: np.ndarray,
    basis: np.ndarray,
    parent_sg: int = 1,
) -> np.ndarray | None:
    """Supercell fractional shift that phases a non-Gamma Bloch mode by ~1/4.

    Integer conventional translations are not enough at zone-boundary k
    (LD g=1/2: one parent cell only flips sign).  Parent centering
    translations (I: ½½½) supply the missing quarter-period.
    """
    k = np.asarray(kvec, dtype=float)
    if np.allclose(k, 0.0, atol=1e-8):
        return None
    try:
        inv_b = np.linalg.inv(np.asarray(basis, dtype=float))
    except np.linalg.LinAlgError:
        return None
    best = None
    best_dist = 99.0
    # Search at least one parent conventional cell past the child supercell,
    # so a primitive lattice with a long lock-in (e.g. 1×1×24) can still
    # hit k·T ≈ 1/4.  A fixed 0..5 window misses that translation.
    span = int(np.ceil(float(np.max(np.abs(np.asarray(basis, dtype=float)))))) + 1
    span = min(max(span, 2), 24)
    for extra in _parent_centering_vectors(parent_sg):
        for nx, ny, nz in product(range(0, span), repeat=3):
            t_parent = extra + np.array([nx, ny, nz], dtype=float)
            if np.allclose(t_parent, 0.0):
                continue
            t_sc = t_parent @ inv_b
            if np.allclose(t_sc, np.round(t_sc), atol=1e-6):
                continue
            phase = float(np.dot(k, t_parent))
            phase -= round(phase)
            away = abs(phase)
            if away < 0.08 or away > 0.92:
                continue
            dist = abs(away - 0.25)
            if dist < best_dist:
                best_dist = dist
                best = t_sc
    return best


def _shift_displacement(
    disp: np.ndarray,
    structure: Structure,
    frac_shift: np.ndarray,
) -> np.ndarray:
    """d'(r) = d(r − t) with t in supercell fractional coordinates."""
    coords = np.asarray(structure.frac_coords, dtype=float)
    species = [site.species_string for site in structure]
    shift = np.asarray(frac_shift, dtype=float).reshape(3)
    out = np.zeros_like(disp)
    for i in range(len(structure)):
        src = (coords[i] - shift) % 1.0
        j = _nearest_index(
            src, coords, species[i], species, structure.lattice.matrix,
        )
        if j < 0:
            continue
        out[i] = disp[j]
    return out


def _quadrature_partner(
    disp: np.ndarray,
    structure: Structure,
    kvec: np.ndarray,
    basis: np.ndarray,
    parent_sg: int = 1,
) -> np.ndarray | None:
    """Second real component of a 2-D IR from a parent-cell phase shift."""
    shift = _quadrature_supercell_shift(kvec, basis, parent_sg=parent_sg)
    if shift is None:
        return None
    moved = _shift_displacement(disp, structure, shift)
    v0 = np.asarray(disp, dtype=float).ravel()
    v1 = np.asarray(moved, dtype=float).ravel()
    n0 = float(np.linalg.norm(v0))
    if n0 < 1e-12:
        return None
    u0 = v0 / n0
    v1 = v1 - float(np.dot(v1, u0)) * u0
    n1 = float(np.linalg.norm(v1))
    if n1 < 0.05 * n0:
        return None
    out = v1.reshape(disp.shape)
    nrm = float(np.max(np.abs(out)))
    if nrm < 1e-12:
        return None
    return out / nrm


def _same_parent_real_k(
    first: np.ndarray,
    second: np.ndarray,
    parent_sg: int,
) -> bool:
    """Whether two arms represent the same real ±k space."""
    centering = _parent_centering_vectors(parent_sg)
    for sign in (1.0, -1.0):
        delta = np.asarray(first, dtype=float) - sign * np.asarray(second, dtype=float)
        nearest = np.rint(delta)
        if not np.allclose(delta, nearest, atol=1e-7):
            continue
        if all(
            np.isclose(
                float(np.dot(nearest, translation)),
                round(float(np.dot(nearest, translation))),
                atol=1e-7,
            )
            for translation in centering
        ):
            return True
    return False


def _is_self_conjugate_parent_k(kvec: np.ndarray, parent_sg: int) -> bool:
    """Whether ``k`` and ``-k`` differ by a parent reciprocal vector."""
    doubled = 2.0 * np.asarray(kvec, dtype=float)
    nearest = np.rint(doubled)
    if not np.allclose(doubled, nearest, atol=1e-7):
        return False
    return all(
        np.isclose(
            float(np.dot(nearest, translation)),
            round(float(np.dot(nearest, translation))),
            atol=1e-7,
        )
        for translation in _parent_centering_vectors(parent_sg)
    )


def _parent_little_group_orbit(
    disp: np.ndarray,
    structure: Structure,
    basis: np.ndarray,
    parent_sg: int,
    kvec: np.ndarray,
    dimension: int,
    actions: Sequence[tuple[np.ndarray, np.ndarray]] | None = None,
) -> list[np.ndarray]:
    """Complete a degenerate smodes seed using actual parent operations.

    ``smodes`` emits one representative vector for each copy of a degenerate
    irrep.  Rotating its atom table geometrically is insufficient for
    nonsymmorphic operations because the translational part controls the site
    permutation.  Transform each full parent ``SymmOp`` into child fractional
    axes, apply it to the mapped supercell field, and retain an independent
    Cartesian orbit of the requested irrep dimension.
    """
    target = max(int(dimension), 1)
    seed = np.asarray(disp, dtype=float)
    child_lattice = np.asarray(structure.lattice.matrix, dtype=float)
    orbit = [seed]
    first_cart = (seed @ child_lattice).ravel()
    first_norm = float(np.linalg.norm(first_cart))
    if first_norm < 1e-12 or target <= 1:
        return orbit
    cart_basis = [first_cart / first_norm]
    if actions is None:
        actions = _parent_little_group_actions(
            structure, basis, parent_sg, kvec,
        )
    for targets, rotation in actions:
        transformed = _apply_displacement_action(seed, targets, rotation)
        cart = (transformed @ child_lattice).ravel()
        residual = cart.copy()
        for previous in cart_basis:
            residual -= float(np.dot(residual, previous)) * previous
        norm = float(np.linalg.norm(residual))
        if norm <= 1e-8 * max(float(np.linalg.norm(cart)), 1.0):
            continue
        cart_basis.append(residual / norm)
        orbit.append(transformed)
        if len(orbit) >= target:
            break
    return orbit


@lru_cache(maxsize=32)
def _parent_space_group_ops(parent_sg: int) -> tuple[SymmOp, ...]:
    try:
        return tuple(SpaceGroup.from_int_number(int(parent_sg)).symmetry_ops)
    except (ValueError, KeyError, TypeError):
        return (
            SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3)),
        )


def _parent_little_group_actions(
    structure: Structure,
    basis: np.ndarray,
    parent_sg: int,
    kvec: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Site permutations for the parent little group in child coordinates."""
    basis_t = np.asarray(basis, dtype=float).T
    try:
        inv_basis_t = np.linalg.inv(basis_t)
    except np.linalg.LinAlgError:
        return []
    result: list[tuple[np.ndarray, np.ndarray]] = []
    for parent_op in _parent_space_group_ops(int(parent_sg)):
        rotation = np.asarray(parent_op.rotation_matrix, dtype=float)
        try:
            transformed_k = np.linalg.inv(rotation).T @ np.asarray(kvec, dtype=float)
        except np.linalg.LinAlgError:
            continue
        if not _same_parent_real_k(
            transformed_k, np.asarray(kvec, dtype=float), parent_sg,
        ):
            continue
        child_rotation = inv_basis_t @ rotation @ basis_t
        child_translation = (
            inv_basis_t @ np.asarray(parent_op.translation_vector, dtype=float)
        )
        child_op = SymmOp.from_rotation_and_translation(
            child_rotation, child_translation,
        )
        action = _displacement_actions(structure, [child_op])
        if action is not None:
            result.append(action[0])
    return result


def _parent_map(parent: Structure, supercell: Structure, basis: np.ndarray
                ) -> tuple[np.ndarray, np.ndarray]:
    """For each supercell atom: parent index and integer translation t."""
    inv_b = np.linalg.inv(np.asarray(basis, dtype=float))
    parent_frac_in_parent = (
        np.asarray(supercell.frac_coords, dtype=float) @ np.asarray(basis, dtype=float)
    )
    pc = np.asarray(parent.frac_coords, dtype=float)
    idx = np.zeros(len(supercell), dtype=int)
    trans = np.zeros((len(supercell), 3), dtype=float)
    for j in range(len(supercell)):
        best_i, best = 0, 1e9
        for i in range(len(parent)):
            if supercell[j].species_string != parent[i].species_string:
                continue
            dist = float(parent.lattice.get_distance_and_image(
                pc[i], parent_frac_in_parent[j],
            )[0])
            if dist < best:
                best, best_i = dist, i
        idx[j] = best_i
        t = parent_frac_in_parent[j] - pc[best_i]
        trans[j] = np.round(t)
    _ = inv_b
    return idx, trans


def _site_point_groups(parent: Structure) -> dict[int, str]:
    """Return spglib's oriented site-symmetry symbols unchanged.

    Dots carry crystallographic orientation information (for example ``m..``
    versus ``.m.``), so normalising them away makes the result unsuitable for
    representation labels.  This helper is retained for diagnostic callers;
    labels themselves come from ISO's oriented microscopic tables.
    """
    try:
        cfg = get_config()
        sga = SpacegroupAnalyzer(
            parent,
            symprec=cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance=cfg.symmetry_angle_tolerance_degrees,
        )
        data = sga.get_symmetry_dataset()
        symbols = list(data.get("site_symmetry_symbols") or [])
    except Exception:  # noqa: BLE001
        return {}
    out: dict[int, str] = {}
    for i, sym in enumerate(symbols):
        out[i] = str(sym)
    return out


def _site_transport_rotations(
    parent: Structure,
    wyckoff_sites: Sequence[dict],
) -> dict[int, np.ndarray]:
    """Return parent operations carrying each orbit representative to a site.

    Site-symmetry irreps are defined at one representative point.  A vector
    printed on an equivalent point must therefore be pulled back through the
    corresponding coset transporter before its direction is classified.  The
    old largest-global-axis shortcut confused, for example, the orthorhombic
    ``B2u`` and ``B3u`` axes after an orbit permutation.
    """
    cfg = get_config()
    try:
        operations = SpacegroupAnalyzer(
            parent,
            symprec=cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance=cfg.symmetry_angle_tolerance_degrees,
        ).get_symmetry_operations(cartesian=False)
    except (ValueError, TypeError, np.linalg.LinAlgError, AttributeError):
        return {}

    rotations: dict[int, np.ndarray] = {}
    for site in wyckoff_sites:
        try:
            representative = int(site["representative_index"])
        except (KeyError, TypeError, ValueError):
            continue
        if representative < 0 or representative >= len(parent):
            continue
        rep_coord = np.asarray(parent[representative].frac_coords, dtype=float)
        indices = [
            int(value)
            for value in (site.get("equivalent_indices") or [representative])
        ]
        for target in indices:
            if target < 0 or target >= len(parent):
                continue
            target_coord = np.asarray(parent[target].frac_coords, dtype=float)
            best_rotation: np.ndarray | None = None
            best_distance = float("inf")
            for operation in operations:
                mapped = operation.operate(rep_coord)
                distance = float(
                    parent.lattice.get_distance_and_image(mapped, target_coord)[0]
                )
                if distance < best_distance:
                    best_distance = distance
                    best_rotation = np.asarray(operation.rotation_matrix, dtype=float)
            if (
                best_rotation is not None
                and best_distance <= cfg.symmetry_cartesian_tolerance_angstrom
            ):
                rotations[target] = best_rotation
    return rotations


def _site_irrep_for_disp(site_pg: str, disp: np.ndarray) -> str:
    """Compatibility shim for callers that lack an audited mode identity.

    One arbitrary displacement vector need not span an irreducible subspace,
    especially after a degenerate basis rotation.  A label therefore cannot
    be inferred from its largest Cartesian component.  Returning the empty
    unresolved value keeps staged callers fail-closed until they consume the
    ISO microscopic :class:`ModeIdentity` directly.
    """

    _ = site_pg, disp
    return ""


def _looks_numeric(token: str) -> bool:
    try:
        parse_k_token(token)
        return True
    except (ValueError, ZeroDivisionError):
        return False


def _harmonic_order(kvec: np.ndarray, primary: np.ndarray | None) -> int | None:
    """Smallest n ≥ 0 such that k ≈ n q + G (positive harmonics, then Gamma)."""
    if primary is None:
        return None
    k = np.asarray(kvec, dtype=float)
    q = np.asarray(primary, dtype=float)
    if np.allclose(k, 0.0, atol=1e-5):
        return 0
    qn = float(np.linalg.norm(q))
    if qn < 1e-8:
        return None
    for n in range(1, 13):
        residual = k - n * q
        residual -= np.round(residual)
        if np.allclose(residual, 0.0, atol=1e-5):
            return n
    return None


def _k_family_rank(k_joined: str, primary: _KCandidate | None) -> tuple:
    """Sort: primary |n|=1, Gamma, |n|=2,3,…, then other folding k (M,…)."""
    try:
        kvec = np.array([parse_k_token(p) for p in k_joined.split(",")], dtype=float)
    except (ValueError, ZeroDivisionError):
        return (50, k_joined)
    if primary is None:
        if np.allclose(kvec, 0.0):
            return (1, k_joined)
        return (40, k_joined)
    order = _harmonic_order(kvec, primary.coords)
    if order == 1:
        return (0, k_joined)
    if order == 0:
        return (1, k_joined)
    if order is not None:
        return (1 + order, k_joined)
    return (40, k_joined)


def _irrep_sort_key(irrep: str) -> tuple:
    m = re.match(r"^([A-Za-z]+)(\d*)([+-]?)$", irrep or "")
    if not m:
        return (irrep or "", 0, "")
    return (m.group(1), int(m.group(2) or 0), m.group(3))


def _site_irrep_sort_key(site_ir: str) -> tuple:
    m = re.match(r"^([A-Za-z]+)(\d*)([ug]?)(?:_(\d+))?$", site_ir or "")
    if not m:
        return (site_ir or "", 0, "", 0)
    return (m.group(1), int(m.group(2) or 0), m.group(3) or "", int(m.group(4) or 0))


def _display_site_label(
    wyckoff_sites: Sequence[dict],
    letter: str,
    orbit_id: str = "",
) -> tuple[str, str]:
    counters: dict[str, int] = {}
    for ordinal, site in enumerate(wyckoff_sites):
        let = str(site.get("wyckoff_letter") or "")
        elem = str(site.get("species") or "X")
        counters[elem] = counters.get(elem, 0) + 1
        label = str(site.get("display_label") or f"{elem}{counters[elem]}")
        if let == letter and (
            not orbit_id or _wyckoff_orbit_id(site, ordinal) == orbit_id
        ):
            return label, elem
    return "X1", "X"


def _iso_kpoint_named(kpoints: Sequence | None, label: str):
    want = (label or "").strip()
    for kp in kpoints or []:
        if str(getattr(kp, "label", "") or "").strip() == want:
            return kp
    return None


def smodes_k_spec(
    parent_sg: int,
    cand: _KCandidate,
    kpoints: Sequence | None,
) -> tuple[str, list[str] | None]:
    """smodes input uses iso/Kovalev parameter scale, not the website ``g``.

    I4/mmm LD ``g=1/6`` is website ``(0,0,1/6)`` but smodes ``LD 1/12``
    because the iso line is ``(0,0,2a)``.
    """
    params = list(cand.params) if cand.params else None
    iso_kp = _iso_kpoint_named(kpoints, cand.label)
    if params and iso_kp is not None:
        try:
            params = official_kparams_to_iso(
                parent_sg, cand.label, params, iso_kp,
            )
        except (ValueError, ZeroDivisionError, TypeError):
            params = list(cand.params)
    return cand.label, SmodesWrapper.smodes_kparams(params)


def _collect_smodes_blocks(
    smodes: SmodesWrapper | None,
    parent: Structure,
    parent_sg: int,
    wyckoff_sites: Sequence[dict],
    specs: Sequence[tuple[_KCandidate, tuple[str, list[str] | None]]],
) -> list[SmodesModeBlock]:
    """One smodes call per k and lift each compatible star arm.

    ``smodes`` emits the representative arm only.  The website, however,
    labels a star by that representative even when an OPD activates a rotated
    arm.  Rotate the representative real-space pattern with the parent point
    operation retained by :func:`folding_k_candidates` before tiling it into
    the child cell.
    """
    if smodes is None or not hasattr(smodes, "compute_mode_blocks"):
        return []
    blocks: list[SmodesModeBlock] = []
    parent_lattice = np.asarray(parent.lattice.matrix, dtype=float)
    inv_parent_lattice = np.linalg.inv(parent_lattice)

    try:
        parent_rotations = [
            np.asarray(op.rotation_matrix, dtype=float)
            for op in SpaceGroup.from_int_number(int(parent_sg)).symmetry_ops
        ]
    except (ValueError, KeyError, TypeError):
        parent_rotations = [np.eye(3)]
    unique_parent_rotations: list[np.ndarray] = []
    rotation_keys: set[tuple[int, ...]] = set()
    for rotation in parent_rotations:
        key = tuple(np.rint(rotation).astype(int).ravel())
        if key in rotation_keys:
            continue
        rotation_keys.add(key)
        unique_parent_rotations.append(rotation)
    def _same_real_k(first: np.ndarray, second: np.ndarray) -> bool:
        """Whether two arms represent the same real ±k subspace."""
        return _same_parent_real_k(first, second, parent_sg)

    def _little_cogroup_rotations(kvec: np.ndarray) -> list[np.ndarray]:
        """Parent rotations preserving the real ±k arm.

        ``smodes`` emits one seed per degenerate little-group mode.  Rotating
        that seed through the little co-group supplies the missing Cartesian
        partner(s), most visibly the second component of GM5-/M5- and of a
        parametric Y/DT irrep.  Reynolds projection later removes duplicates
        and components forbidden by the selected child subgroup.
        """
        out: list[np.ndarray] = []
        for rotation in unique_parent_rotations:
            try:
                transformed = np.linalg.inv(rotation).T @ np.asarray(kvec, dtype=float)
            except np.linalg.LinAlgError:
                continue
            if _same_real_k(transformed, kvec):
                out.append(rotation)
        return out or [np.eye(3)]

    def _rotate_block(block: SmodesModeBlock, rot: np.ndarray | None) -> SmodesModeBlock:
        if rot is None or np.allclose(rot, np.eye(3), atol=1e-8):
            return replace(block)
        # Row-cartesian transform corresponding to x' = R x in fractional
        # column notation: c' = c L^-1 R^T L.
        cart_rot = inv_parent_lattice @ np.asarray(rot, dtype=float).T @ parent_lattice
        atoms = [
            SmodesAtom(atom.index, atom.species, np.asarray(atom.cart) @ cart_rot)
            for atom in block.atoms
        ]
        return replace(
            block,
            superlattice_cart=np.asarray(block.superlattice_cart) @ cart_rot,
            atoms=atoms,
            displacements_cart=np.asarray(block.displacements_cart) @ cart_rot,
        )

    # One process invocation can accept every k specification.  Starting WSL
    # once per folding k dominates a large (e.g. 6x6) export by minutes.
    unique_specs: dict[tuple[str, tuple[str, ...]], tuple[str, list[str] | None]] = {}
    representative: dict[tuple[str, tuple[str, ...]], _KCandidate] = {}
    for cand, spec in specs:
        cache_key = (spec[0], tuple(spec[1] or ()))
        unique_specs.setdefault(cache_key, spec)
        representative.setdefault(cache_key, cand)
    raw_blocks = smodes.compute_mode_blocks(
        parent,
        parent_sg,
        list(wyckoff_sites),
        list(unique_specs.values()),
    ) or []
    cache: dict[tuple[str, tuple[str, ...]], list[SmodesModeBlock]] = {
        key: [] for key in unique_specs
    }
    keys_by_label: dict[str, list[tuple[str, tuple[str, ...]]]] = {}
    for key in unique_specs:
        keys_by_label.setdefault(key[0], []).append(key)
    for raw_block in raw_blocks:
        choices = keys_by_label.get(raw_block.k_label, [])
        if not choices:
            continue
        if len(choices) == 1:
            cache[choices[0]].append(raw_block)
            continue
        try:
            raw_k = np.array(
                [parse_k_token(token) for token in raw_block.k_coords],
                dtype=float,
            )
        except (ValueError, ZeroDivisionError):
            raw_k = np.zeros(3, dtype=float)

        def _distance(
            key: tuple[str, tuple[str, ...]],
            raw_k_value: np.ndarray = raw_k,
        ) -> float:
            expected = np.asarray(representative[key].coords, dtype=float)
            delta = raw_k_value - expected
            delta -= np.round(delta)
            return float(np.linalg.norm(delta))

        cache[min(choices, key=_distance)].append(raw_block)

    for cand, spec in specs:
        cache_key = (spec[0], tuple(spec[1] or ()))
        part = cache.get(cache_key, [])
        if not part:
            continue
        tokens = tuple(format_k_token(x) for x in cand.coords)
        active = cand.active_coords if cand.active_coords is not None else cand.coords
        active_tokens = tuple(format_k_token(x) for x in active)
        for raw_block in part:
            block = _rotate_block(raw_block, cand.arm_rotation)
            block.k_label = cand.label
            block.k_coords = tokens
            block.active_k_coords = active_tokens
            blocks.append(block)
            degeneracy = int(getattr(raw_block, "degeneracy", 1) or 1)
            if degeneracy < 2:
                continue
            # Keep only a linearly independent co-group orbit.  A tetragonal
            # E irrep needs one 90-degree partner, not all 15 symmetry-related
            # copies; avoiding those duplicates is material for 144-atom ZIPs.
            orbit_basis: list[np.ndarray] = []
            seed = np.asarray(block.dense_displacements(), dtype=float).ravel()
            seed_norm = float(np.linalg.norm(seed))
            if seed_norm > 1e-12:
                orbit_basis.append(seed / seed_norm)
            active_k = active
            for rotation in _little_cogroup_rotations(active_k):
                if np.allclose(rotation, np.eye(3), atol=1e-8):
                    continue
                partner = _rotate_block(block, rotation)
                vector = np.asarray(partner.dense_displacements(), dtype=float).ravel()
                for previous in orbit_basis:
                    vector -= float(np.dot(vector, previous)) * previous
                norm = float(np.linalg.norm(vector))
                if norm < 1e-8:
                    continue
                orbit_basis.append(vector / norm)
                partner.k_label = cand.label
                partner.k_coords = tokens
                partner.active_k_coords = active_tokens
                blocks.append(partner)
                if len(orbit_basis) >= degeneracy:
                    break
    return blocks


def _primary_candidate(subgroup: SubgroupInfo, parent_sg: int) -> _KCandidate | None:
    label = (subgroup.k_point_label or "").strip()
    if not label:
        return None
    tokens = official_special_k_coords(
        parent_sg, label, list(subgroup.k_coordinates or []),
        list(subgroup.k_parameters or []),
    )
    if not tokens:
        tokens = [str(c) for c in (subgroup.k_coordinates or [])]
    if len(tokens) != 3:
        return None
    try:
        coords = np.array([parse_k_token(t) for t in tokens], dtype=float)
    except (ValueError, ZeroDivisionError):
        return None
    params = list(subgroup.k_parameters or []) or None
    return _KCandidate(label, coords, params, tuple(tokens))


def _linearly_independent_mode_items(
    items: Sequence[tuple[str, DistortionMode, np.ndarray, str]],
    lattice: np.ndarray,
    *,
    rtol: float = 1e-8,
) -> list[tuple[str, DistortionMode, np.ndarray, str]]:
    """Select an exact independent subset without changing its mode labels.

    A cosine-like threshold is not a rank test: two legitimate modes can be
    almost parallel, while several moderately correlated duplicates can still
    be linearly dependent.  Modified Gram--Schmidt in Cartesian coordinates
    retains the original labelled vectors and uses residual norm only to
    decide whether each adds a new dimension.
    """
    cart_basis: list[np.ndarray] = []
    kept: list[tuple[str, DistortionMode, np.ndarray, str]] = []
    lattice_array = np.asarray(lattice, dtype=float)
    for item in items:
        original = (np.asarray(item[2], dtype=float) @ lattice_array).ravel()
        original_norm = float(np.linalg.norm(original))
        if original_norm < 1e-12:
            continue
        residual = original.copy()
        for previous in cart_basis:
            residual -= float(np.dot(residual, previous)) * previous
        residual_norm = float(np.linalg.norm(residual))
        if residual_norm <= rtol * original_norm:
            continue
        # A second pass stabilises long, highly correlated mode lists.
        residual /= residual_norm
        for previous in cart_basis:
            residual -= float(np.dot(residual, previous)) * previous
        residual_norm = float(np.linalg.norm(residual))
        if residual_norm <= rtol:
            continue
        cart_basis.append(residual / residual_norm)
        kept.append(item)
    return kept


def compute_parametric_modes(
    parent: Structure,
    symmetry_info: dict,
    subgroup: SubgroupInfo,
    wyckoff_letters: Sequence[str],
    smodes: SmodesWrapper,
    *,
    kpoints: Sequence | None = None,
    nmod: int = 0,
    wyckoff_orbit_ids: Sequence[str] | None = None,
) -> ParametricModeResult:
    """Complete displacive modes for a parametric-k isotropy subgroup."""
    parent_sg = int(symmetry_info.get("space_group_number") or subgroup.parent_sg or 0)
    wyckoff_sites = list(symmetry_info.get("wyckoff_sites") or [])
    basis = np.asarray(
        subgroup.basis_vectors or [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        dtype=float,
    )
    supercell = build_supercell(parent, basis)
    parent_idx, parent_translations = _parent_map(parent, supercell, basis)
    child_origins = _child_origin_candidates(
        subgroup.origin or [0.0, 0.0, 0.0],
        basis,
        parent_translations,
    )
    ops: list[SymmOp] = []
    displacement_actions = None
    for child_origin in child_origins:
        trial_ops = _child_ops(
            int(subgroup.space_group_number or 1),
            origin=child_origin,
            structure=supercell,
        )
        trial_actions = _displacement_actions(supercell, trial_ops)
        if trial_actions is None:
            continue
        ops = trial_ops
        displacement_actions = trial_actions
        break
    if displacement_actions is None:
        return ParametricModeResult(
            modes=[], supercell_displacements={}, labels={}, nmod=max(int(nmod), 0),
            note="child space-group operations do not act completely on the supercell",
    )
    allowed_letters = {str(x) for x in wyckoff_letters if str(x)}
    requested_orbits = {
        str(value) for value in (wyckoff_orbit_ids or ()) if str(value)
    }
    site_by_orbit: dict[str, dict] = {}
    orbit_of_parent: dict[int, str] = {}
    letter_of_orbit: dict[str, str] = {}
    for ordinal, site in enumerate(wyckoff_sites):
        orbit_id = _wyckoff_orbit_id(site, ordinal)
        letter = str(site.get("wyckoff_letter") or "")
        site_by_orbit[orbit_id] = site
        letter_of_orbit[orbit_id] = letter
        orbit_of_parent[int(site["representative_index"])] = orbit_id
        for eq in site.get("equivalent_indices") or []:
            orbit_of_parent[int(eq)] = orbit_id
    allowed_orbits = {
        orbit_id for orbit_id, letter in letter_of_orbit.items()
        if (not allowed_letters or letter in allowed_letters)
        and (not requested_orbits or orbit_id in requested_orbits)
    }
    primary = _primary_candidate(subgroup, parent_sg)
    nmod = max(int(nmod), 0)

    if kpoints is None:
        kpoints = []
    child_symbol = str(subgroup.space_group_symbol or "P").strip()
    child_centering = next(
        (char.upper() for char in child_symbol if char.upper() in "PABCIFR"),
        "P",
    )
    candidates = folding_k_candidates(
        parent_sg,
        kpoints,
        basis,
        primary,
        nmod,
        child_centering=child_centering,
    )

    specs: list[tuple[_KCandidate, tuple[str, list[str] | None]]] = []
    seen_spec: set[tuple[str, tuple[str, ...], tuple[str, ...]]] = set()
    for cand in candidates:
        active = cand.active_coords if cand.active_coords is not None else cand.coords
        key = (
            cand.label,
            tuple(cand.params or ()),
            tuple(format_k_token(value) for value in active),
        )
        if key in seen_spec:
            continue
        seen_spec.add(key)
        specs.append((cand, smodes_k_spec(parent_sg, cand, kpoints)))

    blocks = _collect_smodes_blocks(
        smodes, parent, parent_sg, wyckoff_sites, specs,
    )

    parent_sym = hm_symbol(parent_sg) or str(symmetry_info.get("space_group_symbol") or "P1")
    parent_compact = parent_sym.replace(" ", "")

    collected: list[tuple[str, DistortionMode, np.ndarray, str]] = []
    used_keys: set[str] = set()

    def _unique_key(base: str) -> str:
        key = base
        n = 2
        while key in used_keys:
            key = f"{base}_{n}"
            n += 1
        used_keys.add(key)
        return key

    orbit_order = [
        _wyckoff_orbit_id(site, ordinal)
        for ordinal, site in enumerate(wyckoff_sites)
        if _wyckoff_orbit_id(site, ordinal) in allowed_orbits
    ]
    orbit_letter_counts: dict[str, int] = {}
    for orbit_id in orbit_order:
        letter = letter_of_orbit.get(orbit_id, "")
        orbit_letter_counts[letter] = orbit_letter_counts.get(letter, 0) + 1
    def _k_tokens_for(block: SmodesModeBlock) -> tuple[str, ...]:
        raw = tuple(
            format_k_token(parse_k_token(tok)) if _looks_numeric(tok) else str(tok)
            for tok in (block.k_coords or ())
        )
        if raw:
            return raw
        match = next((c for c in candidates if c.label == block.k_label), None)
        if match is not None:
            return tuple(format_k_token(x) for x in match.coords)
        return ("0", "0", "0")

    # Group smodes blocks by (k, irrep, letter) to assign (a)/(b) and _1/_2.
    grouped: dict[tuple, list[tuple[SmodesModeBlock, np.ndarray, str]]] = {}
    parent_action_cache: dict[
        tuple[str, ...], list[tuple[np.ndarray, np.ndarray]]
    ] = {}
    for block in blocks:
        disp = _map_smodes_block_to_supercell(block, supercell)
        if disp is None:
            continue
        by_orbit: dict[str, np.ndarray] = {}
        for j in range(len(supercell)):
            orbit_id = orbit_of_parent.get(int(parent_idx[j]), "")
            if orbit_id not in allowed_orbits:
                continue
            if orbit_id not in by_orbit:
                by_orbit[orbit_id] = np.zeros_like(disp)
            by_orbit[orbit_id][j] = disp[j]
        k_tokens = _k_tokens_for(block)
        for orbit_id, part in by_orbit.items():
            if float(np.max(np.abs(part))) < 1e-10:
                continue
            letter = letter_of_orbit.get(orbit_id, "")
            group_key = (k_tokens, block.irrep, orbit_id, letter)
            grouped.setdefault(group_key, []).append((block, part, letter))
            if int(getattr(block, "degeneracy", 1) or 1) < 2:
                continue
            try:
                phase_tokens = block.active_k_coords or k_tokens
                kvec = np.array([parse_k_token(t) for t in phase_tokens], dtype=float)
            except (ValueError, ZeroDivisionError):
                continue
            # At a self-conjugate special k, a multidimensional real irrep
            # needs its little-group orbit.  At a general/parametric k,
            # smodes' conjugate partner is instead generated by the phase
            # translation below; adding both constructions mixes harmonic
            # labels even though the total invariant-space rank is unchanged.
            if _is_self_conjugate_parent_k(kvec, parent_sg):
                action_key = tuple(str(token) for token in phase_tokens)
                parent_actions = parent_action_cache.get(action_key)
                if parent_actions is None:
                    parent_actions = _parent_little_group_actions(
                        supercell, basis, parent_sg, kvec,
                    )
                    parent_action_cache[action_key] = parent_actions
                for orbit_part in _parent_little_group_orbit(
                    part,
                    supercell,
                    basis,
                    parent_sg,
                    kvec,
                    int(getattr(block, "degeneracy", 1) or 1),
                    actions=parent_actions,
                )[1:]:
                    grouped[group_key].append((block, orbit_part, letter))
            partner = _quadrature_partner(
                part, supercell, kvec, basis, parent_sg=parent_sg,
            )
            if partner is None:
                continue
            grouped[group_key].append((block, partner, letter))

    for (k_tokens, irrep, orbit_id, letter), group_items in grouped.items():
        # A smodes basis is not generally aligned with the selected OPD.  An
        # invariant mode can be a linear combination of several degenerate
        # blocks/star arms, so project the whole seed space and retain an
        # independent Cartesian basis instead of testing blocks one by one.
        projected_items: list[tuple[SmodesModeBlock, np.ndarray, str]] = []
        cart_basis: list[np.ndarray] = []
        child_lattice = np.asarray(supercell.lattice.matrix, dtype=float)
        inv_child_lattice = np.linalg.inv(child_lattice)
        for block, part, item_letter in group_items:
            projected = _project_invariant(
                part,
                supercell,
                ops,
                actions=displacement_actions,
            )
            if projected is None:
                continue
            cart = (projected @ child_lattice).ravel()
            for previous in cart_basis:
                cart -= float(np.dot(cart, previous)) * previous
            norm = float(np.linalg.norm(cart))
            if norm < 1e-8:
                continue
            cart /= norm
            cart_basis.append(cart.copy())
            frac = cart.reshape((-1, 3)) @ inv_child_lattice
            projected_items.append((block, frac, item_letter))
        projected_items_final = projected_items
        if not projected_items_final:
            continue

        def _weight(item: tuple[SmodesModeBlock, np.ndarray, str]) -> tuple:
            _blk, part, _let = item
            cart = part @ np.asarray(supercell.lattice.matrix, dtype=float)
            return (-float(np.sum(cart * cart)),)
        sorted_items = sorted(projected_items_final, key=_weight)
        # SMODES/Gram--Schmidt columns define a valid numerical fixed space,
        # but their individual orientation does not define a local site irrep,
        # repeated-copy order, or official component name.  Keep every
        # independent vector with an explicit unresolved identity.  A later
        # ISO microscopic query may replace the *whole* group by a canonical
        # basis after rank/principal-angle validation; it must never label
        # these columns one by one.
        for basis_index, (_blk, part, _let) in enumerate(sorted_items):
            orbit_site = site_by_orbit.get(orbit_id) or {}
            representative_index = int(orbit_site.get("representative_index", 0))
            representative_arrow = np.zeros(3, dtype=float)
            for child_index in range(len(supercell)):
                if int(parent_idx[child_index]) != representative_index:
                    continue
                arrow = np.asarray(part[child_index], dtype=float) @ basis
                if float(np.linalg.norm(arrow)) > float(np.linalg.norm(representative_arrow)):
                    representative_arrow = np.array(arrow, dtype=float, copy=True)
            component = _component_label(basis_index)
            nrm = float(np.max(np.abs(part)))
            normalized_part = part
            if nrm > 1e-16:
                normalized_part = part / nrm
            site_name, _elem = _display_site_label(
                wyckoff_sites, letter, orbit_id,
            )
            k_joined = ",".join(k_tokens) if k_tokens else "0,0,0"
            stem = irrep[:2] if irrep else (subgroup.k_point_label or "")
            k_label = _blk.k_label or stem
            site_token = letter
            if orbit_letter_counts.get(letter, 0) > 1:
                site_token = f"{letter}@{_orbit_key_token(orbit_id)}"
            key = _unique_key(
                f"{irrep}[{k_joined}]__{site_token}__"
                f"unresolved({component})"
            )
            pretty = (
                f"{parent_compact}[{k_joined}]{irrep}"
                f"[{site_name}:{letter}:dsp]?"
            )
            identity = ModeIdentity.unresolved(
                parent_sg=parent_sg,
                global_irrep=irrep,
                k_coordinates=tuple(k_tokens),
                wyckoff_letter=letter,
                orbit_id=orbit_id,
                reason="smodes_basis_has_no_iso_microscopic_identity",
            )
            mode = DistortionMode(
                irrep_label=irrep,
                dimension=max(len(sorted_items), 1),
                mode_type="displacive",
                wyckoff_site=letter,
                k_point_label=k_label,
                opd_symbol=subgroup.opd_symbol,
                bush_modes=[
                    BushMode(
                        irrep_label=irrep,
                        opd_symbol=subgroup.opd_symbol,
                        wyckoff_letter=letter,
                        point=list(parent[representative_index].frac_coords),
                        displacements=[representative_arrow.tolist()],
                    )
                ],
                amplitude_key=key,
                site_irrep="",
                k_coords_label=k_joined,
                opd_component=component,
                wyckoff_orbit_id=orbit_id,
                mode_identity=identity,
            )
            collected.append((key, mode, normalized_part, pretty))

    if not collected:
        return ParametricModeResult(
            modes=[], supercell_displacements={}, labels={}, nmod=nmod,
            note="parametric k: smodes returned no child-invariant displacive modes",
        )

    # Deduplicate the complete projected space by rank, not cosine similarity.
    kept = _linearly_independent_mode_items(
        collected,
        np.asarray(supercell.lattice.matrix, dtype=float),
    )

    def _kept_sort(item: tuple[str, DistortionMode, np.ndarray, str]) -> tuple:
        mode = item[1]
        orbit_id = str(mode.wyckoff_orbit_id or "")
        site_i = orbit_order.index(orbit_id) if orbit_id in orbit_order else 99
        fam = _k_family_rank(mode.k_coords_label or "", primary)
        return (
            site_i,
            fam,
            _irrep_sort_key(mode.irrep_label or ""),
            _site_irrep_sort_key(mode.site_irrep or ""),
            mode.opd_component or "a",
        )

    kept.sort(key=_kept_sort)
    modes = [m for _k, m, _d, _lab in kept]
    sc_disp = {k: d for k, _m, d, _lab in kept}
    labels = {k: lab for k, _m, _d, lab in kept}
    return ParametricModeResult(
        modes=modes,
        supercell_displacements=sc_disp,
        labels=labels,
        nmod=nmod,
    )
