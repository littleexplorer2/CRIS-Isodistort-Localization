"""Finite-cell Reynolds construction. See docs/reference_library/LOCAL_MODE_DERIVATION.md.

No website results, numerical ISO projection coefficients, or case tables enter here.
The input H is supplied by the official generic ISO database and is checked against G.
"""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from ..structure import SymmetryOperation, expand_structure

IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO = (F(0),) * 3
MAX_ATOMS = 256
MAX_GROUP = 384
MAX_TRANSLATION_TRIALS = 1000000


class ModeError(ValueError):
    def __init__(self, code, message, **details):
        super().__init__(message)
        self.code, self.details = code, details


def rational_matrix(matrix):
    return tuple(tuple(F(x) for x in row) for row in matrix)


def determinant(a):
    return sum(a[0][i] * (a[1][(i+1) % 3] * a[2][(i+2) % 3]
                          - a[1][(i+2) % 3] * a[2][(i+1) % 3]) for i in range(3))


def inverse(a):
    d = determinant(a)
    if not d:
        raise ModeError("singular_basis", "The subgroup basis is singular.")
    return tuple(tuple((a[(j+1) % 3][(i+1) % 3] * a[(j+2) % 3][(i+2) % 3]
                        - a[(j+1) % 3][(i+2) % 3] * a[(j+2) % 3][(i+1) % 3]) / d
                       for j in range(3)) for i in range(3))


def mv(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(3)) for i in range(3))


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
                 for i in range(3))


def close_group(generators):
    """Exact finite Seitz closure modulo this cell's integer translations."""
    generators = tuple(dict.fromkeys(generators))
    identity = SymmetryOperation(IDENTITY, ZERO)
    result, seen = [identity], {identity}
    for op in result:
        for generator in generators:
            new = generator.compose(op)
            if new not in seen:
                if len(result) >= MAX_GROUP:
                    raise ModeError("group_limit", f"Finite group exceeds {MAX_GROUP} operations.")
                seen.add(new)
                result.append(new)
    return tuple(result)


def complete_parent_group(cosets, primitive_basis):
    generators = [SymmetryOperation(op.rotation, op.translation) for op in cosets]
    generators += [SymmetryOperation(IDENTITY, tuple(primitive_basis[i][j] for i in range(3)))
                   for j in range(3)]
    return close_group(generators)


def _supercell_atoms(structure, c, origin):
    atoms = expand_structure(structure)
    ci, det = inverse(c), determinant(c)
    if det <= 0:
        raise ModeError("left_handed_basis", "Only positive-volume subgroup bases are supported.")
    expected = len(atoms) * det
    if expected.denominator != 1 or expected > MAX_ATOMS:
        raise ModeError("atom_limit", "Cell atom count must be integral and at most 256.",
                        expected=str(expected))
    corners = [tuple(origin[i] + v[i] for i in range(3))
               for v in (mv(c, corner) for corner in product((F(0), F(1)), repeat=3))]
    output, seen, trials = [], set(), 0
    for atom in atoms:
        ranges = [range(math.floor(min(v[i] for v in corners) - atom.fractional[i]),
                        math.ceil(max(v[i] for v in corners) - atom.fractional[i]) + 1)
                  for i in range(3)]
        trials += math.prod(len(r) for r in ranges)
        if trials > MAX_TRANSLATION_TRIALS:
            raise ModeError("translation_limit", "Supercell enumeration exceeds the explicit budget.")
        for translation in product(*ranges):
            y = mv(ci, tuple(atom.fractional[i] + translation[i] - origin[i] for i in range(3)))
            if not all(0 <= x < 1 for x in y):
                continue
            if y in seen:
                raise ModeError("duplicate_supercell_atom", "Supercell enumeration generated overlapping atoms.")
            seen.add(y)
            output.append({"fractional_exact": y, "symbol": atom.type_symbol,
                           "occupancy": atom.occupancy, "source_label": atom.source_label,
                           "source_index": atom.source_index,
                           "parent_translation": translation})
    if len(output) != expected:
        raise ModeError("atom_count_mismatch", "Supercell enumeration violates volume/multiplicity.",
                        expected=str(expected), actual=len(output))
    for index, atom in enumerate(output):
        atom["label"] = f"a{index + 1}"
    return output


def _canonical_basis(projector, expected_rank, tol):
    """Deterministic twice-reorthogonalized columns, with a rank/span check."""
    columns = []
    for candidate in projector.T:
        value = candidate.copy()
        for _ in range(2):
            for column in columns:
                value -= column * np.dot(column, value)
        norm = np.linalg.norm(value)
        if norm > tol:
            value /= norm
            pivot = np.flatnonzero(np.abs(value) > tol)[0]
            if value[pivot] < 0:
                value *= -1
            columns.append(value)
    basis = np.column_stack(columns) if columns else np.empty((len(projector), 0))
    if len(columns) != expected_rank:
        raise ModeError("rank_failure", "Projector rank disagrees with the independent character count.",
                        expected=expected_rank, actual=len(columns))
    span_error = float(np.linalg.norm(projector - basis @ basis.T, ord=np.inf))
    if span_error > tol:
        raise ModeError("span_failure", "Basis does not span the projector.", residual=span_error)
    return basis


def _check_projector(p, character_rank, tol):
    rank = round(character_rank)
    residuals = {"symmetry": float(np.linalg.norm(p-p.T, ord=np.inf)),
                 "idempotence": float(np.linalg.norm(p@p-p, ord=np.inf)),
                 "trace_integer": abs(float(character_rank)-rank),
                 "trace_agreement": abs(float(np.trace(p))-character_rank)}
    if max(residuals.values()) > tol:
        raise ModeError("projector_failure", "Reynolds projector invariants failed.", **residuals)
    return rank, residuals


def _kelvin():
    tensors = []
    for i in range(3):
        t = np.zeros((3, 3)); t[i, i] = 1
        tensors.append(t)
    for i, j in ((1, 2), (0, 2), (0, 1)):
        t = np.zeros((3, 3)); t[i, j] = t[j, i] = 1/math.sqrt(2)
        tensors.append(t)
    return np.array(tensors)


def construct_modes(structure, parent_group, subgroup, subgroup_parent):
    """Inputs are typed ISO records; basis matrices have printed vectors as columns."""
    c, origin = rational_matrix(subgroup.basis), tuple(map(F, subgroup.origin))
    ci = inverse(c)
    g = complete_parent_group(parent_group.operations, parent_group.primitive_basis)
    if set(g) != set(structure.symmetry_operations):
        raise ModeError("parent_setting_mismatch",
                        "CIF operations differ from ISO INTERNATIONAL NEW conventional setting.")
    if any(site.occupancy != 1 for site in structure.sites):
        raise ModeError("occupancy_unsupported", "The mode workflow currently requires full occupancy.")
    pure_g = {op.translation for op in g if op.rotation == IDENTITY}
    for j in range(3):
        if tuple(c[i][j] % 1 for i in range(3)) not in pure_g:
            raise ModeError("incompatible_translation", "A subgroup cell vector is not a parent translation.")
    m_g = len(pure_g)
    p_h = rational_matrix(subgroup_parent.primitive_basis)
    if determinant(p_h) <= 0:
        raise ModeError("primitive_basis_orientation", "Unsupported primitive basis orientation.")
    multiplicity_h = 1 / determinant(p_h)
    expected_size = determinant(c) * m_g / multiplicity_h
    if expected_size != subgroup.size:
        raise ModeError("cell_size_mismatch", "ISO Size disagrees with the exact primitive volume ratio.",
                        calculated=str(expected_size), reported=subgroup.size)
    generators = []
    for op in subgroup.operations:
        rp, tp = rational_matrix(op.rotation), tuple(map(F, op.translation))
        if SymmetryOperation(rp, tp) not in set(g):
            raise ModeError("not_parent_subgroup", "An H operation is absent from the CIF parent group.")
        rs = mm(mm(ci, rp), c)
        ts = mv(ci, tuple(x+y-z for x, y, z in zip(mv(rp, origin), tp, origin)))
        generators.append(SymmetryOperation(rs, ts))
    generators += [SymmetryOperation(IDENTITY, tuple(p_h[i][j] for i in range(3)))
                   for j in range(3)]
    h = close_group(generators)
    for op in h:
        rp = mm(mm(c, op.rotation), ci)
        tp = tuple(a+b-d for a, b, d in zip(origin, mv(c, op.translation), mv(rp, origin)))
        if SymmetryOperation(rp, tp) not in set(g):
            raise ModeError("closure_not_parent", "Closed H contains a non-parent operation.")
    # |G:H| = volume(H conventional)/volume(G conventional) * |G/T_Gconv|/|H/T_Hconv|.
    calculated_index = determinant(c) * len(g) / len(h)
    if calculated_index != subgroup.index:
        raise ModeError("group_index_mismatch", "Group orders and volume disagree with the ISO index.",
                        calculated=str(calculated_index), reported=subgroup.index)
    atoms = _supercell_atoms(structure, c, origin)
    n, a = len(atoms), np.asarray(structure.lattice.matrix, dtype=float) @ np.array(c, float)
    condition = float(np.linalg.cond(a))
    # Conservative accumulation bound in binary64, not a tolerance fitted to output.
    tol = 256*np.finfo(float).eps*max(3*n, len(h), 6)*max(1., condition**2)
    if tol > 1e-7:
        raise ModeError("ill_conditioned_cell", "Cell condition exceeds the supported precision bound.",
                        condition=condition, error_bound=tol)
    lookup = {atom["fractional_exact"]: i for i, atom in enumerate(atoms)}
    permutations, rotations = [], []
    projector, trans_projector, strain_projector = np.zeros((3*n, 3*n)), np.zeros((3, 3)), np.zeros((6, 6))
    tensors, character = _kelvin(), 0.
    orthogonal_error = 0.
    for operation in h:
        try:
            perm = [lookup[operation.apply(atom["fractional_exact"])] for atom in atoms]
        except KeyError as exc:
            raise ModeError("site_mapping_failure", "H does not act on the enumerated supercell.") from exc
        if len(set(perm)) != n or any(atoms[i]["symbol"] != atoms[j]["symbol"]
                                     or atoms[i]["source_index"] != atoms[j]["source_index"]
                                     for i, j in enumerate(perm)):
            raise ModeError("site_species_failure", "H does not preserve atomic identities.")
        q = a @ np.asarray(operation.rotation, float) @ np.linalg.inv(a)
        orthogonal_error = max(orthogonal_error, float(np.linalg.norm(q.T@q-np.eye(3), ord=np.inf)))
        for i, j in enumerate(perm):
            projector[3*j:3*j+3, 3*i:3*i+3] += q / len(h)
        trans_projector += q / len(h)
        representation = np.einsum("aij,bij->ab", tensors, q @ tensors @ q.T)
        strain_projector += representation / len(h)
        character += sum(i == j for i, j in enumerate(perm))*np.trace(q)/len(h)
        rotations.append(q); permutations.append(perm)
    if orthogonal_error > tol:
        raise ModeError("metric_precision_failure",
                        "The supplied metric is not compatible with an orthogonal group action at this bound.",
                        residual=orthogonal_error, tolerance=tol)
    # Exact permutation homomorphism independently checks each composition and convention.
    h_lookup = {operation: index for index, operation in enumerate(h)}
    for i, left in enumerate(h):
        for j, right in enumerate(h):
            k = h_lookup[left.compose(right)]
            if [permutations[i][p] for p in permutations[j]] != permutations[k]:
                raise ModeError("permutation_composition_failure", "Site action violates group composition.")
    rank, displacement_residuals = _check_projector(projector, character, tol)
    trank, translation_residuals = _check_projector(trans_projector, np.trace(trans_projector), tol)
    t_basis = _canonical_basis(trans_projector, trank, tol)
    uniform = np.tile(t_basis, (n, 1))/math.sqrt(n)
    internal_p = projector - uniform@uniform.T
    internal = _canonical_basis(internal_p, rank-trank, tol)
    basis = np.column_stack((internal, uniform))
    srank, strain_residuals = _check_projector(strain_projector, np.trace(strain_projector), tol)
    strain_basis = _canonical_basis(strain_projector, srank, tol)
    invariant_error = 0.
    for perm, q in zip(permutations, rotations):
        transformed = np.empty_like(basis.reshape(n, 3, rank))
        transformed[perm] = np.einsum("ij,njk->nik", q, basis.reshape(n, 3, rank))
        invariant_error = max(invariant_error, float(np.max(np.abs(
            transformed.reshape(3*n, rank)-basis), initial=0.)))
    orthonormal_error = float(np.max(np.abs(basis.T@basis-np.eye(rank)), initial=0.))
    if max(invariant_error, orthonormal_error) > tol:
        raise ModeError("basis_validation_failure", "Final basis invariants failed.")
    return {"atoms": atoms, "lattice": a, "operations": h, "permutations": permutations,
            "basis": basis, "strain_basis": strain_basis, "strain_tensors": tensors,
            "rigid_translations": trank, "internal_modes": rank-trank,
            "parent_primitive_cells": int(determinant(c)*m_g),
            "basis_convention": "CRIS_REA_H_invariant_Cartesian_orthonormal_v1",
            "validation": {"tolerance": tol, "cell_condition": condition,
                           "rotation_orthogonality": orthogonal_error,
                           "mode_invariance": invariant_error, "orthonormality": orthonormal_error,
                           "displacement_projector": displacement_residuals,
                           "translation_projector": translation_residuals,
                           "strain_projector": strain_residuals,
                           "character_dimension": rank, "group_order": len(h),
                           "group_index": int(calculated_index), "primitive_size": int(expected_size),
                           "exact_group_closure": True, "exact_permutation_composition": True,
                           "complete_in_declared_H_fixed_space": True}}


def distort(model, amplitudes=None, strain_amplitudes=None):
    u, s = model["basis"], model["strain_basis"]
    def coefficients(values, length, unit):
        result = np.zeros(length) if values is None else np.asarray(values, float)
        if result.shape != (length,) or not np.all(np.isfinite(result)):
            raise ModeError("invalid_amplitudes", f"Supply exactly {length} finite coefficients ({unit}).")
        return result
    amplitudes = coefficients(amplitudes, u.shape[1], "angstrom")
    strain_amplitudes = coefficients(strain_amplitudes, s.shape[1], "dimensionless")
    displacement = (u @ amplitudes).reshape(-1, 3)
    strain = np.einsum("a,aij->ij", s @ strain_amplitudes, model["strain_tensors"])
    deformation = np.eye(3) + strain
    if np.linalg.eigvalsh(deformation).min() <= 0:
        raise ModeError("nonpositive_deformation", "I + strain must be positive definite.")
    coordinates = np.array([atom["fractional_exact"] for atom in model["atoms"]], float)
    coordinates = (coordinates + displacement @ np.linalg.inv(model["lattice"]).T) % 1
    lattice = deformation @ model["lattice"]
    # A small distance threshold detects collisions, not atom matching or symmetry snapping.
    for i in range(len(coordinates)):
        delta = coordinates[i+1:] - coordinates[i]
        delta -= np.rint(delta)
        if len(delta) and np.any(np.linalg.norm(delta @ lattice.T, axis=1) < 1e-9):
            raise ModeError("atom_collision", "Chosen amplitudes cause coincident atomic sites.")
    return {"lattice": lattice, "fractional": coordinates, "displacements_angstrom": displacement,
            "amplitudes_angstrom": amplitudes, "strain_amplitudes": strain_amplitudes,
            "strain_tensor": strain, "deformation": deformation}


def seitz_expression(operation):
    expressions = []
    for row, translation in zip(operation.rotation, operation.translation):
        text = ""
        for coefficient, name in zip(row, "xyz"):
            if coefficient:
                term = name if abs(coefficient) == 1 else f"{abs(coefficient)}*{name}"
                text += ("+" if coefficient > 0 else "-") + term
        if translation:
            text += "+" + str(translation)
        expressions.append(text.lstrip("+") or "0")
    return ",".join(expressions)
