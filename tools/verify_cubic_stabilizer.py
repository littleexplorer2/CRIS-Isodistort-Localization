"""Exact finite-group checks for the generic Gamma polar-vector probe.

This is independent analysis, not a production crystallographic solver.
It neither reads website results nor imports CRIS code.
"""
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import json

def matrix(perm, signs):
    return tuple(tuple(signs[i] if j == perm[i] else 0 for j in range(3)) for i in range(3))

def apply(m, v):
    return tuple(sum(a*b for a,b in zip(row,v)) for row in m)

def multiply(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))

def check():
    group = {matrix(p,s) for p in permutations(range(3)) for s in product((-1,1), repeat=3)}
    assert len(group) == 48
    assert all(multiply(a,b) in group for a in group for b in group)
    records = []
    for direction, expected in [((1,0,0),8),((1,1,0),4),((1,1,1),6),((1,2,3),1),((0,0,0),48)]:
        subgroup = {m for m in group if apply(m,direction) == direction}
        assert len(subgroup) == expected
        assert all(multiply(a,b) in subgroup for a in subgroup for b in subgroup)
        projector = tuple(tuple(sum(Fraction(m[i][j],len(subgroup)) for m in subgroup)
                                for j in range(3)) for i in range(3))
        assert multiply(projector,projector) == projector
        assert apply(projector,direction) == direction
        records.append(dict(direction=direction,group_order=48,stabilizer_order=len(subgroup),
                            index=48//len(subgroup),projector=[[str(x) for x in row] for row in projector]))
    # A group operation sending eta to -eta does not fix a nonzero polar OP.
    inversion = ((-1,0,0),(0,-1,0),(0,0,-1))
    assert apply(inversion,(1,0,0)) != (1,0,0)
    return dict(status="passed",arithmetic="exact integers and rational fractions; no tolerance",
                scope="Gamma polar-vector action of cubic signed-permutation group only",
                limitations=["Does not prove ISO implementation for general irreps, nonzero k, or superspace",
                             "Space-group labels are not derived by this check"],
                directions=records)

if __name__ == "__main__":
    result = check()
    path = Path(__file__).resolve().parents[1] / "evidence/native/cubic-stabilizer-check.json"
    path.write_text(json.dumps(result,indent=2), encoding="utf-8")
    print(json.dumps(result))
