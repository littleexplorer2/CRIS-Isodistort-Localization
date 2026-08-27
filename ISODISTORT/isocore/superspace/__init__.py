"""isocore (3+d) 超空间计算内核。

``nmod``（官网 # of independent incommensurate modulations）即附加维度 d。
``nmod=0`` 为普通三维；``nmod=1`` 为 (3+1)；取位固定 standard (IT-C)。
"""
from .constants import EPS, MAX_NMOD, PHYSICAL_DIM, SERIAL_KIND, SERIAL_KIND_3P1, eps, max_nmod
from .group import (
    SuperspaceGroup,
    SuperspaceOperation,
    lift_operation,
    lift_operation_3p1,
)
from .kvector import KsVector, SuperspaceLattice, assert_ks_compatible, default_lattice_matrix
from .representation import (
    Irrep,
    OpdOperator,
    SuperspaceMode,
    apply_operation_to_mode,
    build_irreps,
    filter_allowed_opds,
    generate_modes,
    generate_opds,
    little_group,
    transform_opd,
)
from .validate import validate_nmod, validate_space_group_number
from .workflow import (
    SuperspaceResult,
    format_superspace_report,
    load_superspace_result,
    run_superspace_cli,
    run_superspace_workflow,
    save_superspace_result,
)

__all__ = [
    "EPS",
    "MAX_NMOD",
    "PHYSICAL_DIM",
    "SERIAL_KIND",
    "SERIAL_KIND_3P1",
    "Irrep",
    "KsVector",
    "OpdOperator",
    "SuperspaceGroup",
    "SuperspaceLattice",
    "SuperspaceMode",
    "SuperspaceOperation",
    "SuperspaceResult",
    "apply_operation_to_mode",
    "assert_ks_compatible",
    "build_irreps",
    "default_lattice_matrix",
    "eps",
    "filter_allowed_opds",
    "format_superspace_report",
    "generate_modes",
    "generate_opds",
    "lift_operation",
    "lift_operation_3p1",
    "little_group",
    "load_superspace_result",
    "max_nmod",
    "run_superspace_cli",
    "run_superspace_workflow",
    "save_superspace_result",
    "transform_opd",
    "validate_nmod",
    "validate_space_group_number",
]
