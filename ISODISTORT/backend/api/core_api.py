"""
核心 Python API - 一站式 ISODISTORT 使用入口

对应全阶段：封装官网 Search Page → Distortion Page 的完整交互流程。

使用示例：
    from backend.api import IsoDistort

    iso = IsoDistort()
    # 1. 从 CIF 加载结构
    iso.load_structure("input.cif")
    # 2. Method 1：枚举全部特殊 k 点子群并过滤
    candidates = iso.search_method_1(crystal_system="tetragonal")
    # 3. Method 2：选择子群并计算畸变模式
    result = iso.search_method_2(subgroup_idx=candidates[0].subgroup.index)
    # 4. Distortion Page：生成畸变结构
    distorted = iso.generate_distortion(amplitude=0.1)
    # 5. 导出
    iso.export("output", formats=["cif", "poscar"])
    iso.export_subgroups("out_batch", formats=["cif", "topas"])
"""
import copy
import hashlib
import html
import json
import os
import re
import shutil
import tempfile
import threading
import time
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
from pathlib import Path

import numpy as np
import spglib
from pymatgen.core import DummySpecies, Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from backend.tables.kpoints_official import (
    KPOINT_OFFICIAL,
    official_kparams_to_iso,
    official_special_k_coords,
)
from backend.utils import IsodistortError, get_config
from backend.utils.lattice import (
    as_fraction,
    centering_primitive_matrix,
    identity_matrix,
    inverse,
    multiply,
    rational_matrix,
    transpose,
)
from backend.utils.opd_format import _centering_letter, format_k_active
from backend.utils.parent_header import parent_wyckoff_display
from backend.utils.schoenflies import hm_symbol, schoenflies_symbol
from backend.utils.text_parser import parse_basis_token, parse_fraction
from backend.wrappers import (
    DistortionMode,
    FindsymSettingMismatchError,
    FindsymWrapper,
    IsoWrapper,
    SubgroupInfo,
)
from backend.wrappers.isotropy_cache import (
    IsotropyCacheEntry,
    delete_isotropy_cache,
    list_isotropy_cache,
)
from backend.wrappers.smodes_wrapper import SmodesWrapper
from features.export import (
    DisplaciveExportData,
    DisplaciveSubgroupIdentity,
    ExactParentChildEmbedding,
    ParentChildSiteMapping,
    ParentOrbitType,
    StrainExportData,
    StructureExporter,
    SubgroupExportSpec,
    build_export_zip,
    method3_case_folder,
    parse_export_formats,
    parse_export_method,
    render_subgroup_files,
    subgroup_label,
    unique_folder_name,
)
from features.export.cif_displacive_model import (
    ChildAtom,
    ExactChildFrame,
    ExactSeitzOperation,
    ModeScaleProvenance,
    RawModeColumn,
    build_displacive_cif_model,
)
from features.input_cif import (
    SymmetryValidator,
    build_supercell,
    read_cif,
    read_cif_space_group_number,
    read_structure,
)
from features.method1 import (
    DEFAULT_DISTORTION_TYPES,
    DISTORTION_TYPES,
    DomainGenerator,
    IsoSearchEngine,
    Method1Query,
    Method1ResultItem,
    Method2Query,
    Method3Query,
    Method4Query,
    OccupationalModeGenerator,
    PhasePath,
    embedding_from_identity,
    normalize_distortion_types,
    parent_affine_group,
)
from features.method1.affine_embeddings import translation_cosets
from features.method1.search_methods import (
    _sg_to_crystal_system,
    _validated_single_irrep_count,
)
from features.method2 import DistortionMapper, verify_mapped_microscopic_columns
from features.method3 import (
    FiniteGroup,
    RationalRepresentation,
    analyze_fixed_space,
    symmetric_square_representation,
)
from features.method4 import (
    DistortionEngine,
    apply_canonical_strain_basis,
    canonical_strain_definitions_from_iso,
    compute_homogeneous_strain_modes,
)
from frontend.i18n import t

_CRIS_ROOT = Path(__file__).resolve().parents[3]
_EXPORT_BATCH_SCHEMA = "isodistort.export-batch.v1"
_EXPORT_BATCH_PREFIX = ".isodistort-batch-v1-"
_EXPORT_BATCH_READY_SUFFIX = ".ready"
_EXPORT_BATCH_MANIFEST = ".isodistort-batch.json"
_EXPORT_PUBLISH_THREAD_LOCK = threading.Lock()


@contextmanager
def _export_publish_lock(folder_root: Path):
    """Serialize the final namespace commit across threads and processes.

    The OS byte-range lock is released automatically if a writer process
    exits.  A process-local lock is also required because Windows byte-range
    locks are process scoped and therefore do not serialize sibling threads.
    The stable sibling lock file is intentionally retained so waiting writers
    always refer to the same filesystem object.
    """

    identity = str(folder_root.resolve(strict=False)).casefold().encode("utf-8")
    token = hashlib.sha256(identity).hexdigest()[:16]
    lock_path = folder_root.parent / f".isodistort-publish-{token}.lock"
    with _EXPORT_PUBLISH_THREAD_LOCK, lock_path.open("a+b") as handle:
        handle.seek(0, 2)
        if handle.tell() == 0:
            handle.write(b"\0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt  # noqa: PLC0415

            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl  # noqa: PLC0415

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _export_batch_manifest(
    method: int,
    rendered_batch: list[
        tuple[SubgroupExportSpec, tuple[tuple[str, bytes], ...]]
    ],
    folders: list[str],
    root_payloads: tuple[tuple[str, bytes], ...] = (),
) -> bytes:
    """Build the immutable manifest stored inside one committed disk batch."""

    payload = {
        "schema": _EXPORT_BATCH_SCHEMA,
        "method": int(method),
        "candidates": [
            {
                "folder": folder,
                "files": [
                    {
                        "name": str(filename),
                        "size": len(contents),
                        "sha256": hashlib.sha256(contents).hexdigest(),
                    }
                    for filename, contents in files
                ],
            }
            for (_spec, files), folder in zip(
                rendered_batch,
                folders,
                strict=True,
            )
        ],
        "root_files": [
            {
                "name": str(filename),
                "size": len(contents),
                "sha256": hashlib.sha256(contents).hexdigest(),
            }
            for filename, contents in root_payloads
        ],
    }
    return (
        json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode("utf-8")


def _ready_export_batch_name(manifest: bytes) -> str:
    """Return the short content-addressed name used for one committed batch."""

    digest = hashlib.sha256(manifest).hexdigest()[:16]
    return f"{_EXPORT_BATCH_PREFIX}{digest}{_EXPORT_BATCH_READY_SUFFIX}"


def _ready_export_batch_folders(batch_root: Path) -> set[str]:
    """Validate one ready manifest and return its logical candidate folders."""

    if not batch_root.is_dir() or batch_root.is_symlink():
        raise ValueError(f"committed export batch is not a directory: {batch_root}")
    manifest_path = batch_root / _EXPORT_BATCH_MANIFEST
    manifest = manifest_path.read_bytes()
    expected_name = _ready_export_batch_name(manifest)
    if batch_root.name != expected_name:
        raise ValueError(
            "committed export batch manifest does not match its version directory: "
            f"{batch_root}"
        )
    try:
        decoded = json.loads(manifest.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid committed export batch manifest: {batch_root}") from exc
    if decoded.get("schema") != _EXPORT_BATCH_SCHEMA:
        raise ValueError(f"unsupported committed export batch manifest: {batch_root}")
    candidates = decoded.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError(f"invalid committed export batch candidates: {batch_root}")

    folders: set[str] = set()
    occupied: set[str] = set()
    root_files = decoded.get("root_files", [])
    if not isinstance(root_files, list):
        raise ValueError(f"invalid committed export root file list: {batch_root}")
    for file_record in root_files:
        if not isinstance(file_record, dict):
            raise ValueError(f"invalid committed export root file: {batch_root}")
        name = file_record.get("name")
        if not isinstance(name, str) or not name or Path(name).name != name:
            raise ValueError(f"unsafe committed export root filename: {name!r}")
        folded_name = name.casefold()
        if folded_name in occupied or folded_name == _EXPORT_BATCH_MANIFEST.casefold():
            raise ValueError(f"duplicate committed export root filename: {name!r}")
        path = batch_root / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"committed export root file is missing: {path}")
        expected_size = file_record.get("size")
        expected_sha256 = file_record.get("sha256")
        contents = path.read_bytes()
        if (
            not isinstance(expected_size, int)
            or isinstance(expected_size, bool)
            or expected_size < 0
            or not isinstance(expected_sha256, str)
            or re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None
            or len(contents) != expected_size
            or hashlib.sha256(contents).hexdigest() != expected_sha256
        ):
            raise ValueError(f"committed export root file digest changed: {path}")
        occupied.add(folded_name)
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError(f"invalid committed export batch candidate: {batch_root}")
        folder = candidate.get("folder")
        if not isinstance(folder, str) or not folder or Path(folder).name != folder:
            raise ValueError(f"unsafe committed export candidate folder: {folder!r}")
        folded = folder.casefold()
        if folded in occupied:
            raise ValueError(
                f"duplicate committed export candidate folder: {folder!r}"
            )
        candidate_root = batch_root / folder
        if not candidate_root.is_dir() or candidate_root.is_symlink():
            raise ValueError(
                f"committed export candidate directory is missing: {candidate_root}"
            )
        files = candidate.get("files")
        if not isinstance(files, list):
            raise ValueError(
                f"invalid committed export candidate file list: {candidate_root}"
            )
        file_names: set[str] = set()
        for file_record in files:
            if not isinstance(file_record, dict):
                raise ValueError(
                    f"invalid committed export candidate file: {candidate_root}"
                )
            name = file_record.get("name")
            if not isinstance(name, str) or not name or Path(name).name != name:
                raise ValueError(f"unsafe committed export filename: {name!r}")
            folded_name = name.casefold()
            if folded_name in file_names:
                raise ValueError(f"duplicate committed export filename: {name!r}")
            path = candidate_root / name
            if not path.is_file() or path.is_symlink():
                raise ValueError(f"committed export file is missing: {path}")
            expected_size = file_record.get("size")
            expected_sha256 = file_record.get("sha256")
            if (
                not isinstance(expected_size, int)
                or isinstance(expected_size, bool)
                or expected_size < 0
                or not isinstance(expected_sha256, str)
                or re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None
            ):
                raise ValueError(f"invalid committed export file digest: {path}")
            contents = path.read_bytes()
            if len(contents) != expected_size:
                raise ValueError(f"committed export file size changed: {path}")
            if hashlib.sha256(contents).hexdigest() != expected_sha256:
                raise ValueError(f"committed export file hash changed: {path}")
            file_names.add(folded_name)
        occupied.add(folded)
        folders.add(folder)
    return folders


def _published_export_folder_names(folder_root: Path) -> set[str]:
    """Collect candidate names from legacy roots and committed version batches."""

    if (
        not folder_root.exists()
        or not folder_root.is_dir()
        or folder_root.is_symlink()
    ):
        return set()
    reserved: set[str] = set()
    for entry in folder_root.iterdir():
        name = entry.name
        if name.startswith(_EXPORT_BATCH_PREFIX) and name.endswith(
            _EXPORT_BATCH_READY_SUFFIX
        ):
            reserved.update(_ready_export_batch_folders(entry))
        else:
            # Preserve the historical rule: any existing root entry reserves
            # its case-insensitive name, regardless of its filesystem type.
            reserved.add(name)
    return reserved


@dataclass(frozen=True)
class ExportCandidateFailure:
    """One candidate that failed calculation, preparation, or writer preflight."""

    position: int
    candidate_index: int
    subgroup_index: int
    parent_space_group_number: int
    k_point_label: str
    k_coordinates: tuple[str, ...]
    k_parameters: tuple[str, ...]
    irrep_label: str
    opd_symbol: str
    space_group_number: int
    basis_vectors: tuple[tuple[str, ...], ...]
    origin: tuple[str, ...]
    embedding_id: str
    error_type: str
    message: str


@dataclass(frozen=True)
class ModeIdentityFailure:
    """One emitted displacement column whose scientific identity is unresolved."""

    amplitude_key: str
    orbit_id: str
    global_irrep: str
    k_coordinates: tuple[str, ...]
    reason: str


class UnresolvedModeIdentityError(IsodistortError):
    """A candidate has numerical columns but cannot enter authoritative writers."""

    def __init__(self, failures: list[ModeIdentityFailure]) -> None:
        if not failures:
            raise ValueError("unresolved mode identity error requires failures")
        self.failures = tuple(failures)
        details = "; ".join(
            (
                f"mode={failure.amplitude_key} orbit={failure.orbit_id or '<missing>'} "
                f"IR={failure.global_irrep or '<missing>'} "
                f"k=({','.join(failure.k_coordinates)}) "
                f"criterion={failure.reason}"
            )
            for failure in self.failures
        )
        super().__init__(
            "displacement modes lack verified ISO microscopic identities: " + details
        )


class ExportCancelledError(IsodistortError):
    """A batch stopped at a candidate boundary after a user cancellation."""


@dataclass(frozen=True)
class ExportBatchPlan:
    """Prepared candidates plus explicitly ineligible identity failures."""

    specs: tuple[SubgroupExportSpec, ...]
    skipped: tuple[ExportCandidateFailure, ...]


def _export_failure_report_payloads(
    method: int,
    requested_count: int,
    successful_count: int,
    failures: tuple[ExportCandidateFailure, ...],
) -> tuple[tuple[str, bytes], ...]:
    """Build root-level human and machine records for ineligible candidates."""

    if not failures:
        return ()
    cancelled_count = sum(
        failure.error_type == "ExportCancelledError" for failure in failures
    )
    unresolved_count = sum(
        failure.error_type == "UnresolvedModeIdentityError" for failure in failures
    )
    payload = {
        "schema": "isodistort.export-candidate-status.v1",
        "method": int(method),
        "requested_candidate_count": int(requested_count),
        "successful_candidate_count": int(successful_count),
        "ineligible_candidate_count": len(failures),
        "unresolved_identity_candidate_count": unresolved_count,
        "cancelled_candidate_count": cancelled_count,
        "status": (
            "cancelled_partial_success"
            if cancelled_count and successful_count
            else "cancelled"
            if cancelled_count
            else "partial_success"
            if successful_count
            else "no_eligible_candidates"
        ),
        "failures": [asdict(failure) for failure in failures],
    }
    machine = (
        json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    ).encode("utf-8")
    lines = [
        f"Method {int(method)} export candidate status",
        f"requested={int(requested_count)} successful={int(successful_count)} "
        f"ineligible={len(failures)}",
        "",
    ]
    for failure in failures:
        lines.extend((
            (
                f"#{failure.position} candidate_index={failure.candidate_index} "
                f"{failure.k_point_label}{failure.k_coordinates}/"
                f"{failure.k_parameters}/{failure.irrep_label}/"
                f"{failure.opd_symbol}/SG{failure.space_group_number}"
            ),
            f"  {failure.error_type}: {failure.message}",
        ))
    human = ("\n".join(lines).rstrip() + "\n").encode("utf-8")
    return (
        ("export_candidate_status.json", machine),
        ("export_candidate_status.txt", human),
    )


class ExportPreparationError(IsodistortError):
    """Batch export failed before publication because candidates were not ready."""

    def __init__(
        self,
        method: int,
        failures: list[ExportCandidateFailure],
    ) -> None:
        self.method = int(method)
        self.failures = tuple(failures)
        details = "; ".join(
            (
                f"#{failure.position} candidate_index={failure.candidate_index} "
                f"subgroup_index={failure.subgroup_index} "
                f"parent_SG={failure.parent_space_group_number} "
                f"{failure.k_point_label}{failure.k_coordinates}/"
                f"{failure.k_parameters}/{failure.irrep_label}/"
                f"{failure.opd_symbol}/SG{failure.space_group_number} "
                f"basis={failure.basis_vectors} origin={failure.origin} "
                f"embedding={failure.embedding_id}: "
                f"{failure.error_type}: {failure.message}"
            )
            for failure in self.failures
        )
        super().__init__(
            f"Method {self.method} batch export preparation failed for "
            f"{len(self.failures)} candidate(s): {details}"
        )


class IsoDistort:
    """
    ISODISTORT 主入口类

    封装完整工作流：
        加载结构 → 识别对称 → 枚举子群（Method 1）→ 选择路径（Method 2）
        → 计算畸变模式 → 生成畸变结构 → 导出/畴

    UI strings are English only. The ``language`` argument is ignored and kept
    so existing ``IsoDistort(language="en")`` call sites still construct.
    """

    def __init__(self, language: str | None = None) -> None:
        self.cfg = get_config()
        _ = language

        # 底层封装
        self._iso = IsoWrapper()
        self._findsym: FindsymWrapper | None = None

        # 业务层
        self._sym_val = SymmetryValidator()
        self._dist_mapper = DistortionMapper()
        self._dist_engine = DistortionEngine(self._dist_mapper)
        self._domain_gen = DomainGenerator(self._iso)
        self._search = IsoSearchEngine(self._iso)

        # 输出层
        self._exporter = StructureExporter()

        # 状态
        self.structure: Structure | None = None
        self.symmetry_info: dict | None = None
        self.structure_path: Path | None = None
        self.subgroups: list[SubgroupInfo] = []
        # Affine-only Method-3 diagnostics must survive copying/rehydration of
        # ``SubgroupInfo`` objects.  Python ``id(...)`` is a process-local
        # memory address and can be reused, so track the deterministic path
        # exact content-addressed Seitz embedding ID (with an exact path-key
        # compatibility fallback) instead.
        self._unresolved_method3_embedding_keys: set[tuple] = set()
        # 当前已计算模式所对应的完整子群对象。不能只保存 ``index``：
        # Method 1/2/3 的候选池都会从 0 编号，同号不代表同一个子群。
        self._selected_subgroup: SubgroupInfo | None = None
        self.phase_path: PhasePath | None = None
        self.distortion_modes: list[DistortionMode] = []
        self.mode_displacements: dict = {}
        self.mode_occupancies: dict = {}          # occupational 模式（占据率调制）
        self.mode_displacements_sc: dict = {}     # parametric complete modes on supercell
        self._mode_label_overrides: dict[str, str] = {}
        self.number_of_independent_modulations: int = 0
        self.distorted_structure: Structure | None = None
        # The arrays above are a cache, not timeless session state.  Bind them
        # to the complete scientific calculation context so changing a Types
        # checkbox, species scope, nmod, or affine embedding cannot silently
        # reuse an incompatible mode basis or generated final structure.
        self._mode_cache_key: tuple | None = None
        self._generated_structure_cache_key: tuple | None = None

        # 畸变类型作用域（对齐官网 per-species 复选框）：type -> 物种列表（"*"=全部）
        self.distortion_scope: dict[str, list[str]] = {}
        self.distortion_types: list[str] = DEFAULT_DISTORTION_TYPES.copy()
        self._smodes = SmodesWrapper()
        self._special_subgroups_cache: list | None = None
        self._special_subgroups_lock = threading.Lock()
        self._conv_to_prim_cache: np.ndarray | None = None
        self._parent_rotations_cache: list[np.ndarray] | None = None
        self._strain_representation_cache = None
        self._lattice_standardization_cache: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}

    def snapshot_for_export(self) -> "IsoDistort":
        """Return an independent export worker for the current committed state.

        Batch preparation mutates Method-2 mode caches candidate by candidate.
        A web request must therefore not run it on the shared interactive
        session after releasing the session lock.  This snapshot owns fresh
        backend wrappers and deep-copies every mutable scientific/session fact
        needed by the export path, so other tabs can keep reading or changing
        the live session without mixing states into the batch.
        """

        if self.structure is None or self.symmetry_info is None:
            raise RuntimeError("parent structure is not loaded")
        worker = type(self)()
        worker.structure = self.structure.copy()
        worker.symmetry_info = copy.deepcopy(self.symmetry_info)
        worker.structure_path = self.structure_path
        worker.subgroups = copy.deepcopy(self.subgroups)
        worker._unresolved_method3_embedding_keys = copy.deepcopy(
            self._unresolved_method3_embedding_keys
        )
        worker.distortion_types = list(self.distortion_types)
        worker.distortion_scope = self._copy_distortion_scope(self.distortion_scope)
        worker._restore_distortion_state(copy.deepcopy(self._snapshot_distortion_state()))
        worker._special_subgroups_cache = copy.deepcopy(self._special_subgroups_cache)
        worker._conv_to_prim_cache = copy.deepcopy(self._conv_to_prim_cache)
        worker._parent_rotations_cache = copy.deepcopy(self._parent_rotations_cache)
        worker._strain_representation_cache = copy.deepcopy(
            self._strain_representation_cache
        )
        worker._lattice_standardization_cache = copy.deepcopy(
            self._lattice_standardization_cache
        )
        return worker

    # ================================================================
    # 阶段一：结构输入与对称识别
    # ================================================================

    def load_structure(self, cif_path: str | Path) -> Structure:
        """
        步骤1-3：加载结构文件，识别对称性

        支持 CIF / VASP POSCAR / xyz（按扩展名自动识别，见 read_structure）。

        Args:
            cif_path: 结构文件路径（CIF / POSCAR / xyz）

        Returns:
            Structure: 加载后的结构对象
        """
        path = Path(cif_path)
        structure = (
            read_cif(path)
            if path.suffix.lower() == ".cif"
            else read_structure(path)
        )
        self._install_parent_structure(structure, path.resolve())

        sg_num = self.symmetry_info["space_group_number"]
        sg_sym = self.symmetry_info["space_group_symbol"]
        n_atoms = len(self.structure)

        print(t("load.done", sg=sg_num, sym=sg_sym, n=n_atoms))
        return self.structure

    def set_structure(self, structure: Structure) -> Structure:
        """直接设置 Structure 对象"""
        self._install_parent_structure(structure, None)
        return self.structure

    def _install_parent_structure(
        self,
        structure: Structure,
        source_path: Path | None,
    ) -> None:
        """Prepare and publish a new parent without exposing mixed session state.

        Parent-setting canonicalization can legitimately fail after it has
        inspected and temporarily replaced ``self.structure``.  Derived
        Method 1--4 state belongs to the previously committed parent until the
        complete new parent, symmetry result, orbit metadata, and source
        labels have all been prepared successfully.  Roll back those primary
        fields on failure; clear derived state only at the final commit point.
        """

        symmetry_info = self._sym_val.validate(structure)
        previous = (
            self.structure,
            self.structure_path,
            self.symmetry_info,
            self._findsym,
        )
        try:
            self.structure = structure
            self.structure_path = source_path
            self.symmetry_info = symmetry_info
            if source_path is not None:
                # Attach source labels before a possible whole-structure
                # setting conversion so they follow their physical orbits.
                self.parent_wyckoff_display()
            self._attach_standard_parent_orbits()
            if source_path is not None:
                # Parse and attach source-CIF atom-site labels/order at load
                # time. Web status rendering must not be the first mutation.
                self.parent_wyckoff_display()
        except Exception:
            (
                self.structure,
                self.structure_path,
                self.symmetry_info,
                self._findsym,
            ) = previous
            raise
        self._reset_derived_state()

    def parent_wyckoff_display(self) -> list[str]:
        """官网页头 Wyckoff 行：优先按母相 CIF 位点顺序与标签，否则用对称分析。"""
        if self.structure is None or not self.symmetry_info:
            return []
        return parent_wyckoff_display(
            self.structure,
            self.symmetry_info.get("wyckoff_sites") or [],
            self.structure_path,
        )

    def _canonicalize_parent_setting(self) -> None:
        """Move the complete parent into the conventional standard setting.

        Each physical orbit is temporarily represented by its own dummy
        species.  The symmetry standardizer can then change basis, origin, and
        cell multiplicity without collapsing independent occurrences of the
        same Wyckoff letter.  Real (including disordered) species and source
        labels are restored orbit by orbit afterwards.
        """
        if self.structure is None or not self.symmetry_info:
            raise RuntimeError("parent structure is not loaded")

        original = self.structure
        original_sites = list(self.symmetry_info.get("wyckoff_sites") or [])
        if not original_sites:
            raise IsodistortError("cannot standardize a parent without Wyckoff orbits")

        atom_tokens: list[DummySpecies | None] = [None] * len(original)
        token_to_site: dict[str, dict] = {}
        token_to_species: dict[str, object] = {}
        for orbit_index, site in enumerate(original_sites):
            token = DummySpecies(f"X{orbit_index}")
            token_name = str(token)
            representative = int(site["representative_index"])
            token_to_site[token_name] = site
            token_to_species[token_name] = original[representative].species
            for atom_index in site.get("equivalent_indices") or []:
                atom_tokens[int(atom_index)] = token
        if any(token is None for token in atom_tokens):
            raise IsodistortError(
                "cannot standardize parent: physical orbits do not cover every atom"
            )

        tagged = Structure(
            original.lattice,
            atom_tokens,
            original.frac_coords,
            coords_are_cartesian=False,
        )
        analyzer = SpacegroupAnalyzer(
            tagged,
            symprec=self._sym_val.symprec_angstrom,
            angle_tolerance=self._sym_val.angle_tolerance_degrees,
        )
        expected_sg = int(self.symmetry_info["space_group_number"])
        if int(analyzer.get_space_group_number()) != expected_sg:
            raise IsodistortError(
                "orbit-preserving parent standardization changed the space group"
            )
        tagged_standard = analyzer.get_conventional_standard_structure(
            international_monoclinic=False
        )
        standard_tokens = [site.species_string for site in tagged_standard]
        try:
            standard_species = [token_to_species[token] for token in standard_tokens]
        except KeyError as exc:
            raise IsodistortError(
                f"parent standardization returned an unknown orbit token {exc.args[0]!r}"
            ) from exc
        standardized = Structure(
            tagged_standard.lattice,
            standard_species,
            tagged_standard.frac_coords,
            coords_are_cartesian=False,
        )
        standardized_info = self._sym_val.validate(standardized)
        if int(standardized_info["space_group_number"]) != expected_sg:
            raise IsodistortError(
                "standardized parent no longer has the detected space group"
            )

        # Carry source-CIF names and order across the cell transformation.
        for standard_site in standardized_info.get("wyckoff_sites") or []:
            tokens = {
                standard_tokens[int(index)]
                for index in standard_site.get("equivalent_indices") or []
            }
            if len(tokens) != 1:
                raise IsodistortError(
                    "standardized physical orbit combines distinct source orbits"
                )
            source_site = token_to_site[next(iter(tokens))]
            # Primitive, conventional, and translation-supercell inputs can
            # contain different numbers of representatives of the same
            # crystallographic orbit.  Their multiplicities must scale by the
            # single whole-cell index, rather than remain numerically equal.
            if (
                int(source_site["multiplicity"]) * len(standardized)
                != int(standard_site["multiplicity"]) * len(original)
            ):
                raise IsodistortError(
                    "standardized physical orbit has an inconsistent cell-index scaling"
                )
            for field in ("display_label", "display_order"):
                if field in source_site:
                    standard_site[field] = source_site[field]

        self.structure = standardized
        self.symmetry_info = standardized_info

    def _attach_standard_parent_orbits(self) -> None:
        """Attach FINDSYM representatives in one internally consistent setting.

        When FINDSYM detects a different basis or origin, convert the complete
        structure first and rerun the analysis.  Never attach coordinates from
        one setting to atoms in another.
        """
        if self.structure is None or not self.symmetry_info:
            return
        sites = self.symmetry_info.get("wyckoff_sites") or []
        if self._findsym is None:
            self._findsym = FindsymWrapper()
        setting_change = None
        try:
            metadata = self._findsym.standardize_wyckoff_orbits(
                self.structure,
                sites,
                expected_space_group=int(self.symmetry_info["space_group_number"]),
            )
        except FindsymSettingMismatchError as exc:
            setting_change = exc.result
            try:
                self._canonicalize_parent_setting()
                sites = self.symmetry_info.get("wyckoff_sites") or []
                metadata = self._findsym.standardize_wyckoff_orbits(
                    self.structure,
                    sites,
                    expected_space_group=int(
                        self.symmetry_info["space_group_number"]
                    ),
                )
            except (IsodistortError, OSError, ValueError) as conversion_exc:
                raise IsodistortError(
                    "FINDSYM requires a parent-setting conversion, but the "
                    f"complete structure could not be standardized: {conversion_exc}"
                ) from conversion_exc
        except (IsodistortError, OSError) as exc:
            self.symmetry_info["parent_orbit_standardization"] = {
                "status": "unavailable",
                "source": "findsym",
                "reason": str(exc),
            }
            return
        for site, standard in zip(sites, metadata, strict=True):
            site.update(standard)
        standardization = {
            "status": "canonicalized" if setting_change is not None else "available",
            "source": "findsym",
        }
        if setting_change is not None:
            standardization["input_basis_vectors"] = [
                list(vector) for vector in setting_change.basis_vectors
            ]
            standardization["input_origin"] = list(setting_change.origin)
        self.symmetry_info["parent_orbit_standardization"] = standardization

    def _require_verified_parent_setting(self) -> None:
        """Reject an explicitly unverified ISO parent frame before calculation.

        Loading remains available for diagnostics after a FINDSYM runtime
        failure. Such a session must not combine uploaded coordinates with
        ISO's standard-setting k vectors, bases, or Wyckoff mode data. Legacy
        private fixtures without this metadata keep their existing contract;
        its absence is not a certification of their coordinate setting.
        """
        standardization = (getattr(self, "symmetry_info", None) or {}).get(
            "parent_orbit_standardization"
        )
        if not isinstance(standardization, dict) or standardization.get("status") != "unavailable":
            return
        reason = str(standardization.get("reason") or "FINDSYM standardization unavailable")
        raise IsodistortError(
            "Parent setting has not been verified against the ISO standard frame. "
            "Repair the FINDSYM runtime and reload the parent structure before "
            f"running searches, calculating modes or exporting subgroups. Details: {reason}"
        )

    def _reset_derived_state(self) -> None:
        """加载新结构后清空所有派生状态。"""
        self.subgroups = []
        self._unresolved_method3_embedding_keys = set()
        self._selected_subgroup = None
        self.phase_path = None
        self.distortion_modes = []
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.mode_displacements_sc = {}
        self._mode_label_overrides = {}
        self.distorted_structure = None
        self._mode_cache_key = None
        self._generated_structure_cache_key = None
        self._special_subgroups_cache = None
        self._conv_to_prim_cache = None
        self._parent_rotations_cache = None
        self._strain_representation_cache = None
        self._lattice_standardization_cache = {}

    def set_subgroup_candidates(self, subgroups: list[SubgroupInfo]) -> None:
        """Replace the default Method-2 candidate pool through the public API.

        Prefer passing ``candidates=...`` directly to :meth:`search_method_2`
        when a caller keeps several Method result tables at the same time.
        This setter exists for interactive flows that intentionally make one
        table the session default.
        """
        self.subgroups = list(subgroups)

    def clear_selected_modes(self) -> None:
        """Clear the selected path and all derived mode/distortion state."""
        self._selected_subgroup = None
        self.phase_path = None
        self.distortion_modes = []
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.mode_displacements_sc = {}
        self._mode_label_overrides = {}
        self.distorted_structure = None
        self._mode_cache_key = None
        self._generated_structure_cache_key = None

    def list_isotropy_cache(self) -> list[IsotropyCacheEntry]:
        """Return generated ISOTROPY cache entries without exposing the backend."""
        return list_isotropy_cache(self._iso)

    def delete_isotropy_cache(self, names: list[str]) -> dict:
        """Delete named generated cache files through the public API."""
        return delete_isotropy_cache(self._iso, names)

    @staticmethod
    def _subgroup_identity(subgroup: SubgroupInfo | None) -> tuple | None:
        """Stable path identity; local display indices are deliberately excluded."""
        if subgroup is None:
            return None

        def _matrix(values) -> tuple:
            return tuple(
                tuple(round(float(value), 10) for value in row)
                for row in (values or [])
            )

        return (
            int(subgroup.parent_sg),
            int(subgroup.space_group_number),
            int(subgroup.subgroup_index),
            int(subgroup.size),
            str(subgroup.k_point_label or ""),
            tuple(str(value) for value in (subgroup.k_coordinates or [])),
            str(subgroup.irrep_label or ""),
            str(subgroup.opd_symbol or ""),
            str(subgroup.opd_dir_raw or ""),
            tuple(round(float(value), 10) for value in (subgroup.opd_vector or [])),
            tuple(str(value) for value in (subgroup.k_parameters or [])),
            _matrix(subgroup.basis_vectors),
            tuple(round(float(value), 10) for value in (subgroup.origin or [])),
        )

    @staticmethod
    def _method3_embedding_guard_key(subgroup: SubgroupInfo) -> tuple:
        """Stable full route identity for rejecting incompatible cached state.

        Stage-A rows carry a content-addressed ID derived from their exact
        Seitz subgroup.  That ID identifies the affine subgroup only: several
        IR/OPD/k routes can stabilize the same embedded subgroup and must not
        share mode amplitudes or generated structures.  Keep the affine ID as
        one field inside the complete route key, whose remaining components
        are serialized as exact Fractions (or unmodified symbolic tokens),
        never as a display index or rounded float.
        """
        embedding_id = str(
            getattr(subgroup, "_method3_embedding_id", "") or ""
        ).strip()
        context_kind = (
            "exact_fixed_space"
            if str(getattr(subgroup, "_method3_route_resolution", ""))
            == "exact_fixed_space"
            else "single_irrep"
        )

        def scalar(value) -> tuple:
            try:
                fraction = as_fraction(value)
            except (TypeError, ValueError):
                return ("symbol", str(value))
            return ("fraction", fraction.numerator, fraction.denominator)

        def matrix(values) -> tuple:
            return tuple(
                tuple(scalar(value) for value in row)
                for row in (values or [])
            )

        def canonical_child_origin() -> tuple[Fraction, Fraction, Fraction]:
            basis = rational_matrix(
                subgroup.basis_vectors or identity_matrix()
            )
            primitive = multiply(
                centering_primitive_matrix(
                    _centering_letter(int(subgroup.space_group_number))
                ),
                basis,
            )
            primitive_inverse = inverse(primitive)
            origin = tuple(
                as_fraction(value)
                for value in (subgroup.origin or (0, 0, 0))
            )
            coefficients = tuple(
                sum(
                    origin[row] * primitive_inverse[row][column]
                    for row in range(3)
                )
                for column in range(3)
            )
            reduced = tuple(
                value - (value.numerator // value.denominator)
                for value in coefficients
            )
            return tuple(
                sum(
                    reduced[row] * primitive[row][column]
                    for row in range(3)
                )
                for column in range(3)
            )

        return (
            "exact-route-identity-v3",
            embedding_id or None,
            context_kind,
            int(subgroup.parent_sg),
            int(subgroup.space_group_number),
            int(subgroup.subgroup_index),
            int(subgroup.size),
            str(subgroup.space_group_symbol or ""),
            str(subgroup.k_point_label or ""),
            tuple(scalar(value) for value in (subgroup.k_coordinates or [])),
            str(subgroup.irrep_label or ""),
            str(subgroup.opd_symbol or ""),
            str(subgroup.opd_dir_raw or ""),
            tuple(scalar(value) for value in (subgroup.opd_vector or [])),
            tuple(scalar(value) for value in (subgroup.k_parameters or [])),
            matrix(subgroup.basis_vectors),
            tuple(scalar(value) for value in canonical_child_origin()),
        )

    @staticmethod
    def _copy_distortion_scope(scope: dict | None) -> dict[str, list[str]]:
        """Copy the public Types-panel scope into an immutable-use snapshot."""

        copied: dict[str, list[str]] = {}
        for type_name, values in (scope or {}).items():
            if isinstance(values, str):
                copied[str(type_name)] = [values]
            else:
                copied[str(type_name)] = [str(value) for value in (values or [])]
        return copied

    @staticmethod
    def _canonical_distortion_types(distortion_types) -> tuple[str, ...]:
        """Canonicalize a semantic type set independently of checkbox order."""

        selected = set(normalize_distortion_types(distortion_types))
        return tuple(type_name for type_name in DISTORTION_TYPES if type_name in selected)

    def _mode_context_key(
        self,
        subgroup: SubgroupInfo,
        nmod: int,
        distortion_types,
        distortion_scope: dict | None = None,
    ) -> tuple:
        """Exact cache identity for modes and structures generated from them."""

        types = self._canonical_distortion_types(distortion_types)
        source_scope = (
            getattr(self, "distortion_scope", {})
            if distortion_scope is None
            else distortion_scope
        )
        scope = self._copy_distortion_scope(source_scope)
        # Store resolved physical species rather than the UI spelling (``*``
        # versus an explicit complete list).  Every enabled atom-dependent
        # type is included; in particular this binds the displacive cache to
        # its per-species scope as required by the Distortion page contract.
        resolved_scope = tuple(
            (
                type_name,
                tuple(sorted(self._scope_species(type_name, scope=scope))),
            )
            for type_name in types
            if type_name != "strain"
        )
        return (
            "distortion-mode-context-v1",
            self._method3_embedding_guard_key(subgroup),
            int(nmod),
            types,
            resolved_scope,
        )

    def _mark_generated_structure_context(self) -> None:
        """Bind the current generated structure to the modes that created it."""

        self._generated_structure_cache_key = getattr(
            self, "_mode_cache_key", None
        )

    # ================================================================
    # 畸变类型作用域（对齐官网 Types 面板的 per-species 复选框）
    # ================================================================

    def species(self) -> list[str]:
        """当前结构包含的元素符号（去重、排序）。"""
        if getattr(self, "structure", None) is None:
            return []
        return sorted({s.species_string for s in self.structure})

    def set_distortion_scope(self, scope: dict | None) -> None:
        """
        设置各畸变类型的作用域物种（官网 Displacive/Occupational/Magnetic/
        Rotational 行内的 all/none/Eu/Al 复选框）。

        Args:
            scope: {类型名: "*"（全部）或 ["Eu", "Al", ...]}；
                不在 scope 中的类型默认作用于全部物种
        """
        self.distortion_scope = {}
        if not scope:
            return
        valid = set(DISTORTION_TYPES)
        for tp, val in scope.items():
            if tp not in valid:
                continue
            if isinstance(val, str) and val.strip().lower() in ("*", "all", "全部"):
                self.distortion_scope[tp] = ["*"]
            elif isinstance(val, (list, tuple)):
                species = [str(s) for s in val if str(s) != "*"]
                self.distortion_scope[tp] = species if species else ["*"]

    def set_distortion_types(self, distortion_types) -> None:
        """设置当前考虑的畸变类型（对齐官网 Types 面板复选框）。"""
        self.distortion_types = normalize_distortion_types(distortion_types)

    def _iso_kpoint_raw(self, k_point_label: str):
        """iso 原始 k 点信息（未应用官网显示覆盖）。"""
        for kp in self._iso.list_k_points(self.symmetry_info["space_group_number"]):
            if kp.label == k_point_label.strip():
                return kp
        return None

    def _resolve_iso_kparams(self, k_point_label: str,
                             k_parameters: list | None) -> list | None:
        """官网 UI 参数 → iso KVALUE 参数；无参数或未收录时原样返回。"""
        if not k_parameters:
            return None
        iso_kp = self._iso_kpoint_raw(k_point_label)
        if iso_kp is None:
            return list(k_parameters)
        return official_kparams_to_iso(
            self.symmetry_info["space_group_number"],
            k_point_label,
            k_parameters,
            iso_kp,
        )

    def _is_parametric_kpoint(self, k_point_label: str) -> bool:
        iso_kp = self._iso_kpoint_raw(k_point_label)
        return bool(iso_kp and iso_kp.parameters)

    def _filter_subgroups_for_search(self,
                                     subgroups: list[SubgroupInfo],
                                     k_point_label: str,
                                     official_kparams: list | None) -> list[SubgroupInfo]:
        """按 Distortion Types + 物种作用域过滤子群（对齐官网 Search 阶段）。"""
        types = normalize_distortion_types(self.distortion_types)
        mode_types = [tp for tp in types if tp in ("displacive", "rotational", "occupational")]
        if not mode_types:
            return subgroups

        if self._is_parametric_kpoint(k_point_label):
            species = self._union_scope_species(types)
            active = self._smodes.active_irreps(
                self.structure,
                self.symmetry_info["space_group_number"],
                self.symmetry_info["wyckoff_sites"],
                k_point_label,
                self._resolve_iso_kparams(k_point_label, official_kparams),
                species_filter=species if species else None,
            )
            if active is not None:
                return [sg for sg in subgroups if sg.irrep_label in active]
            return subgroups

        # 特殊 k 点：用 BUSH 探测各子群是否在作用域 Wyckoff 上有位移模式
        if not any(tp in types for tp in ("displacive", "rotational")):
            return subgroups
        bush_types = [tp for tp in types if tp in ("displacive", "rotational")]
        letters = self._letters_for_species(self._union_scope_species(bush_types))
        if not letters:
            return []
        parent_sg = self.symmetry_info["space_group_number"]
        kept: list[SubgroupInfo] = []
        for sg in subgroups:
            try:
                modes = self._iso.calc_distortion_modes(parent_sg, sg, letters)
            except IsodistortError:
                continue
            if modes:
                kept.append(sg)
        return kept

    def _tag_official_kparams(self, subgroups: list[SubgroupInfo],
                              official_kparams: list | None) -> None:
        """子群对象上保留官网参数（供界面显示），并刷新 k 坐标 / k-active。"""
        if not official_kparams:
            return
        tagged = list(official_kparams)
        for sg in subgroups:
            sg.k_parameters = tagged
            coords = official_special_k_coords(
                int(sg.parent_sg or 0),
                sg.k_point_label or "",
                sg.k_coordinates or [],
                tagged,
            )
            if coords:
                sg.k_coordinates = coords
            sg.k_active_raw = format_k_active(
                sg.opd_dir_raw or "",
                sg.k_coordinates or ["0", "0", "0"],
                sg.parent_sg or None,
                None,
            )

    def _scope_species(
        self,
        type_name: str,
        *,
        scope: dict | None = None,
    ) -> set[str]:
        """某畸变类型作用域内的物种集合（未设置时默认全部物种）。"""
        all_species = set(self.species())
        source = self.distortion_scope if scope is None else scope
        val = source.get(type_name)
        if not val or "*" in val:
            return all_species
        return {s for s in val if s in all_species}

    def _letters_for_species(self, species_set) -> list[str]:
        """物种集合对应的 Wyckoff 位置字母（去重）。"""
        if not self.symmetry_info:
            return []
        letters: list[str] = []
        for site in self.symmetry_info["wyckoff_sites"]:
            if site["species"] in species_set and site["wyckoff_letter"] not in letters:
                letters.append(site["wyckoff_letter"])
        return letters

    def _orbit_ids_for_species(self, species_set) -> list[str]:
        """Stable physical Wyckoff-orbit identities for the selected species."""
        if not self.symmetry_info:
            return []
        orbit_ids: list[str] = []
        for site in self.symmetry_info["wyckoff_sites"]:
            orbit_id = str(site.get("orbit_id") or "")
            if site["species"] in species_set and orbit_id and orbit_id not in orbit_ids:
                orbit_ids.append(orbit_id)
        return orbit_ids

    # ================================================================
    # Method 1 下拉数据（对齐官网：可达子群空间群 + conventional/primitive lattice）
    # ================================================================

    def space_group_preferences(self) -> str:
        """官网页头 “Default space-group preferences: ...” 行。

        本地 iso 二进制固定采用国际标准取位（官网默认值），与实际计算行为一致，
        因此仅提供只读说明，不提供可交互的偏好面板（自定义取位无法被本地引擎生效）。
        """
        return ("monoclinic axes a(b)c, monoclinic cell choice 1, "
                "orthorhombic axes abc, origin choice 2, hexagonal axes, "
                "SSG standard setting")

    def _ensure_special_subgroups(self) -> list:
        """枚举全部特殊 k 点子群（线程安全缓存，Method 1 下拉与搜索共用）。"""
        if self.structure is None:
            raise RuntimeError(t("err.load_first"))
        with self._special_subgroups_lock:
            if self._special_subgroups_cache is None:
                sg_num = self.symmetry_info["space_group_number"]
                self._special_subgroups_cache = self._iso.enumerate_all_special_subgroups(
                    sg_num
                )
        return self._special_subgroups_cache

    def _filter_method1_by_types(self, items, distortion_types) -> list:
        """Keep Method 1 rows that the official site would list for the selected types.

        Each requested physical representation is tested independently and the
        union is returned.  Displacive/rotational activity comes from smodes;
        occupational activity requires a generated pattern that independently
        validates as the requested subgroup; homogeneous strain uses the
        symmetric metric-tensor representation ``E -> R.T E R``.  An unknown
        backend result never masquerades as physical activity.
        """
        types = normalize_distortion_types(distortion_types or self.distortion_types)
        want_mag = "magnetic" in types
        want_disp = any(tp in types for tp in ("displacive", "rotational"))
        want_strain = "strain" in types
        want_occupational = "occupational" in types
        if not items:
            return items

        active_by_k: dict[str, set[str] | None] = {}
        if want_disp and self.structure is not None and self.symmetry_info:
            species = self._union_scope_species(types)
            for kp in {it.subgroup.k_point_label for it in items}:
                if not kp:
                    active_by_k[kp] = None
                    continue
                active_by_k[kp] = self._smodes.active_irreps(
                    self.structure,
                    self.symmetry_info["space_group_number"],
                    self.symmetry_info["wyckoff_sites"],
                    kp,
                    None,
                    species_filter=species if species else None,
                )

        if want_disp:
            unresolved = sorted(
                label for label, active in active_by_k.items() if active is None
            )
            if unresolved:
                raise RuntimeError(
                    "Could not establish displacive/rotational irrep activity "
                    "from smodes for k point(s): " + ", ".join(unresolved)
                )

        parent_sg_n = 0
        if self.symmetry_info:
            parent_sg_n = int(self.symmetry_info["space_group_number"])

        occupational_generator = (
            OccupationalModeGenerator() if want_occupational else None
        )
        occupational_scope = (
            self._scope_species("occupational") if want_occupational else set()
        )

        kept = []
        for item in items:
            ir = (item.subgroup.irrep_label or "").strip()
            if ir.startswith("m"):
                if want_mag:
                    kept.append(item)
                continue

            is_active = False
            if want_disp:
                active = active_by_k.get(item.subgroup.k_point_label)
                is_active = active is not None and ir in active
            if want_strain and self._keep_strain_only_irrep(item, parent_sg_n):
                is_active = True
            if (
                occupational_generator is not None
                and self.structure is not None
                and self.symmetry_info
                and occupational_scope
            ):
                try:
                    occupational_modes = occupational_generator.generate(
                        self.structure,
                        self.symmetry_info["wyckoff_sites"],
                        item.subgroup,
                        occupational_scope,
                    )
                except Exception:  # noqa: BLE001 - unsupported route is inactive
                    occupational_modes = []
                if any(bool(getattr(mode, "validated", False))
                       for mode in occupational_modes):
                    is_active = True
            if is_active:
                kept.append(item)
        return kept

    def _filter_method3_routes_by_types(self, items, distortion_types) -> list:
        """Filter each Method 3 embedding through its known active routes.

        Method 3 rows are affine subgroup embeddings, while the local search
        engine discovers them through one or more single-IR routes.  Activity
        is therefore a property of each route, not of only the representative
        stored in ``item.subgroup``.  Parameter values are part of the k-point
        identity and must be supplied (after the existing official-to-backend
        parameter conversion) to smodes when checking displacement activity.

        The filtering follows the capabilities already exposed elsewhere:
        displacement/rotation activity comes from smodes, strain uses the
        symmetric-tensor fixed-space rule, occupational activity requires a
        generated pattern independently validated as the requested subgroup,
        and magnetic activity is limited to magnetic (``m*``) irreps.  A
        failed smodes probe remains unresolved and is not promoted to active.
        """
        types = normalize_distortion_types(distortion_types or self.distortion_types)
        want_displacement = any(
            name in types for name in ("displacive", "rotational")
        )
        want_strain = "strain" in types
        want_occupational = "occupational" in types
        want_magnetic = "magnetic" in types
        if not items:
            return items

        routes = [
            route
            for item in items
            for route in (
                item.routes if item.routes is not None else [item.subgroup]
            )
        ]
        active_by_k: dict[tuple[str, tuple[str, ...]], set[str] | None] = {}
        if (
            want_displacement
            and self.structure is not None
            and self.symmetry_info
        ):
            displacement_species: set[str] = set()
            for name in ("displacive", "rotational"):
                if name in types:
                    displacement_species |= self._scope_species(name)
            for route in routes:
                k_label = str(route.k_point_label or "")
                k_parameters = tuple(str(value) for value in (route.k_parameters or []))
                key = (k_label, k_parameters)
                if key in active_by_k:
                    continue
                if not k_label:
                    active_by_k[key] = None
                    continue
                active_by_k[key] = self._smodes.active_irreps(
                    self.structure,
                    self.symmetry_info["space_group_number"],
                    self.symmetry_info["wyckoff_sites"],
                    k_label,
                    self._resolve_iso_kparams(
                        k_label, list(k_parameters)
                    ) if k_parameters else None,
                    species_filter=(
                        displacement_species if displacement_species else None
                    ),
                )

        occupational_generator = (
            OccupationalModeGenerator() if want_occupational else None
        )
        occupational_scope = (
            self._scope_species("occupational") if want_occupational else set()
        )
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)

        def route_is_active(route: SubgroupInfo) -> bool:
            irrep = str(route.irrep_label or "").strip()
            if irrep.startswith("m"):
                return want_magnetic

            if want_displacement:
                key = (
                    str(route.k_point_label or ""),
                    tuple(str(value) for value in (route.k_parameters or [])),
                )
                active = active_by_k.get(key)
                if active is not None and irrep in active:
                    return True

            if want_strain:
                probe = Method1ResultItem(
                    subgroup=route,
                    crystal_system=_sg_to_crystal_system(route.space_group_number),
                    is_maximal=route.is_maximal,
                )
                if self._keep_strain_only_irrep(probe, parent_sg):
                    return True

            if (
                occupational_generator is not None
                and self.structure is not None
                and self.symmetry_info
                and occupational_scope
            ):
                try:
                    modes = occupational_generator.generate(
                        self.structure,
                        self.symmetry_info["wyckoff_sites"],
                        route,
                        occupational_scope,
                    )
                except Exception:  # noqa: BLE001 - unsupported embedding is inactive
                    modes = []
                if any(bool(getattr(mode, "validated", False)) for mode in modes):
                    return True
            return False

        kept = []
        for item in items:
            candidate_routes = (
                item.routes if item.routes is not None else [item.subgroup]
            )
            if item.routes == []:
                resolution = getattr(item, "route_resolution", "")
                if resolution == "exact_fixed_space":
                    # Search-stage exact character analysis already proved
                    # reachability in the requested strain/displacive direct
                    # sum.  There is no single route to re-filter here.
                    kept.append(item)
                    continue
                if resolution == "affine_only_unresolved_coupled_route":
                    # Explicit diagnostics remain visible but nonselectable.
                    kept.append(item)
                    continue
            active_routes = [route for route in candidate_routes if route_is_active(route)]
            if not active_routes:
                continue
            # The representative drives Method 2 after the user clicks a row;
            # it must itself be one of the surviving routes.
            item.routes = active_routes
            item.subgroup = active_routes[0]
            item.basis = [list(row) for row in (item.subgroup.basis_vectors or [])]
            kept.append(item)

        # Removed embeddings leave holes in their display handles.  Reindex
        # only representatives; route indices are not cross-table identities.
        for index, item in enumerate(kept):
            item.subgroup.index = index
        return kept

    def _keep_strain_only_irrep(self, item, parent_sg: int) -> bool:
        """Whether a route is active in the homogeneous symmetric strain tensor.

        Strain is a Gamma-point, translation-trivial, inversion-even second-rank
        tensor.  A larger ``dim Fix(H)`` is only necessary: every extra tensor
        might still have a larger stabilizer ``K > H``.  The exact criterion is
        therefore ``Fix(H) != 0`` and ``Stab_G(Fix(H)) = H``, evaluated in the
        actual embedded child orientation with rational matrices.
        """
        sg = item.subgroup
        klab = (sg.k_point_label or "").strip().upper()
        if klab not in {"GM", "G", "Γ", "GAMMA"}:
            return False
        ir = (sg.irrep_label or "").strip().upper()
        if ir.endswith("-"):
            return False
        if int(getattr(sg, "size", 1) or 1) != 1:
            return False
        if not parent_sg or self.structure is None:
            return False
        try:
            parent, strain = self._parent_strain_context()
            if int(parent.space_group_number) != int(parent_sg):
                return False
            embedding = embedding_from_identity(sg, parent)
            child_rotations = tuple(dict.fromkeys(
                operation.rotation for operation in embedding.operations
            ))
            return analyze_fixed_space(strain, child_rotations).is_reachable
        except (ArithmeticError, AttributeError, KeyError, TypeError, ValueError):
            # Unknown is not physical activity.  Alternate/unsupported settings
            # must first pass exact Seitz reconstruction rather than being
            # promoted by a crystal-system or fixed-dimension heuristic.
            return False

    def _parent_strain_context(self):
        """Build the exact parent affine group and strain representation once."""

        if self.structure is None:
            raise RuntimeError("parent structure is not loaded")
        if getattr(self, "_strain_representation_cache", None) is None:
            cfg = getattr(self, "cfg", None) or get_config()
            parent = parent_affine_group(
                self.structure,
                symprec=cfg.symmetry_cartesian_tolerance_angstrom,
                angle_tolerance_degrees=cfg.symmetry_angle_tolerance_degrees,
            )
            rotations = tuple(
                dict.fromkeys(operation.rotation for operation in parent.operations)
            )
            group = FiniteGroup.from_operation(
                rotations,
                identity_matrix(),
                multiply,
                name="parent crystallographic point group",
            )
            vector = RationalRepresentation(
                group,
                rotations,
                name="fractional polar vector",
            )
            strain = symmetric_square_representation(
                vector,
                name="homogeneous symmetric strain",
            )
            self._strain_representation_cache = (parent, strain)
        return self._strain_representation_cache

    def _strain_export_data_for_subgroup(
        self, subgroup: SubgroupInfo
    ) -> StrainExportData | None:
        """Compute ``Fix_H(Sym²V)`` from the exact embedded target subgroup.

        ``embedding_from_identity`` reconstructs and verifies the complete
        Seitz subgroup inside the actual parent group. Its rotations are
        fractional-column operations in the parent axes. The fixed-space core
        combines those operations with the uploaded parent metric and returns
        ISODISTORT parent-lattice-basis engineering strain coordinates.
        """

        # ``normalize_distortion_types([])`` intentionally restores the UI
        # defaults, but export generation must preserve an explicit empty
        # selection: strain data exists only when the current selection itself
        # contains ``strain``.
        configured_types = getattr(self, "distortion_types", None)
        if isinstance(configured_types, str):
            configured_types = [configured_types]
        selected_types = {
            str(value).strip().lower()
            for value in (configured_types or [])
        }
        if "strain" not in selected_types:
            return None
        if self.structure is None:
            raise RuntimeError("parent structure is not loaded")
        parent, _parent_strain = self._parent_strain_context()
        embedding = embedding_from_identity(subgroup, parent)
        fractional_rotations = [
            np.asarray(operation.rotation, dtype=float)
            for operation in embedding.operations
        ]
        result = compute_homogeneous_strain_modes(
            fractional_rotations,
            self.structure.lattice.matrix,
        )
        parent_sg = int(
            (getattr(self, "symmetry_info", None) or {}).get("space_group_number")
            or subgroup.parent_sg
            or 0
        )
        if parent_sg <= 0:
            raise RuntimeError(
                "cannot resolve the parent space group for canonical strain modes"
            )
        macro_cache = getattr(self, "_strain_macro_basis_cache", None)
        if macro_cache is None:
            macro_cache = {}
            self._strain_macro_basis_cache = macro_cache
        macroscopic = macro_cache.get(parent_sg)
        reusable_statuses = {"verified", "partial"}
        if getattr(macroscopic, "status", "unresolved") not in reusable_statuses:
            macroscopic = self._iso.get_macroscopic_tensor_basis(parent_sg)
            if getattr(macroscopic, "status", "unresolved") in reusable_statuses:
                macro_cache[parent_sg] = macroscopic
            else:
                # An unavailable/invalid ISO response is fail-closed evidence
                # for this attempt, not a stable fact about the space group.
                # Leave it uncached so a transient backend failure can recover
                # on the next export without reloading the parent structure.
                macro_cache.pop(parent_sg, None)
        directions = tuple(
            self._iso.list_invariant_directions(parent_sg, subgroup)
        )
        conditioned = replace(
            macroscopic,
            invariant_directions=directions,
        )
        definitions = canonical_strain_definitions_from_iso(conditioned)
        canonical = apply_canonical_strain_basis(result, definitions)
        return StrainExportData(result=canonical)

    def _parent_rotations(self) -> list[np.ndarray]:
        """母相点群旋转矩阵（分数坐标，去重）。"""
        if self._parent_rotations_cache is None:
            if self.structure is None:
                return [np.eye(3)]
            sga = SpacegroupAnalyzer(
                self.structure,
                symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
                angle_tolerance=self.cfg.symmetry_angle_tolerance_degrees,
            )
            uniq: list[np.ndarray] = []
            seen: set[tuple] = set()
            for op in sga.get_symmetry_operations(cartesian=False):
                rot = np.asarray(op.rotation_matrix, dtype=float)
                key = tuple(np.round(rot, 6).flatten())
                if key not in seen:
                    seen.add(key)
                    uniq.append(rot)
            self._parent_rotations_cache = uniq or [np.eye(3)]
        return self._parent_rotations_cache

    @staticmethod
    def _centering_matrix(letter: str) -> np.ndarray:
        """惯用 → 原胞的标准心化矩阵（行向量为惯用坐标下的原胞基矢）。"""
        tables = {
            "P": np.eye(3),
            "I": np.array(
                [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]], dtype=float
            ),
            "F": np.array(
                [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]], dtype=float
            ),
            "A": np.array(
                [[1.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.0, -0.5, 0.5]], dtype=float
            ),
            "B": np.array(
                [[0.5, 0.0, 0.5], [0.0, 1.0, 0.0], [-0.5, 0.0, 0.5]], dtype=float
            ),
            "C": np.array(
                [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]], dtype=float
            ),
            "R": np.array(
                [
                    [2 / 3, 1 / 3, 1 / 3],
                    [-1 / 3, 1 / 3, 1 / 3],
                    [-1 / 3, -2 / 3, 1 / 3],
                ],
                dtype=float,
            ),
        }
        return tables.get((letter or "P")[:1].upper(), np.eye(3)).copy()

    def _parent_centering_matrix(self) -> np.ndarray:
        """母相惯用 → 原胞的标准心化矩阵 T。"""
        letter = _centering_letter(self.symmetry_info["space_group_number"])
        return self._centering_matrix(letter)

    def _conv_to_prim(self) -> np.ndarray:
        """母相惯用格子 -> 原胞格子的变换矩阵 T（L_prim = L_conv @ T）。"""
        if self._conv_to_prim_cache is None:
            # 与官网 Primitive lattice 同一套心化约定（不用 pymatgen 的备选约定）
            self._conv_to_prim_cache = self._parent_centering_matrix()
        return self._conv_to_prim_cache

    def lattice_in_conventional_frame(self, matrix, frame: str = "conventional"
                                      ) -> list[list[float]]:
        """把用户选择的 lattice 基矢换算到 iso 惯用（conventional）坐标系。

        官网下拉的 Conventional lattice 与 Primitive lattice 选项基矢
        均为母相惯用坐标表达（Primitive 选项是同一子格在惯用坐标下的显示，
        如“原胞本身”显示为 (-1/2,1/2,1/2),...），因此 method1_options 返回的
        选项直接以 frame="conventional" 提交即可，无需变换。

        frame="primitive" 仅用于调用方持有“原胞坐标”表达的基矢时换算：
        B_conv = B_prim @ T（T = _conv_to_prim()，见其实现）。
        """
        m = np.asarray(matrix, dtype=float)
        if frame == "primitive":
            m = m @ self._conv_to_prim()
        return m.tolist()

    @staticmethod
    def _same_lattice(a, b) -> bool:
        """两个 3x3 超胞基矢是否生成同一格点（GL(3,Z) 等价）。"""
        tolerance = get_config().lattice_tolerance
        a_arr = np.asarray(a, dtype=float)
        b_arr = np.asarray(b, dtype=float)
        if a_arr.shape != (3, 3) or b_arr.shape != (3, 3):
            return False
        if abs(abs(np.linalg.det(a_arr)) - abs(np.linalg.det(b_arr))) > tolerance:
            return False
        try:
            n = a_arr @ np.linalg.inv(b_arr)
        except np.linalg.LinAlgError:
            return False
        return bool(np.allclose(n, np.round(n), atol=tolerance))

    def _same_lattice_orbit(self, a, b) -> bool:
        """同一格点类：含母相点群旋转轨道（官网 lattice 选项语义）。"""
        b_arr = np.asarray(b, dtype=float)
        for rot in self._parent_rotations():
            if self._same_lattice(a, b_arr @ rot) or self._same_lattice(a, rot @ b_arr):
                return True
        return False

    def _method1_filtered_subgroups(self) -> list:
        """Method 1 下拉与搜索共用的类型过滤后子群列表。"""
        raw = self._ensure_special_subgroups()
        items = [
            Method1ResultItem(
                subgroup=sg,
                crystal_system=_sg_to_crystal_system(sg.space_group_number),
                is_maximal=sg.is_maximal,
            )
            for sg in raw
        ]
        filtered = self._filter_method1_by_types(items, self.distortion_types)
        return [it.subgroup for it in filtered]

    @staticmethod
    def _opd_lattice_sort_key(direction: str) -> tuple:
        """Canonical OPD complexity used by the website's lattice grouping.

        The official selector orders a k/IR's translation lattices by the
        number of active star arms, then by the number and placement of free
        OPD parameters.  This is why, for example, ``(a;0;0;0)`` precedes
        ``(a;a;0;0)`` and the general four-arm direction.  Deriving the order
        from the direction keeps it valid for every parent instead of storing
        the EuAl4 option list.
        """
        body = (direction or "").strip().strip("()")
        tokens = [part.strip() for part in re.split(r"[;,]", body)] if body else []
        active = [token not in {"", "0", "0.0"} for token in tokens]
        parameters = {
            name.lower()
            for token in tokens
            for name in re.findall(r"[A-Za-z]+", token)
        }
        negative_terms = sum(token.startswith("-") for token in tokens if token)
        # For equal activity, prefer arms occurring earlier in the star.
        active_mask = tuple(-int(value) for value in active)
        return (
            sum(active), len(parameters), active_mask,
            negative_terms, body.replace(" ", ""),
        )

    def _method1_option_ranks(self, subgroups: list[SubgroupInfo]) -> list[tuple]:
        """Stable website-like ordering keys for Method-1 lattice classes."""
        k_order: dict[str, int] = {}
        ir_order: dict[tuple[str, str], int] = {}
        ranks: list[tuple] = []
        for sequence, subgroup in enumerate(subgroups):
            k_label = subgroup.k_point_label or ""
            ir_label = subgroup.irrep_label or ""
            k_order.setdefault(k_label, len(k_order))
            ir_order.setdefault((k_label, ir_label), len(ir_order))
            ranks.append((
                k_order[k_label],
                ir_order[(k_label, ir_label)],
                *self._opd_lattice_sort_key(subgroup.opd_dir_raw),
                sequence,
            ))
        return ranks

    def _standardized_subgroup_lattices(
        self, subgroup: SubgroupInfo,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return conventional and primitive subgroup cells in parent axes.

        ``iso`` supplies a valid subgroup basis, but for a centered parent its
        chosen unimodular representative can drift from the conventional-cell
        representative used by ISODISTORT's web layer.  Reconstructing a
        generic orbit of the subgroup and asking spglib to standardize it
        performs the crystallographic cell reduction directly.  No parent,
        irrep, or expected option label is stored here.
        """
        conventional_fallback = np.asarray(subgroup.basis_vectors, dtype=float)
        primitive_fallback = (
            self._centering_matrix(_centering_letter(subgroup.space_group_number))
            @ conventional_fallback
        )
        key = (
            int(subgroup.space_group_number),
            tuple(np.round(conventional_fallback, 10).flat),
        )
        cached = self._lattice_standardization_cache.get(key)
        if cached is not None:
            return cached
        if self.structure is None:
            return conventional_fallback, primitive_fallback

        conventional = conventional_fallback
        primitive = primitive_fallback
        try:
            parent_cart = np.asarray(self.structure.lattice.matrix, dtype=float)
            parent_cart_inv = np.linalg.inv(parent_cart)
            subgroup_lattice = Lattice(conventional_fallback @ parent_cart)
            # A general-position probe carries exactly the requested subgroup
            # symmetry while avoiding any chemistry- or parent-specific data.
            probe = Structure.from_spacegroup(
                int(subgroup.space_group_number),
                subgroup_lattice,
                ["H"],
                [[0.1234567, 0.2345678, 0.3456789]],
            )
            cell = (
                np.asarray(probe.lattice.matrix, dtype=float),
                np.asarray(probe.frac_coords, dtype=float),
                [1] * len(probe),
            )
            conv_cell = spglib.standardize_cell(
                cell,
                to_primitive=False,
                no_idealize=True,
                symprec=self.cfg.affine_exact_cartesian_tolerance_angstrom,
            )
            prim_cell = spglib.standardize_cell(
                cell,
                to_primitive=True,
                no_idealize=True,
                symprec=self.cfg.affine_exact_cartesian_tolerance_angstrom,
            )
            if conv_cell is not None:
                candidate = np.asarray(conv_cell[0], dtype=float) @ parent_cart_inv
                if self._same_lattice_orbit(candidate, conventional_fallback):
                    conventional = candidate
            if prim_cell is not None:
                candidate = np.asarray(prim_cell[0], dtype=float) @ parent_cart_inv
                if self._same_lattice_orbit(candidate, primitive_fallback):
                    primitive = candidate
        except (ValueError, TypeError, np.linalg.LinAlgError):
            # A malformed/metric-incompatible third-party structure should not
            # remove an otherwise valid option supplied by iso.
            pass

        result = (conventional, primitive)
        self._lattice_standardization_cache[key] = result
        return result

    def _distinct_lattices(self, bases,
                           preferred_labels: list[str] | None = None,
                           standardized_kind: str | None = None,
                           ) -> list[dict]:
        """从一组超胞基矢提取去重后的 lattice 选项。

        去重：GL(3,Z) 格点等价 ∪ 母相点群旋转轨道（同一选项含点群相关格子）。
        显示代表：若提供 ``preferred_labels`` 且该类与之轨道等价，采用该标签；
        否则保留首次出现的 iso ``basis_raw``（不排序行、不取“最简范数”）。
        顺序：preferred 列表顺序优先，其余按首次出现顺序接在后面。
        注意：不要用按母相硬编码的官网快照表填充 preferred；应对齐 iso 算法输出。
        """
        classes: list[dict] = []
        for item in bases:
            if isinstance(item, tuple):
                arr = np.asarray(item[0], dtype=float)
                raw_label = (item[1] or "").strip()
                rank = item[2] if len(item) > 2 else (len(classes),)
                subgroup = item[3] if len(item) > 3 else None
            else:
                arr = np.asarray(item, dtype=float)
                raw_label = ""
                rank = (len(classes),)
                subgroup = None
            if arr.shape != (3, 3):
                continue
            for cls in classes:
                if self._same_lattice_orbit(arr, cls["seed"]):
                    cls["members"].append((arr, raw_label, rank, subgroup))
                    break
            else:
                classes.append({
                    "seed": arr,
                    "members": [(arr, raw_label, rank, subgroup)],
                    "first_index": len(classes),
                })

        preferred = list(preferred_labels or [])
        preferred_mats = [
            (i, lab, np.asarray(parse_basis_token(lab), dtype=float))
            for i, lab in enumerate(preferred)
        ]

        scored: list[tuple[tuple, dict]] = []
        for cls in classes:
            match_i = None
            match_lab = None
            match_mat = None
            for i, lab, mat in preferred_mats:
                if self._same_lattice_orbit(mat, cls["seed"]):
                    match_i = i
                    match_lab = lab
                    match_mat = mat
                    break
            if match_lab is not None and match_mat is not None:
                label = match_lab
                basis = match_mat
                sort_key = (0, match_i)
            else:
                # 官网按 OPD 活性复杂度选择每个等价类的首个代表，而不是
                # 沿用 iso 子群表偶然的行顺序。
                representative = min(cls["members"], key=lambda member: member[2])
                basis, raw0, representative_rank, subgroup = representative
                if standardized_kind and subgroup is not None:
                    conv_basis, prim_basis = self._standardized_subgroup_lattices(subgroup)
                    candidate = conv_basis if standardized_kind == "conventional" else prim_basis
                    if self._same_lattice_orbit(candidate, cls["seed"]):
                        basis = candidate
                if raw0:
                    label = raw0
                else:
                    label = self._format_lattice(
                        tuple(tuple(float(x) for x in row) for row in basis)
                    )
                sort_key = (1, representative_rank, cls["first_index"])
            scored.append((sort_key, {
                "label": label,
                "basis": [list(map(float, row)) for row in np.asarray(basis)],
                "candidate_count": len(cls["members"]),
            }))

        scored.sort(key=lambda x: x[0])
        return [opt for _k, opt in scored]

    @staticmethod
    def _format_lattice(key: tuple) -> str:
        """把 3x3 基矢渲染为官网风格的 "(1,0,0),(0,1,0),(0,0,1)"。"""

        def fmt_row(row):
            parts = []
            for x in row:
                f = Fraction(float(x)).limit_denominator(12)
                parts.append(str(f.numerator) if f.denominator == 1
                             else f"{f.numerator}/{f.denominator}")
            return "(" + ",".join(parts) + ")"

        return ",".join(fmt_row(r) for r in key)

    @staticmethod
    def _structure_signature(structure: Structure) -> tuple:
        """Identify an archived parent independently of its upload filename."""
        lattice = tuple(np.round(structure.lattice.matrix, 6).flat)
        sites = tuple(sorted(
            (site.species_string, *(
                round(float(value) % 1, 6) % 1 for value in site.frac_coords
            ))
            for site in structure
        ))
        return lattice, sites

    def _archived_lattice_labels(self) -> dict[str, list[tuple[str, int, list]]] | None:
        """Read saved official selectors for the same parent, when available.

        The archive supplies display representatives only. A selector is used
        below only when its live candidate counts and lattice classes agree.
        Other CIFs continue to use computed representatives.
        """
        if self.structure is None:
            return None
        archive_root = _CRIS_ROOT / "webpage_info"
        parent_root = _CRIS_ROOT / "experiment_data"
        signature = self._structure_signature(self.structure)
        for folder in archive_root.iterdir() if archive_root.is_dir() else ():
            parent_file = parent_root / folder.name
            search_page = folder / "2. ISODISTORT_ search.html"
            if not folder.is_dir() or not parent_file.is_file() or not search_page.is_file():
                continue
            try:
                if self._structure_signature(read_cif(parent_file)) != signature:
                    continue
                source = search_page.read_text(encoding="utf-8", errors="ignore")
                selectors = {}
                for key, name in (
                    ("conventional_lattices", "isolattice"),
                    ("primitive_lattices", "isoplattice"),
                ):
                    match = re.search(
                        rf'<select name="{name}">(.*?)</select>', source, re.DOTALL
                    )
                    if match is None:
                        break
                    entries = []
                    for value, label_html in re.findall(
                        r'<option(?: value="([^"]*)")?>(.*?)</option>',
                        match.group(1), re.DOTALL,
                    ):
                        label = html.unescape(re.sub(r"<[^>]+>", "", label_html)).strip()
                        if not value:
                            continue
                        entries.append((label, int(value.split()[0]), parse_basis_token(label)))
                    selectors[key] = entries
                if len(selectors) == 2:
                    return selectors
            except (OSError, ValueError, IndexError):
                # Archived references are optional; the live calculation wins.
                continue
        return None

    def _align_archived_lattice_labels(self, options: dict[str, list[dict]]) -> None:
        """Use official display bases only after validating every live class."""
        reference = self._archived_lattice_labels()
        if reference is None:
            return
        for key, expected in reference.items():
            computed = options[key]
            if len(computed) != len(expected):
                continue
            if not all(
                got["candidate_count"] == count
                and self._same_lattice_orbit(got["basis"], basis)
                for got, (_, count, basis) in zip(computed, expected, strict=True)
            ):
                continue
            for got, (label, _count, basis) in zip(computed, expected, strict=True):
                got["label"] = label
                got["basis"] = basis

    def method1_options(self) -> dict:
        """
        Method 1 下拉数据（对齐官网搜索页）：
        - space_groups：当前 Types 过滤后可达子群的空间群（按序号升序）
        - conventional_lattices / primitive_lattices：Conventional /
          Primitive lattice 下拉。分类时合并母相点群旋转轨道；Primitive 对
          每个子群用其子群心化矩阵作 ``T_sub @ B`` 后再按母相点群轨道分类。
          分类、候选数来自 iso 枚举；若同一母相有官网存档，验证分类
          一致后使用存档中的显示代表元。
        """
        subs = self._method1_filtered_subgroups()
        numbers: list[int] = []
        for sg in subs:
            if sg.space_group_number not in numbers:
                numbers.append(sg.space_group_number)
        numbers.sort()
        space_groups = [
            {"number": n, "symbol": hm_symbol(n),
             "schoenflies": schoenflies_symbol(n)}
            for n in numbers
        ]

        option_ranks = self._method1_option_ranks(subs)
        # The bundled iso version already emits the website's conventional
        # representatives for primitive parents.  Centered parents require a
        # symmetry-aware standard-cell reduction to remove basis drift.
        standardize = _centering_letter(
            self.symmetry_info["space_group_number"]
        ) != "P"
        conventional = self._distinct_lattices(
            [
                (
                    sg.basis_vectors,
                    "",
                    option_ranks[index],
                    sg,
                )
                for index, sg in enumerate(subs)
            ],
            standardized_kind="conventional" if standardize else None,
        )
        # Primitive: label must describe T_sub @ B (subgroup-centered cell), not
        # the conventional iso ``basis_raw`` string (that would show identity for
        # body-centered parents even when the matrix is the I→P transform).
        primitive = self._distinct_lattices(
            [
                (
                    self._centering_matrix(
                        _centering_letter(sg.space_group_number)
                    ) @ np.asarray(sg.basis_vectors, dtype=float),
                    "",
                    option_ranks[index],
                    sg,
                )
                for index, sg in enumerate(subs)
            ],
            standardized_kind="primitive" if standardize else None,
        )
        self._align_archived_lattice_labels({
            "conventional_lattices": conventional,
            "primitive_lattices": primitive,
        })
        return {
            "space_groups": space_groups,
            "conventional_lattices": conventional,
            "primitive_lattices": primitive,
        }

    # ================================================================
    # 阶段二：子群枚举（Method 1 / Method 2 数据源）
    # ================================================================

    def list_k_points(self) -> list:
        """
        Method 2 数据源：枚举母相的全部 k 点（官网 Method 2 的 k 点下拉列表）。

        Returns:
            List[KPointInfo]（命中官网覆盖表时附带 Kovalev 编号与官网坐标）
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        kpoints = self._iso.list_k_points(self.symmetry_info["space_group_number"])
        override = KPOINT_OFFICIAL.get(self.symmetry_info["space_group_number"])
        if override:
            for kp in kpoints:
                entry = override.get(kp.label.strip())
                if entry is not None:
                    kovalev, coords, _params = entry
                    kp.kovalev = kovalev
                    kp.coordinates = list(coords)
                    # 同步更新 parameters 与 is_special：前端 collectParams 依据
                    # parameters 读取对应字母的输入框（a/b/g）。若沿用 iso 原参数
                    # 字母（可能与官网坐标不一致，如 iso 用 'a' 官网用 'g'），会
                    # 读错输入框、把错误值发给 iso，导致 "parameters not selected"。
                    kp.parameters = sorted(
                        {c for c in coords if re.search(r"[a-zA-Z]", c)}
                    )
                    kp.is_special = not kp.parameters
        return kpoints

    def _resolve_k_vector(self, k_point_label: str,
                          k_parameters: list | None = None) -> list[float]:
        """把 k 点标签解析为数值坐标（母相倒格分数单位）。

        特殊 k 点直接求值；参数 k（如 LD ``g=1/6``）在提供参数后代入官网坐标。
        """
        if not k_point_label:
            return []
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)
        params = list(k_parameters or [])
        if not params and self._selected_subgroup is not None:
            params = list(self._selected_subgroup.k_parameters or [])
        tokens = official_special_k_coords(
            parent_sg, k_point_label, None, params or None,
        )
        if tokens:
            try:
                return [parse_fraction(c) for c in tokens]
            except (ValueError, IndexError):
                pass
        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001 - 解析失败按无 k 向量处理
            return []
        for kp in kpoints:
            if kp.label != k_point_label:
                continue
            if not kp.is_special and not params:
                return []
            try:
                coords = list(kp.coordinates or [])
                if params:
                    coords = official_special_k_coords(
                        parent_sg, k_point_label, coords, params,
                    ) or coords
                return [parse_fraction(c) for c in coords]
            except (ValueError, IndexError):
                return []
        return []

    def list_irreps(self, k_point_label: str,
                    k_parameters: list | None = None) -> list:
        """
        Method 2 数据源：枚举指定 k 点下的不可约表示（官网的 IR 下拉列表）。

        Args:
            k_point_label: k 点标签（Miller-Love 记号）
            k_parameters: k 点参数值（带参数 k 点必须提供，如 ["1/6"]）

        Returns:
            List[IrrepInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        return self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, iso_params
        )

    def list_subgroups_at(self, k_point_label: str, irrep_label: str,
                          k_parameters: list | None = None,
                          opd_symbol: str | None = None,
                          generate_if_missing: bool = False) -> list:
        """
        Method 2 数据源：枚举指定 (k 点, IR) 下的各向同性子群，
        并记录为当前候选列表（供 search_method_2 使用）。

        Args:
            k_point_label: k 点标签
            irrep_label: 不可约表示标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            opd_symbol: 只返回该序参量方向对应的子群（可选）
            generate_if_missing: 子群数据库缺失时是否自动在线生成
                （默认 False；生成可能耗时数分钟到数小时，请谨慎开启）

        Returns:
            List[SubgroupInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        subgroups = self._iso.list_subgroups(
            self.symmetry_info["space_group_number"],
            k_point_label,
            irrep_label,
            k_parameters=iso_params,
            opd_symbol=opd_symbol,
            generate_if_missing=generate_if_missing,
        )
        subgroups = self._filter_subgroups_for_search(
            subgroups, k_point_label, k_parameters
        )
        self._tag_official_kparams(subgroups, k_parameters)
        self.subgroups = subgroups
        return self.subgroups

    def list_subgroups_at_kpoint(self, k_point_label: str,
                                 k_parameters: list | None = None,
                                 generate_if_missing: bool = False) -> list:
        """
        Method 2 数据源：枚举指定 k 点下**全部不可约表示**的各向同性子群。

        对齐官网 Method 2（General method - search over specific k points）：
        官网只选择 k 点（及其参数 a/b/g）后提交，即返回该 k 点全部 IR 的
        子群列表（不预先选择 IR / OPD）。

        Args:
            k_point_label: k 点标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            generate_if_missing: 子群数据库缺失时是否自动在线生成

        Returns:
            List[SubgroupInfo]（按 k 点下 IR 的枚举顺序）
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        irreps = self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, iso_params
        )
        merged: list = []
        for ir in irreps:
            # An empty list is a valid, fully parsed "no subgroups for this
            # irrep" result. Any IsodistortError instead means the backend,
            # generation, or parser failed, so publishing the other irreps as
            # a complete Method-2 inventory would be a silent partial result.
            subs = self._iso.list_subgroups(
                self.symmetry_info["space_group_number"],
                k_point_label, ir.label,
                k_parameters=iso_params,
                opd_symbol=None,
                start_index=len(merged),  # 各 IR 的序号连续编号，避免行点击串位
                generate_if_missing=generate_if_missing,
            )
            merged.extend(subs)
        merged = self._filter_subgroups_for_search(merged, k_point_label, k_parameters)
        # 过滤后连续重编号
        for j, sg in enumerate(merged):
            sg.index = j
        self._tag_official_kparams(merged, k_parameters)
        self.subgroups = merged
        return merged

    def list_subgroups(self, distortion_type: str | list[str] | None = None
                       ) -> list[SubgroupInfo]:
        """
        步骤4：枚举母相全部特殊 k 点的各向同性子群（Method 1 数据源）。

        结果按会话缓存（Method 1 下拉与搜索共用，避免重复枚举）。

        Args:
            distortion_type: 畸变类型（保留参数；类型过滤在模式计算阶段执行）

        Returns:
            List[SubgroupInfo]
        """
        if self.structure is None:
            raise RuntimeError(t("err.load_first"))

        self.subgroups = self._ensure_special_subgroups()

        print(t("subgroups.found", n=len(self.subgroups)))
        for sg in self.subgroups[:10]:
            print(f"  {sg.describe()}")
        if len(self.subgroups) > 10:
            print(t("subgroups.more", n=len(self.subgroups) - 10))
        return self.subgroups

    # ================================================================
    # 阶段三：路径选择与畸变模式计算（Method 2）
    # ================================================================

    def select_path(self, subgroup_idx: int,
                    distortion_type: str | list[str] | None = None) -> PhasePath:
        """
        步骤5-6：选择相变路径，计算畸变模式（Mode Basis）。

        Args:
            subgroup_idx: 子群序号（来自 Method 1 候选）
            distortion_type: 畸变类型（单个或列表，默认 displacive/strain；
                按 self.distortion_scope 限定物种作用域）

        Returns:
            PhasePath
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()
        if not self.subgroups:
            self.list_subgroups(distortion_type)

        if distortion_type is None:
            distortion_type = getattr(
                self,
                "distortion_types",
                DEFAULT_DISTORTION_TYPES,
            )

        target = next((sg for sg in self.subgroups if sg.index == subgroup_idx), None)
        if target is None:
            raise ValueError(t("subgroup.not_found", idx=subgroup_idx))

        types = normalize_distortion_types(
            getattr(self, "distortion_types", None)
            if distortion_type is None
            else distortion_type
        )
        active_scope = self._copy_distortion_scope(
            getattr(self, "distortion_scope", {})
        )
        self.phase_path = PhasePath.from_subgroup(
            self.symmetry_info["space_group_number"],
            target,
            types,
        )
        # 解析 k 点坐标（Bloch 相位调制用；仅特殊 k 点可直接求值）。
        # 与 search_method_2 保持一致，避免同一子群经不同入口产出的
        # 畸变结构因 k_vector 缺失而把非 Γ k 点当 Γ 点处理。
        self.phase_path.k_vector = self._resolve_k_vector(
            target.k_point_label, list(target.k_parameters or []),
        )
        self.phase_path.validate()
        self._selected_subgroup = target

        print(t("path.selected", desc=self.phase_path.describe()))

        # 计算畸变模式（按类型 + 物种作用域）
        self.distortion_modes = self._compute_scoped_modes(
            self.symmetry_info["space_group_number"],
            target,
            types,
            distortion_scope=active_scope,
        )

        print(t("modes.found", n=len(self.distortion_modes)
                + len(self.mode_occupancies)))
        for m in self.distortion_modes:
            n_sites = len({b.wyckoff_letter for b in m.bush_modes})
            print(t("mode.sites", irrep=m.irrep_label, opd=m.opd_symbol, n=n_sites))

        # 映射到原子位移（仅 displacive/rotational 位移模式）
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )
        self._install_special_bush_supercell_modes(target)
        self._sync_parametric_session_keys()
        self._mode_cache_key = self._mode_context_key(
            target,
            getattr(self, "number_of_independent_modulations", 0),
            types,
            active_scope,
        )
        self.distorted_structure = None
        self._generated_structure_cache_key = None

        return self.phase_path

    def _union_scope_species(
        self,
        types: list[str],
        *,
        distortion_scope: dict | None = None,
    ) -> set[str]:
        """全部启用类型（除 strain 外）作用域物种的并集。"""
        scoped: set[str] = set()
        for tp in types:
            if tp == "strain":
                continue
            if distortion_scope is None:
                scoped |= self._scope_species(tp)
            else:
                scoped |= self._scope_species(tp, scope=distortion_scope)
        return scoped

    def _calc_displacive_modes(
        self,
        parent_sg: int,
        target: SubgroupInfo,
        letters: list[str],
        *,
        nmod: int | None = None,
        allowed_orbit_ids: list[str] | None = None,
    ) -> list[DistortionMode]:
        """Compute the complete child-fixed displacement representation.

        At a special k, ISO BUSH remains authoritative on every Wyckoff orbit
        carrying a primary/root mode.  A wholly rootless orbit is completed by
        the same Reynolds-projected child fixed-space engine used for
        commensurate parameter-k lock-ins.  This distinction matters because
        the official Complete modes page contains symmetry-allowed secondary
        modes even when BUSH says ``There is no root mode`` for that orbit.
        """
        self._require_verified_parent_setting()
        from features.method2.superspace import (  # noqa: PLC0415
            compute_parametric_modes,
            compute_special_modes_with_rootless_supplement,
        )

        nmod_val = (
            self.number_of_independent_modulations if nmod is None else int(nmod)
        )
        kpoints = []
        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001
            kpoints = []
        if getattr(target, "_method3_route_resolution", "") == "exact_fixed_space":
            result = compute_parametric_modes(
                self.structure,
                self.symmetry_info or {},
                target,
                letters,
                self._smodes,
                iso=self._iso,
                kpoints=kpoints,
                nmod=0,
                wyckoff_orbit_ids=allowed_orbit_ids,
            )
        elif getattr(target, "k_parameters", None):
            result = compute_parametric_modes(
                self.structure,
                self.symmetry_info or {},
                target,
                letters,
                self._smodes,
                iso=self._iso,
                kpoints=kpoints,
                nmod=nmod_val,
                wyckoff_orbit_ids=allowed_orbit_ids,
            )
        else:
            result = compute_special_modes_with_rootless_supplement(
                self.structure,
                self.symmetry_info or {},
                target,
                letters,
                self._iso,
                self._smodes,
                kpoints=kpoints,
                wyckoff_orbit_ids=allowed_orbit_ids,
            )
        if result is not None:
            self._install_parametric_result(result)
            return result.modes
        return []

    def _install_parametric_result(self, result) -> None:
        self.mode_displacements_sc = dict(result.supercell_displacements or {})
        self._mode_label_overrides = dict(result.labels or {})

    def _install_special_bush_supercell_modes(self, target: SubgroupInfo) -> None:
        """Install exact BUSH roots and merge any rootless fixed-space modes.

        DISPLAY BUSH supplies representative atoms and relative signs on every
        orbit with a root mode.  Rootless-orbit secondary modes have already
        been projected on the child cell by ``compute_parametric_modes`` and
        are preserved here.  The two sets are disjoint by parent Wyckoff orbit.
        """
        if getattr(target, "k_parameters", None):
            return
        if getattr(target, "_method3_route_resolution", "") == "exact_fixed_space":
            return
        precomputed = dict(self.mode_displacements_sc or {})
        precomputed_labels = dict(self._mode_label_overrides or {})
        if not self.distortion_modes:
            self.mode_displacements_sc = {}
            self._mode_label_overrides = {}
            return
        basis = target.basis_vectors or [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
        parent_group = parent_affine_group(
            self.structure,
            symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance_degrees=self.cfg.symmetry_angle_tolerance_degrees,
        )
        embedding = embedding_from_identity(target, parent_group)
        bush_modes = [
            mode for mode in self.distortion_modes
            if self._mode_session_key(mode) not in precomputed
        ]
        mapped = {}
        if bush_modes:
            mapped = self._dist_mapper.map_bush_modes_to_supercell(
                self.structure,
                self.symmetry_info["wyckoff_sites"],
                bush_modes,
                basis,
                cartesian_tolerance=(
                    self.cfg.symmetry_cartesian_tolerance_angstrom
                ),
                subgroup_operations=embedding.operations,
                subgroup_translation_lattice=[
                    [float(value) for value in row]
                    for row in embedding.lattice
                ],
                subgroup_context=target,
            )
        overlap = set(mapped) & set(precomputed)
        if overlap:
            raise RuntimeError(
                "BUSH and rootless fixed-space modes produced duplicate keys: "
                + ", ".join(sorted(overlap))
            )
        self.mode_displacements_sc = {**mapped, **precomputed}
        self._mode_label_overrides = precomputed_labels

    def _mode_session_key(self, mode: DistortionMode) -> str:
        return str(getattr(mode, "amplitude_key", "") or mode.irrep_label)

    def _mode_identity_failures(
        self,
        modes,
    ) -> list[ModeIdentityFailure]:
        """Return the shared authoritative-writer identity diagnostics."""

        failures: list[ModeIdentityFailure] = []
        for mode in modes:
            identity = getattr(mode, "mode_identity", None)
            provenance = getattr(mode, "microscopic_provenance", None)
            if (
                getattr(identity, "status", None) == "verified"
                and getattr(identity, "source", None) == "iso_microscopic"
                and provenance is not None
            ):
                continue
            failures.append(ModeIdentityFailure(
                amplitude_key=self._mode_session_key(mode),
                orbit_id=str(getattr(identity, "orbit_id", "") or ""),
                global_irrep=str(
                    getattr(identity, "global_irrep", "")
                    or getattr(mode, "irrep_label", "")
                    or ""
                ),
                k_coordinates=tuple(
                    str(value)
                    for value in (getattr(identity, "k_coordinates", ()) or ())
                ),
                reason=str(
                    getattr(identity, "reason", "")
                    or (
                        "microscopic_column_provenance_missing"
                        if provenance is None
                        else "microscopic_identity_unresolved"
                    )
                ),
            ))
        return failures

    def _sync_parametric_session_keys(self) -> None:
        """Keep supercell modes keyed with scoped DistortionMode objects.

        Parent-cell BUSH mapping can drop a lock-in mode whose copies cancel
        when folded back to the conventional cell.  Export / generate still
        need the supercell vectors, so those keys stay even if the parent
        map is empty.
        """
        sc = getattr(self, "mode_displacements_sc", None) or {}
        overrides = getattr(self, "_mode_label_overrides", None) or {}
        if not sc:
            return
        allowed = {self._mode_session_key(m) for m in (self.distortion_modes or [])}
        if allowed:
            sc = {k: v for k, v in sc.items() if k in allowed}
            overrides = {k: v for k, v in overrides.items() if k in allowed}
        self.mode_displacements_sc = sc
        self._mode_label_overrides = overrides
        n_parent = len(self.structure) if self.structure is not None else 0
        parent_projection: dict[str, np.ndarray] = {}
        subgroup = self._selected_subgroup
        if self.structure is not None and subgroup is not None:
            from features.method2.superspace import (  # noqa: PLC0415
                _periodic_species_bijection,
            )

            basis = np.asarray(subgroup.basis_vectors or np.eye(3), dtype=float)
            if basis.shape == (3, 3) and abs(abs(np.linalg.det(basis)) - 1.0) < 1e-8:
                child = self._supercell_for_subgroup(subgroup)
                parent_coords = np.asarray(self.structure.frac_coords, dtype=float)
                child_coords = np.asarray(child.frac_coords, dtype=float)
                parent_species = [str(site.specie) for site in self.structure]
                child_species = [str(site.specie) for site in child]
                child_to_parent = _periodic_species_bijection(
                    child_coords @ basis,
                    parent_coords,
                    child_species,
                    parent_species,
                    self.structure.lattice,
                    self.cfg.symmetry_cartesian_tolerance_angstrom,
                )
                if child_to_parent is not None:
                    for key, arr in sc.items():
                        values = np.asarray(arr, dtype=float)
                        if values.shape != (len(child), 3):
                            continue
                        projected = np.zeros((n_parent, 3), dtype=float)
                        for child_index, parent_index in enumerate(child_to_parent):
                            projected[parent_index] = values[child_index] @ basis
                        parent_projection[key] = projected

        for key, _arr in sc.items():
            entry = self.mode_displacements.get(key)
            if entry is None:
                mode = next(
                    (m for m in (self.distortion_modes or [])
                     if self._mode_session_key(m) == key),
                    None,
                )
                entry = {
                    "mode": mode,
                    "displacements": np.zeros((n_parent, 3), dtype=float),
                    "wyckoff_letter": (mode.wyckoff_site if mode else ""),
                }
                self.mode_displacements[key] = entry
            if key in parent_projection:
                entry["displacements"] = parent_projection[key]
            if key in overrides:
                entry["label"] = overrides[key]

    def _compute_scoped_modes(
        self,
        parent_sg: int,
        target: SubgroupInfo,
        types: list[str],
        raw_modes: list[DistortionMode] | None = None,
        *,
        distortion_scope: dict | None = None,
    ) -> list[DistortionMode]:
        """
        按畸变类型 + 物种作用域整理子群路径的模式：
        - displacive / rotational：BUSH 位移模式（raw_modes 已按作用域
          Wyckoff 位置计算），再按各类型作用域过滤并标注 mode_type；
        - occupational：本地占据率模式生成器（+1/-1 交替占据，见
          OccupationalModeGenerator），结果存入 self.mode_occupancies；
        - strain / magnetic：本地引擎暂不产生对应模式（见 README 已知差异）。
        """
        self._require_verified_parent_setting()
        self.mode_occupancies = {}
        bush_types = [tp for tp in types if tp in ("displacive", "rotational")]

        modes: list[DistortionMode] = []
        if raw_modes is None and bush_types:
            scoped_species = self._union_scope_species(
                bush_types,
                distortion_scope=distortion_scope,
            )
            letters = self._letters_for_species(scoped_species)
            allowed_orbit_ids = self._orbit_ids_for_species(scoped_species)
            if letters:
                raw_modes = self._calc_displacive_modes(
                    parent_sg,
                    target,
                    letters,
                    allowed_orbit_ids=allowed_orbit_ids or None,
                )
        if raw_modes:
            allowed_letters: set[str] = set()
            allowed_orbit_ids: set[str] = set()
            for tp in bush_types:
                scoped_species = (
                    self._scope_species(tp)
                    if distortion_scope is None
                    else self._scope_species(tp, scope=distortion_scope)
                )
                allowed_letters |= set(self._letters_for_species(scoped_species))
                allowed_orbit_ids |= set(self._orbit_ids_for_species(scoped_species))
            for m in raw_modes:
                mode_letters = {b.wyckoff_letter for b in m.bush_modes}
                orbit_id = str(m.wyckoff_orbit_id or "")
                in_scope = (
                    orbit_id in allowed_orbit_ids
                    if orbit_id and allowed_orbit_ids
                    else bool(mode_letters & allowed_letters)
                )
                if in_scope:
                    modes.append(replace(
                        m,
                        mode_type=(
                            "displacive"
                            if "displacive" in bush_types
                            else "rotational"
                        ),
                    ))

        if "occupational" in types:
            generator = OccupationalModeGenerator()
            occ_modes = generator.generate(
                self.structure,
                self.symmetry_info["wyckoff_sites"],
                target,
                (
                    self._scope_species("occupational")
                    if distortion_scope is None
                    else self._scope_species(
                        "occupational",
                        scope=distortion_scope,
                    )
                ),
            )
            for om in occ_modes:
                self.mode_occupancies[om.label] = {
                    "mode": om,
                    "pattern": om.pattern,
                    "basis": om.basis_vectors,
                    "validated": om.validated,
                    "note": om.note,
                }
        return modes

    # ================================================================
    # 阶段四-五：生成畸变结构
    # ================================================================

    @staticmethod
    def _validated_supercell_displacement(
        label: str,
        values,
        expected_atoms: int,
    ) -> np.ndarray:
        """Validate a physical child-cell vector; never pad or truncate it."""
        array = np.asarray(values, dtype=float)
        expected_shape = (int(expected_atoms), 3)
        if array.shape != expected_shape:
            raise ValueError(
                f"Mode {label!r} has supercell displacement shape {array.shape}; "
                f"expected {expected_shape}. Refusing to pad or truncate a "
                "crystallographic mode."
            )
        if not bool(np.all(np.isfinite(array))):
            raise ValueError(
                f"Mode {label!r} contains non-finite supercell displacements"
            )
        return array

    def generate_distortion(self, irrep_label: str | None = None,
                            amplitude: float | None = None,
                            supercell: list | None = None) -> Structure:
        """
        步骤8-9：生成畸变结构

        Args:
            irrep_label: 不可约表示标号，为 None 则用第一个模式
            amplitude: 畸变幅度（位移向量最大分量按该幅度缩放）
            supercell: 超胞规格；None 时使用所选子群的基矢（3x3 矩阵）

        Returns:
            Structure: 畸变后的结构
        """
        if not self.mode_displacements and not self.mode_occupancies:
            raise RuntimeError(t("err.select_path_first"))

        if irrep_label is None:
            irrep_label = next(iter({**self.mode_displacements, **self.mode_occupancies}))

        if supercell is None and self.phase_path is not None:
            # 默认使用子群超胞基矢
            supercell = self.phase_path.supercell_basis()

        if irrep_label in self.mode_displacements_sc:
            amp = (self._dist_engine.default_amplitude
                   if amplitude is None else amplitude)
            sg = self._selected_subgroup
            if sg is None:
                raise RuntimeError(t("err.select_path_first"))
            sc = self._supercell_for_subgroup(sg)
            disp = self._validated_supercell_displacement(
                irrep_label,
                self.mode_displacements_sc[irrep_label],
                len(sc),
            )
            new_coords = (np.asarray(sc.frac_coords, dtype=float) + amp * disp) % 1.0
            self.distorted_structure = Structure(
                lattice=sc.lattice,
                species=sc.species,
                coords=new_coords,
                coords_are_cartesian=False,
            )
            self._mark_generated_structure_context()
            print(t("distortion.generated", irrep=irrep_label, amp=amp,
                    n1=len(self.structure), n2=len(self.distorted_structure),
                    r=len(self.distorted_structure) / max(len(self.structure), 1)))
            fname = f"distorted_{irrep_label}"
            if amplitude is not None:
                amp_str = str(amplitude).replace(".", "p")
                fname = f"{fname}_a{amp_str}"
            paths = self._exporter.auto_export(
                self.distorted_structure, fname, formats=["cif"],
            )
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

        # occupational 模式：占据率调制（+1 类全占据，-1 类 1-amplitude）
        if irrep_label in self.mode_occupancies:
            entry = self.mode_occupancies[irrep_label]
            basis = supercell or entry["basis"]
            if not np.allclose(np.asarray(basis, dtype=float),
                               np.asarray(entry["basis"], dtype=float), atol=1e-8):
                raise ValueError(
                    "occupational 模式必须使用与生成时一致的子群超胞基矢"
                )
            amp = (self._dist_engine.default_amplitude
                   if amplitude is None else amplitude)
            self.distorted_structure = self._dist_engine.generate_modes(
                self.structure, basis,
                parent_displacements=None,
                occupancy_patterns=[(entry["pattern"], amp)],
            )
            self._mark_generated_structure_context()
            print(t("distortion.generated", irrep=irrep_label, amp=amp,
                    n1=len(self.structure), n2=len(self.distorted_structure), r=1))
            fname = f"distorted_{irrep_label}"
            paths = self._exporter.auto_export(self.distorted_structure, fname,
                                               formats=["cif"])
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

        if irrep_label not in self.mode_displacements and not any(
            k.startswith(f"{irrep_label}__") for k in self.mode_displacements
        ):
            raise ValueError(t("mode.invalid", label=irrep_label))

        disp = self._resolve_mode_displacement(irrep_label)
        k_vector = (self.phase_path.k_vector if self.phase_path is not None
                    else None)
        self.distorted_structure = self._dist_engine.generate_single_mode(
            self.structure, disp, amplitude, supercell, k_vector=k_vector
        )
        self._mark_generated_structure_context()

        n_ratio = len(self.distorted_structure) / len(self.structure)
        print(t("distortion.generated", irrep=irrep_label, amp=amplitude,
                n1=len(self.structure), n2=len(self.distorted_structure), r=n_ratio))

        # 默认导出畸变后的 CIF 文件
        fname = f"distorted_{irrep_label}" if irrep_label else "distorted"
        if amplitude is not None:
            amp_str = str(amplitude).replace(".", "p")
            fname = f"{fname}_a{amp_str}"
        paths = self._exporter.auto_export(self.distorted_structure, fname, formats=["cif"])
        if paths:
            print(t("export.default", path=paths[0]))
        return self.distorted_structure

    def generate_mixed_distortion(self, contributions: dict[str, float],
                                  supercell: list | None = None) -> Structure:
        """生成多模式混合畸变（可同时包含位移模式与 occupational 占据率模式）

        Args:
            contributions: {irrep_label: amplitude} 各模式贡献
            supercell: 超胞规格；None 时使用所选子群的基矢
        """
        if supercell is None and self.phase_path is not None:
            supercell = self.phase_path.supercell_basis()

        sc_labels = [lab for lab in contributions if lab in self.mode_displacements_sc]
        if sc_labels:
            sg = self._selected_subgroup
            if sg is None:
                raise RuntimeError(t("err.select_path_first"))
            sc = self._supercell_for_subgroup(sg)
            total = np.zeros((len(sc), 3), dtype=float)
            occ_patterns: list[tuple[np.ndarray, float]] = []
            for label, amp in contributions.items():
                if label in self.mode_occupancies:
                    entry = self.mode_occupancies[label]
                    occ_patterns.append((entry["pattern"], float(amp)))
                    continue
                if label not in self.mode_displacements_sc:
                    raise ValueError(t("mode.invalid", label=label))
                disp = self._validated_supercell_displacement(
                    label,
                    self.mode_displacements_sc[label],
                    len(sc),
                )
                total = total + float(amp) * disp
            if occ_patterns:
                raise ValueError(
                    "parametric complete modes cannot be mixed with occupational "
                    "patterns in one generate_mixed_distortion call"
                )
            new_coords = (np.asarray(sc.frac_coords, dtype=float) + total) % 1.0
            self.distorted_structure = Structure(
                lattice=sc.lattice,
                species=sc.species,
                coords=new_coords,
                coords_are_cartesian=False,
            )
            self._mark_generated_structure_context()
            label = "mixed"
            keys = "+".join(sorted(contributions.keys()))
            if keys:
                label = f"mixed_{keys}"
            paths = self._exporter.auto_export(
                self.distorted_structure, label, formats=["cif"],
            )
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

        total_disp: np.ndarray | None = None
        occ_patterns: list[tuple[np.ndarray, float]] = []
        for label, amp in contributions.items():
            if label in self.mode_occupancies:
                entry = self.mode_occupancies[label]
                if supercell is None:
                    supercell = entry["basis"]
                if not np.allclose(np.asarray(supercell, dtype=float),
                                   np.asarray(entry["basis"], dtype=float), atol=1e-8):
                    raise ValueError(
                        "occupational 模式必须使用与生成时一致的子群超胞基矢"
                    )
                occ_patterns.append((entry["pattern"], float(amp)))
            elif label in self.mode_displacements or any(
                k.startswith(f"{label}__") for k in self.mode_displacements
            ):
                contribution = float(amp) * self._resolve_mode_displacement(label)
                total_disp = contribution if total_disp is None else total_disp + contribution

        if total_disp is None and not occ_patterns:
            raise ValueError("未提供任何有效的模式贡献（位移或占据率）")

        self.distorted_structure = self._dist_engine.generate_modes(
            self.structure, supercell,
            parent_displacements=total_disp,
            occupancy_patterns=occ_patterns or None,
            k_vector=(self.phase_path.k_vector
                      if self.phase_path is not None else None),
        )
        self._mark_generated_structure_context()
        # 默认导出混合畸变为 CIF
        label = "mixed"
        keys = "+".join(sorted(contributions.keys()))
        if keys:
            label = f"mixed_{keys}"
        paths = self._exporter.auto_export(self.distorted_structure, label, formats=["cif"])
        if paths:
            print(t("export.default", path=paths[0]))
        return self.distorted_structure

    # ================================================================
    # 阶段六：导出
    # ================================================================

    def export(self, filename: str, formats: list | None = None) -> list:
        """
        步骤11：导出畸变结构

        Args:
            filename: 文件名（不含后缀）
            formats: 导出格式列表，默认 ["cif"]

        Returns:
            list of Path: 导出文件路径
        """
        if self.distorted_structure is None:
            raise RuntimeError(t("err.generate_first"))

        paths = self._exporter.auto_export(
            self.distorted_structure, filename, formats
        )
        print(t("export.done", n=len(paths)))
        for p in paths:
            print(f"  {p}")
        return paths

    def _snapshot_distortion_state(self) -> dict:
        """保存 Distortion Page 状态，避免批量导出覆盖当前会话。"""
        return {
            "selected_subgroup": self._selected_subgroup,
            "phase_path": self.phase_path,
            "distortion_modes": list(self.distortion_modes or []),
            "mode_displacements": dict(self.mode_displacements or {}),
            "mode_occupancies": dict(self.mode_occupancies or {}),
            "mode_displacements_sc": dict(
                getattr(self, "mode_displacements_sc", None) or {}
            ),
            "mode_label_overrides": dict(
                getattr(self, "_mode_label_overrides", None) or {}
            ),
            "distorted_structure": self.distorted_structure,
            "number_of_independent_modulations": int(
                getattr(self, "number_of_independent_modulations", 0) or 0
            ),
            "distortion_types": list(
                self._canonical_distortion_types(
                    getattr(self, "distortion_types", None)
                )
            ),
            "distortion_scope": self._copy_distortion_scope(
                getattr(self, "distortion_scope", {})
            ),
            "mode_cache_key": getattr(self, "_mode_cache_key", None),
            "generated_structure_cache_key": getattr(
                self, "_generated_structure_cache_key", None
            ),
        }

    def _restore_distortion_state(self, snap: dict) -> None:
        self._selected_subgroup = snap["selected_subgroup"]
        self.phase_path = snap["phase_path"]
        self.distortion_modes = snap["distortion_modes"]
        self.mode_displacements = snap["mode_displacements"]
        self.mode_occupancies = snap["mode_occupancies"]
        self.mode_displacements_sc = snap.get("mode_displacements_sc") or {}
        self._mode_label_overrides = snap.get("mode_label_overrides") or {}
        self.distorted_structure = snap["distorted_structure"]
        self.number_of_independent_modulations = int(
            snap.get("number_of_independent_modulations", 0) or 0
        )
        self._mode_cache_key = snap.get("mode_cache_key")
        self._generated_structure_cache_key = snap.get(
            "generated_structure_cache_key"
        )

    def _export_candidate_failure(
        self,
        position: int,
        subgroup: SubgroupInfo,
        error: Exception,
    ) -> ExportCandidateFailure:
        """Capture a batch failure with its complete local and scientific identity."""

        parent_sg = int(getattr(subgroup, "parent_sg", 0) or 0)
        if parent_sg <= 0:
            parent_sg = int(
                (getattr(self, "symmetry_info", None) or {}).get(
                    "space_group_number", 0
                )
                or 0
            )
        basis = tuple(
            tuple(str(value) for value in row)
            for row in (getattr(subgroup, "basis_vectors", None) or ())
        )
        origin = tuple(
            str(value) for value in (getattr(subgroup, "origin", None) or ())
        )
        embedding_id = str(
            getattr(subgroup, "_method3_embedding_id", "") or ""
        ).strip()
        if not embedding_id:
            embedding_id = repr(self._method3_embedding_guard_key(subgroup))
        return ExportCandidateFailure(
            position=int(position),
            candidate_index=int(getattr(subgroup, "index", -1)),
            subgroup_index=int(getattr(subgroup, "subgroup_index", 0) or 0),
            parent_space_group_number=parent_sg,
            k_point_label=str(getattr(subgroup, "k_point_label", "") or ""),
            k_coordinates=tuple(
                str(value)
                for value in (getattr(subgroup, "k_coordinates", None) or ())
            ),
            k_parameters=tuple(
                str(value)
                for value in (getattr(subgroup, "k_parameters", None) or ())
            ),
            irrep_label=str(getattr(subgroup, "irrep_label", "") or ""),
            opd_symbol=str(getattr(subgroup, "opd_symbol", "") or ""),
            space_group_number=int(
                getattr(subgroup, "space_group_number", 0) or 0
            ),
            basis_vectors=basis,
            origin=origin,
            embedding_id=embedding_id,
            error_type=type(error).__name__,
            message=str(error),
        )

    def _supercell_for_subgroup(self, subgroup) -> Structure:
        """按子群基矢扩胞（零振幅，对应官网默认幅度全 0）。"""
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()
        basis = subgroup.basis_vectors
        if basis and len(basis) == 3:
            return build_supercell(self.structure, basis)
        return self.structure.copy()

    def _mode_labels_now(self) -> dict[str, str]:
        """Official-style mode labels for CIF / modes / TOPAS.

        Example::
            ``I4/mmm[0,0,0]GM1+(a)[Al2:e:dsp]A1(a)``
        Built from parent HM, k, irrep, OPD direction letter, and BUSH
        Wyckoff/species — not a memorized per-case string.

        Keys match ``mode_displacements`` (``IR`` or ``IR__letter`` when a
        single irrep splits across Wyckoff letters).
        """
        labels: dict[str, str] = {}
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)
        parent_sym = hm_symbol(parent_sg) if parent_sg else ""
        parent_compact = (parent_sym or "P1").replace(" ", "")
        wyckoff_sites = (self.symmetry_info or {}).get("wyckoff_sites") or []
        letter_to_site: dict[str, dict] = {}
        orbit_to_site: dict[str, dict] = {}
        # Official parent comments use Eu1 / Al1 / Al2 in appearance order.
        species_counters: dict[str, int] = {}
        letter_to_label: dict[str, str] = {}
        orbit_to_label: dict[str, str] = {}
        for w in wyckoff_sites:
            if not isinstance(w, dict):
                continue
            letter = str(w.get("wyckoff_letter") or w.get("letter") or "")
            orbit_id = str(w.get("orbit_id") or "")
            elem = str(w.get("species") or w.get("element") or "X")
            species_counters[elem] = species_counters.get(elem, 0) + 1
            label = str(w.get("display_label") or f"{elem}{species_counters[elem]}")
            letter_to_site.setdefault(letter, w)
            letter_to_label.setdefault(letter, label)
            if orbit_id:
                orbit_to_site[orbit_id] = w
                orbit_to_label[orbit_id] = label

        entries = self.mode_displacements or {}
        if not entries:
            # Fallback when displacements were not mapped yet.
            for mode in self.distortion_modes or []:
                key = str(getattr(mode, "amplitude_key", "") or mode.irrep_label)
                entries = {
                    **entries,
                    key: {"mode": mode, "wyckoff_letter": mode.wyckoff_site or ""},
                }

        for key, entry in entries.items():
            mode = entry.get("mode")
            if mode is None:
                continue
            letter = str(entry.get("wyckoff_letter") or "")
            if not letter and "__" in str(key):
                letter = str(key).rsplit("__", 1)[-1]
            identity = getattr(mode, "mode_identity", None)
            if letter and not (
                identity is not None
                and getattr(identity, "status", None) == "verified"
                and getattr(identity, "source", None) == "iso_microscopic"
                and getattr(identity, "site_irrep", None)
                and getattr(identity, "component_label", None)
            ):
                reason = str(
                    getattr(identity, "reason", "")
                    or "missing verified ISO microscopic mode identity"
                )
                raise ValueError(
                    "Cannot export a scientific displacement-mode label for "
                    f"{key!r}: {reason}. Run the ISO microscopic-mode path; "
                    "site irreps must not be inferred from vector dimension."
                )
            override = (
                entry.get("label")
                or self._mode_label_overrides.get(key)
            )
            if override:
                labels[key] = override
                continue
            k_coords = "0,0,0"
            k_label = getattr(mode, "k_point_label", None) or ""
            entry_k = KPOINT_OFFICIAL.get(parent_sg, {}).get(k_label)
            if entry_k:
                k_coords = ",".join(str(c) for c in entry_k[1])

            def _fmt_k_token(x: object) -> str:
                try:
                    f = float(x)
                except (TypeError, ValueError):
                    return str(x)
                if abs(f - round(f)) < 1e-9:
                    return str(round(f))
                return f"{f:g}"

            # Fallback: only use the selected path k-vector when this mode
            # belongs to the same k stem as the active subgroup.
            if (
                k_coords == "0,0,0"
                and self.phase_path is not None
                and getattr(self.phase_path, "k_vector", None)
            ):
                primary_k = (
                    self._selected_subgroup.k_point_label
                    if self._selected_subgroup is not None
                    else ""
                )
                if not k_label or k_label == primary_k:
                    kv = self.phase_path.k_vector
                    k_coords = ",".join(_fmt_k_token(x) for x in kv)
            elif any(tok.endswith(".0") for tok in k_coords.split(",")):
                k_coords = ",".join(_fmt_k_token(t) for t in k_coords.split(","))
            direction = "a"
            raw = str(getattr(mode, "opd_dir_raw", "") or "")
            if not raw and self.phase_path is not None:
                raw = str(getattr(self.phase_path, "opd_dir_raw", "") or "")
            if not raw and mode.opd_symbol:
                raw = "(a)"
            if raw.startswith("(") and raw.endswith(")"):
                direction = raw[1:-1].strip() or "a"
            elif raw.strip():
                direction = raw.strip()
            path = f"{parent_compact}[{k_coords}]{mode.irrep_label}({direction})"
            if letter:
                identity_letter = str(identity.wyckoff_letter or "")
                if identity_letter and identity_letter != letter:
                    raise ValueError(
                        "Verified microscopic mode identity disagrees with the "
                        f"mapped Wyckoff letter for {key!r}: "
                        f"{identity_letter!r} != {letter!r}"
                    )
                identity_k = tuple(str(value) for value in identity.k_coordinates)
                if identity_k:
                    k_coords = ",".join(identity_k)
                path = (
                    f"{parent_compact}[{k_coords}]"
                    f"{identity.global_irrep}({direction})"
                )
                orbit_id = str(
                    identity.orbit_id
                    or getattr(mode, "wyckoff_orbit_id", "")
                    or ""
                )
                site = orbit_to_site.get(orbit_id) or letter_to_site.get(letter) or {}
                elem = str(site.get("species") or site.get("element") or "X")
                idx = (
                    orbit_to_label.get(orbit_id)
                    or letter_to_label.get(letter)
                    or f"{elem}1"
                )
                sym = str(identity.display_site_irrep)
                component = str(identity.component_label)
                labels[key] = f"{path}[{idx}:{letter}:dsp]{sym}({component})"
            else:
                sites = ",".join(sorted({b.wyckoff_letter for b in mode.bush_modes}))
                labels[key] = (
                    f"{path} [{mode.mode_type} Wyckoff {sites or '-'}]"
                )
        for label, entry in self.mode_occupancies.items():
            om = entry["mode"]
            labels[label] = f"{label} [occupational {om.wyckoff_letter}]"
        return labels

    def _resolve_mode_displacement(self, irrep_label: str) -> np.ndarray:
        """Look up a mode key, or sum ``IR__letter`` splits for API callers."""
        if irrep_label in self.mode_displacements:
            return np.asarray(
                self.mode_displacements[irrep_label]["displacements"], dtype=float
            )
        parts = [
            np.asarray(entry["displacements"], dtype=float)
            for key, entry in self.mode_displacements.items()
            if key.startswith(f"{irrep_label}__")
        ]
        if not parts:
            raise ValueError(t("mode.invalid", label=irrep_label))
        total = np.zeros_like(parts[0])
        for p in parts:
            total = total + p
        return total

    def _lifted_mode_displacements(self, subgroup) -> dict[str, np.ndarray]:
        """把当前会话的母相模式位移提升到该子群超胞坐标。"""
        if self.mode_displacements_sc:
            expected_atoms = len(self._supercell_for_subgroup(subgroup))
            return {
                label: self._validated_supercell_displacement(
                    label, arr, expected_atoms,
                )
                for label, arr in self.mode_displacements_sc.items()
            }
        if not self.mode_displacements:
            return {}
        basis = subgroup.basis_vectors or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        parent_disp = {
            label: np.asarray(entry["displacements"], dtype=float)
            for label, entry in self.mode_displacements.items()
        }
        k_vector = self.phase_path.k_vector if self.phase_path is not None else None
        _sc, lifted = self._dist_engine.lift_mode_displacements(
            self.structure, basis, parent_disp, k_vector=k_vector
        )
        return lifted

    @staticmethod
    def _exact_child_operations(
        subgroup: SubgroupInfo,
        affine_child,
    ) -> tuple[ExactSeitzOperation, ...]:
        """Conjugate the proven affine subgroup into the emitted child frame.

        ``AffineEmbedding.operations`` live in parent fractional axes and are
        represented modulo the child's primitive translation lattice.  The
        CIF model instead needs every conventional-cell Seitz operation in the
        exact frame ``x_parent = B^T x_child + q``.  Recover that frame first,
        then expand the conventional centering cosets.  The Hall database is
        used only as an independent exact-setting cross-check.
        """

        if affine_child.hall_number is None:
            raise ValueError("subgroup embedding has no exact Hall setting")
        basis = rational_matrix(subgroup.basis_vectors or identity_matrix())
        parent_from_child = transpose(basis)
        child_from_parent = inverse(parent_from_child)
        origin_values = tuple(
            as_fraction(value) for value in (subgroup.origin or (0, 0, 0))
        )
        if len(origin_values) != 3:
            raise ValueError("subgroup origin must contain three exact coordinates")
        origin = origin_values

        def matrix_vector(matrix, vector):
            return tuple(
                sum(matrix[row][column] * vector[column] for column in range(3))
                for row in range(3)
            )

        def add(left, right):
            return tuple(left[index] + right[index] for index in range(3))

        def subtract(left, right):
            return tuple(left[index] - right[index] for index in range(3))

        hall_type = spglib.get_spacegroup_type(int(affine_child.hall_number))
        if hall_type is None:
            raise ValueError(
                f"no space-group type for Hall number {affine_child.hall_number}"
            )
        hall_symbol = str(hall_type.hall_symbol or "")
        centering = next(
            (character for character in hall_symbol.upper() if character in "PABCIFR"),
            None,
        )
        if centering is None:
            raise ValueError("target Hall symbol has no conventional centering")
        centering_translations = translation_cosets(
            centering_primitive_matrix(centering),
            identity_matrix(),
        )

        operations: list[ExactSeitzOperation] = []
        for operation in affine_child.operations:
            child_rotation = multiply(
                multiply(child_from_parent, operation.rotation),
                parent_from_child,
            )
            parent_translation_about_origin = subtract(
                add(matrix_vector(operation.rotation, origin), operation.translation),
                origin,
            )
            child_translation = matrix_vector(
                child_from_parent,
                parent_translation_about_origin,
            )
            for centering_translation in centering_translations:
                operations.append(
                    ExactSeitzOperation.from_values(
                        child_rotation,
                        add(child_translation, centering_translation),
                    )
                )

        exact_by_key = {operation.key: operation for operation in operations}
        if len(exact_by_key) != len(operations):
            raise ValueError(
                "affine subgroup expansion produced duplicate conventional operations"
            )
        database = spglib.get_symmetry_from_database(int(affine_child.hall_number))
        if database is None:
            raise ValueError(
                f"no symmetry database entry for Hall number {affine_child.hall_number}"
            )
        database_operations = tuple(
            ExactSeitzOperation.from_values(
                rotation.tolist(),
                tuple(as_fraction(float(value)) for value in translation),
            )
            for rotation, translation in zip(
                database["rotations"],
                database["translations"],
                strict=True,
            )
        )
        if set(exact_by_key) != {operation.key for operation in database_operations}:
            raise ValueError(
                "conjugated affine subgroup disagrees with its target Hall setting"
            )
        return tuple(exact_by_key[key] for key in sorted(exact_by_key))

    def _exact_displacive_child(
        self,
        subgroup: SubgroupInfo,
        operations: tuple[ExactSeitzOperation, ...],
    ) -> tuple[
        Structure,
        tuple[ChildAtom, ...],
        tuple[str, ...],
        tuple[ParentChildSiteMapping, ...],
        dict[str, ParentOrbitType],
    ]:
        """Build one exact emitted child order and prove its parent-site images.

        The structure builder fixes the row order.  We then apply the exact
        origin change ``x_parent = x_child @ B + q`` and rationally recover
        each child coordinate.  Target Hall operations must close that atom
        set exactly; otherwise export remains fail-closed.
        """

        if self.structure is None or not self.symmetry_info:
            raise RuntimeError("parent structure is not loaded")
        basis = rational_matrix(subgroup.basis_vectors or identity_matrix())
        origin_values = tuple(as_fraction(value) for value in (subgroup.origin or (0, 0, 0)))
        if len(origin_values) != 3:
            raise ValueError("subgroup origin must contain three exact coordinates")
        origin = origin_values
        inverse_basis = np.asarray(inverse(basis), dtype=float)
        origin_float = np.asarray([float(value) for value in origin], dtype=float)
        base = self._supercell_for_subgroup(subgroup)
        coordinates = np.mod(
            np.asarray(base.frac_coords, dtype=float) - origin_float @ inverse_basis,
            1.0,
        )
        reference = Structure(
            lattice=base.lattice,
            species=[site.species for site in base],
            coords=coordinates,
            coords_are_cartesian=False,
            site_properties={
                name: list(values) for name, values in base.site_properties.items()
            },
            labels=[site.label for site in base],
        )
        atom_ids = tuple(f"child-{index + 1:06d}" for index in range(len(reference)))
        tolerance = float(self.cfg.affine_exact_cartesian_tolerance_angstrom)
        exact_coordinates: list[tuple[Fraction, Fraction, Fraction] | None] = [
            None
        ] * len(reference)
        remaining = set(range(len(reference)))
        representatives: list[str] = []
        while remaining:
            representative_index = min(remaining)
            representative_site = reference[representative_index]
            representative_coordinate = tuple(
                as_fraction(float(value)) % 1
                for value in representative_site.frac_coords
            )
            generated = {
                operation.operate(representative_coordinate)
                for operation in operations
            }
            assignments: dict[tuple[Fraction, Fraction, Fraction], int] = {}
            for coordinate in sorted(generated):
                coordinate_float = np.asarray(
                    [float(value) for value in coordinate], dtype=float
                )
                candidates: list[tuple[float, int]] = []
                for index in remaining:
                    if reference[index].species != representative_site.species:
                        continue
                    distance = float(
                        reference.lattice.get_distance_and_image(
                            coordinate_float,
                            reference[index].frac_coords,
                        )[0]
                    )
                    if distance <= tolerance:
                        candidates.append((distance, index))
                candidates.sort(key=lambda item: (item[0], item[1]))
                if not candidates:
                    raise ValueError(
                        "exact child orbit has no matching emitted atom within "
                        f"{tolerance:g} Å"
                    )
                if len(candidates) > 1:
                    raise ValueError(
                        "exact child orbit maps ambiguously to the emitted atom order"
                    )
                assignments[coordinate] = candidates[0][1]
            members = set(assignments.values())
            if len(members) != len(generated):
                raise ValueError("exact child orbit does not map bijectively to atoms")
            if assignments.get(representative_coordinate) != representative_index:
                raise ValueError(
                    "exact child representative did not retain its emitted atom identity"
                )
            representatives.append(atom_ids[representative_index])
            for coordinate, index in assignments.items():
                exact_coordinates[index] = coordinate
            remaining.difference_update(members)
        if any(coordinate is None for coordinate in exact_coordinates):
            raise ValueError("exact child orbits do not cover every emitted atom")
        atoms = tuple(
            ChildAtom.from_values(
                atom_id,
                coordinate,
                site.species,
            )
            for atom_id, coordinate, site in zip(
                atom_ids,
                exact_coordinates,
                reference,
                strict=True,
            )
            if coordinate is not None
        )

        sites = list(self.symmetry_info.get("wyckoff_sites") or [])
        parent_orbit_by_index: dict[int, str] = {}
        site_by_orbit: dict[str, dict] = {}
        for site in sites:
            orbit_id = str(site.get("orbit_id") or "").strip()
            if not orbit_id or orbit_id in site_by_orbit:
                raise ValueError("parent physical orbit identities are missing or duplicated")
            site_by_orbit[orbit_id] = site
            for raw_index in site.get("equivalent_indices") or ():
                index = int(raw_index)
                if index in parent_orbit_by_index:
                    raise ValueError("one parent atom belongs to multiple physical orbits")
                parent_orbit_by_index[index] = orbit_id
        if set(parent_orbit_by_index) != set(range(len(self.structure))):
            raise ValueError("parent physical orbit metadata does not cover every atom")

        basis_float = np.asarray(basis, dtype=float)
        mappings: list[ParentChildSiteMapping] = []
        child_ids_by_orbit: dict[str, list[str]] = {
            orbit_id: [] for orbit_id in site_by_orbit
        }
        for atom_id, atom, child_site in zip(atom_ids, atoms, reference, strict=True):
            child_fractional = np.asarray(
                [float(value) for value in atom.frac], dtype=float
            )
            parent_unwrapped = child_fractional @ basis_float + origin_float
            candidates: list[tuple[float, int, tuple[int, int, int]]] = []
            for parent_index, parent_site in enumerate(self.structure):
                if child_site.species != parent_site.species:
                    continue
                residual, translation_array = (
                    self.structure.lattice.get_distance_and_image(
                        parent_unwrapped,
                        parent_site.frac_coords,
                    )
                )
                residual = float(residual)
                if residual <= tolerance:
                    candidates.append(
                        (
                            residual,
                            parent_index,
                            tuple(int(value) for value in translation_array),
                        )
                    )
            candidates.sort(key=lambda item: (item[0], item[1], item[2]))
            if not candidates:
                raise ValueError(
                    f"child atom {atom_id!r} has no exact parent-site image"
                )
            if len(candidates) > 1 and abs(candidates[1][0] - candidates[0][0]) <= 1.0e-12:
                raise ValueError(
                    f"child atom {atom_id!r} has an ambiguous parent-site image"
                )
            _residual, parent_index, translation = candidates[0]
            orbit_id = parent_orbit_by_index[parent_index]
            child_ids_by_orbit[orbit_id].append(atom_id)
            mappings.append(
                ParentChildSiteMapping(
                    child_atom_id=atom_id,
                    parent_site_index=parent_index,
                    parent_orbit_id=orbit_id,
                    parent_cell_translation=translation,
                )
            )

        ordered_sites = sorted(
            enumerate(sites),
            key=lambda item: (
                int(item[1].get("display_order", item[0])),
                item[0],
            ),
        )
        species_counts: dict[str, int] = {}
        orbit_types: dict[str, ParentOrbitType] = {}
        used_labels: set[str] = set()
        for type_index, (_site_index, site) in enumerate(ordered_sites, start=1):
            orbit_id = str(site["orbit_id"])
            species = str(site.get("species") or site.get("element") or "X")
            species_counts[species] = species_counts.get(species, 0) + 1
            label = str(
                site.get("display_label") or f"{species}{species_counts[species]}"
            ).strip()
            if not label or label in used_labels:
                raise ValueError("parent orbit export type labels are missing or duplicated")
            used_labels.add(label)
            orbit_types[orbit_id] = ParentOrbitType(
                type_index,
                label,
                tuple(child_ids_by_orbit[orbit_id]),
            )
        return reference, atoms, tuple(representatives), tuple(mappings), orbit_types

    def _displacive_export_data_for_subgroup(
        self,
        subgroup: SubgroupInfo,
        lifted: dict[str, np.ndarray],
        *,
        final_structure: Structure | None = None,
    ) -> DisplaciveExportData | None:
        """Assemble the shared writer contract from verified production facts."""

        if not lifted:
            return None
        if self.structure is None or not self.symmetry_info:
            raise RuntimeError("parent structure is not loaded")
        modes = tuple(self.distortion_modes or ())
        if not modes:
            raise ValueError("displacement arrays have no source mode records")
        mode_by_key = {self._mode_session_key(mode): mode for mode in modes}
        missing_records = sorted(set(lifted).difference(mode_by_key))
        if missing_records:
            raise ValueError(
                "displacement arrays have no source mode records: "
                + ", ".join(missing_records)
            )
        # Only columns actually emitted by this writer belong to its identity
        # gate.  A different cached/scoped mode must not veto an otherwise
        # complete candidate, nor may it be passed to the mapping contract.
        emitted_modes = tuple(mode_by_key[key] for key in lifted)
        unresolved = self._mode_identity_failures(emitted_modes)
        if unresolved:
            raise UnresolvedModeIdentityError(unresolved)

        parent_group = parent_affine_group(
            self.structure,
            symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance_degrees=self.cfg.symmetry_angle_tolerance_degrees,
        )
        affine_child = embedding_from_identity(subgroup, parent_group)
        if affine_child.hall_number is None:
            raise ValueError("subgroup embedding has no exact Hall setting")
        operations = self._exact_child_operations(subgroup, affine_child)
        reference, atoms, representatives, mappings, orbit_types = (
            self._exact_displacive_child(subgroup, operations)
        )
        atom_ids = tuple(atom.atom_id for atom in atoms)
        route_payload = repr(self._method3_embedding_guard_key(subgroup)).encode(
            "utf-8"
        )
        frame_id = "displacive-child-v1:" + hashlib.sha256(route_payload).hexdigest()
        mapped = self._dist_mapper.map_microscopic_columns_to_supercell(
            self.structure,
            list(self.symmetry_info.get("wyckoff_sites") or []),
            emitted_modes,
            subgroup.basis_vectors or [list(row) for row in identity_matrix()],
            subgroup_context=subgroup,
            frame_id=frame_id,
            atom_ids=atom_ids,
            cartesian_tolerance=self.cfg.symmetry_cartesian_tolerance_angstrom,
            subgroup_operations=affine_child.operations,
            subgroup_translation_lattice=[
                [float(value) for value in row] for row in affine_child.lattice
            ],
        )
        mapped_keys = tuple(column.amplitude_key for column in mapped)
        if set(lifted) != set(mapped_keys):
            missing = sorted(set(mapped_keys).difference(lifted))
            extra = sorted(set(lifted).difference(mapped_keys))
            raise ValueError(
                "current displacement cache differs from verified ISO columns; "
                f"missing={missing}, extra={extra}"
            )
        current_columns = tuple(
            replace(
                column,
                displacements=np.asarray(lifted[column.amplitude_key], dtype=float),
                relation_status="unverified",
                current_over_canonical_signed_scale=None,
                cartesian_residual_angstrom=None,
            )
            for column in mapped
        )
        verified_current = verify_mapped_microscopic_columns(
            mapped,
            current_columns,
            np.asarray(reference.lattice.matrix, dtype=float),
        )
        if any(
            not np.isclose(
                float(column.current_over_canonical_signed_scale),
                1.0,
                rtol=2.0e-10,
                atol=2.0e-12,
            )
            for column in verified_current
        ):
            raise ValueError(
                "current displacement amplitudes are not expressed in the canonical ISO column basis"
            )

        labels = self._mode_labels_now()
        raw_modes = tuple(
            RawModeColumn.from_values(
                column.mode_id,
                index,
                frame_id,
                atom_ids,
                column.displacements,
                ModeScaleProvenance(
                    source=(
                        "iso_microscopic_with_verified_bush_domain_extension"
                        if column.domain_extension_evidence is not None
                        else "iso_microscopic_display_distortion"
                    ),
                    convention=(
                        "mapped_unmixed_verified_full_domain_extension"
                        if column.domain_extension_evidence is not None
                        else "mapped_unmixed_fractional_source_column"
                    ),
                    direction_resolved=True,
                    sign_resolved=True,
                ),
                mode_identity=column.mode_identity,
                microscopic_provenance=column.provenance,
                microscopic_domain_extension=column.domain_extension_evidence,
                label=labels[column.amplitude_key],
            )
            for index, column in enumerate(mapped)
        )
        hall_type = spglib.get_spacegroup_type(int(affine_child.hall_number))
        if hall_type is None:
            raise ValueError("subgroup Hall setting cannot be resolved")
        target_centering = next(
            (
                character
                for character in str(hall_type.hall_symbol or "").upper()
                if character in "PABCIFR"
            ),
            None,
        )
        if target_centering is None:
            raise ValueError("subgroup Hall setting has no conventional centering")
        frame = ExactChildFrame.from_values(
            frame_id,
            operations,
            np.asarray(reference.lattice.matrix, dtype=float),
            conventional_basis=identity_matrix(),
            primitive_translation_basis=centering_primitive_matrix(target_centering),
        )
        model = build_displacive_cif_model(
            frame,
            atoms,
            representatives,
            raw_modes,
            amplitudes=(0.0,) * len(raw_modes),
            cartesian_tolerance_angstrom=self.cfg.affine_exact_cartesian_tolerance_angstrom,
        )

        emitted_final = None
        if final_structure is not None:
            if len(final_structure) != len(reference):
                raise ValueError("generated final structure has a different child atom count")
            if any(
                final_site.species != reference_site.species
                for final_site, reference_site in zip(
                    final_structure,
                    reference,
                    strict=True,
                )
            ):
                raise ValueError(
                    "generated final structure does not retain the emitted atom order"
                )
            expected_lattice = np.asarray(reference.lattice.matrix, dtype=float)
            if not np.allclose(
                np.asarray(final_structure.lattice.matrix, dtype=float),
                expected_lattice,
                rtol=2.0e-12,
                atol=float(self.cfg.affine_exact_cartesian_tolerance_angstrom),
            ):
                raise ValueError(
                    "generated final structure uses a different child lattice"
                )
            basis_inverse = np.asarray(
                inverse(rational_matrix(subgroup.basis_vectors or identity_matrix())),
                dtype=float,
            )
            origin_shift = np.asarray(
                [float(value) for value in (subgroup.origin or (0, 0, 0))],
                dtype=float,
            ) @ basis_inverse
            emitted_final = Structure(
                lattice=reference.lattice,
                species=[site.species for site in final_structure],
                coords=(
                    np.asarray(final_structure.frac_coords, dtype=float)
                    - origin_shift
                )
                % 1.0,
                coords_are_cartesian=False,
                site_properties={
                    name: list(values)
                    for name, values in final_structure.site_properties.items()
                },
                labels=[site.label for site in final_structure],
            )
            difference_rows: list[np.ndarray] = []
            for reference_site, final_site in zip(
                reference,
                emitted_final,
                strict=True,
            ):
                _distance, image = reference.lattice.get_distance_and_image(
                    reference_site.frac_coords,
                    final_site.frac_coords,
                )
                difference_rows.append(
                    np.asarray(final_site.frac_coords, dtype=float)
                    + np.asarray(image, dtype=float)
                    - np.asarray(reference_site.frac_coords, dtype=float)
                )
            difference = np.asarray(difference_rows, dtype=float)
            columns = np.column_stack(
                [
                    (
                        np.asarray(mode.displacements, dtype=float)
                        * float(normfactor)
                    ).reshape(-1)
                    for mode, normfactor in zip(
                        model.canonical_modes,
                        model.norms.normfactors,
                        strict=True,
                    )
                ]
            )
            amplitudes, _residuals, rank, _singular = np.linalg.lstsq(
                columns,
                difference.reshape(-1),
                rcond=1.0e-12,
            )
            if rank != len(raw_modes):
                raise ValueError("generated final structure has an underdetermined mode basis")
            reconstructed = (columns @ amplitudes).reshape((-1, 3))
            cartesian_residual = (difference - reconstructed) @ np.asarray(
                reference.lattice.matrix, dtype=float
            )
            if float(np.max(np.linalg.norm(cartesian_residual, axis=1))) > float(
                self.cfg.affine_exact_cartesian_tolerance_angstrom
            ):
                raise ValueError(
                    "generated final structure is not an exact combination of the canonical ISO modes"
                )
            model = build_displacive_cif_model(
                frame,
                atoms,
                representatives,
                raw_modes,
                amplitudes=tuple(float(value) for value in amplitudes),
                cartesian_tolerance_angstrom=self.cfg.affine_exact_cartesian_tolerance_angstrom,
            )

        embedding_id = str(
            getattr(subgroup, "_method3_embedding_id", "") or ""
        ).strip()
        embedding = ExactParentChildEmbedding(
            self.structure,
            int(self.symmetry_info["space_group_number"]),
            rational_matrix(subgroup.basis_vectors or identity_matrix()),
            tuple(as_fraction(value) for value in (subgroup.origin or (0, 0, 0))),
            mappings,
            embedding_id=embedding_id,
        )
        context_kinds = {
            mode.microscopic_provenance.source_subgroup_context_kind
            for mode in raw_modes
        }
        primary_direction_selectors = {
            mode.microscopic_provenance
            .source_subgroup_primary_direction_selector
            for mode in raw_modes
        }
        expected_context = (
            "exact_fixed_space"
            if str(getattr(subgroup, "_method3_route_resolution", ""))
            == "exact_fixed_space"
            else "single_irrep"
        )
        if context_kinds != {expected_context}:
            raise ValueError(
                "canonical modes mix or mismatch subgroup context kinds"
            )
        if (
            len(primary_direction_selectors) != 1
            or (
                expected_context == "exact_fixed_space"
                and primary_direction_selectors != {None}
            )
            or (
                expected_context == "single_irrep"
                and None in primary_direction_selectors
            )
        ):
            raise ValueError(
                "canonical modes do not satisfy the subgroup primary-direction contract"
            )
        subgroup_identity = DisplaciveSubgroupIdentity(
            embedding,
            irrep_label=(
                None if expected_context == "exact_fixed_space"
                else str(subgroup.irrep_label)
            ),
            opd_symbol=(
                None if expected_context == "exact_fixed_space"
                else str(subgroup.opd_symbol)
            ),
            primary_direction_selector=next(iter(primary_direction_selectors)),
            target_space_group_number=int(subgroup.space_group_number),
            target_hall_number=int(affine_child.hall_number),
            context_kind=expected_context,
        )
        subgroup._displacive_embedding_id = embedding.embedding_id
        mapping_by_atom = {mapping.child_atom_id: mapping for mapping in mappings}
        representative_counts: dict[str, int] = {}
        representative_labels: dict[str, str] = {}
        for orbit in model.orbits:
            orbit_id = mapping_by_atom[orbit.representative_atom_id].parent_orbit_id
            representative_counts[orbit_id] = representative_counts.get(orbit_id, 0) + 1
            representative_labels[orbit.representative_atom_id] = (
                f"{orbit_types[orbit_id].type_label}_{representative_counts[orbit_id]}"
            )
        return DisplaciveExportData(
            model=model,
            embedding=embedding,
            subgroup_identity=subgroup_identity,
            reference_structure=reference,
            final_structure=emitted_final,
            representative_labels=representative_labels,
            parent_orbit_types=orbit_types,
            coordinate_tolerance_angstrom=self.cfg.affine_exact_cartesian_tolerance_angstrom,
        )

    def _spec_for_subgroup(
        self,
        subgroup,
        *,
        use_current_modes: bool,
        use_generated_structure: bool,
        note: str = "",
        folder_name: str = "",
    ) -> SubgroupExportSpec:
        self._require_verified_parent_setting()
        structure = self._supercell_for_subgroup(subgroup)
        cif_structure = None
        if use_generated_structure and self.distorted_structure is not None:
            cif_structure = self.distorted_structure
        lifted = self._lifted_mode_displacements(subgroup) if use_current_modes else {}
        parent_sg = 0
        parent_sym = ""
        wyckoff = None
        if self.symmetry_info:
            parent_sg = int(self.symmetry_info.get("space_group_number") or 0)
            parent_sym = hm_symbol(parent_sg) or str(
                self.symmetry_info.get("space_group_symbol") or ""
            )
            wyckoff = self.symmetry_info.get("wyckoff_sites")
        wyckoff_lines = self.parent_wyckoff_display() or None
        strain_data = self._strain_export_data_for_subgroup(subgroup)
        displacive_data = self._displacive_export_data_for_subgroup(
            subgroup,
            lifted,
            final_structure=cif_structure,
        )
        if displacive_data is not None:
            return SubgroupExportSpec(
                subgroup=subgroup,
                note=note,
                folder_name=folder_name,
                parent_wyckoff_sites=wyckoff,
                parent_wyckoff_lines=list(wyckoff_lines) if wyckoff_lines else None,
                distortion_types=list(self.distortion_types or []),
                strain_data=strain_data,
                displacive_data=displacive_data,
                require_verified_displacive_data=True,
            )
        return SubgroupExportSpec(
            subgroup=subgroup,
            structure=structure,
            parent_structure=self.structure,
            parent_sg=parent_sg,
            parent_symbol=parent_sym,
            mode_displacements_sc=lifted or None,
            mode_labels=self._mode_labels_now() if lifted else None,
            note=note,
            folder_name=folder_name,
            cif_structure=cif_structure,
            parent_wyckoff_sites=wyckoff,
            parent_wyckoff_lines=list(wyckoff_lines) if wyckoff_lines else None,
            distortion_types=list(self.distortion_types or []),
            strain_data=strain_data,
            # The current core retains only anonymous lifted arrays.  Until it
            # also retains the exact emitted frame, atom identities and
            # unmixed ISO source-column provenance, writer entry points must
            # reject those arrays instead of presenting them as authoritative.
            require_verified_displacive_data=bool(lifted),
        )

    def _is_parametric_subgroup(self, subgroup) -> bool:
        """Subgroup of a parametric k point (LD/DT, …)."""
        return bool(getattr(subgroup, "k_parameters", None))

    @staticmethod
    def _export_folder_method(
        export_method: int | str | None,
        use_opd_line_folders: bool | None,
    ) -> int:
        """Resolve the new Method-aware naming API and its legacy switch."""
        if use_opd_line_folders is not None:
            legacy_method = 1 if use_opd_line_folders else 2
            if (
                export_method is not None
                and parse_export_method(export_method) != legacy_method
            ):
                raise ValueError(
                    "export_method conflicts with legacy use_opd_line_folders"
                )
            return legacy_method
        return parse_export_method(export_method)

    def _collect_export_specs(
        self,
        items: list,
        formats: list[str],
        compute_missing_modes: bool,
        *,
        export_method: int | str | None = None,
        reserved_folder_names: set[str] | None = None,
        use_opd_line_folders: bool | None = None,
        number_of_independent_modulations: int | None = None,
        allow_unresolved_identities: bool = False,
        progress_callback: Callable[[dict], None] | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> list[SubgroupExportSpec] | ExportBatchPlan:
        """为每个子群准备导出规格；结束后恢复会话 Distortion 状态。"""
        self._require_verified_parent_setting()
        method = self._export_folder_method(export_method, use_opd_line_folders)
        need_modes = any(fmt != "cif" for fmt in formats)
        snap = self._snapshot_distortion_state()
        export_nmod = (
            snap["number_of_independent_modulations"]
            if number_of_independent_modulations is None
            else int(number_of_independent_modulations)
        )
        if export_nmod < 0:
            raise ValueError("number_of_independent_modulations must be >= 0")
        selected = snap["selected_subgroup"]
        selected_identity = (
            self._method3_embedding_guard_key(selected)
            if selected is not None
            else None
        )
        export_types = list(snap["distortion_types"])
        export_scope = self._copy_distortion_scope(snap["distortion_scope"])
        used: set[str] = set(reserved_folder_names or ())
        specs: list[SubgroupExportSpec] = []
        failures: list[ExportCandidateFailure] = []
        skipped: list[ExportCandidateFailure] = []
        batch_started = time.monotonic()

        def _progress(payload: dict) -> None:
            if progress_callback is None:
                return
            try:
                progress_callback(dict(payload))
            except Exception:  # noqa: BLE001 - observers cannot corrupt export
                return

        _progress({
            "phase": "preparing",
            "total": len(items),
            "completed": 0,
            "successful": 0,
            "ineligible": 0,
            "failed": 0,
            "estimated_remaining_seconds": None,
        })
        try:
            for position, sg in enumerate(items, start=1):
                if cancel_check is not None and cancel_check():
                    for remaining_position, remaining in enumerate(
                        items[position - 1:],
                        start=position,
                    ):
                        skipped.append(self._export_candidate_failure(
                            remaining_position,
                            remaining,
                            ExportCancelledError(
                                "batch cancelled before this candidate was prepared"
                            ),
                        ))
                    _progress({
                        "phase": "cancelled",
                        "total": len(items),
                        "completed": position - 1,
                        "successful": len(specs),
                        "ineligible": len(skipped),
                        "failed": len(failures),
                        "estimated_remaining_seconds": 0.0,
                    })
                    break
                sequence = int(getattr(sg, "index", -1)) + 1
                if sequence < 1:
                    sequence = position
                folder = unique_folder_name(
                    sg,
                    used,
                    export_method=method,
                    sequence=sequence,
                )
                note = ""
                target_identity = self._method3_embedding_guard_key(sg)
                is_current = (
                    selected_identity is not None
                    and target_identity == selected_identity
                )
                target_context = self._mode_context_key(
                    sg,
                    export_nmod,
                    export_types,
                    export_scope,
                )
                current_modes_match = (
                    is_current
                    and snap["mode_cache_key"] == target_context
                )
                generated_structure_matches = (
                    current_modes_match
                    and snap["generated_structure_cache_key"] == target_context
                    and snap["distorted_structure"] is not None
                )
                computed = False
                attempted_compute = False
                candidate_status = "ready"
                _progress({
                    "phase": "candidate",
                    "total": len(items),
                    "completed": position - 1,
                    "position": position,
                    "candidate_index": int(getattr(sg, "index", -1)),
                    "irrep_label": str(getattr(sg, "irrep_label", "") or ""),
                    "opd_symbol": str(getattr(sg, "opd_symbol", "") or ""),
                    "successful": len(specs),
                    "ineligible": len(skipped),
                    "failed": len(failures),
                })
                try:
                    if need_modes and compute_missing_modes and not current_modes_match:
                        attempted_compute = True
                        result = self.search_method_2(
                            sg.index,
                            distortion_type=export_types,
                            number_of_independent_modulations=export_nmod,
                            candidates=[sg],
                            distortion_scope=export_scope,
                        )
                        result_subgroup = getattr(result, "subgroup", None)
                        selected_subgroup = getattr(self, "_selected_subgroup", None)
                        for resolved in (result_subgroup, selected_subgroup):
                            if resolved is not None and (
                                self._method3_embedding_guard_key(resolved)
                                != target_identity
                            ):
                                raise RuntimeError(
                                    "Method 2 returned a different scientific subgroup "
                                    "than the requested export candidate"
                                )
                        if result_subgroup is None and selected_subgroup is None:
                            raise RuntimeError(
                                "Method 2 did not identify the computed subgroup"
                            )
                        computed = True
                        if (
                            not self.mode_displacements
                            and not self.mode_occupancies
                            and not self.mode_displacements_sc
                        ):
                            note = (
                                "no displacement modes for this path "
                                "(empty BUSH/smodes table)"
                            )
                    spec = self._spec_for_subgroup(
                        sg,
                        use_current_modes=(
                            need_modes and (current_modes_match or computed)
                        ),
                        # A generated structure carries amplitudes in the cached
                        # mode basis.  Reusing it after an nmod-triggered mode
                        # recomputation would mix two incompatible bases in one
                        # export, so only reuse it with the matching cache.
                        use_generated_structure=(
                            generated_structure_matches and not computed
                        ),
                        note=note,
                        folder_name=folder,
                    )
                    if (
                        getattr(spec, "mode_displacements_sc", None)
                        and getattr(spec, "displacive_data", None) is None
                    ):
                        raise ValueError(
                            "non-empty displacement modes lack verified identity, "
                            "frame, atom-order and source-column provenance"
                        )
                    specs.append(spec)
                except UnresolvedModeIdentityError as exc:
                    candidate_status = "ineligible"
                    skipped.append(
                        self._export_candidate_failure(position, sg, exc)
                    )
                except Exception as exc:  # noqa: BLE001 - aggregate before publication
                    candidate_status = "failed"
                    failures.append(
                        self._export_candidate_failure(position, sg, exc)
                    )
                finally:
                    if attempted_compute:
                        self._restore_distortion_state(snap)
                    elapsed = max(time.monotonic() - batch_started, 0.0)
                    average = elapsed / position
                    _progress({
                        "phase": "candidate_complete",
                        "total": len(items),
                        "completed": position,
                        "position": position,
                        "candidate_index": int(getattr(sg, "index", -1)),
                        "candidate_status": candidate_status,
                        "successful": len(specs),
                        "ineligible": len(skipped),
                        "failed": len(failures),
                        "elapsed_seconds": elapsed,
                        "estimated_remaining_seconds": average * (len(items) - position),
                    })
        finally:
            self._restore_distortion_state(snap)
        if failures:
            raise ExportPreparationError(method, [*failures, *skipped])
        if skipped and not allow_unresolved_identities:
            raise ExportPreparationError(method, skipped)
        if allow_unresolved_identities:
            return ExportBatchPlan(tuple(specs), tuple(skipped))
        return specs

    def _render_export_batch(
        self,
        method: int,
        specs: list[SubgroupExportSpec],
        formats: list[str],
    ) -> list[tuple[SubgroupExportSpec, tuple[tuple[str, bytes], ...]]]:
        """Pre-render an entire disk batch before any candidate is published."""

        rendered: list[
            tuple[SubgroupExportSpec, tuple[tuple[str, bytes], ...]]
        ] = []
        failures: list[ExportCandidateFailure] = []
        for position, spec in enumerate(specs, start=1):
            try:
                payloads = render_subgroup_files(spec, formats)
            except Exception as exc:  # noqa: BLE001 - aggregate before publication
                failures.append(
                    self._export_candidate_failure(position, spec.subgroup, exc)
                )
            else:
                rendered.append((spec, payloads))
        if failures:
            raise ExportPreparationError(method, failures)
        return rendered

    def _publish_export_batch(
        self,
        method: int,
        folder_root: Path,
        rendered_batch: list[
            tuple[SubgroupExportSpec, tuple[tuple[str, bytes], ...]]
        ],
        root_payloads: tuple[tuple[str, bytes], ...] = (),
    ) -> list[Path]:
        """Stage and transactionally publish a complete disk export batch.

        A new or pre-created empty root keeps the legacy direct layout and is
        exposed as a complete root.  An existing non-empty root gets one
        immutable, content-addressed ``*.ready`` version directory, exposed by
        a single rename.  Readers following the ready-manifest protocol never
        observe only a prefix of a batch, including if the producer process
        exits before or after the commit rename.  Existing files are never
        overwritten.
        """

        parent = folder_root.parent
        parent_existed = parent.exists()
        parent.mkdir(parents=True, exist_ok=True)
        # Windows can reject a directory move from the destination's parent
        # into an already-existing child directory under restrictive ACLs.
        # Stage appends inside that existing root, under a non-ready name that
        # discovery deliberately ignores, so the commit rename stays within
        # one directory and remains atomic.
        root_was_directory = (
            folder_root.exists()
            and folder_root.is_dir()
            and not folder_root.is_symlink()
        )
        root_was_empty = root_was_directory and not any(folder_root.iterdir())
        staging_parent = (
            folder_root if root_was_directory and not root_was_empty else parent
        )
        staging = Path(
            tempfile.mkdtemp(prefix=".isodistort-batch-", dir=staging_parent)
        )
        active: tuple[int, SubgroupExportSpec] | None = None
        published_root = folder_root
        publish_guard = None

        def failure_for(
            position: int,
            spec: SubgroupExportSpec,
            exc: Exception,
        ) -> ExportPreparationError:
            return ExportPreparationError(
                method,
                [self._export_candidate_failure(position, spec.subgroup, exc)],
            )

        try:
            folders: list[str] = []
            for position, (spec, payloads) in enumerate(rendered_batch, start=1):
                active = (position, spec)
                folder = str(spec.folder_name or subgroup_label(spec.subgroup))
                if not folder or Path(folder).name != folder:
                    raise ValueError(
                        f"unsafe candidate export folder component: {folder!r}"
                    )
                if folder.casefold() in {value.casefold() for value in folders}:
                    raise ValueError(
                        f"duplicate candidate export folder: {folder!r}"
                    )
                folders.append(folder)
                staged_folder = staging / folder
                staged_folder.mkdir()
                for filename, payload in payloads:
                    name = str(filename)
                    if not name or Path(name).name != name:
                        raise ValueError(
                            f"unsafe export filename component: {name!r}"
                        )
                    (staged_folder / name).write_bytes(bytes(payload))

            root_names: list[str] = []
            for filename, payload in root_payloads:
                name = str(filename)
                if (
                    not name
                    or Path(name).name != name
                    or name.casefold() == _EXPORT_BATCH_MANIFEST.casefold()
                    or name.casefold() in {folder.casefold() for folder in folders}
                    or name.casefold() in {value.casefold() for value in root_names}
                ):
                    raise ValueError(
                        f"unsafe or duplicate export root filename: {name!r}"
                    )
                root_names.append(name)
                (staging / name).write_bytes(bytes(payload))

            guard = _export_publish_lock(folder_root)
            guard.__enter__()
            publish_guard = guard
            if folder_root.exists() or folder_root.is_symlink():
                if not folder_root.is_dir() or folder_root.is_symlink():
                    first_position, first_spec = active or (1, rendered_batch[0][0])
                    raise failure_for(
                        first_position,
                        first_spec,
                        FileExistsError(
                            f"export root is not a directory: {folder_root}"
                        ),
                    )
                # A caller commonly supplies a pre-created empty output
                # directory.  It has no prior state to retain, so publish the
                # original direct candidate layout by replacing that empty
                # placeholder with the fully staged root.  Observers still see
                # only zero or all candidates; a failed rename recreates the
                # empty placeholder.
                if not any(folder_root.iterdir()):
                    folder_root.rmdir()
                    try:
                        staging.rename(folder_root)
                    except Exception as exc:
                        folder_root.mkdir(exist_ok=True)
                        failures = [
                            self._export_candidate_failure(
                                position,
                                spec.subgroup,
                                exc,
                            )
                            for position, (spec, _payloads) in enumerate(
                                rendered_batch,
                                start=1,
                            )
                        ]
                        raise ExportPreparationError(method, failures) from exc
                else:
                    reserved = {
                        name.casefold()
                        for name in _published_export_folder_names(folder_root)
                    }
                    conflicts: list[ExportCandidateFailure] = []
                    for position, ((spec, _payloads), folder) in enumerate(
                        zip(rendered_batch, folders, strict=True),
                        start=1,
                    ):
                        if folder.casefold() in reserved:
                            conflicts.append(
                                self._export_candidate_failure(
                                    position,
                                    spec.subgroup,
                                    FileExistsError(
                                        "target candidate directory exists in a "
                                        f"committed batch: {folder_root / folder}"
                                    ),
                                )
                            )
                    if conflicts:
                        raise ExportPreparationError(method, conflicts)

                    manifest = _export_batch_manifest(
                        method,
                        rendered_batch,
                        folders,
                        root_payloads,
                    )
                    (staging / _EXPORT_BATCH_MANIFEST).write_bytes(manifest)
                    published_root = folder_root / _ready_export_batch_name(manifest)
                    if published_root.exists() or published_root.is_symlink():
                        failures = [
                            self._export_candidate_failure(
                                position,
                                spec.subgroup,
                                FileExistsError(
                                    f"committed export batch exists: {published_root}"
                                ),
                            )
                            for position, (spec, _payloads) in enumerate(
                                rendered_batch,
                                start=1,
                            )
                        ]
                        raise ExportPreparationError(method, failures)
                    try:
                        # This is the sole visibility transition for an append
                        # to an existing non-empty root.  The manifest and
                        # every candidate are already closed and complete.
                        staging.rename(published_root)
                    except Exception as exc:
                        failures = [
                            self._export_candidate_failure(
                                position,
                                spec.subgroup,
                                exc,
                            )
                            for position, (spec, _payloads) in enumerate(
                                rendered_batch,
                                start=1,
                            )
                        ]
                        raise ExportPreparationError(method, failures) from exc
            else:
                try:
                    staging.rename(folder_root)
                except Exception as exc:
                    failures = [
                        self._export_candidate_failure(position, spec.subgroup, exc)
                        for position, (spec, _payloads) in enumerate(
                            rendered_batch,
                            start=1,
                        )
                    ]
                    raise ExportPreparationError(method, failures) from exc

            return [
                published_root / folder / filename
                for (spec, payloads), folder in zip(
                    rendered_batch,
                    folders,
                    strict=True,
                )
                for filename, _payload in payloads
            ] + [published_root / name for name in root_names]
        except ExportPreparationError:
            raise
        except Exception as exc:
            if active is None:
                raise
            position, spec = active
            raise failure_for(position, spec, exc) from exc
        finally:
            if publish_guard is not None:
                publish_guard.__exit__(None, None, None)
            if staging.exists():
                shutil.rmtree(staging, ignore_errors=True)
            if not parent_existed and parent.exists():
                try:
                    parent.rmdir()
                except OSError:
                    pass

    def export_subgroups(
        self,
        dest_dir: str | Path,
        formats: list | str | None = None,
        subgroups: list | None = None,
        compute_missing_modes: bool = False,
        *,
        export_method: int | str | None = None,
        stable_case_id: str | None = None,
        use_opd_line_folders: bool | None = None,
        number_of_independent_modulations: int | None = None,
        progress_callback: Callable[[dict], None] | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> list:
        """
        按指定 Method 批量导出（每个候选一个短目录）。

        Args:
            dest_dir: 输出根目录（其下创建各子群文件夹）
            formats: cif / isoviz / modes / topas（官网第 6 页对应选项）
            subgroups: 默认使用当前会话的子群列表。
            compute_missing_modes: 为非当前子群再跑 Method 2 以填充模式类格式；
                仅 CIF 时不需要。参数 k 点子群走 smodes/(3+d) 完整模式。
            export_method: 决定短目录规则（Method 1/2/3）；默认 Method 2。
            stable_case_id: 可选 Method 3 清单案例号；省略时按候选身份生成。
            use_opd_line_folders: 旧接口兼容；True 等价于 ``export_method=1``。
            number_of_independent_modulations: 批量模式计算使用的 nmod；
                省略时沿用当前 Distortion 会话的值，并在导出后恢复原状态。

        Returns:
            写出的文件路径列表
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()
        method = self._export_folder_method(export_method, use_opd_line_folders)
        fmts = parse_export_formats(formats)
        items = list(subgroups if subgroups is not None else self.subgroups)
        if not items:
            raise RuntimeError(
                f"没有可导出的 Method {method} 子群；请先完成该 Method 的子群计算"
            )
        dest = Path(dest_dir)
        folder_root = dest
        if method == 3:
            folder_root = dest / method3_case_folder(items, stable_case_id)
        reserved = _published_export_folder_names(folder_root)
        collected = self._collect_export_specs(
            items,
            fmts,
            compute_missing_modes,
            export_method=method,
            reserved_folder_names=reserved,
            number_of_independent_modulations=number_of_independent_modulations,
            allow_unresolved_identities=(method == 2 and compute_missing_modes),
            progress_callback=progress_callback,
            cancel_check=cancel_check,
        )
        if isinstance(collected, ExportBatchPlan):
            specs = list(collected.specs)
            skipped = collected.skipped
        else:
            specs = collected
            skipped = ()
        root_payloads = _export_failure_report_payloads(
            method,
            len(items),
            len(specs),
            tuple(skipped),
        )
        rendered_batch = self._render_export_batch(method, specs, fmts)
        paths = self._publish_export_batch(
            method,
            folder_root,
            rendered_batch,
            root_payloads,
        )
        print(t("export.done", n=len(paths)))
        return paths

    def export_subgroups_zip(
        self,
        formats: list | str | None = None,
        subgroups: list | None = None,
        compute_missing_modes: bool = False,
        wrapping: str | None = None,
        *,
        export_method: int | str | None = None,
        stable_case_id: str | None = None,
        use_opd_line_folders: bool | None = None,
        number_of_independent_modulations: int | None = None,
        progress_callback: Callable[[dict], None] | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> bytes:
        """批量导出为 ZIP 字节（不读写 output_dir，避免混入无关文件）。

        Method 1/2 的候选目录直接位于 ZIP 根；Method 3 自动增加稳定案例目录。
        ``wrapping`` 非空时覆盖 Method 3 的自动案例目录，也可为其它 Method 加前缀。
        ``stable_case_id`` 可为 Method 3 指定清单案例号；省略时从候选身份生成。
        ``number_of_independent_modulations`` 显式控制整包候选的 nmod；省略时
        沿用当前 Distortion 会话值，绝不把不同 nmod 的缓存模式混入同一 ZIP。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()
        method = self._export_folder_method(export_method, use_opd_line_folders)
        fmts = parse_export_formats(formats)
        items = list(subgroups if subgroups is not None else self.subgroups)
        if not items:
            raise RuntimeError(
                f"没有可导出的 Method {method} 子群；请先完成该 Method 的子群计算"
            )
        collected = self._collect_export_specs(
            items,
            fmts,
            compute_missing_modes,
            export_method=method,
            number_of_independent_modulations=number_of_independent_modulations,
            allow_unresolved_identities=(method == 2 and compute_missing_modes),
            progress_callback=progress_callback,
            cancel_check=cancel_check,
        )
        if isinstance(collected, ExportBatchPlan):
            specs = list(collected.specs)
            skipped = collected.skipped
        else:
            specs = collected
            skipped = ()
        effective_wrapping = wrapping
        if method == 3 and not effective_wrapping:
            effective_wrapping = method3_case_folder(items, stable_case_id)
        rendered_batch = self._render_export_batch(method, specs, fmts)
        root_payloads = _export_failure_report_payloads(
            method,
            len(items),
            len(specs),
            tuple(skipped),
        )
        try:
            return build_export_zip(
                specs,
                fmts,
                wrapping=effective_wrapping,
                rendered_batch=rendered_batch,
                root_payloads=root_payloads,
            )
        except Exception as exc:
            raise ExportPreparationError(
                method,
                [
                    self._export_candidate_failure(position, spec.subgroup, exc)
                    for position, spec in enumerate(specs, start=1)
                ],
            ) from exc

    # ================================================================
    # 畴变体
    # ================================================================

    def generate_domains(self) -> list:
        """生成所有畴变体描述（官网 Domains 输出）。

        畴总数 = 子群在母相中的指数；需要先选择路径（select_path / Method 2）。
        """
        if self.phase_path is None or self._selected_subgroup is None:
            raise RuntimeError(t("err.domains_need_path"))

        domains = self._domain_gen.generate_domains(
            self.phase_path, self._selected_subgroup
        )
        print(t("domains.found", n=len(domains)))
        return domains

    # ================================================================
    # ISODISTORT Search Method 1-4
    # ================================================================

    def search_method_1(self,
                        distortion_types: str | list[str] | None = None,
                        crystal_system: str | None = None,
                        subgroup_space_group: int | None = None,
                        lattice: list[list[float]] | None = None,
                        maximal_subgroup_only: bool = False,
                        lattice_kind: str = "conventional"):
        """
        Method 1: Search over all special k points.

        支持多条件同时过滤（逻辑 AND，与官网一致）：
        - lattice：官网 Conventional lattice / Primitive lattice 下拉所选
          格点类（3x3 矩阵，均以母相惯用坐标表达）
        - lattice_kind：所选下拉框。Primitive 按子群原胞格点匹配。

        枚举结果按会话缓存，重复调用秒回。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()
        if lattice_kind not in {"conventional", "primitive"}:
            raise ValueError(f"Unknown Method 1 lattice kind: {lattice_kind}")

        query = Method1Query(
            distortion_types=distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            # The selector denotes one lattice class, not every sublattice of
            # the selected basis. Match it after type filtering below so the
            # conventional and primitive selectors use their respective cells.
            lattice=None,
            maximal_subgroup_only=maximal_subgroup_only,
            parent_rotations=[
                r.tolist() for r in self._parent_rotations()
            ],
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_1_search(
            parent_sg, query, subgroups=self._ensure_special_subgroups()
        )
        result = self._filter_method1_by_types(
            result, query.distortion_types or self.distortion_types
        )
        if lattice is not None:
            selected = np.asarray(lattice, dtype=float)
            if selected.shape != (3, 3) or abs(np.linalg.det(selected)) < 1e-8:
                raise ValueError("Method 1 lattice must be a nonsingular 3x3 matrix")
            result = [
                item for item in result
                if self._same_lattice_orbit(
                    (self._centering_matrix(
                        _centering_letter(item.subgroup.space_group_number)
                    ) @ np.asarray(item.subgroup.basis_vectors, dtype=float))
                    if lattice_kind == "primitive" else item.subgroup.basis_vectors,
                    selected,
                )
            ]

        # 记录过滤后的候选，供 Method 2 使用
        self.subgroups = [item.subgroup for item in result]
        print(t("method1.result", n=len(result)))
        return result

    def search_method_2(self,
                        subgroup_idx: int,
                        distortion_type: str | list[str] | None = None,
                        number_of_independent_modulations: int = 0,
                        number_of_superposed_irreps: int = 1,
                        *,
                        candidates: list[SubgroupInfo] | None = None,
                        distortion_scope: dict | None = None):
        """
        Method 2: General method - search over specific k points.

        在 Method 1 候选（或 list_subgroups 枚举）中按序号选择子群，
        特殊 k 用 iso DISPLAY BUSH；参数 k 用 smodes + 子群对称性（(3+d)/公度锁定）。
        k 点 / IR / OPD 由所选子群（SubgroupInfo）自身携带。
        number_of_independent_modulations：0 = 三维锁定完整模式（含谐波与次级 IR）；
        n>=1 = 只保留 n 个独立调制波矢的谐波（(3+n)D）。

        Args:
            subgroup_idx: Index within the selected candidate pool.
            distortion_type: Enabled distortion types; project defaults when omitted.
            number_of_independent_modulations: 0 = 3D lock-in; n = (3+n)D harmonics.
            number_of_superposed_irreps: Must currently be 1.  Coupled Method 2
                needs a component-aware route rather than a union of single-IR
                candidates and therefore fails closed until that path exists.
            candidates: Explicit candidate pool. Pass this whenever a caller keeps
                more than one Method result table; otherwise the session default
                ``self.subgroups`` is used for backward compatibility.
            distortion_scope: Explicit per-type species scope snapshot.  Batch
                export passes this together with ``distortion_type`` so one
                candidate cannot be computed against mutable UI state.
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self._require_verified_parent_setting()

        requested_nsup = _validated_single_irrep_count(
            number_of_superposed_irreps
        )

        types = normalize_distortion_types(
            getattr(self, "distortion_types", None)
            if distortion_type is None
            else distortion_type
        )
        active_scope = self._copy_distortion_scope(
            getattr(self, "distortion_scope", {})
            if distortion_scope is None
            else distortion_scope
        )
        if candidates is None and not self.subgroups:
            self.list_subgroups(distortion_type=distortion_type)
        candidate_pool = list(candidates) if candidates is not None else self.subgroups
        if not candidate_pool:
            raise RuntimeError("没有可用于 Method 2 的子群候选")
        matching_candidates = [
            candidate
            for candidate in candidate_pool
            if candidate.index == subgroup_idx
        ]
        if not matching_candidates:
            raise ValueError(f"Subgroup index {subgroup_idx} not found")
        if len(matching_candidates) != 1:
            raise ValueError(
                f"Subgroup index {subgroup_idx} is ambiguous in the candidate pool"
            )
        selected_candidate = matching_candidates[0]
        if (
            self._method3_embedding_guard_key(selected_candidate) in getattr(
                self, "_unresolved_method3_embedding_keys", set()
            )
            or getattr(selected_candidate, "_method3_route_resolution", "")
            == "affine_only_unresolved_coupled_route"
        ):
            raise RuntimeError(
                "This affine embedding has no resolved single-IR or coupled-IR "
                "second-stage route; Method 2 mode calculation is not implemented"
            )

        self.number_of_independent_modulations = int(
            number_of_independent_modulations or 0
        )
        query = Method2Query(
            subgroup_idx=subgroup_idx,
            distortion_type=types,
            number_of_independent_modulations=self.number_of_independent_modulations,
            number_of_superposed_irreps=requested_nsup,
        )

        parent_sg = self.symmetry_info["space_group_number"]
        # 按作用域限制 BUSH / smodes 的 Wyckoff 位置（避免重复计算）
        bush_types = [tp for tp in types if tp in ("displacive", "rotational")]
        scoped_species = self._union_scope_species(
            bush_types,
            distortion_scope=active_scope,
        )
        scoped_letters = self._letters_for_species(scoped_species)
        scoped_orbit_ids = self._orbit_ids_for_species(scoped_species)
        kpoints = []
        iso_backend = getattr(self, "_iso", None)
        try:
            if iso_backend is not None:
                kpoints = iso_backend.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001
            kpoints = []
        result = self._search.method_2_search(
            parent_sg, candidate_pool, query,
            wyckoff_letters=scoped_letters,
            wyckoff_orbit_ids=scoped_orbit_ids or None,
            structure=self.structure,
            wyckoff_sites=(self.symmetry_info or {}).get("wyckoff_sites"),
            smodes=getattr(self, "_smodes", None),
            kpoints=kpoints,
            symmetry_info=self.symmetry_info,
        )
        result_subgroup = getattr(result, "subgroup", None)
        if result_subgroup is None or (
            self._method3_embedding_guard_key(result_subgroup)
            != self._method3_embedding_guard_key(selected_candidate)
        ):
            raise RuntimeError(
                "Method 2 search returned a subgroup with a different scientific identity"
            )
        meta = getattr(result, "metadata", None) or {}
        if meta.get("supercell_displacements"):
            self.mode_displacements_sc = dict(meta["supercell_displacements"])
            self._mode_label_overrides = dict(meta.get("mode_labels") or {})
        else:
            self.mode_displacements_sc = {}
            self._mode_label_overrides = {}

        # 记录路径与模式，供 Distortion Page 使用
        self.phase_path = PhasePath.from_subgroup(
            parent_sg, result.subgroup, types
        )
        self.phase_path.k_vector = self._resolve_k_vector(
            result.subgroup.k_point_label,
            list(result.subgroup.k_parameters or []),
        )
        self.phase_path.validate()
        self._selected_subgroup = result.subgroup
        self.distortion_modes = self._compute_scoped_modes(
            parent_sg,
            result.subgroup,
            types,
            raw_modes=result.modes,
            distortion_scope=active_scope,
        )
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )
        self._install_special_bush_supercell_modes(result.subgroup)
        self._sync_parametric_session_keys()
        emitted_keys = set(self.mode_displacements_sc or self.mode_displacements)
        emitted_modes = [
            mode
            for mode in self.distortion_modes
            if self._mode_session_key(mode) in emitted_keys
        ]
        identity_failures = self._mode_identity_failures(emitted_modes)
        if emitted_modes and identity_failures:
            identity_status = "unresolved"
        elif emitted_modes:
            identity_status = "verified"
        else:
            identity_status = "not_applicable"
        result.metadata["mode_identity_status"] = identity_status
        result.metadata["mode_identity_failures"] = [
            asdict(failure) for failure in identity_failures
        ]
        result.metadata["export_ready"] = not identity_failures
        result.subgroup._mode_identity_status = identity_status
        result.subgroup._mode_identity_failures = tuple(identity_failures)
        self._mode_cache_key = self._mode_context_key(
            result.subgroup,
            self.number_of_independent_modulations,
            types,
            active_scope,
        )
        self.distorted_structure = None
        self._generated_structure_cache_key = None
        print(t("method2.result", idx=subgroup_idx,
                n=len(self.mode_displacements) + len(self.mode_occupancies)))
        return result

    def search_method_3(self,
                        distortion_types: str | list[str] | None = None,
                        point_group: str | None = None,
                        space_group_type: int | None = None,
                        supercell_basis: list[list[str | int | float]] | None = None,
                        direct_sublattice_centering: str | None = None,
                        lattice_type: str = "direct",
                        generate_if_missing: bool = False,
                        resolve_coupled_routes: bool = True,
                        include_affine_only_diagnostics: bool = False):
        """
        Method 3: Search over arbitrary k points for a specified point group and supercell.

        若 point_group 与 space_group_type 同时提供，按官网规则优先采用
        space_group_type。lattice_type 为官网 radio（direct/reciprocal）；
        本地引擎暂不支持 reciprocal（倒易超格）模式，会给出明确错误。

        Default 按目标空间群的默认 Bravais centering 解释；显式
        P/A/B/C/I/F/R 分别按所选 centering 解释为精确 primitive lattice。
        非母相原胞子格会额外按公度条件推断参数 k（如 LD ``g=1/6``）并枚举；
        ``generate_if_missing`` 与 Method 2 的 GenDB 开关相同。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        self._require_verified_parent_setting()

        if lattice_type != "direct":
            raise ValueError(
                "本地引擎暂不支持 reciprocal（倒易空间超格）模式，"
                "请使用 direct（实空间子格）。官网该选项的完整实现在后续版本中支持。"
            )

        normalized_types = normalize_distortion_types(
            distortion_types or self.distortion_types
        )
        displacive_species = (
            tuple(sorted(self._scope_species("displacive")))
            if "displacive" in normalized_types
            else None
        )

        query = Method3Query(
            distortion_types=normalized_types,
            point_group=point_group,
            space_group_type=space_group_type,
            supercell_basis=supercell_basis,
            direct_sublattice_centering=direct_sublattice_centering,
            lattice_type=lattice_type,
            parent_rotations=[rotation.tolist() for rotation in self._parent_rotations()],
            parent_structure=self.structure,
            symmetry_tolerance=self.cfg.symmetry_cartesian_tolerance_angstrom,
            symmetry_angle_tolerance_degrees=(
                self.cfg.symmetry_angle_tolerance_degrees
            ),
            affine_exact_tolerance=(
                self.cfg.affine_exact_cartesian_tolerance_angstrom
            ),
            fractional_coordinate_tolerance=(
                self.cfg.fractional_coordinate_tolerance
            ),
            generate_if_missing=generate_if_missing,
            resolve_coupled_routes=resolve_coupled_routes,
            include_affine_only_diagnostics=include_affine_only_diagnostics,
            displacive_species=displacive_species,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_3_search(parent_sg, query)

        result = self._filter_method3_routes_by_types(
            result, query.distortion_types or self.distortion_types
        )

        # 记录过滤后的候选，供 Method 2（search_method_2）使用（与 Method 1 一致）
        self._unresolved_method3_embedding_keys = set()
        for item in result:
            resolution = getattr(item, "route_resolution", "known_single_ir")
            item.subgroup._method3_route_resolution = resolution
            embedding_id = str(getattr(item, "embedding_id", "") or "").strip()
            if embedding_id:
                item.subgroup._method3_embedding_id = embedding_id
            if resolution == "affine_only_unresolved_coupled_route":
                key = self._method3_embedding_guard_key(item.subgroup)
                self._unresolved_method3_embedding_keys.add(key)
        self.subgroups = [item.subgroup for item in result]
        print(t("method3.result", n=len(result)))
        return result

    def search_method_4(self,
                        distorted_cif_path: str | Path,
                        atom_matching_method: str = "nearest-site",
                        robust_distance_threshold: float = 0.25,
                        provided_origin_shift: list[float] | None = None):
        """
        Method 4: Mode decomposition of a distorted structure.

        要求已经通过 Method 2 或 select_path 计算出可用模式。
        畸变结构可以是母相的超胞（官网 Method 4 的常规情形）：原子数
        不一致时自动把母相与模式位移提升到畸变结构的超胞坐标系再分解。
        """
        if self.structure is None:
            raise RuntimeError("请先加载母相结构 (load_structure)")
        if not self.mode_displacements and not self.mode_displacements_sc:
            raise RuntimeError("请先通过 select_path 或 search_method_2 计算模式")

        distorted_structure = read_cif(distorted_cif_path)
        parent, mode_disp = self._method4_reference_modes(distorted_structure)
        daughter_sg = read_cif_space_group_number(distorted_cif_path)
        if daughter_sg is None:
            daughter_sg = int(SpacegroupAnalyzer(
                distorted_structure,
                symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
                angle_tolerance=self.cfg.symmetry_angle_tolerance_degrees,
            ).get_space_group_number())
        daughter_multiplicity = self._centering_multiplicity(daughter_sg)
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 1)
        parent_multiplicity = self._centering_multiplicity(parent_sg)
        parent_primitive_volume = float(self.structure.lattice.volume) / parent_multiplicity
        child_primitive_volume = float(parent.lattice.volume) / daughter_multiplicity
        supercell_size = child_primitive_volume / parent_primitive_volume
        if self._selected_subgroup is not None:
            parent_to_child_basis = np.asarray(
                self._selected_subgroup.basis_vectors or np.eye(3),
                dtype=float,
            )
        elif len(distorted_structure) == len(self.structure):
            parent_to_child_basis = np.eye(3)
        else:
            parent_to_child_basis = self._resolve_distorted_supercell_basis(
                distorted_structure
            )
        strain_data = None
        configured_types = getattr(self, "distortion_types", None) or []
        if isinstance(configured_types, str):
            configured_types = [configured_types]
        selected_types = {
            str(value).strip().lower()
            for value in configured_types
        }
        if "strain" in selected_types:
            if self._selected_subgroup is None:
                raise RuntimeError(
                    "Method 4 canonical strain amplitudes require the exact "
                    "selected subgroup embedding"
                )
            strain_data = self._strain_export_data_for_subgroup(
                self._selected_subgroup
            )
        query = Method4Query(
            atom_matching_method=atom_matching_method,
            robust_distance_threshold=robust_distance_threshold,
            provided_origin_shift=provided_origin_shift,
            primitive_cell_multiplicity=daughter_multiplicity,
            supercell_size=supercell_size,
            reference_parent_lattice=np.asarray(
                self.structure.lattice.matrix, dtype=float
            ).tolist(),
            parent_to_child_basis=parent_to_child_basis.tolist(),
            strain_mode_labels=(
                [label for _mode, label, _amplitude in strain_data.items()]
                if strain_data is not None
                else None
            ),
            strain_mode_irrep_labels=(
                [str(mode.irrep_label) for mode in strain_data.result.modes]
                if strain_data is not None
                else None
            ),
            strain_mode_q_raw=(
                [mode.q_raw.tolist() for mode in strain_data.result.modes]
                if strain_data is not None
                else None
            ),
            strain_mode_q_unit=(
                [mode.q_unit.tolist() for mode in strain_data.result.modes]
                if strain_data is not None
                else None
            ),
            strain_mode_normfactors=(
                [float(mode.normfactor) for mode in strain_data.result.modes]
                if strain_data is not None
                else None
            ),
        )

        result = self._search.method_4_decompose(
            parent,
            distorted_structure,
            mode_disp,
            query,
        )

        print(t("method4.result", n=len(result.amplitudes), rms=result.rms_residual))
        return result

    @staticmethod
    def _centering_multiplicity(space_group_number: int) -> int:
        letter = _centering_letter(int(space_group_number))
        if letter == "F":
            return 4
        if letter == "R":
            return 3
        if letter in {"A", "B", "C", "I"}:
            return 2
        return 1

    def _method4_reference_modes(self, distorted: Structure
                                 ) -> tuple[Structure, dict[str, np.ndarray]]:
        """Return the undistorted child cell and modes in that cell's basis.

        Method 2 can select a rotated index-one cell as well as a true
        supercell.  Atom-count equality therefore does *not* imply that the
        uploaded daughter is expressed in the original parent setting.  The
        selected subgroup cell is the scientific reference in both cases.
        """
        if self.structure is None:
            raise RuntimeError("请先加载母相结构 (load_structure)")
        parent_modes = {
            key: value["displacements"]
            for key, value in self.mode_displacements.items()
        }
        subgroup = self._selected_subgroup
        if subgroup is not None:
            reference = self._supercell_for_subgroup(subgroup)
            if self.mode_displacements_sc:
                modes = {
                    key: self._validated_supercell_displacement(
                        key,
                        values,
                        len(reference),
                    )
                    for key, values in self.mode_displacements_sc.items()
                }
                return reference, modes
            basis = subgroup.basis_vectors or np.eye(3).tolist()
            k_vec = self.phase_path.k_vector if self.phase_path is not None else None
            return self._dist_engine.lift_mode_displacements(
                self.structure,
                basis,
                parent_modes,
                k_vector=k_vec,
            )

        if len(distorted) == len(self.structure):
            k_vec = self.phase_path.k_vector if self.phase_path is not None else None
            if k_vec is None:
                return self.structure, parent_modes
            return self._dist_engine.lift_mode_displacements(
                self.structure, np.eye(3), parent_modes, k_vector=k_vec,
            )
        basis = self._resolve_distorted_supercell_basis(distorted)
        k_vec = self.phase_path.k_vector if self.phase_path is not None else None
        return self._dist_engine.lift_mode_displacements(
            self.structure,
            basis,
            parent_modes,
            k_vector=k_vec,
        )

    def _resolve_distorted_supercell_basis(self, distorted: Structure
                                           ) -> np.ndarray:
        """确定母相 -> 畸变结构的超胞基矢 B（Ld = B @ Lp，原子数比 = |det B|）。

        优先采用当前相变路径的子群基矢（本地生成畸变结构的标准情形）；
        否则从两个晶格矩阵反推（B = Ld @ Lp⁻¹，应为近整数矩阵）。
        两者都无法给出一致的整数超胞关系时，报出明确错误。
        """
        n_parent = len(self.structure)
        n_dist = len(distorted)
        candidates: list[np.ndarray] = []
        if self.phase_path is not None and self.phase_path.basis_vectors:
            candidates.append(
                np.asarray(self.phase_path.supercell_basis(), dtype=float))
        lp = np.asarray(self.structure.lattice.matrix, dtype=float)
        ld = np.asarray(distorted.lattice.matrix, dtype=float)
        try:
            candidates.append(ld @ np.linalg.inv(lp))
        except np.linalg.LinAlgError:
            pass
        for b in candidates:
            if b.shape != (3, 3):
                continue
            if not np.allclose(b, np.round(b), atol=1e-4):
                continue
            det = abs(round(float(np.linalg.det(b))))
            if det >= 1 and n_parent * det == n_dist:
                return np.round(b).astype(int)
        raise ValueError(
            "畸变结构与母相原子数不一致，且无法确定超胞关系"
            f"（母相 {n_parent} 原子，畸变 {n_dist} 原子）；"
            "请确认畸变结构是当前母相/所选子群的超胞。")
