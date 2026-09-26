"""Read-only scientific audit of saved Method 1/2/3 output pairs.

Candidates are paired by crystallographic header identity rather than directory
order or the printed basis. Point-group-equivalent bases may use different rows
while generating the same lattice.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import platform
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.io.cif import CifParser

WORKSPACE = Path(__file__).resolve().parents[3]
VALIDATOR_ROOT = WORKSPACE / "ISODISTORT_VALIDATE"
if str(VALIDATOR_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATOR_ROOT))

from isodistort_validate.compare_cif import compare_cif  # noqa: E402

DEFAULT_OUTPUT_ROOT = WORKSPACE / "output_compare"
OFFICIAL_DIR = "官网"
LOCAL_DIR = "现有网页版交互"
CORE_LOCAL_FILES = (
    "subgroup.cif",
    "data.isoviz",
    "topas.str",
    "Complete modes details.txt",
)
LIVE12_CHECKPOINT_SCHEMA = 2
ISO_DATA_BASENAMES = (
    "const.dat",
    "data_diperiodic.txt",
    "data_images.txt",
    "data_irreps.txt",
    "data_isotropy.txt",
    "data_little.txt",
    "data_magnetic.txt",
    "data_space.txt",
    "data_ssg.txt",
    "data_ssgmag.txt",
    "data_wyckoff.txt",
)
RUNTIME_DISTRIBUTIONS = ("numpy", "pymatgen", "spglib", "PyYAML")
LIVE12_GROUP_COUNT = 4


@dataclass(frozen=True)
class CandidateHeader:
    ir: str
    identity: str
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    origin: tuple[Fraction, Fraction, Fraction]
    parametric_k: bool


@dataclass
class AuditIssue:
    severity: str
    code: str
    case: str
    candidate: str
    message: str


@dataclass
class AuditSummary:
    output_root: str
    groups: dict[str, dict[str, int]] = field(default_factory=dict)
    paired_candidates: int = 0
    direct_cif_matches: int = 0
    equivalent_setting_matches: int = 0
    basis_equivalent_pairs: int = 0
    basis_exact_pairs: int = 0
    bases_normalized_to_official: int = 0
    format_files_checked: int = 0
    mode_files_compared: int = 0
    mode_count_matches: int = 0
    issues: list[AuditIssue] = field(default_factory=list)

    @property
    def failures(self) -> int:
        return sum(issue.severity == "error" for issue in self.issues)

    @property
    def warnings(self) -> int:
        return sum(issue.severity == "warning" for issue in self.issues)


_MATRIX_RE = re.compile(r"\(([^()]*)\)")
_IR_RE = re.compile(r"^#\s*IR:\s*([^,\s]+)", re.MULTILINE)
_PARAMETER_RE = re.compile(
    r"^#\s*k point:.*\b([a-z])\s*=\s*[-+]?\d+(?:/\d+)?\s*$",
    re.MULTILINE,
)


def _fraction(value: str) -> Fraction:
    return Fraction(value.strip())


def _canonical_number(value: Any) -> str:
    """Format an integer/float/string as the exact short rational used by the site."""
    number = Fraction(str(value)).limit_denominator(10_000)
    return str(number.numerator) if number.denominator == 1 else f"{number.numerator}/{number.denominator}"


def _canonical_origin(origin: tuple[Fraction, Fraction, Fraction]) -> str:
    """Return an origin modulo integer parent-lattice translations."""
    return ",".join(_canonical_number(value % 1) for value in origin)


def _parse_vectors(value: str, expected: int) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    vectors: list[tuple[Fraction, Fraction, Fraction]] = []
    for raw in _MATRIX_RE.findall(value):
        parts = tuple(_fraction(item) for item in raw.split(","))
        if len(parts) == 3:
            vectors.append(parts)
    if len(vectors) != expected:
        raise ValueError(f"expected {expected} three-component vectors, got {len(vectors)}: {value}")
    return tuple(vectors)


def _candidate_line(text: str) -> str:
    for line in text.splitlines():
        if line.startswith("#") and "basis={" in line and ", origin=" in line:
            return line.removeprefix("#").strip()
    raise ValueError("CIF header has no candidate line containing basis and origin")


def parse_candidate_header(path: Path) -> CandidateHeader:
    text = path.read_text(encoding="utf-8-sig", errors="strict")
    line = _candidate_line(text)
    ir_match = _IR_RE.search(text)
    if ir_match is None:
        raise ValueError("CIF header has no IR line")
    basis_match = re.search(r"basis=\{(.*?)\},\s*origin=", line)
    origin_match = re.search(r"origin=(\([^)]*\))", line)
    if basis_match is None or origin_match is None:
        raise ValueError("candidate line has an invalid basis/origin")
    basis = _parse_vectors(basis_match.group(1), 3)
    origin = _parse_vectors(origin_match.group(1), 1)[0]
    identity = re.sub(r",?\s*basis=\{.*?\}(?=,\s*origin=)", "", line)
    identity = re.sub(r"origin=\([^)]*\)", f"origin=({_canonical_origin(origin)})", identity)
    identity = " ".join(identity.split())
    return CandidateHeader(
        ir=ir_match.group(1),
        identity=f"{ir_match.group(1)} | {identity}",
        basis=basis,
        origin=origin,
        parametric_k=_PARAMETER_RE.search(text) is not None,
    )


def _determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _inverse(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    det = _determinant(matrix)
    if det == 0:
        raise ValueError("candidate basis is singular")
    a, b, c = matrix
    cofactors = (
        (b[1] * c[2] - b[2] * c[1], -(b[0] * c[2] - b[2] * c[0]), b[0] * c[1] - b[1] * c[0]),
        (-(a[1] * c[2] - a[2] * c[1]), a[0] * c[2] - a[2] * c[0], -(a[0] * c[1] - a[1] * c[0])),
        (a[1] * b[2] - a[2] * b[1], -(a[0] * b[2] - a[2] * b[0]), a[0] * b[1] - a[1] * b[0]),
    )
    return tuple(tuple(cofactors[col][row] / det for col in range(3)) for row in range(3))


def _multiply(
    first: tuple[tuple[Fraction, ...], ...],
    second: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(sum(first[row][k] * second[k][col] for k in range(3)) for col in range(3))
        for row in range(3)
    )


def bases_generate_same_lattice(
    first: tuple[tuple[Fraction, ...], ...],
    second: tuple[tuple[Fraction, ...], ...],
) -> bool:
    """Return whether two row bases differ by an integer unimodular matrix."""
    transform = _multiply(first, _inverse(second))
    return all(value.denominator == 1 for row in transform for value in row) and abs(_determinant(transform)) == 1


def _load_structure(path: Path):
    structures = CifParser(str(path), occupancy_tolerance=100).parse_structures(primitive=False)
    if len(structures) != 1:
        raise ValueError(f"expected one structure, got {len(structures)}")
    return structures[0]


_MATCHER = StructureMatcher(
    ltol=1e-5,
    stol=2e-5,
    angle_tol=1e-3,
    # Equivalent conventional settings can have beta and 180-beta, or a
    # different centered-cell representative.  Primitive reduction removes
    # that presentation choice before the strict site-distance comparison.
    primitive_cell=True,
    scale=False,
    attempt_supercell=False,
    comparator=ElementComparator(),
)


def _equivalent_structures(local: Path, official: Path) -> bool:
    local_structure = _load_structure(local)
    official_structure = _load_structure(official)
    if local_structure.composition != official_structure.composition:
        return False
    # Six-decimal CIF lattice values accumulate into a larger absolute volume
    # error for large supercells.  A 1e-6 relative gate remains tighter than
    # the precision carried by those lattice fields without rejecting the same
    # cell solely because the two writers rounded intermediate lengths.
    if not np.isclose(local_structure.volume, official_structure.volume, rtol=1e-6, atol=1e-5):
        return False
    return bool(_MATCHER.fit(local_structure, official_structure))


def _index_candidates(root: Path) -> tuple[dict[str, Path], list[tuple[str, list[Path]]]]:
    grouped: dict[str, list[Path]] = defaultdict(list)
    for cif in sorted(root.rglob("subgroup.cif")):
        grouped[parse_candidate_header(cif).identity].append(cif.parent)
    unique: dict[str, Path] = {}
    for key, values in grouped.items():
        if len(values) == 1:
            unique[key] = values[0]
            continue
        # A saved reference set can contain a copied/misnamed directory. Keep
        # the one whose directory label agrees with its internal IR/OPD header,
        # and report the remaining aliases as a reference conflict.
        ir, body = (part.strip() for part in key.split("|", 1))
        opd = body.split(None, 1)[0]
        matching = [path for path in values if path.name.startswith(f"{ir} {opd}")]
        if len(matching) == 1:
            unique[key] = matching[0]
    duplicates = [(key, values) for key, values in grouped.items() if len(values) > 1]
    return unique, duplicates


def _add_issue(
    summary: AuditSummary,
    severity: str,
    code: str,
    case: str,
    candidate: str,
    message: str,
) -> None:
    summary.issues.append(AuditIssue(severity, code, case, candidate, message))


def _check_formats(summary: AuditSummary, case: str, candidate: str, folder: Path, parametric: bool) -> None:
    for name in CORE_LOCAL_FILES:
        path = folder / name
        if not path.is_file():
            _add_issue(summary, "error", "missing-local-format", case, candidate, f"missing {name}")
            continue
        summary.format_files_checked += 1
        if path.stat().st_size == 0:
            _add_issue(summary, "error", "empty-local-format", case, candidate, f"empty {name}")
    markers = {
        "data.isoviz": ("!isoversion", "!atomcoordlist"),
        "topas.str": ("space_group", "'{{{distorted parameters"),
        "Complete modes details.txt": ("Complete modes details", "Subgroup:"),
    }
    for name, required in markers.items():
        path = folder / name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8-sig", errors="replace")
        missing = [marker for marker in required if marker not in text]
        if missing:
            _add_issue(
                summary,
                "error",
                "invalid-local-format",
                case,
                candidate,
                f"{name} missing markers: {', '.join(missing)}",
            )
    if not parametric:
        return
    details = folder / "Complete modes details.txt"
    if details.is_file():
        text = details.read_text(encoding="utf-8-sig", errors="replace")
        if _count_local_displacive_modes(text) == 0:
            _add_issue(
                summary,
                "error",
                "empty-parametric-modes",
                case,
                candidate,
                "Complete modes details.txt contains no displacive parameter-k modes",
            )
    isoviz = folder / "data.isoviz"
    if isoviz.is_file():
        text = isoviz.read_text(encoding="utf-8-sig", errors="replace")
        marker = text.lower().find("!displacivemodelist")
        tail = text[marker:] if marker >= 0 else ""
        if marker < 0 or re.search(r"(?m)^\s*\d+\s+\d+\s+[-+0-9.]", tail) is None:
            _add_issue(
                summary,
                "error",
                "empty-parametric-isoviz-modes",
                case,
                candidate,
                "data.isoviz contains no displacive parameter-k mode rows",
            )
    topas = folder / "topas.str"
    if topas.is_file():
        text = topas.read_text(encoding="utf-8-sig", errors="replace")
        if re.search(r"(?m)^\s*prm\s+!a\d+\b", text) is None:
            _add_issue(
                summary,
                "error",
                "empty-parametric-topas-modes",
                case,
                candidate,
                "topas.str contains no active displacive mode parameters",
            )


def _count_local_displacive_modes(text: str) -> int:
    """Count actual mode definitions, excluding prose/table headings."""
    return len(
        re.findall(
            r"(?m)^Mode\s+(?!definitions\b|vectors\b|amplitudes\b)\S+",
            text,
        )
    )


def _count_official_displacive_modes(folder: Path) -> int | None:
    html_files = sorted(folder.glob("*complete modes details*.html"))
    if not html_files:
        return None
    text = html_files[0].read_text(encoding="utf-8-sig", errors="replace")
    return len(re.findall(r":dsp\].*?\bnormfactor\s*=", text, flags=re.IGNORECASE))


def _compare_mode_counts(
    summary: AuditSummary,
    case: str,
    candidate: str,
    official_folder: Path,
    local_folder: Path,
) -> None:
    official_count = _count_official_displacive_modes(official_folder)
    details = local_folder / "Complete modes details.txt"
    if official_count is None or not details.is_file():
        return
    local_text = details.read_text(encoding="utf-8-sig", errors="replace")
    local_count = _count_local_displacive_modes(local_text)
    summary.mode_files_compared += 1
    if local_count == official_count:
        summary.mode_count_matches += 1
        return
    _add_issue(
        summary,
        "error",
        "displacive-mode-count-mismatch",
        case,
        candidate,
        f"official={official_count}; local={local_count}",
    )


def audit_output_root(output_root: Path) -> AuditSummary:
    output_root = output_root.expanduser().resolve()
    if not output_root.is_dir():
        raise FileNotFoundError(f"output comparison root not found: {output_root}")
    summary = AuditSummary(output_root=str(output_root))
    for case_root in sorted(path for path in output_root.iterdir() if path.is_dir()):
        # Method 3 rows are affine embeddings and can legitimately have no
        # single-IR header.  They are audited by audit_method3_downloads.py and
        # compare_method3_official_local.py, which bind the saved result page to
        # the manifest and compare exact affine subgroup operators.  Reusing the
        # Method 1/2 IR-path parser here produced false "CIF header has no IR"
        # failures for scientifically valid coupled/route-less embeddings.
        for method in ("Method1", "Method2"):
            official_root = case_root / OFFICIAL_DIR / method
            local_root = case_root / LOCAL_DIR / method
            if not official_root.exists() and not local_root.exists():
                continue
            case = f"{case_root.name}/{method}"
            if not official_root.is_dir() or not local_root.is_dir():
                _add_issue(summary, "error", "missing-side", case, "", "official or local method directory is missing")
                continue
            try:
                official, official_duplicates = _index_candidates(official_root)
                local, local_duplicates = _index_candidates(local_root)
            except (OSError, ValueError) as exc:
                _add_issue(summary, "error", "inventory-parse", case, "", str(exc))
                continue
            summary.groups[case] = {"official": len(official), "local": len(local)}
            official_conflict_names = {
                path.name for _key, paths in official_duplicates for path in paths
            }
            for key, paths in official_duplicates:
                _add_issue(
                    summary,
                    "warning",
                    "reference-conflict",
                    case,
                    key,
                    "duplicate internal identity: " + "; ".join(map(str, paths)),
                )
            for key, paths in local_duplicates:
                _add_issue(summary, "error", "duplicate-identity", case, key, "; ".join(map(str, paths)))
            for key in sorted(official.keys() - local.keys()):
                _add_issue(summary, "error", "missing-local-candidate", case, key, official[key].name)
            for key in sorted(local.keys() - official.keys()):
                if local[key].name in official_conflict_names:
                    _add_issue(
                        summary,
                        "warning",
                        "reference-conflict",
                        case,
                        key,
                        f"official directory {local[key].name!r} contains another candidate's CIF",
                    )
                else:
                    _add_issue(summary, "error", "extra-local-candidate", case, key, local[key].name)
            for key in sorted(official.keys() & local.keys()):
                summary.paired_candidates += 1
                official_folder = official[key]
                local_folder = local[key]
                official_cif = official_folder / "subgroup.cif"
                local_cif = local_folder / "subgroup.cif"
                official_header = parse_candidate_header(official_cif)
                local_header = parse_candidate_header(local_cif)
                try:
                    same_lattice = bases_generate_same_lattice(local_header.basis, official_header.basis)
                except ValueError as exc:
                    _add_issue(summary, "error", "invalid-basis", case, key, str(exc))
                    same_lattice = False
                if same_lattice:
                    summary.basis_equivalent_pairs += 1
                    if local_header.basis == official_header.basis:
                        summary.basis_exact_pairs += 1
                    else:
                        # All subsequent semantic checks use the official file
                        # as the target setting.  This is an explicit
                        # local->official unimodular normalization, not a
                        # tolerance-based waiver of the representation.
                        summary.bases_normalized_to_official += 1
                else:
                    _add_issue(
                        summary,
                        "error",
                        "non-equivalent-basis",
                        case,
                        key,
                        f"official={official_folder.name}; local={local_folder.name}",
                    )
                try:
                    direct = compare_cif(local_cif, official_cif, ignore_atom_order=True)
                    if direct.structure_equal:
                        summary.direct_cif_matches += 1
                    elif same_lattice and _equivalent_structures(local_cif, official_cif):
                        summary.equivalent_setting_matches += 1
                    else:
                        _add_issue(
                            summary,
                            "error",
                            "cif-semantic-mismatch",
                            case,
                            key,
                            "; ".join(direct.details["issues"]),
                        )
                except (OSError, TypeError, ValueError) as exc:
                    _add_issue(summary, "error", "cif-compare-error", case, key, str(exc))
                _check_formats(summary, case, key, local_folder, local_header.parametric_k)
                _compare_mode_counts(summary, case, key, official_folder, local_folder)
    return summary


def _payload(summary: AuditSummary) -> dict[str, Any]:
    result = asdict(summary)
    result.update(failures=summary.failures, warnings=summary.warnings, passed=summary.failures == 0)
    return result


def _print_summary(summary: AuditSummary) -> None:
    print("Method output scientific audit")
    for case, counts in summary.groups.items():
        print(f"  {case}: official={counts['official']} local={counts['local']}")
    print(f"Paired candidates: {summary.paired_candidates}")
    print(f"Direct CIF semantic matches: {summary.direct_cif_matches}")
    print(f"Equivalent-setting CIF matches: {summary.equivalent_setting_matches}")
    print(f"Equivalent basis pairs: {summary.basis_equivalent_pairs}")
    print(f"Exact official basis representations: {summary.basis_exact_pairs}")
    print(f"Bases normalized to official representation: {summary.bases_normalized_to_official}")
    print(f"Local core files checked: {summary.format_files_checked}")
    print(
        "Displacive mode-count matches: "
        f"{summary.mode_count_matches}/{summary.mode_files_compared}"
    )
    for code, count in sorted(Counter(issue.code for issue in summary.issues).items()):
        print(f"  {code}: {count}")
    for issue in summary.issues[:40]:
        candidate = f" [{issue.candidate}]" if issue.candidate else ""
        print(
            f"{issue.severity.upper()} {issue.case} {issue.code}{candidate}: "
            f"{issue.message}"
        )
    if len(summary.issues) > 40:
        print(f"... {len(summary.issues) - 40} additional issues (use --json for full details)")
    print(f"Conclusion: {'PASS' if summary.failures == 0 else 'FAIL'}")


def run_live_method3_smoke() -> list[dict[str, Any]]:
    """Run the two fixed Method 3 cases through the real iso/WSL backend.

    The expected fields are a test manifest captured from the official site,
    never consulted by production code.  ZIPs stay in memory so the check does
    not leave downloads or extracted trees behind.
    """
    project_root = WORKSPACE / "ISODISTORT"
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from isocore.api import IsoDistort  # noqa: PLC0415

    cases = (
        {
            "name": "EuAl4 Parent.cif",
            "target": 99,
            "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 6]],
            "expected": {("LD1", "C1", "(0,0,1/6)", "0,0,0")},
        },
        {
            "name": "NdNiO2 own.cif",
            "target": 47,
            "basis": [[-3, 0, 0], [0, 0, 1], [0, 2, 0]],
            "expected": {
                ("Y1", "P1", "(1/3,1/2,0)", "0,0,0"),
                ("Y4", "P1", "(1/3,1/2,0)", "0,1/2,0"),
            },
        },
    )
    reports: list[dict[str, Any]] = []
    for case in cases:
        api = IsoDistort()
        api.load_structure(WORKSPACE / "experiment_data" / case["name"])
        items = api.search_method_3(
            distortion_types=["strain", "displacive"],
            space_group_type=case["target"],
            supercell_basis=case["basis"],
            direct_sublattice_centering="d",
            lattice_type="direct",
            generate_if_missing=True,
        )
        actual = {
            (
                item.subgroup.irrep_label,
                item.subgroup.opd_symbol,
                item.subgroup.official_fields()["k_active"],
                ",".join(_canonical_number(value) for value in item.subgroup.origin),
            )
            for item in items
        }
        if actual != case["expected"]:
            raise AssertionError(
                f"{case['name']} Method 3 candidates differ from the official manifest: "
                f"expected={sorted(case['expected'])}, actual={sorted(actual)}"
            )
        requested_basis = tuple(
            tuple(Fraction(value) for value in row) for row in case["basis"]
        )
        for item in items:
            subgroup = item.subgroup
            if subgroup.space_group_number != case["target"]:
                raise AssertionError(
                    f"{case['name']}: expected SG {case['target']}, got {subgroup.space_group_number}"
                )
            candidate_basis = tuple(
                tuple(Fraction(str(value)) for value in row)
                for row in subgroup.basis_vectors
            )
            if not bases_generate_same_lattice(candidate_basis, requested_basis):
                raise AssertionError(f"{case['name']}: candidate basis does not generate the requested lattice")

        body = api.export_subgroups_zip(
            formats=["cif"],
            subgroups=[item.subgroup for item in items],
            compute_missing_modes=False,
        )
        with zipfile.ZipFile(io.BytesIO(body)) as archive:
            cif_names = [name for name in archive.namelist() if name.endswith("/subgroup.cif")]
            if len(cif_names) != len(items):
                raise AssertionError(
                    f"{case['name']}: ZIP contains {len(cif_names)} CIFs for {len(items)} candidates"
                )
            for name in cif_names:
                cif_text = archive.read(name).decode("utf-8-sig")
                match = re.search(r"(?im)^_symmetry_Int_Tables_number\s+(\d+)\s*$", cif_text)
                if match is None or int(match.group(1)) != case["target"]:
                    raise AssertionError(f"{case['name']}: {name} has the wrong declared space group")
                parsed = CifParser.from_str(cif_text, occupancy_tolerance=100).parse_structures(
                    primitive=False
                )
                if len(parsed) != 1:
                    raise AssertionError(f"{case['name']}: {name} is not a single parseable structure")
        reports.append(
            {
                "case": case["name"],
                "target_space_group": case["target"],
                "candidates": len(items),
                "zip_cifs": len(cif_names),
                "identities": sorted(" | ".join(fields) for fields in actual),
            }
        )
    return reports


def _live_subgroup_identity(subgroup: Any) -> str:
    """Build the same basis-independent identity used for saved CIF headers."""
    line = subgroup.opd_line_body()
    line = re.sub(r",?\s*basis=\{.*?\}(?=,\s*origin=)", "", line)
    origin = tuple(
        Fraction(str(value)).limit_denominator(10_000)
        for value in subgroup.origin
    )
    line = re.sub(
        r"origin=\([^)]*\)",
        f"origin=({_canonical_origin(origin)})",
        line,
    )
    return f"{subgroup.irrep_label} | {' '.join(line.split())}"


def _live12_candidate_fingerprint(subgroup: Any) -> str:
    """Sign all subgroup fields that can steer the downstream mode calculation."""
    fields = (
        "index",
        "space_group_number",
        "space_group_symbol",
        "subgroup_index",
        "size",
        "is_maximal",
        "opd_symbol",
        "opd_vector",
        "basis_vectors",
        "origin",
        "k_point_label",
        "irrep_label",
        "k_parameters",
        "k_coordinates",
        "parent_sg",
        "opd_dir_raw",
        "basis_raw",
        "origin_raw",
        "k_active_raw",
    )
    payload = {
        "opd_line_body": subgroup.opd_line_body(),
        "fields": {name: getattr(subgroup, name, None) for name in fields},
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(encoded, usedforsecurity=False).hexdigest()


def _live12_now() -> str:
    return datetime.now(UTC).isoformat()


def _live12_runtime_versions() -> dict[str, Any]:
    """Describe the interpreter and scientific packages that affect results."""
    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - signatures must survive broken package metadata
            error_type = f"{type(exc).__module__}.{type(exc).__qualname__}"
            distributions[name] = f"<unavailable:{error_type}>"
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
            "build": sys.version,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _live12_configured_path(config_path: Path, value: object, fallback: str) -> Path:
    raw = value if isinstance(value, str) and value.strip() else fallback
    path = Path(raw)
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def _live12_backend_dependencies(package_root: Path) -> list[tuple[str, Path]]:
    """Return immutable ISO/smodes inputs; generated WSL ``i*.iso`` is excluded."""
    config_path = package_root / "config" / "settings.yaml"
    isobyu: dict[str, Any] = {}
    try:
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict) and isinstance(loaded.get("isobyu"), dict):
            isobyu = loaded["isobyu"]
    except (OSError, UnicodeError, yaml.YAMLError):
        # The configuration file is itself content-signed below.  Deterministic
        # fallbacks preserve an explicit missing-file marker for a broken config.
        pass

    bin_dir = _live12_configured_path(config_path, isobyu.get("bin_dir"), "../isobyu")
    data_dir = _live12_configured_path(config_path, isobyu.get("data_dir"), "../isobyu")
    iso_name = isobyu.get("iso_bin")
    if not isinstance(iso_name, str) or not iso_name.strip():
        iso_name = "iso"
    smodes_name = isobyu.get("smodes_bin")
    if not isinstance(smodes_name, str) or not smodes_name.strip():
        smodes_name = "smodes"

    dependencies = [
        ("iso-binary", bin_dir / iso_name),
        ("smodes-binary", bin_dir / smodes_name),
    ]
    data_names = set(ISO_DATA_BASENAMES)
    try:
        data_names.update(path.name for path in data_dir.glob("data_*.txt"))
    except OSError:
        pass
    dependencies.extend(
        (f"iso-data:{name}", data_dir / name)
        for name in sorted(data_names, key=str.casefold)
    )
    return dependencies


def _live12_checkpoint_policy() -> dict[str, Any]:
    """Document why mutable generated subgroup caches are not batch-signed."""
    return {
        "scientific_inputs": "content_signed",
        "candidate_mode_counts": (
            "reuse_only_after_the_same_identity_and_full_candidate_fingerprint_are_rediscovered_"
            "by_current_live_enumeration"
        ),
        "generated_wsl_i_star_iso_cache": {
            "included_in_global_signature": False,
            "reason": (
                "generate_if_missing may create or update ~/.id/tmp/i*.iso during this audit; "
                "signing that mutable global directory would invalidate the audit's own checkpoint"
            ),
            "safety_guard": (
                "the live candidate set is enumerated on every invocation and a saved count is "
                "used only when identity plus basis/origin/k/OPD/subgroup fields match the current "
                "live object and the identity is also present in the official set"
            ),
        },
    }


def _live12_hash_path(digest: Any, path: Path, label: str) -> None:
    digest.update(label.encode("utf-8"))
    digest.update(b"\0")
    if path.is_file():
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
    else:
        digest.update(b"<missing>")
    digest.update(b"\0")


def _live12_signature_state(
    output_root: Path,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> dict[str, Any]:
    """Fingerprint every immutable scientific input needed for safe resume.

    This deliberately does not walk WSL's mutable ``~/.id/tmp/i*.iso`` cache.
    The accompanying policy explains the identity-level guard that makes
    interrupted Method 1/2 count reuse safe without self-invalidating a run.
    """
    workspace = (workspace or WORKSPACE).resolve()
    package_root = (package_root or workspace / "ISODISTORT").resolve()
    output_root = output_root.resolve()
    runtime = _live12_runtime_versions()
    policy = _live12_checkpoint_policy()
    digest = hashlib.sha256(usedforsecurity=False)
    digest.update(
        json.dumps(
            {"runtime_dependencies": runtime, "checkpoint_reuse_policy": policy},
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    digest.update(b"\0")

    inputs: list[tuple[str, Path]] = [
        ("validator", Path(__file__).resolve()),
        ("config", package_root / "config" / "settings.yaml"),
    ]
    inputs.extend(
        (f"isocore:{path.relative_to(package_root).as_posix()}", path)
        for path in sorted(
            (package_root / "isocore").rglob("*.py"),
            key=lambda item: item.as_posix().casefold(),
        )
    )
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        inputs.append((f"parent:{case_name}", workspace / "experiment_data" / case_name))

    backend_dependencies = _live12_backend_dependencies(package_root)
    inputs.extend(backend_dependencies)
    official_reference_count = 0
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        for method in ("Method1", "Method2"):
            official_root = output_root / case_name / OFFICIAL_DIR / method
            prefix = f"official:{case_name}:{method}"
            if not official_root.is_dir():
                inputs.append((f"{prefix}:<missing-directory>", official_root))
                continue
            references = [
                path
                for path in official_root.rglob("*")
                if path.is_file()
                and (
                    path.name.casefold() == "subgroup.cif"
                    or (
                        "complete modes details" in path.name.casefold()
                        and path.suffix.casefold() in {".html", ".htm"}
                    )
                )
            ]
            if not references:
                inputs.append((f"{prefix}:<no-reference-files>", official_root / ".missing"))
            for path in sorted(
                references,
                key=lambda item: item.relative_to(official_root).as_posix().casefold(),
            ):
                relative = path.relative_to(official_root).as_posix()
                inputs.append((f"{prefix}:{relative}", path))
                official_reference_count += 1

    for label, path in inputs:
        _live12_hash_path(digest, path, label)
    return {
        "signature": digest.hexdigest(),
        "runtime_dependencies": runtime,
        "checkpoint_reuse_policy": policy,
        "signature_inputs": {
            "config": str((package_root / "config" / "settings.yaml").resolve()),
            "python_source_file_count": 1 + sum(label.startswith("isocore:") for label, _ in inputs),
            "official_reference_file_count": official_reference_count,
            "backend_dependencies": [
                {
                    "label": label,
                    "path": str(path.resolve()),
                    "exists": path.is_file(),
                    "content_signed": True,
                }
                for label, path in backend_dependencies
            ],
        },
    }


def _live12_signature(
    output_root: Path,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> str:
    """Return the content/runtime signature used by the Method 1/2 checkpoint."""
    return str(
        _live12_signature_state(
            output_root,
            workspace=workspace,
            package_root=package_root,
        )["signature"]
    )


def _save_live12_checkpoint(path: Path, payload: dict[str, Any]) -> None:
    """Atomically persist audit progress so terminal interruption is harmless."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(path)


def _live12_report_payload(checkpoint: dict[str, Any]) -> dict[str, Any]:
    reports = list(checkpoint.get("reports", {}).values())
    return {
        "schema": LIVE12_CHECKPOINT_SCHEMA,
        "signature": checkpoint["signature"],
        "status": checkpoint["status"],
        "execution_complete": checkpoint["execution_complete"],
        "validation_complete": checkpoint.get("validation_complete", False),
        "passed": checkpoint.get("passed", False),
        "created_at": checkpoint["created_at"],
        "updated_at": checkpoint["updated_at"],
        "completed_at": checkpoint.get("completed_at"),
        "last_invocation_started_at": checkpoint["last_invocation_started_at"],
        "last_invocation_elapsed_seconds": checkpoint["last_invocation_elapsed_seconds"],
        "cumulative_elapsed_seconds": checkpoint["cumulative_elapsed_seconds"],
        "resumed_from_checkpoint": checkpoint["resumed_from_checkpoint"],
        "checkpoint_reuse_policy": checkpoint["signature_state"]["checkpoint_reuse_policy"],
        "runtime_dependencies": checkpoint["signature_state"]["runtime_dependencies"],
        "signature_inputs": checkpoint["signature_state"]["signature_inputs"],
        "progress": checkpoint["progress"],
        "last_error": checkpoint.get("last_error"),
        "reports": reports,
    }


def _update_live12_timing(
    checkpoint: dict[str, Any],
    invocation_started: float,
    previous_elapsed: float,
) -> None:
    elapsed = max(0.0, time.perf_counter() - invocation_started)
    checkpoint["updated_at"] = _live12_now()
    checkpoint["last_invocation_elapsed_seconds"] = round(elapsed, 6)
    checkpoint["cumulative_elapsed_seconds"] = round(previous_elapsed + elapsed, 6)


def _save_live12_state(
    checkpoint_path: Path,
    report_path: Path | None,
    checkpoint: dict[str, Any],
) -> None:
    _save_live12_checkpoint(checkpoint_path, checkpoint)
    if report_path is not None:
        _save_live12_checkpoint(report_path, _live12_report_payload(checkpoint))


def run_live_method12_audit(
    output_root: Path,
    checkpoint_path: Path,
    report_path: Path | None = None,
    *,
    resume: bool = True,
) -> dict[str, Any]:
    """Recompute every saved Method 1/2 candidate and compare mode counts.

    The official folders are read only for candidate identities and expected
    displacive-mode counts.  All local results stay in memory, which separates
    stale saved downloads from defects in the current computation engine.
    Successful candidate counts are atomically checkpointed after every item.
    Reuse requires an exact signature over source, configuration, interpreter,
    scientific package versions, parent/reference files, ISO/smodes binaries,
    and bundled data tables.  Mutable generated WSL ``i*.iso`` files are not
    globally signed; instead, identities are always enumerated afresh and only
    a current identity with the same complete subgroup fingerprint can reuse
    its saved mode count.
    """
    project_root = WORKSPACE / "ISODISTORT"
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from isocore.api import IsoDistort  # noqa: PLC0415

    cases = (
        ("EuAl4 Parent.cif", "LD", ["1/6"]),
        ("NdNiO2 own.cif", "Y", ["1/3"]),
    )
    invocation_started = time.perf_counter()
    invocation_started_at = _live12_now()
    signature_state = _live12_signature_state(output_root)
    signature = str(signature_state["signature"])
    checkpoint: dict[str, Any] = {
        "schema": LIVE12_CHECKPOINT_SCHEMA,
        "signature": signature,
        "signature_state": signature_state,
        "status": "incomplete",
        "execution_complete": False,
        "validation_complete": False,
        "passed": False,
        "created_at": invocation_started_at,
        "updated_at": invocation_started_at,
        "last_invocation_started_at": invocation_started_at,
        "last_invocation_elapsed_seconds": 0.0,
        "cumulative_elapsed_seconds": 0.0,
        "resumed_from_checkpoint": False,
        "counts": {},
        "reports": {},
        "progress": {
            "completed_groups": 0,
            "total_groups": LIVE12_GROUP_COUNT,
            "current_group": None,
            "current_candidate": 0,
            "current_group_candidates": 0,
            "checkpointed_candidate_counts": 0,
        },
    }
    if resume and checkpoint_path.is_file():
        try:
            saved = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            saved = None
        if (
            isinstance(saved, dict)
            and saved.get("schema") == LIVE12_CHECKPOINT_SCHEMA
            and saved.get("signature") == signature
            and isinstance(saved.get("counts"), dict)
            and isinstance(saved.get("reports"), dict)
        ):
            checkpoint = saved
            checkpoint["signature_state"] = signature_state
            checkpoint["resumed_from_checkpoint"] = any(
                isinstance(values, dict) and bool(values)
                for values in checkpoint["counts"].values()
            )
            print(f"Resuming live Method 1/2 audit from {checkpoint_path}")
        else:
            print("Ignoring incompatible live Method 1/2 checkpoint")
    checkpoint.setdefault("created_at", invocation_started_at)
    previous_elapsed_raw = checkpoint.get("cumulative_elapsed_seconds", 0.0)
    previous_elapsed = (
        float(previous_elapsed_raw)
        if isinstance(previous_elapsed_raw, (int, float)) and not isinstance(previous_elapsed_raw, bool)
        else 0.0
    )
    checkpoint.update(
        {
            "signature_state": signature_state,
            "status": "incomplete",
            "execution_complete": False,
            "validation_complete": False,
            "passed": False,
            "last_invocation_started_at": invocation_started_at,
            "last_invocation_elapsed_seconds": 0.0,
            "reports": {},
            "progress": {
                "completed_groups": 0,
                "total_groups": LIVE12_GROUP_COUNT,
                "current_group": None,
                "current_candidate": 0,
                "current_group_candidates": 0,
                "checkpointed_candidate_counts": sum(
                    len(values)
                    for values in checkpoint["counts"].values()
                    if isinstance(values, dict)
                ),
            },
        }
    )
    checkpoint.pop("completed_at", None)
    checkpoint.pop("last_error", None)
    _save_live12_state(checkpoint_path, report_path, checkpoint)

    try:
        for case_name, k_label, k_parameters in cases:
            for method in ("Method1", "Method2"):
                group_started = time.perf_counter()
                group_key = f"{case_name}|{method}"
                official_root = output_root / case_name / OFFICIAL_DIR / method
                official, duplicates = _index_candidates(official_root)
                expected = {
                    identity: count
                    for identity, folder in official.items()
                    if (count := _count_official_displacive_modes(folder)) is not None
                }

                api = IsoDistort()
                with contextlib.redirect_stdout(io.StringIO()):
                    api.load_structure(WORKSPACE / "experiment_data" / case_name)
                    if method == "Method1":
                        rows = api.search_method_1(
                            distortion_types=["strain", "displacive"],
                        )
                        candidates = [row.subgroup for row in rows]
                    else:
                        candidates = api.list_subgroups_at_kpoint(
                            k_label,
                            k_parameters=k_parameters,
                            generate_if_missing=True,
                        )

                live: dict[str, Any] = {}
                duplicate_live: set[str] = set()
                for subgroup in candidates:
                    identity = _live_subgroup_identity(subgroup)
                    if identity in live:
                        duplicate_live.add(identity)
                    live[identity] = subgroup

                missing = sorted(official.keys() - live.keys())
                extra = sorted(live.keys() - official.keys())
                mismatches: list[dict[str, Any]] = []
                matched = 0
                identities = sorted(expected.keys() & live.keys())
                saved_counts = checkpoint["counts"].get(group_key, {})
                counts = {
                    identity: record
                    for identity, record in saved_counts.items()
                    if identity in identities
                    and isinstance(record, dict)
                    and isinstance(record.get("count"), int)
                    and not isinstance(record.get("count"), bool)
                    and record.get("candidate_fingerprint")
                    == _live12_candidate_fingerprint(live[identity])
                } if isinstance(saved_counts, dict) else {}
                checkpoint["counts"][group_key] = counts
                computed_count = 0
                reused_count = 0
                checkpoint["progress"].update(
                    {
                        "current_group": group_key,
                        "current_candidate": 0,
                        "current_group_candidates": len(identities),
                    }
                )
                _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                _save_live12_state(checkpoint_path, report_path, checkpoint)

                for position, identity in enumerate(identities, start=1):
                    subgroup = live[identity]
                    candidate_fingerprint = _live12_candidate_fingerprint(subgroup)
                    saved_record = counts.get(identity)
                    if (
                        isinstance(saved_record, dict)
                        and isinstance(saved_record.get("count"), int)
                        and not isinstance(saved_record.get("count"), bool)
                        and saved_record.get("candidate_fingerprint") == candidate_fingerprint
                    ):
                        actual_count = saved_record["count"]
                        reused_count += 1
                    else:
                        with contextlib.redirect_stdout(io.StringIO()):
                            api.search_method_2(
                                subgroup_idx=subgroup.index,
                                distortion_type=["displacive"],
                                number_of_independent_modulations=0,
                                candidates=candidates,
                            )
                        actual_count = len(api.mode_displacements)
                        counts[identity] = {
                            "count": actual_count,
                            "candidate_fingerprint": candidate_fingerprint,
                        }
                        computed_count += 1
                    checkpoint["progress"].update(
                        {
                            "current_candidate": position,
                            "checkpointed_candidate_counts": sum(
                                len(values)
                                for values in checkpoint["counts"].values()
                                if isinstance(values, dict)
                            ),
                        }
                    )
                    _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                    _save_live12_state(checkpoint_path, report_path, checkpoint)
                    if position % 10 == 0 or position == len(identities):
                        print(
                            f"Live audit progress {case_name} {method}: "
                            f"{position}/{len(identities)}"
                        )
                    expected_count = expected[identity]
                    if actual_count == expected_count:
                        matched += 1
                    else:
                        mismatches.append(
                            {
                                "candidate": identity,
                                "official": expected_count,
                                "live": actual_count,
                            }
                        )

                duplicate_official = [
                    {
                        "candidate": identity,
                        "directories": [str(path) for path in paths],
                    }
                    for identity, paths in duplicates
                ]
                reference_complete = bool(official) and not duplicate_official and len(expected) == len(official)
                validation_passed = bool(
                    reference_complete
                    and not missing
                    and not extra
                    and not duplicate_live
                    and not mismatches
                )
                report = {
                    "case": case_name,
                    "method": method,
                    "status": (
                        "passed"
                        if validation_passed
                        else "failed"
                        if reference_complete
                        else "incomplete_reference"
                    ),
                    "execution_complete": True,
                    "validation_complete": reference_complete,
                    "validation_passed": validation_passed,
                    "official_candidates": len(official),
                    "official_candidates_with_mode_reference": len(expected),
                    "live_candidates": len(live),
                    "mode_files_compared": len(expected.keys() & live.keys()),
                    "mode_count_matches": matched,
                    "checkpoint_counts_reused": reused_count,
                    "mode_counts_computed": computed_count,
                    "elapsed_seconds": round(time.perf_counter() - group_started, 6),
                    "missing": missing,
                    "extra": extra,
                    "duplicate_official": duplicate_official,
                    "duplicate_live": sorted(duplicate_live),
                    "mismatches": mismatches,
                }
                checkpoint["reports"][group_key] = report
                checkpoint["progress"].update(
                    {
                        "completed_groups": len(checkpoint["reports"]),
                        "current_group": None,
                        "current_candidate": 0,
                        "current_group_candidates": 0,
                    }
                )
                _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
                _save_live12_state(checkpoint_path, report_path, checkpoint)

        reports = list(checkpoint["reports"].values())
        validation_complete = (
            len(reports) == LIVE12_GROUP_COUNT
            and all(report["validation_complete"] for report in reports)
        )
        passed = validation_complete and all(report["validation_passed"] for report in reports)
        checkpoint.update(
            {
                "status": (
                    "complete_passed"
                    if passed
                    else "complete_failed"
                    if validation_complete
                    else "complete_validation_incomplete"
                ),
                "execution_complete": True,
                "validation_complete": validation_complete,
                "passed": passed,
                "completed_at": _live12_now(),
            }
        )
        _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
        _save_live12_state(checkpoint_path, report_path, checkpoint)
        return _live12_report_payload(checkpoint)
    except BaseException as exc:
        checkpoint.update(
            {
                "status": "incomplete_error",
                "execution_complete": False,
                "validation_complete": False,
                "passed": False,
                "last_error": f"{type(exc).__name__}: {exc}",
            }
        )
        _update_live12_timing(checkpoint, invocation_started, previous_elapsed)
        _save_live12_state(checkpoint_path, report_path, checkpoint)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--json", action="store_true", help="write the complete audit result to stdout as JSON")
    parser.add_argument(
        "--live-method3",
        action="store_true",
        help="also run the fixed EuAl4/NdNiO2 Method 3 cases through iso/WSL",
    )
    parser.add_argument(
        "--live-method12",
        action="store_true",
        help="recompute every saved Method 1/2 candidate and compare current mode counts",
    )
    parser.add_argument(
        "--live-method12-checkpoint",
        type=Path,
        default=WORKSPACE / "ISODISTORT" / "output" / "validation" / "live_method12_checkpoint.json",
        help="atomic resume/checkpoint JSON for --live-method12",
    )
    parser.add_argument(
        "--live-method12-report",
        type=Path,
        default=WORKSPACE / "ISODISTORT" / "output" / "validation" / "live_method12_report.json",
        help="atomic complete/incomplete report JSON for --live-method12",
    )
    parser.add_argument(
        "--live-method12-restart",
        action="store_true",
        help="ignore an otherwise compatible Method 1/2 checkpoint and recompute all mode counts",
    )
    args = parser.parse_args(argv)
    try:
        summary = audit_output_root(args.output_root)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    try:
        live_method3 = run_live_method3_smoke() if args.live_method3 else []
        live_method12 = (
            run_live_method12_audit(
                args.output_root.expanduser().resolve(),
                args.live_method12_checkpoint.expanduser().resolve(),
                args.live_method12_report.expanduser().resolve(),
                resume=not args.live_method12_restart,
            )
            if args.live_method12 else {}
        )
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: live validation failed: {exc}", file=sys.stderr)
        return 1
    if args.json:
        payload = _payload(summary)
        payload["live_method3"] = live_method3
        payload["live_method12"] = live_method12
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        _print_summary(summary)
        for report in live_method3:
            print(
                f"Live Method 3 PASS {report['case']}: "
                f"SG={report['target_space_group']} candidates={report['candidates']} "
                f"ZIP_CIFs={report['zip_cifs']}"
            )
        if live_method12:
            print(
                "Live Method 1/2 audit: "
                f"status={live_method12['status']} "
                f"execution_complete={live_method12['execution_complete']} "
                f"validation_complete={live_method12['validation_complete']} "
                f"elapsed={live_method12['last_invocation_elapsed_seconds']:.3f}s"
            )
            print(
                "  Checkpoint reuse: "
                f"resumed={live_method12['resumed_from_checkpoint']}; "
                "generated WSL i*.iso cache is excluded; live identities and candidate "
                "fingerprints are re-evaluated"
            )
        for report in live_method12.get("reports", []):
            print(
                f"Live {report['method']} {report['case']}: "
                f"status={report['status']} "
                f"candidates={report['live_candidates']}/"
                f"{report['official_candidates']} mode-counts="
                f"{report['mode_count_matches']}/"
                f"{report['mode_files_compared']} missing={len(report['missing'])} "
                f"extra={len(report['extra'])} mismatches={len(report['mismatches'])}"
            )
            for mismatch in report["mismatches"][:10]:
                print(
                    "  MODE MISMATCH "
                    f"{mismatch['candidate']}: official={mismatch['official']} "
                    f"live={mismatch['live']}"
                )
    live12_failed = bool(live_method12) and not bool(live_method12.get("passed"))
    return 0 if summary.failures == 0 and not live12_failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
