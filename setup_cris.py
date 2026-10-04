"""Install, download, and diagnose the CRIS runtime.

This is the single deployment entry point for all three CRIS subprojects.
It intentionally uses only the Python standard library so that it can create
the shared ``.venv`` before third-party packages are available.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import venv
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERSION_FILE = ROOT / "VERSION"
VENV_DIR = ROOT / ".venv"
MIN_PYTHON = (3, 10)


@dataclass(frozen=True)
class Project:
    key: str
    directory: Path
    requirement_file: Path
    dev_requirement_file: Path


PROJECTS = {
    "isodistort": Project(
        "isodistort",
        ROOT / "ISODISTORT",
        ROOT / "ISODISTORT" / "requirements.txt",
        ROOT / "ISODISTORT" / "requirements-dev.txt",
    ),
    "validate": Project(
        "validate",
        ROOT / "ISODISTORT_VALIDATE",
        ROOT / "ISODISTORT_VALIDATE" / "requirements.txt",
        ROOT / "ISODISTORT_VALIDATE" / "requirements-dev.txt",
    ),
    "isoviz": Project(
        "isoviz",
        ROOT / "ISOVIZ_INPUT",
        ROOT / "ISOVIZ_INPUT" / "requirements.txt",
        ROOT / "ISOVIZ_INPUT" / "requirements-dev.txt",
    ),
}


class SetupError(RuntimeError):
    """An actionable deployment error."""


@dataclass(frozen=True)
class Check:
    component: str
    name: str
    status: str
    detail: str
    fix: str = ""


class DoctorReport:
    def __init__(self, selected: Sequence[str]) -> None:
        self.selected = list(selected)
        self.checks: list[Check] = []

    def add(
        self,
        component: str,
        name: str,
        status: str,
        detail: str,
        fix: str = "",
    ) -> None:
        self.checks.append(Check(component, name, status, detail, fix))

    def summary(self) -> dict[str, int]:
        return {
            status: sum(check.status == status for check in self.checks)
            for status in ("pass", "warn", "fail")
        }

    def as_payload(self) -> dict[str, object]:
        return {
            "cris_version": cris_version(),
            "repository_root": str(ROOT),
            "selected_projects": self.selected,
            "checks": [asdict(check) for check in self.checks],
            "summary": self.summary(),
        }

    def print_text(self) -> None:
        print(f"CRIS {cris_version()} environment report")
        print(f"Repository: {ROOT}")
        for check in self.checks:
            print(f"[{check.status.upper():4s}] {check.component}/{check.name}: {check.detail}")
            if check.fix:
                print(f"       Fix: {check.fix}")
        summary = self.summary()
        print(
            "Summary: "
            f"pass={summary['pass']}, warn={summary['warn']}, fail={summary['fail']}"
        )


def cris_version() -> str:
    try:
        return VERSION_FILE.read_text(encoding="utf-8").strip()
    except OSError:
        return "unknown"


def _is_windows() -> bool:
    return sys.platform.startswith("win")


def _venv_python() -> Path:
    if _is_windows():
        return VENV_DIR / "Scripts" / "python.exe"
    return VENV_DIR / "bin" / "python"


def _run(
    command: Sequence[str | Path],
    *,
    cwd: Path | None = None,
    capture: bool = False,
    timeout: int | None = None,
    check: bool = False,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    child_env = os.environ.copy()
    if env:
        child_env.update(env)
    if _is_windows():
        # Requirements files are UTF-8, while a normal Windows Python can
        # otherwise make pip decode them with the active ANSI code page.
        child_env["PYTHONUTF8"] = "1"
    return subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd) if cwd else None,
        env=child_env,
        check=check,
        text=True,
        errors="replace",
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        timeout=timeout,
    )


def _run_wsl(command: Sequence[str], *, timeout: int = 30) -> tuple[int, str]:
    """Run wsl.exe while keeping Linux stdout and Windows diagnostics separate."""
    completed = subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        timeout=timeout,
    )

    output_parts = [
        _decode_process_output(completed.stdout or b""),
        _decode_process_output(completed.stderr or b""),
    ]
    output = "\n".join(part for part in output_parts if part)
    return completed.returncode, output.strip()


def _decode_process_output(raw: bytes) -> str:
    """Decode UTF-8 program output or UTF-16LE Windows diagnostics."""
    if not raw:
        return ""
    if b"\x00" in raw:
        return raw.decode("utf-16-le", errors="replace").strip("\x00\r\n ")
    return raw.decode("utf-8", errors="replace").strip()


def _selected_projects(name: str) -> list[Project]:
    if name == "all":
        return list(PROJECTS.values())
    return [PROJECTS[name]]


def _requirement_files(projects: Sequence[Project], *, dev: bool) -> list[Path]:
    files: list[Path] = []
    for project in projects:
        candidates = [project.requirement_file]
        if dev:
            candidates.append(project.dev_requirement_file)
        for path in candidates:
            if path not in files:
                files.append(path)
    return files


def _require_python_version() -> None:
    if sys.version_info < MIN_PYTHON:
        required = ".".join(map(str, MIN_PYTHON))
        actual = f"{sys.version_info.major}.{sys.version_info.minor}"
        raise SetupError(f"Python >= {required} is required; current interpreter is {actual}.")


def _verify_repository() -> None:
    required = [VERSION_FILE, *(project.directory for project in PROJECTS.values())]
    missing = [path for path in required if not path.exists()]
    if missing:
        joined = ", ".join(str(path) for path in missing)
        raise SetupError(f"setup_cris.py must remain in the CRIS repository root; missing: {joined}")


def _ensure_venv(*, recreate: bool) -> Path:
    python = _venv_python()
    if recreate and VENV_DIR.exists():
        resolved = VENV_DIR.resolve()
        if resolved.parent != ROOT.resolve() or resolved.name != ".venv":
            raise SetupError(f"Refusing to remove unexpected virtualenv path: {resolved}")
        print(f"[venv] Removing requested environment: {resolved}")
        shutil.rmtree(resolved)

    if not VENV_DIR.exists():
        print(f"[venv] Creating shared environment: {VENV_DIR}")
        venv.EnvBuilder(with_pip=True).create(VENV_DIR)
    else:
        print(f"[venv] Reusing shared environment: {VENV_DIR}")

    if not python.is_file():
        raise SetupError(f"Virtualenv interpreter is missing: {python}")
    probe = _run(
        [python, "-c", "import sys; print('.'.join(map(str, sys.version_info[:3])))"],
        capture=True,
        timeout=30,
    )
    if probe.returncode != 0:
        raise SetupError(f"Virtualenv interpreter cannot run: {python}")
    version = probe.stdout.strip()
    major_minor = tuple(int(part) for part in version.split(".")[:2])
    if major_minor < MIN_PYTHON:
        raise SetupError(
            f"Existing .venv uses Python {version}; Python >= 3.10 is required. "
            "Rerun with --recreate using a supported bootstrap interpreter."
        )
    print(f"[venv] Python {version} ({python})")
    return python


def _ensure_runtime_directories(projects: Sequence[Project]) -> None:
    selected = {project.key for project in projects}
    paths: list[Path] = []
    if "isodistort" in selected:
        paths.extend([ROOT / "ISODISTORT" / "output", ROOT / "ISODISTORT" / "output" / "tmp"])
    if "validate" in selected:
        compare = ROOT / "ISODISTORT_VALIDATE" / "compare"
        paths.extend([compare, compare / "item", compare / "true"])
    if "isoviz" in selected:
        inputs = ROOT / "ISOVIZ_INPUT" / "input_content"
        paths.extend([inputs, inputs / "data.csv", inputs / "subgroup.isoviz"])
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)
        print(f"[path] {path}")


def _pip_requirement_args(files: Sequence[Path]) -> list[str]:
    args: list[str] = []
    for path in files:
        if not path.is_file():
            raise SetupError(f"Requirement file is missing: {path}")
        args.extend(["-r", str(path)])
    return args


def _install(args: argparse.Namespace) -> int:
    _require_python_version()
    _verify_repository()
    projects = _selected_projects(args.project)
    requirement_files = _requirement_files(projects, dev=args.dev)
    python = _ensure_venv(recreate=args.recreate)
    _ensure_runtime_directories(projects)

    source_args: list[str] = []
    if args.wheelhouse:
        wheelhouse = Path(args.wheelhouse).expanduser().resolve()
        if not wheelhouse.is_dir():
            raise SetupError(f"Offline wheelhouse does not exist: {wheelhouse}")
        source_args = ["--no-index", "--find-links", str(wheelhouse)]

    if args.upgrade_tools:
        print("[pip] Updating pip/setuptools/wheel")
        result = _run(
            [
                python,
                "-m",
                "pip",
                "install",
                "--disable-pip-version-check",
                *source_args,
                "--upgrade",
                "pip",
                "setuptools",
                "wheel",
            ]
        )
        if result.returncode != 0:
            raise SetupError("Unable to update pip build tools.")

    print("[pip] Resolving declared dependency versions")
    result = _run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            *source_args,
            *_pip_requirement_args(requirement_files),
        ]
    )
    if result.returncode != 0:
        hint = " Check the network connection and package index."
        if args.wheelhouse:
            hint = " The wheelhouse does not contain every required compatible wheel."
        raise SetupError(f"Dependency installation failed.{hint}")

    if args.no_doctor:
        print("[done] Dependencies installed; diagnostic was skipped by request.")
        return 0

    command = _post_install_doctor_command(python, project=args.project, dev=args.dev)
    print("[doctor] Verifying the installed environment")
    return _run(command, cwd=ROOT).returncode


def _post_install_doctor_command(
    python: Path,
    *,
    project: str,
    dev: bool,
) -> list[str | Path]:
    """Build the doctor command used after dependency installation.

    The physical Windows venv remains the pip installation target.  On this
    OneDrive checkout its Python image cannot create WSL processes, so the
    post-install diagnostic must use the same external-base launcher as normal
    ISODISTORT runs.  Other platforms continue to execute the venv directly.
    """
    arguments: list[str | Path] = [str(Path(__file__).resolve()), "doctor", "--project", project]
    if dev:
        arguments.append("--dev")
    runner = ROOT / "run_cris.ps1"
    powershell = shutil.which("powershell.exe") if _is_windows() else None
    if powershell and runner.is_file():
        return [
            powershell,
            "-NoProfile",
            "-ExecutionPolicy",
            "Bypass",
            "-File",
            runner,
            *arguments,
        ]
    return [python, *arguments]


def _download(args: argparse.Namespace) -> int:
    _require_python_version()
    _verify_repository()
    projects = _selected_projects(args.project)
    files = _requirement_files(projects, dev=args.dev)
    wheelhouse = Path(args.wheelhouse).expanduser().resolve()
    wheelhouse.mkdir(parents=True, exist_ok=True)
    binary_args = [] if args.allow_source else ["--only-binary=:all:"]
    print(f"[download] Target: {wheelhouse}")
    result = _run(
        [
            sys.executable,
            "-m",
            "pip",
            "download",
            "--disable-pip-version-check",
            "--dest",
            wheelhouse,
            *binary_args,
            "pip",
            "setuptools",
            "wheel",
            *_pip_requirement_args(files),
        ]
    )
    if result.returncode != 0:
        raise SetupError(
            "Dependency download failed. The wheelhouse is platform/Python-specific; "
            "rerun on the target platform, or use --allow-source if source archives are acceptable."
        )
    print("[done] Offline dependency bundle is ready.")
    print(f"Install with: py -3.10 setup_cris.py install --wheelhouse \"{wheelhouse}\"")
    return 0


_REQUIREMENTS_PROBE = r"""
import json
import sys
from importlib import metadata
try:
    from packaging.requirements import Requirement
except ImportError:
    from pip._vendor.packaging.requirements import Requirement

rows = []
for filename in sys.argv[1:]:
    for raw in open(filename, encoding="utf-8"):
        text = raw.split("#", 1)[0].strip()
        if not text:
            continue
        req = Requirement(text)
        if req.marker is not None and not req.marker.evaluate():
            continue
        try:
            installed = metadata.version(req.name)
        except metadata.PackageNotFoundError:
            installed = None
        compatible = installed is not None and (
            not req.specifier or req.specifier.contains(installed, prereleases=True)
        )
        rows.append({
            "requirement": text,
            "name": req.name,
            "installed": installed,
            "compatible": compatible,
        })
print(json.dumps(rows))
"""


def _python_version(python: Path) -> tuple[bool, str]:
    probe = _run(
        [python, "-c", "import sys; print('.'.join(map(str, sys.version_info[:3])))"],
        capture=True,
        timeout=30,
    )
    if probe.returncode != 0:
        return False, (probe.stdout or "interpreter failed").strip()
    return True, probe.stdout.strip()


def _probe_json(
    python: Path,
    code: str,
    *,
    cwd: Path,
    arguments: Sequence[Path] = (),
    timeout: int = 60,
) -> tuple[dict[str, object] | list[object] | None, str]:
    probe = _run(
        [python, "-c", code, *(str(path) for path in arguments)],
        cwd=cwd,
        capture=True,
        timeout=timeout,
    )
    output = (probe.stdout or "").strip()
    if probe.returncode != 0:
        return None, output or f"process exited with {probe.returncode}"
    try:
        return json.loads(output), ""
    except json.JSONDecodeError:
        return None, f"invalid diagnostic output: {output[:500]}"


def _check_repository(report: DoctorReport, projects: Sequence[Project]) -> None:
    expected_version = cris_version()
    if re.fullmatch(r"\d+\.\d+\.\d+", expected_version):
        report.add("repository", "version", "pass", expected_version)
    else:
        report.add("repository", "version", "fail", f"invalid VERSION value: {expected_version!r}")

    for project in projects:
        pyproject = project.directory / "pyproject.toml"
        if not pyproject.is_file():
            report.add(project.key, "layout", "fail", f"missing {pyproject}")
            continue
        match = re.search(
            r"(?m)^version\s*=\s*[\"']([^\"']+)[\"']",
            pyproject.read_text(encoding="utf-8"),
        )
        actual = match.group(1) if match else "missing"
        if actual == expected_version:
            report.add(project.key, "version", "pass", actual)
        else:
            report.add(
                project.key,
                "version",
                "fail",
                f"pyproject={actual}, repository={expected_version}",
                "Keep VERSION, pyproject.toml, and package __version__ synchronized.",
            )


def _check_venv(
    report: DoctorReport,
    projects: Sequence[Project],
    *,
    dev: bool,
) -> Path | None:
    python = _venv_python()
    if not python.is_file():
        report.add(
            "python",
            "shared-venv",
            "fail",
            f"missing {python}",
            "Run: py -3.10 setup_cris.py install",
        )
        return None
    usable, version = _python_version(python)
    if not usable:
        report.add(
            "python",
            "shared-venv",
            "fail",
            version,
            "Run: py -3.10 setup_cris.py install --recreate",
        )
        return None
    major_minor = tuple(int(part) for part in version.split(".")[:2])
    if major_minor < MIN_PYTHON:
        report.add(
            "python",
            "version",
            "fail",
            version,
            "Recreate .venv with Python >= 3.10.",
        )
    else:
        report.add("python", "shared-venv", "pass", f"Python {version} at {python}")

    requirement_files = _requirement_files(projects, dev=dev)
    payload, error = _probe_json(
        python,
        _REQUIREMENTS_PROBE,
        cwd=ROOT,
        arguments=requirement_files,
    )
    if not isinstance(payload, list):
        report.add("python", "requirements", "fail", error, "Rerun setup_cris.py install.")
    else:
        incompatible = [row for row in payload if not row.get("compatible")]
        if incompatible:
            detail = "; ".join(
                f"{row['requirement']} (installed={row['installed']})"
                for row in incompatible
            )
            report.add(
                "python",
                "requirements",
                "fail",
                detail,
                "Rerun setup_cris.py install for the selected project(s).",
            )
        else:
            report.add("python", "requirements", "pass", f"{len(payload)} declarations satisfied")

    pip_check = _run([python, "-m", "pip", "check"], capture=True, timeout=60)
    pip_detail = (pip_check.stdout or "").strip()
    if pip_check.returncode == 0:
        report.add("python", "pip-check", "pass", pip_detail or "no broken requirements")
    else:
        report.add(
            "python",
            "pip-check",
            "fail",
            pip_detail,
            "Rerun setup_cris.py install; use --recreate only if the environment remains inconsistent.",
        )
    return python


def _import_probe(report: DoctorReport, python: Path, project: Project, code: str) -> None:
    probe = _run(
        [python, "-c", code],
        cwd=project.directory,
        capture=True,
        timeout=60,
    )
    output = (probe.stdout or "").strip()
    reported_version = output.splitlines()[-1] if output else ""
    if probe.returncode == 0 and reported_version == cris_version():
        report.add(project.key, "imports", "pass", f"package version {reported_version}")
    elif probe.returncode == 0:
        report.add(
            project.key,
            "imports",
            "fail",
            f"package version {reported_version or 'missing'} != repository {cris_version()}",
            "Synchronize VERSION, pyproject.toml, and package __version__.",
        )
    else:
        report.add(
            project.key,
            "imports",
            "fail",
            output or f"process exited with {probe.returncode}",
            f"Run: py -3.10 setup_cris.py install --project {project.key}",
        )


def _check_wsl(report: DoctorReport) -> bool:
    if not _is_windows():
        report.add("isodistort", "linux-runtime", "pass", "native Linux does not require WSL")
        return True
    executable = shutil.which("wsl.exe")
    if not executable:
        report.add(
            "isodistort",
            "wsl",
            "fail",
            "wsl.exe is not in PATH",
            "Install WSL and a default Linux distribution, then run wsl -e sh -c 'echo ok'.",
        )
        return False
    try:
        returncode, output = _run_wsl(
            [executable, "-e", "sh", "-c", "printf CRIS_WSL_OK"]
        )
    except subprocess.TimeoutExpired:
        report.add(
            "isodistort",
            "wsl",
            "fail",
            "default WSL distribution did not respond within 30 seconds",
            "Start the default WSL distribution once and retry.",
        )
        return False
    if returncode == 0 and "CRIS_WSL_OK" in output:
        report.add("isodistort", "wsl", "pass", f"shell available via {executable}")
        return True
    report.add(
        "isodistort",
        "wsl",
        "fail",
        output or "WSL shell check failed",
        "Set a working default Linux distribution and retry.",
    )
    return False


def _check_isodistort(report: DoctorReport, python: Path) -> None:
    project = PROJECTS["isodistort"]
    _import_probe(
        report,
        python,
        project,
        "from isocore import __version__; from isocore.api import IsoDistort; "
        "from web.server import main; import main_terminal; print(__version__)",
    )
    payload, error = _probe_json(
        python,
        "import json, os; from isocore.utils import get_config; c=get_config(); "
        "print(json.dumps({'iso': str(c.iso_bin), 'findsym': str(c.findsym_bin), "
        "'smodes': str(c.resolve_path(c._cfg['isobyu']['bin_dir']) / c._cfg['isobyu']['smodes_bin']), "
        "'data': os.environ.get('ISODATA', ''), "
        "'output': str(c.resolve_path(c._cfg['runtime']['output_dir'])), "
        "'temp': str(c.resolve_path(c._cfg['runtime']['temp_dir']))}))",
        cwd=project.directory,
    )
    if not isinstance(payload, dict):
        report.add("isodistort", "configuration", "fail", error)
        return
    report.add("isodistort", "configuration", "pass", str(project.directory / "config" / "settings.yaml"))
    wsl_ok = _check_wsl(report)
    missing_bins: list[str] = []
    for key in ("iso", "findsym", "smodes"):
        path = Path(str(payload[key]))
        if not path.is_file():
            missing_bins.append(str(path))
    if missing_bins:
        report.add(
            "isodistort",
            "isotropy-binaries",
            "fail",
            "missing " + ", ".join(missing_bins),
            "Download the Linux ISOTROPY Suite manually and place iso/findsym/smodes in ISODISTORT/isobyu/.",
        )
    else:
        report.add(
            "isodistort",
            "isotropy-binaries",
            "pass",
            "iso, findsym and smodes found",
        )
        if wsl_ok and _is_windows():
            for key in ("iso", "findsym", "smodes"):
                executable = shutil.which("wsl.exe") or "wsl.exe"
                probe = _run(
                    [
                        executable,
                        "-e",
                        "sh",
                        "-c",
                        'path=$(wslpath -a "$1") && test -x "$path"',
                        "sh",
                        str(payload[key]),
                    ],
                    capture=True,
                    timeout=30,
                )
                if probe.returncode != 0:
                    report.add(
                        "isodistort",
                        f"{key}-executable",
                        "fail",
                        f"WSL cannot execute {payload[key]}",
                        f"Inside WSL run chmod +x on the {key} binary and retry.",
                    )
                else:
                    report.add("isodistort", f"{key}-executable", "pass", str(payload[key]))
    data_dir = Path(str(payload["data"]))
    data_files = list(data_dir.glob("data_*.txt")) if data_dir.is_dir() else []
    if data_files:
        report.add("isodistort", "isodata", "pass", f"{len(data_files)} data_*.txt files in {data_dir}")
    else:
        report.add(
            "isodistort",
            "isodata",
            "fail",
            f"no data_*.txt files in {data_dir}",
            "Copy the complete ISOTROPY data_*.txt set into the configured data directory.",
        )
    for name in ("output", "temp"):
        path = Path(str(payload[name]))
        status = "pass" if path.is_dir() else "warn"
        fix = "Run setup_cris.py install --project isodistort." if status == "warn" else ""
        report.add("isodistort", f"{name}-directory", status, str(path), fix)


def _check_validate(report: DoctorReport, python: Path) -> None:
    project = PROJECTS["validate"]
    _import_probe(
        report,
        python,
        project,
        "from isodistort_validate import __version__, compare_cif, run_batch; print(__version__)",
    )
    payload, error = _probe_json(
        python,
        "import json; from isodistort_validate.config_loader import get_config; c=get_config(); "
        "print(json.dumps({'root': str(c.compare_root), 'item': str(c.item_dir), "
        "'true': str(c.true_dir)}))",
        cwd=project.directory,
    )
    if not isinstance(payload, dict):
        report.add("validate", "configuration", "fail", error)
        return
    for name in ("root", "item", "true"):
        path = Path(str(payload[name]))
        status = "pass" if path.is_dir() else "warn"
        fix = "Run setup_cris.py install --project validate." if status == "warn" else ""
        report.add("validate", f"compare-{name}", status, str(path), fix)
    item_path = Path(str(payload["item"]))
    true_path = Path(str(payload["true"]))
    item_files = list(item_path.rglob("*.cif")) if item_path.is_dir() else []
    true_files = list(true_path.rglob("*.cif")) if true_path.is_dir() else []
    if item_files and true_files:
        report.add("validate", "input-data", "pass", f"item={len(item_files)}, true={len(true_files)} CIF files")
    else:
        report.add(
            "validate",
            "input-data",
            "warn",
            f"item={len(item_files)}, true={len(true_files)} CIF files; no comparison was run",
            "Place matching relative CIF paths in compare/item and compare/true before use.",
        )


def _check_isoviz(report: DoctorReport, python: Path) -> None:
    project = PROJECTS["isoviz"]
    _import_probe(
        report,
        python,
        project,
        "from isoviz_input import __version__, apply_amplitudes, read_amplitude_csv; print(__version__)",
    )
    payload, error = _probe_json(
        python,
        "import json; from isoviz_input.launcher import find_isoviz_launcher; "
        "from isoviz_input.paths import best_model_root; launcher=find_isoviz_launcher(); "
        "print(json.dumps({'launcher': str(launcher) if launcher else '', "
        "'best_model_root': str(best_model_root())}))",
        cwd=project.directory,
    )
    if not isinstance(payload, dict):
        report.add("isoviz", "configuration", "fail", error)
        return
    launcher_text = str(payload["launcher"])
    if launcher_text:
        launcher = Path(launcher_text)
        report.add("isoviz", "launcher", "pass", str(launcher))
    else:
        launcher = None
        report.add(
            "isoviz",
            "launcher",
            "fail",
            "no configured .lnk, .exe, or .jar was found",
            "Place ISOViz.lnk in the CRIS root, or set ISOVIZ/ISOVIZ_JAR.",
        )

    java = shutil.which("java") or shutil.which("javaw")
    java_required = launcher is None or launcher.suffix.lower() == ".jar"
    if java:
        probe = _run([java, "-version"], capture=True, timeout=30)
        detail = (probe.stdout or java).strip().splitlines()[0]
        status = "pass" if probe.returncode == 0 else "fail"
        report.add("isoviz", "java", status, detail)
    elif java_required:
        report.add(
            "isoviz",
            "java",
            "fail",
            "java/javaw is not in PATH",
            "Install a JRE/JDK or configure an IsoVIZ executable/shortcut that supplies its own runtime.",
        )
    else:
        report.add("isoviz", "java", "warn", "java is not in PATH; configured non-JAR launcher may still work")

    best_model = Path(str(payload["best_model_root"]))
    csv_files = list(best_model.rglob("*.csv")) if best_model.is_dir() else []
    if csv_files:
        report.add("isoviz", "amplitude-input", "pass", f"{len(csv_files)} CSV files under {best_model}")
    else:
        report.add(
            "isoviz",
            "amplitude-input",
            "warn",
            f"no CSV files under {best_model}; GUI recognition was not tested",
            "Run the GD notebook save step or pass --data to ISOVIZ_INPUT/main.py.",
        )


def _doctor(args: argparse.Namespace) -> int:
    _verify_repository()
    projects = _selected_projects(args.project)
    report = DoctorReport([project.key for project in projects])
    _check_repository(report, projects)
    python = _check_venv(report, projects, dev=args.dev)
    if python is not None:
        selected = {project.key for project in projects}
        if "isodistort" in selected:
            _check_isodistort(report, python)
        if "validate" in selected:
            _check_validate(report, python)
        if "isoviz" in selected:
            _check_isoviz(report, python)
    if args.json:
        print(json.dumps(report.as_payload(), ensure_ascii=False, indent=2))
    else:
        report.print_text()
    return 1 if report.summary()["fail"] else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Install, download, or diagnose the shared CRIS environment.",
    )
    parser.add_argument("--version", action="version", version=f"CRIS {cris_version()}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    install = subparsers.add_parser("install", help="Create/reuse .venv and install declared dependencies.")
    install.add_argument("--project", choices=["all", *PROJECTS], default="all")
    install.add_argument("--dev", action="store_true", help="Also install test/lint dependencies.")
    install.add_argument("--recreate", action="store_true", help="Delete and rebuild only CRIS/.venv.")
    install.add_argument("--wheelhouse", help="Install offline from a directory created by the download command.")
    install.add_argument("--upgrade-tools", action="store_true", help="Upgrade pip/setuptools/wheel first.")
    install.add_argument("--no-doctor", action="store_true", help="Skip the post-install diagnostic.")
    install.set_defaults(handler=_install)

    download = subparsers.add_parser("download", help="Download a reusable Python dependency wheelhouse.")
    download.add_argument("--project", choices=["all", *PROJECTS], default="all")
    download.add_argument("--dev", action="store_true", help="Also download test/lint dependencies.")
    download.add_argument("--wheelhouse", required=True, help="Destination directory for downloaded packages.")
    download.add_argument(
        "--allow-source",
        action="store_true",
        help="Allow source archives when no compatible wheel exists (offline builds may need extra tools).",
    )
    download.set_defaults(handler=_download)

    doctor = subparsers.add_parser("doctor", help="Run read-only environment and deployment checks.")
    doctor.add_argument("--project", choices=["all", *PROJECTS], default="all")
    doctor.add_argument("--dev", action="store_true", help="Also require declared test/lint dependencies.")
    doctor.add_argument("--json", action="store_true", help="Emit a machine-readable report.")
    doctor.set_defaults(handler=_doctor)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (OSError, subprocess.SubprocessError, SetupError) as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
