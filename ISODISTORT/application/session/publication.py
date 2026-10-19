"""Render, stage and atomically publish immutable export batches.

This module accepts prepared specs and a failure factory; it never calculates
scientific modes or changes a session.
"""
from __future__ import annotations
import hashlib
import json
import os
import re
import shutil
import tempfile
import threading
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import asdict
from pathlib import Path
from backend.export import SubgroupExportSpec, render_subgroup_files, subgroup_label
from backend.wrappers import SubgroupInfo
from .contracts import ExportCandidateFailure, ExportPreparationError

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


def _export_failure_report_payloads(
    method: int,
    requested_count: int,
    successful_count: int,
    failures: tuple[ExportCandidateFailure, ...],
) -> tuple[tuple[str, bytes], ...]:
    """Build root-level human and machine status records for every batch."""
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
            "success"
            if not failures
            else "cancelled_partial_success"
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


def render_export_batch(
    method: int,
    specs: list[SubgroupExportSpec],
    formats: list[str],
    failure_factory: Callable[[int, SubgroupInfo, Exception], ExportCandidateFailure],
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
                failure_factory(position, spec.subgroup, exc)
            )
        else:
            rendered.append((spec, payloads))
    if failures:
        raise ExportPreparationError(method, failures)
    return rendered


def publish_export_batch(
    method: int,
    folder_root: Path,
    rendered_batch: list[
        tuple[SubgroupExportSpec, tuple[tuple[str, bytes], ...]]
    ],
    failure_factory: Callable[[int, SubgroupInfo, Exception], ExportCandidateFailure],
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
            [failure_factory(position, spec.subgroup, exc)],
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
                        failure_factory(
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
                            failure_factory(
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
                        failure_factory(
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
                        failure_factory(
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
                    failure_factory(position, spec.subgroup, exc)
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
