"""Single-worker export ownership, progress and cancellation.

One lock commits the task's terminal state and releases ownership together.
Callbacks and cancellation always carry task identity, including late arrivals.
"""
from __future__ import annotations

import copy
import threading
from dataclasses import dataclass, field


class ExportBusyError(RuntimeError):
    """A previous export still owns the single export worker."""


def idle_progress() -> dict:
    return {
        "job_id": None, "source_revision": None,
        "active": False, "phase": "idle", "total": 0, "completed": 0,
        "successful": 0, "ineligible": 0, "failed": 0,
        "estimated_remaining_seconds": None, "cancellation_requested": False,
    }


@dataclass
class ExportJob:
    job_id: str
    progress: dict
    cancellation: threading.Event = field(default_factory=threading.Event)


class ExportJobs:
    """Own at most one running export and retain its detached latest status."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._active: ExportJob | None = None
        self._last: dict = idle_progress()
        self._used_ids: set[str] = set()

    def start(self, job_id: str, *, revision: int, method: int, total: int) -> ExportJob:
        with self._lock:
            if self._active is not None:
                raise ExportBusyError("Another subgroup export is already running")
            if job_id in self._used_ids:
                raise ValueError("job_id has already been used; start a new export")
            progress = idle_progress()
            progress.update(job_id=job_id, source_revision=revision, active=True,
                            phase="starting", method=method, total=total)
            job = ExportJob(job_id, progress)
            self._used_ids.add(job_id)
            self._active = job
            return job

    def snapshot(self) -> dict:
        with self._lock:
            return copy.deepcopy(self._active.progress if self._active else self._last)

    def update(self, job: ExportJob, values: dict) -> bool:
        with self._lock:
            if self._active is not job:
                return False
            # Scientific workers supply counts/phase, not ownership metadata.
            protected = {"job_id", "source_revision", "active", "cancellation_requested"}
            job.progress.update({key: value for key, value in values.items()
                                 if key not in protected})
            return True

    def cancel(self, job_id: str | None) -> bool:
        with self._lock:
            job = self._active
            if job is None or job_id != job.job_id:
                return False
            job.cancellation.set()
            job.progress["cancellation_requested"] = True
            if job.progress["phase"] != "cancelled":
                job.progress["phase"] = "cancellation_requested"
            return True

    def finish(self, job: ExportJob, *, error: str | None = None) -> dict:
        with self._lock:
            if self._active is not job:
                raise RuntimeError("export job no longer owns the worker")
            # Only the candidate-boundary callback proves cancellation happened.
            phase = "failed" if error is not None else (
                "cancelled" if job.progress["phase"] == "cancelled" else "complete"
            )
            job.progress.update(active=False, phase=phase, estimated_remaining_seconds=0.0)
            if error is not None:
                job.progress["error"] = error
            self._last = copy.deepcopy(job.progress)
            # The response retains this task's terminal state even if another
            # request acquires the worker before the first download is sent.
            result = copy.deepcopy(self._last)
            self._active = None
            return result

    def reset(self) -> None:
        with self._lock:
            if self._active is not None:
                raise ExportBusyError("Cannot reset progress while an export is running")
            self._last = idle_progress()
