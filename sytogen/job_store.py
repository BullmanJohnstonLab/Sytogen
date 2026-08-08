"""
Shared, cross-process job store for the async /api/sytogen job pipeline.

The previous implementation kept job state in a plain in-process dict
(JOBS = {}). That works fine with a single worker process, but silently
breaks under a multi-worker deployment (e.g. `gunicorn -w 4`): a job
submitted to worker A is invisible to workers B, C, and D, so a client
polling /api/status/<job_id> or downloading /api/sytogen/result/<job_id>
will get a 404 "unknown job" any time their request lands on a different
worker than the one that ran it.

This module replaces that dict with a small SQLite database on disk.
SQLite is used (rather than adding a Redis/Postgres dependency) because
all workers in this deployment model run on the same host and already
share a filesystem - WAL mode gives good enough concurrent read/write
behavior for job-status-sized traffic without any extra infrastructure.
If this ever needs to scale across multiple hosts, swap this module out
for a real shared store (Redis is the natural choice) - the functions
below (create_job/update_job/get_job/delete_job/sweep_expired_jobs) are
the whole interface the rest of the app depends on.
"""

import os
import shutil
import sqlite3
import tempfile
import time

DB_PATH = os.environ.get(
    "SYTOGEN_JOB_DB_PATH",
    os.path.join(tempfile.gettempdir(), "sytogen_jobs.db"),
)

# The only fields callers are allowed to write, to keep update_job() from
# ever building a query out of an arbitrary/unexpected key.
_WRITABLE_COLUMNS = (
    "status",
    "error",
    "traceback",
    "tmpdir",
    "result",
    "finished_at",
)


def _connect():
    # timeout=30: if another worker briefly holds the write lock, wait
    # for it instead of raising "database is locked".
    conn = sqlite3.connect(DB_PATH, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL")
    return conn


def init_db():
    with _connect() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                error TEXT,
                traceback TEXT,
                tmpdir TEXT,
                result TEXT,
                created_at REAL,
                finished_at REAL
            )
            """
        )


def create_job(job_id):
    with _connect() as conn:
        conn.execute(
            "INSERT INTO jobs (job_id, status, created_at) VALUES (?, 'queued', ?)",
            (job_id, time.time()),
        )


def update_job(job_id, **fields):
    """Update one or more columns for a job. Pass tmpdir=None to clear a
    field (e.g. after a temp directory has been removed)."""
    if not fields:
        return

    unknown = set(fields) - set(_WRITABLE_COLUMNS)
    if unknown:
        raise ValueError(f"Unknown job field(s): {sorted(unknown)}")

    set_clause = ", ".join(f"{column} = ?" for column in fields)
    values = [*fields.values(), job_id]

    with _connect() as conn:
        conn.execute(f"UPDATE jobs SET {set_clause} WHERE job_id = ?", values)


def get_job(job_id):
    with _connect() as conn:
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
        ).fetchone()
    return dict(row) if row else None


def delete_job(job_id):
    with _connect() as conn:
        conn.execute("DELETE FROM jobs WHERE job_id = ?", (job_id,))


def sweep_expired_jobs(ttl_seconds):
    """Delete (and clean up the temp directory for) every job that
    finished more than ttl_seconds ago. Safe to call from any worker -
    whichever worker's sweep timer fires first does the cleanup for all
    of them, since they all share this same database."""
    cutoff = time.time() - ttl_seconds

    with _connect() as conn:
        conn.row_factory = sqlite3.Row
        expired = conn.execute(
            "SELECT job_id, tmpdir FROM jobs "
            "WHERE finished_at IS NOT NULL AND finished_at < ?",
            (cutoff,),
        ).fetchall()

        if expired:
            conn.executemany(
                "DELETE FROM jobs WHERE job_id = ?",
                [(row["job_id"],) for row in expired],
            )

    for row in expired:
        if row["tmpdir"]:
            shutil.rmtree(row["tmpdir"], ignore_errors=True)


init_db()
