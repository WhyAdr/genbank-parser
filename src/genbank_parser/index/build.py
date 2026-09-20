"""Crash-safe index construction and update lifecycle."""

from __future__ import annotations

import json
import os
import sqlite3
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .. import __version__
from ..cli_io import (
    InputError,
    OutputError,
    reject_input_output_collision,
    write_text,
)
from ..discovery import (
    DiscoveredInput,
    assign_sample_keys,
    discover_inputs,
    fingerprint_file,
)
from ..io import GenBankInputError, read_genbank
from ..records import RECORD_SCHEMA_VERSION, record_rows
from ..serializers import FEATURE_SCHEMA_VERSION, FeatureRow
from .schema import (
    INDEX_SCHEMA_REVISION,
    INDEX_SCHEMA_VERSION,
    initialize_schema,
    utc_now,
    validate_schema,
)


@dataclass(frozen=True)
class IndexBuildResult:
    report: dict[str, object]
    skipped_sources: tuple[dict[str, str], ...] = ()

    @property
    def exit_code(self) -> int:
        return 3 if self.skipped_sources else 0


def _metadata() -> dict[str, str]:
    import platform
    import sys

    import Bio

    return {
        "schema_version": INDEX_SCHEMA_VERSION,
        "schema_revision": str(INDEX_SCHEMA_REVISION),
        "gbparse_version": __version__,
        "python_version": platform.python_version(),
        "biopython_version": getattr(Bio, "__version__", "unknown"),
        "feature_schema_version": FEATURE_SCHEMA_VERSION,
        "record_schema_version": RECORD_SCHEMA_VERSION,
        "query_semantics_version": "gbparse.query.v1",
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "platform": sys.platform,
    }


def _is_compressed(path: Path) -> int:
    with path.open("rb") as handle:
        return int(handle.read(2) == b"\x1f\x8b")


def _insert_document(
    connection: sqlite3.Connection,
    source: DiscoveredInput,
    sample_key: str,
    fingerprint,
) -> tuple[int, int]:
    document = read_genbank(source.resolved_path)
    if not document.records:
        raise GenBankInputError(f"no GenBank records were parsed from {source.display_path}")
    source_pk = int(
        connection.execute(
            "INSERT INTO sources(display_path, resolved_path_at_index, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                source.display_path,
                str(source.resolved_path),
                sample_key,
                fingerprint.sha256,
                fingerprint.size_bytes,
                fingerprint.mtime_ns,
                _is_compressed(source.resolved_path),
                len(document.records),
                document.total_features,
            ),
        ).lastrowid
    )
    for record_index, record in enumerate(document.records, 1):
        row = record_rows(document)[record_index - 1]
        record_pk = int(
            connection.execute(
                "INSERT INTO records(source_pk, record_index, record_id, record_name, description, length, topology, molecule_type, organism, strain, gc_percent) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    source_pk,
                    record_index,
                    record.id,
                    record.name,
                    record.description,
                    record.length,
                    record.topology,
                    record.molecule_type,
                    row.organism,
                    row.strain,
                    row.gc_percent,
                ),
            ).lastrowid
        )
        for feature in record.features:
            feature_row = FeatureRow.from_feature(feature, source=source.display_path)
            feature_pk = int(
                connection.execute(
                    "INSERT INTO features(record_pk, feature_index, type, locus_tag, gene, product, protein_id, start, end, strand, biological_length, partial, pseudo) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        record_pk,
                        feature.feature_index,
                        feature.type,
                        feature_row.locus_tag,
                        feature_row.gene,
                        feature_row.product,
                        feature_row.protein_id,
                        feature_row.start,
                        feature_row.end,
                        feature_row.strand,
                        feature_row.length,
                        int(feature_row.partial),
                        int(feature_row.pseudo),
                    ),
                ).lastrowid
            )
            for ordinal, (start, end) in enumerate(feature_row.segments):
                connection.execute(
                    "INSERT INTO segments(feature_pk, ordinal, start, end) VALUES (?, ?, ?, ?)",
                    (feature_pk, ordinal, start, end),
                )
            for key, values in sorted(feature_row.qualifiers.items()):
                for ordinal, value in enumerate(values):
                    connection.execute(
                        "INSERT INTO qualifiers(feature_pk, key, ordinal, value) VALUES (?, ?, ?, ?)",
                        (feature_pk, key, ordinal, value),
                    )
            for namespace, values in sorted(feature_row.xref_sources.items()):
                for value, source_field in values:
                    connection.execute(
                        "INSERT OR IGNORE INTO xrefs(feature_pk, namespace, value, source_field) VALUES (?, ?, ?, ?)",
                        (feature_pk, namespace, value, source_field),
                    )
    return len(document.records), document.total_features


def _validate_database(connection: sqlite3.Connection) -> None:
    foreign = connection.execute("PRAGMA foreign_key_check").fetchall()
    if foreign:
        raise ValueError(f"index foreign-key check failed: {foreign[:3]}")
    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise ValueError(f"index integrity check failed: {integrity}")
    connection.commit()


def _new_temp_path(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent, delete=False
    ) as handle:
        return Path(handle.name)


def _discover_for_index(inputs: Iterable[str | Path], destination: Path) -> list[DiscoveredInput]:
    inputs = tuple(inputs)
    reject_input_output_collision(inputs, destination)
    staging_siblings = tuple(destination.parent.glob(f".{destination.name}.*")) if destination.parent.exists() else ()
    discovered = discover_inputs(inputs, exclude=(destination, *staging_siblings))
    if not discovered:
        raise InputError("no GenBank inputs were discovered")
    return discovered


def _write_report(report: dict[str, object], path: str | Path | None, *, inputs: Iterable[str | Path]) -> None:
    if path is not None:
        write_text(
            json.dumps(report, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n",
            path,
            inputs=inputs,
        )


def build_index(
    inputs: Iterable[str | Path],
    destination: str | Path,
    *,
    jobs: int = 1,
    on_error: str = "fail",
    force: bool = False,
    report_path: str | Path | None = None,
) -> IndexBuildResult:
    """Build a validated SQLite index and atomically publish it."""

    if jobs <= 0:
        raise ValueError("jobs must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("on_error must be fail or skip")
    destination_path = Path(destination)
    if destination_path.exists() and not force:
        raise OutputError(f"index already exists; use --force: {destination_path}")
    discovered = _discover_for_index(inputs, destination_path)
    keys = assign_sample_keys(discovered)
    temporary = _new_temp_path(destination_path)
    skipped: list[dict[str, str]] = []
    try:
        connection = sqlite3.connect(temporary)
        initialize_schema(connection, _metadata())
        for source in discovered:
            try:
                fp = fingerprint_file(source.resolved_path)
                with connection:
                    _insert_document(connection, source, keys[source], fp)
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                if on_error == "fail":
                    raise InputError(f"could not index {source.display_path}: {exc}") from exc
                skipped.append({"source": source.display_path, "error": str(exc)})
        _validate_database(connection)
        metadata = {"updated_at": utc_now()}
        connection.executemany(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            sorted(metadata.items()),
        )
        connection.commit()
        connection.close()
        os.replace(temporary, destination_path)
        temporary = None
    finally:
        if connection if "connection" in locals() else False:
            try:
                connection.close()
            except sqlite3.Error:
                pass
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    report = _status_report(destination_path, skipped=skipped, operation="build")
    _write_report(report, report_path, inputs=(destination_path,))
    return IndexBuildResult(report=report, skipped_sources=tuple(skipped))


def _status_report(path: Path, *, skipped: list[dict[str, str]], operation: str, removed: list[str] | None = None) -> dict[str, object]:
    connection = sqlite3.connect(path)
    try:
        metadata = validate_schema(connection)
        counts = {
            "sources": connection.execute("SELECT COUNT(*) FROM sources").fetchone()[0],
            "records": connection.execute("SELECT COUNT(*) FROM records").fetchone()[0],
            "features": connection.execute("SELECT COUNT(*) FROM features").fetchone()[0],
        }
    finally:
        connection.close()
    return {
        "schema_version": "gbparse.index-report.v1",
        "operation": operation,
        "database": str(path),
        "metadata": metadata,
        "counts": counts,
        "skipped_sources": skipped,
        "removed_sources": removed or [],
    }


def update_index(
    database: str | Path,
    inputs: Iterable[str | Path],
    *,
    prune: bool = False,
    jobs: int = 1,
    on_error: str = "fail",
    report_path: str | Path | None = None,
) -> IndexBuildResult:
    """Update an existing index through SQLite backup and atomic replacement."""

    inputs = tuple(inputs)
    if jobs <= 0:
        raise ValueError("jobs must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("on_error must be fail or skip")
    if prune and not any(Path(item).is_dir() for item in inputs):
        raise ValueError("--prune requires at least one directory input root")
    destination = Path(database)
    if not destination.is_file():
        raise InputError(f"index does not exist: {destination}")
    discovered = _discover_for_index(inputs, destination)
    keys = assign_sample_keys(discovered)
    temporary = _new_temp_path(destination)
    skipped: list[dict[str, str]] = []
    removed: list[str] = []
    source_conn = sqlite3.connect(destination)
    temp_conn = sqlite3.connect(temporary)
    try:
        validate_schema(source_conn)
        source_conn.backup(temp_conn)
        temp_conn.close()
        source_conn.close()
        connection = sqlite3.connect(temporary)
        connection.execute("PRAGMA foreign_keys = ON")
        existing = {
            row[0]: row[1]
            for row in connection.execute("SELECT display_path, sha256 FROM sources")
        }
        discovered_paths = {source.display_path for source in discovered}
        if prune:
            for display_path in sorted(set(existing) - discovered_paths):
                connection.execute("DELETE FROM sources WHERE display_path = ?", (display_path,))
                removed.append(display_path)
        for source in discovered:
            try:
                fp = fingerprint_file(source.resolved_path)
                if existing.get(source.display_path) == fp.sha256:
                    continue
                with connection:
                    if source.display_path in existing:
                        connection.execute("DELETE FROM sources WHERE display_path = ?", (source.display_path,))
                    _insert_document(connection, source, keys[source], fp)
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                if on_error == "fail":
                    raise InputError(f"could not update {source.display_path}: {exc}") from exc
                skipped.append({"source": source.display_path, "error": str(exc)})
        _validate_database(connection)
        connection.execute("INSERT OR REPLACE INTO metadata(key, value) VALUES ('updated_at', ?)", (utc_now(),))
        connection.commit()
        connection.close()
        os.replace(temporary, destination)
        temporary = None
    finally:
        for handle in (locals().get("connection"), temp_conn, source_conn):
            try:
                if handle is not None:
                    handle.close()
            except sqlite3.Error:
                pass
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    report = _status_report(destination, skipped=skipped, operation="update", removed=removed)
    _write_report(report, report_path, inputs=(destination,))
    return IndexBuildResult(report=report, skipped_sources=tuple(skipped))


__all__ = ["IndexBuildResult", "build_index", "update_index"]
