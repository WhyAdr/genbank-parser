"""Crash-safe index construction and update lifecycle."""

from __future__ import annotations

import io
import json
import os
import sqlite3
import tempfile
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from .. import __version__
from ..cli_io import (
    InputError,
    OutputError,
    paths_same,
    publish_staged_files,
    reject_input_output_collision,
)
from ..discovery import (
    DiscoveredInput,
    assign_sample_keys,
    discover_inputs,
    snapshot_file,
)
from ..io import GenBankInputError, read_genbank
from ..records import RECORD_SCHEMA_VERSION, record_rows
from ..serializers import FEATURE_SCHEMA_VERSION, FeatureRow
from .schema import (
    INDEX_SCHEMA_REVISION,
    INDEX_SCHEMA_VERSION,
    initialize_schema,
    migrate_schema,
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
        "query_semantics_version": "gbparse.query.v2",
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "platform": sys.platform,
    }


def _prepare_document(source: DiscoveredInput):
    """Read, fingerprint, and parse one immutable source snapshot."""

    raw, fingerprint = snapshot_file(source.resolved_path)
    document = read_genbank(io.BytesIO(raw))
    if not document.records:
        raise GenBankInputError(f"no GenBank records were parsed from {source.display_path}")
    return document, fingerprint, int(raw[:2] == b"\x1f\x8b")


def _prepared_sources(
    discovered: list[DiscoveredInput],
    jobs: int,
) -> dict[str, tuple[object, object] | Exception]:
    """Parse snapshots concurrently while retaining deterministic source order."""

    if jobs == 1 or len(discovered) == 1:
        result: dict[str, tuple[object, object] | Exception] = {}
        for source in discovered:
            try:
                result[source.source_identity] = _prepare_document(source)
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                result[source.source_identity] = exc
        return result
    result = {}
    with ThreadPoolExecutor(max_workers=min(jobs, len(discovered))) as executor:
        futures = {executor.submit(_prepare_document, source): source for source in discovered}
        for future in as_completed(futures):
            source = futures[future]
            try:
                result[source.source_identity] = future.result()
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                result[source.source_identity] = exc
    return result


def _insert_document(
    connection: sqlite3.Connection,
    source: DiscoveredInput,
    sample_key: str,
    fingerprint,
    document,
    compressed: int,
) -> tuple[int, int]:
    source_pk = int(
        connection.execute(
            "INSERT INTO sources(source_identity, display_path, resolved_path_at_index, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                source.source_identity,
                source.display_path,
                str(source.resolved_path),
                sample_key,
                fingerprint.sha256,
                fingerprint.size_bytes,
                fingerprint.mtime_ns,
                compressed,
                len(document.records),
                document.total_features,
            ),
        ).lastrowid
    )
    rows = record_rows(document)
    for record_index, record in enumerate(document.records, 1):
        row = rows[record_index - 1]
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
                for ordinal, (value, source_field) in enumerate(values):
                    connection.execute(
                        "INSERT INTO xrefs(feature_pk, namespace, ordinal, value, source_field) VALUES (?, ?, ?, ?, ?)",
                        (feature_pk, namespace, ordinal, value, source_field),
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


def _discover_for_index(
    inputs: Iterable[str | Path],
    destination: Path,
    *,
    report_path: Path | None = None,
) -> list[DiscoveredInput]:
    inputs = tuple(inputs)
    reject_input_output_collision(inputs, destination)
    if report_path is not None:
        reject_input_output_collision(inputs, report_path)
        if paths_same(destination, report_path):
            raise OutputError("index database and report must use different paths")
    staging_siblings = tuple(destination.parent.glob(f".{destination.name}.*")) if destination.parent.exists() else ()
    report_siblings = (
        tuple(report_path.parent.glob(f".{report_path.name}.*"))
        if report_path is not None and report_path.parent.exists()
        else ()
    )
    discovered = discover_inputs(
        inputs,
        exclude=(destination, report_path, *staging_siblings, *report_siblings),
    )
    if not discovered:
        raise InputError("no GenBank inputs were discovered")
    return discovered


def _preflight_report(
    report_path: Path | None,
    destination: Path,
    inputs: Iterable[str | Path],
    *,
    report_force: bool,
) -> None:
    if report_path is None:
        return
    if paths_same(report_path, destination):
        raise OutputError("index database and report must use different paths")
    reject_input_output_collision(inputs, report_path)
    if report_path.exists() and report_path.is_dir():
        raise OutputError(f"report path is a directory: {report_path}")
    if report_path.exists() and not report_force:
        raise OutputError(f"report already exists; use --report-force: {report_path}")


def _stage_report(report: dict[str, object], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    ) as handle:
        staged = Path(handle.name)
        handle.write(json.dumps(report, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    return staged


def build_index(
    inputs: Iterable[str | Path],
    destination: str | Path,
    *,
    jobs: int = 1,
    on_error: str = "fail",
    force: bool = False,
    report_path: str | Path | None = None,
    report_force: bool = False,
) -> IndexBuildResult:
    """Build a validated SQLite index and atomically publish it."""

    inputs = tuple(inputs)
    if jobs <= 0:
        raise ValueError("jobs must be positive")
    if on_error not in {"fail", "skip"}:
        raise ValueError("on_error must be fail or skip")
    destination_path = Path(destination)
    if destination_path.exists() and not force:
        raise OutputError(f"index already exists; use --force: {destination_path}")
    report_destination = Path(report_path) if report_path is not None else None
    _preflight_report(report_destination, destination_path, inputs, report_force=report_force)
    discovered = _discover_for_index(inputs, destination_path, report_path=report_destination)
    keys = assign_sample_keys(discovered)
    temporary = _new_temp_path(destination_path)
    staged_report: Path | None = None
    skipped: list[dict[str, str]] = []
    successful_sources = 0
    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(temporary)
        initialize_schema(connection, _metadata())
        prepared = _prepared_sources(discovered, jobs)
        for source in discovered:
            try:
                prepared_item = prepared[source.source_identity]
                if isinstance(prepared_item, Exception):
                    raise prepared_item
                document, fp, compressed = prepared_item
                with connection:
                    _insert_document(connection, source, keys[source], fp, document, compressed)
                successful_sources += 1
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                if on_error == "fail":
                    raise InputError(f"could not index {source.display_path}: {exc}") from exc
                skipped.append({"source": source.display_path, "error": str(exc)})
        if successful_sources == 0:
            raise InputError("no usable GenBank sources were indexed; existing index preserved")
        _validate_database(connection)
        metadata = {"updated_at": utc_now()}
        connection.executemany(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            sorted(metadata.items()),
        )
        connection.commit()
        connection.close()
        connection = None
        report = _status_report(
            temporary,
            skipped=skipped,
            operation="build",
            discovered_sources=len(discovered),
            successful_sources=successful_sources,
            database_label=str(destination_path),
        )
        replacements: list[tuple[Path, Path]] = [(temporary, destination_path)]
        if report_destination is not None:
            staged_report = _stage_report(report, report_destination)
            replacements.append((staged_report, report_destination))
        publish_staged_files(replacements, force=True, inputs=inputs)
        temporary = None
        staged_report = None
    finally:
        if connection is not None:
            try:
                connection.close()
            except sqlite3.Error:
                pass
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        if staged_report is not None:
            staged_report.unlink(missing_ok=True)
    return IndexBuildResult(report=report, skipped_sources=tuple(skipped))


def _status_report(
    path: Path,
    *,
    skipped: list[dict[str, str]],
    operation: str,
    removed: list[str] | None = None,
    discovered_sources: int | None = None,
    successful_sources: int | None = None,
    database_label: str | None = None,
) -> dict[str, object]:
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
        "database": database_label or str(path),
        "metadata": metadata,
        "counts": counts,
        "discovered_sources": discovered_sources if discovered_sources is not None else counts["sources"],
        "successful_sources": successful_sources if successful_sources is not None else counts["sources"],
        "skipped_sources": skipped,
        "skipped_count": len(skipped),
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
    report_force: bool = False,
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
    report_destination = Path(report_path) if report_path is not None else None
    _preflight_report(report_destination, destination, inputs, report_force=report_force)
    discovered = _discover_for_index(inputs, destination, report_path=report_destination)
    keys = assign_sample_keys(discovered)
    temporary = _new_temp_path(destination)
    staged_report: Path | None = None
    skipped: list[dict[str, str]] = []
    removed: list[str] = []
    successful_sources = 0
    source_conn: sqlite3.Connection | None = sqlite3.connect(destination)
    temp_conn: sqlite3.Connection | None = sqlite3.connect(temporary)
    connection: sqlite3.Connection | None = None
    try:
        assert source_conn is not None and temp_conn is not None
        source_conn.backup(temp_conn)
        temp_conn.close()
        temp_conn = None
        source_conn.close()
        source_conn = None
        connection = sqlite3.connect(temporary)
        migrate_schema(connection)
        validate_schema(connection)
        existing = {
            row[0]: row[1]
            for row in connection.execute("SELECT source_identity, sha256 FROM sources")
        }
        discovered_paths = {source.source_identity for source in discovered}
        if prune:
            for source_identity in sorted(set(existing) - discovered_paths):
                row = connection.execute(
                    "SELECT display_path FROM sources WHERE source_identity = ?",
                    (source_identity,),
                ).fetchone()
                connection.execute("DELETE FROM sources WHERE source_identity = ?", (source_identity,))
                removed.append(str(row[0]) if row else source_identity)
        prepared = _prepared_sources(discovered, jobs)
        for source in discovered:
            try:
                prepared_item = prepared[source.source_identity]
                if isinstance(prepared_item, Exception):
                    raise prepared_item
                document, fp, compressed = prepared_item
                if existing.get(source.source_identity) == fp.sha256:
                    connection.execute(
                        "UPDATE sources SET display_path = ?, resolved_path_at_index = ?, sample_key = ?, size_bytes = ?, mtime_ns = ?, compressed = ? WHERE source_identity = ?",
                        (
                            source.display_path,
                            str(source.resolved_path),
                            keys[source],
                            fp.size_bytes,
                            fp.mtime_ns,
                            compressed,
                            source.source_identity,
                        ),
                    )
                    successful_sources += 1
                    continue
                with connection:
                    if source.source_identity in existing:
                        connection.execute("DELETE FROM sources WHERE source_identity = ?", (source.source_identity,))
                    _insert_document(connection, source, keys[source], fp, document, compressed)
                successful_sources += 1
            except (GenBankInputError, OSError, ValueError, sqlite3.DatabaseError) as exc:
                if on_error == "fail":
                    raise InputError(f"could not update {source.display_path}: {exc}") from exc
                skipped.append({"source": source.display_path, "error": str(exc)})
        if successful_sources == 0:
            raise InputError("no usable GenBank sources were updated; existing index preserved")
        _validate_database(connection)
        connection.execute("INSERT OR REPLACE INTO metadata(key, value) VALUES ('updated_at', ?)", (utc_now(),))
        connection.commit()
        connection.close()
        connection = None
        report = _status_report(
            temporary,
            skipped=skipped,
            operation="update",
            removed=removed,
            discovered_sources=len(discovered),
            successful_sources=successful_sources,
            database_label=str(destination),
        )
        replacements: list[tuple[Path, Path]] = [(temporary, destination)]
        if report_destination is not None:
            staged_report = _stage_report(report, report_destination)
            replacements.append((staged_report, report_destination))
        publish_staged_files(replacements, force=True, inputs=inputs)
        temporary = None
        staged_report = None
    finally:
        for handle in (connection, temp_conn, source_conn):
            try:
                if handle is not None:
                    handle.close()
            except sqlite3.Error:
                pass
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        if staged_report is not None:
            staged_report.unlink(missing_ok=True)
    return IndexBuildResult(report=report, skipped_sources=tuple(skipped))


def migrate_index(database: str | Path) -> dict[str, str]:
    """Migrate a revision-1 index in place using a transactional schema change."""

    path = Path(database)
    if not path.is_file():
        raise InputError(f"index does not exist: {path}")
    connection = sqlite3.connect(path)
    try:
        metadata = migrate_schema(connection)
        validate_schema(connection)
        return metadata
    finally:
        connection.close()


__all__ = ["IndexBuildResult", "build_index", "migrate_index", "update_index"]
