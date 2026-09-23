"""SQLite schema and compatibility checks for ``gbparse.index.v1``."""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from ..discovery import DiscoveredInput, assign_sample_keys

INDEX_SCHEMA_VERSION = "gbparse.index.v1"
INDEX_SCHEMA_REVISION = 3
SQLITE_USER_VERSION = 3

SCHEMA_SQL = """
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE sources (
    source_pk INTEGER PRIMARY KEY,
    source_identity TEXT NOT NULL UNIQUE,
    display_path TEXT NOT NULL,
    resolved_path_at_index TEXT NOT NULL,
    sample_key TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    mtime_ns INTEGER NOT NULL,
    compressed INTEGER NOT NULL CHECK (compressed IN (0, 1)),
    record_count INTEGER NOT NULL,
    feature_count INTEGER NOT NULL
);
CREATE TABLE records (
    record_pk INTEGER PRIMARY KEY,
    source_pk INTEGER NOT NULL REFERENCES sources(source_pk) ON DELETE CASCADE,
    record_index INTEGER NOT NULL,
    record_id TEXT NOT NULL,
    record_name TEXT NOT NULL,
    description TEXT NOT NULL,
    length INTEGER NOT NULL,
    topology TEXT,
    molecule_type TEXT,
    organism TEXT,
    strain TEXT,
    gc_percent REAL,
    UNIQUE(source_pk, record_index)
);
CREATE TABLE features (
    feature_pk INTEGER PRIMARY KEY,
    record_pk INTEGER NOT NULL REFERENCES records(record_pk) ON DELETE CASCADE,
    feature_index INTEGER NOT NULL,
    type TEXT NOT NULL,
    locus_tag TEXT,
    gene TEXT,
    product TEXT,
    protein_id TEXT,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    strand INTEGER,
    biological_length INTEGER NOT NULL,
    partial INTEGER NOT NULL CHECK (partial IN (0, 1)),
    pseudo INTEGER NOT NULL CHECK (pseudo IN (0, 1)),
    UNIQUE(record_pk, feature_index)
);
CREATE TABLE segments (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    ordinal INTEGER NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    PRIMARY KEY(feature_pk, ordinal)
);
CREATE TABLE qualifiers (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    key TEXT NOT NULL,
    ordinal INTEGER NOT NULL,
    value TEXT NOT NULL,
    PRIMARY KEY(feature_pk, key, ordinal)
);
CREATE TABLE xrefs (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    namespace TEXT NOT NULL,
    ordinal INTEGER NOT NULL,
    value TEXT NOT NULL,
    source_field TEXT NOT NULL,
    PRIMARY KEY(feature_pk, namespace, ordinal),
    UNIQUE(feature_pk, namespace, value, source_field)
);
CREATE INDEX idx_sources_sha256 ON sources(sha256);
CREATE UNIQUE INDEX uq_sources_source_identity ON sources(source_identity);
CREATE UNIQUE INDEX uq_sources_sample_key_nocase ON sources(sample_key COLLATE NOCASE);
CREATE INDEX idx_records_id ON records(record_id);
CREATE INDEX idx_records_name ON records(record_name);
CREATE INDEX idx_records_organism ON records(organism);
CREATE INDEX idx_features_type ON features(type);
CREATE INDEX idx_features_locus_tag ON features(locus_tag);
CREATE INDEX idx_features_gene ON features(gene);
CREATE INDEX idx_features_protein_id ON features(protein_id);
CREATE INDEX idx_features_coordinates ON features(record_pk, start, end);
CREATE INDEX idx_xrefs_namespace_value ON xrefs(namespace, value);
CREATE INDEX idx_qualifiers_key_value ON qualifiers(key, value);
"""

_SOURCES_V3_SQL = """
CREATE TABLE sources (
    source_pk INTEGER PRIMARY KEY,
    source_identity TEXT NOT NULL UNIQUE,
    display_path TEXT NOT NULL,
    resolved_path_at_index TEXT NOT NULL,
    sample_key TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    mtime_ns INTEGER NOT NULL,
    compressed INTEGER NOT NULL CHECK (compressed IN (0, 1)),
    record_count INTEGER NOT NULL,
    feature_count INTEGER NOT NULL
)
"""

_XREFS_V3_SQL = """
CREATE TABLE xrefs (
    feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
    namespace TEXT NOT NULL,
    ordinal INTEGER NOT NULL,
    value TEXT NOT NULL,
    source_field TEXT NOT NULL,
    PRIMARY KEY(feature_pk, namespace, ordinal),
    UNIQUE(feature_pk, namespace, value, source_field)
)
"""

_CURRENT_INDEX_STATEMENTS = (
    "CREATE INDEX idx_sources_sha256 ON sources(sha256)",
    "CREATE UNIQUE INDEX uq_sources_source_identity ON sources(source_identity)",
    "CREATE UNIQUE INDEX uq_sources_sample_key_nocase ON sources(sample_key COLLATE NOCASE)",
    "CREATE INDEX idx_records_id ON records(record_id)",
    "CREATE INDEX idx_records_name ON records(record_name)",
    "CREATE INDEX idx_records_organism ON records(organism)",
    "CREATE INDEX idx_features_type ON features(type)",
    "CREATE INDEX idx_features_locus_tag ON features(locus_tag)",
    "CREATE INDEX idx_features_gene ON features(gene)",
    "CREATE INDEX idx_features_protein_id ON features(protein_id)",
    "CREATE INDEX idx_features_coordinates ON features(record_pk, start, end)",
    "CREATE INDEX idx_xrefs_namespace_value ON xrefs(namespace, value)",
    "CREATE INDEX idx_qualifiers_key_value ON qualifiers(key, value)",
)
_CURRENT_INDEX_SQL = "\n".join(f"{statement};" for statement in _CURRENT_INDEX_STATEMENTS)

_EXPLICIT_INDEX_DEFINITIONS: dict[str, tuple[str, bool, bool, tuple[tuple[str, str], ...]]] = {
    "idx_sources_sha256": ("sources", False, False, (("sha256", "BINARY"),)),
    "uq_sources_source_identity": ("sources", True, False, (("source_identity", "BINARY"),)),
    "uq_sources_sample_key_nocase": ("sources", True, False, (("sample_key", "NOCASE"),)),
    "idx_records_id": ("records", False, False, (("record_id", "BINARY"),)),
    "idx_records_name": ("records", False, False, (("record_name", "BINARY"),)),
    "idx_records_organism": ("records", False, False, (("organism", "BINARY"),)),
    "idx_features_type": ("features", False, False, (("type", "BINARY"),)),
    "idx_features_locus_tag": ("features", False, False, (("locus_tag", "BINARY"),)),
    "idx_features_gene": ("features", False, False, (("gene", "BINARY"),)),
    "idx_features_protein_id": ("features", False, False, (("protein_id", "BINARY"),)),
    "idx_features_coordinates": (
        "features",
        False,
        False,
        (("record_pk", "BINARY"), ("start", "BINARY"), ("end", "BINARY")),
    ),
    "idx_xrefs_namespace_value": (
        "xrefs",
        False,
        False,
        (("namespace", "BINARY"), ("value", "BINARY")),
    ),
    "idx_qualifiers_key_value": (
        "qualifiers",
        False,
        False,
        (("key", "BINARY"), ("value", "BINARY")),
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def configure_new_writer_connection(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = DELETE")
    connection.execute("PRAGMA synchronous = FULL")


def configure_existing_writer_connection(connection: sqlite3.Connection) -> None:
    """Configure a writer without changing an existing journal mode."""

    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA synchronous = FULL")


def configure_connection(connection: sqlite3.Connection) -> None:
    """Backward-compatible alias for creating a new writer database."""

    configure_new_writer_connection(connection)


def configure_reader_connection(connection: sqlite3.Connection) -> None:
    """Configure a reader without changing the database journal mode."""

    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA query_only = ON")


def initialize_schema(connection: sqlite3.Connection, metadata: Mapping[str, str]) -> None:
    configure_new_writer_connection(connection)
    connection.executescript(SCHEMA_SQL)
    connection.execute(f"PRAGMA user_version = {SQLITE_USER_VERSION}")
    connection.executemany(
        "INSERT INTO metadata(key, value) VALUES (?, ?)",
        sorted(metadata.items()),
    )
    connection.commit()


def _legacy_source_identity(value: str) -> str:
    """Normalize a revision-1 stored path without requiring the file to exist."""

    try:
        resolved = Path(value).resolve(strict=False)
    except OSError:
        resolved = Path(os.path.abspath(value))
    return os.path.normcase(os.fspath(resolved))


def _rekey_stored_sources(connection: sqlite3.Connection) -> None:
    """Allocate deterministic keys using only persisted source metadata."""

    rows = connection.execute(
        "SELECT source_pk, source_identity, display_path, resolved_path_at_index, sha256 "
        "FROM sources ORDER BY source_identity, source_pk"
    ).fetchall()
    items = tuple(
        DiscoveredInput(
            requested_path=str(display_path),
            source_identity=str(source_identity),
            resolved_path=Path(str(resolved_path)),
            display_path=str(display_path),
            discovery_root=None,
        )
        for _source_pk, source_identity, display_path, resolved_path, _sha256 in rows
    )
    sha256_by_identity = {
        str(source_identity): str(sha256)
        for _source_pk, source_identity, _display_path, _resolved_path, sha256 in rows
    }
    keys = assign_sample_keys(items, sha256_by_identity=sha256_by_identity)
    connection.executemany(
        "UPDATE sources SET sample_key = ? WHERE source_pk = ?",
        [
            (f"__gbparse-rekey__{int(source_pk)}", int(source_pk))
            for source_pk, _source_identity, _display_path, _resolved_path, _sha256 in rows
        ],
    )
    connection.executemany(
        "UPDATE sources SET sample_key = ? WHERE source_pk = ?",
        [
            (keys[item], int(source_pk))
            for (source_pk, _source_identity, _display_path, _resolved_path, _sha256), item in zip(
                rows, items, strict=True
            )
        ],
    )


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")}


def _assert_existing_index_definition(
    connection: sqlite3.Connection,
    name: str,
    table: str,
    unique: bool,
    partial: bool,
    columns: tuple[tuple[str, str], ...],
) -> None:
    """Reject a same-name index with a definition other than the contract."""

    table_row = connection.execute(
        "SELECT tbl_name FROM sqlite_master WHERE type = 'index' AND name = ?",
        (name,),
    ).fetchone()
    if table_row is None:
        return
    if str(table_row[0]) != table:
        raise ValueError(f"invalid existing index {name!r}; expected table {table!r}")
    index_row = next(
        (row for row in connection.execute(f"PRAGMA index_list({table})") if str(row[1]) == name),
        None,
    )
    if index_row is None:
        raise ValueError(f"invalid existing index {name!r}; index inventory is inconsistent")
    actual_unique = int(index_row[2]) == 1
    actual_partial = bool(int(index_row[4])) if len(index_row) > 4 else False
    actual_columns = tuple(
        (str(row[2]), str(row[4] or "BINARY").upper())
        for row in connection.execute(f"PRAGMA index_xinfo('{name}')")
        if len(row) > 5 and int(row[5]) == 1
    )
    if (
        actual_unique != unique
        or actual_partial != partial
        or actual_columns != columns
    ):
        raise ValueError(
            f"invalid existing index {name!r}; expected "
            f"unique={unique}, partial={partial}, columns={columns}, "
            f"got unique={actual_unique}, partial={actual_partial}, columns={actual_columns}"
        )


def _validate_existing_index_definitions(connection: sqlite3.Connection) -> None:
    for name, (table, unique, partial, columns) in _EXPLICIT_INDEX_DEFINITIONS.items():
        _assert_existing_index_definition(connection, name, table, unique, partial, columns)


def _validate_current_index_definitions(connection: sqlite3.Connection) -> None:
    for name, (table, unique, partial, columns) in _EXPLICIT_INDEX_DEFINITIONS.items():
        if connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'index' AND name = ?",
            (name,),
        ).fetchone() is None:
            raise ValueError(f"invalid gbparse index; required index {name!r} is missing")
        _assert_existing_index_definition(connection, name, table, unique, partial, columns)


def _validate_legacy_contract(connection: sqlite3.Connection, revision: int) -> None:
    required_tables = {"metadata", "sources", "records", "features", "segments", "qualifiers", "xrefs"}
    tables = {
        str(row[0])
        for row in connection.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
    }
    missing_tables = sorted(required_tables - tables)
    if missing_tables:
        raise ValueError("incomplete legacy gbparse index; missing tables: " + ", ".join(missing_tables))
    source_columns = _table_columns(connection, "sources")
    required_source_columns = {
        "source_pk",
        "display_path",
        "resolved_path_at_index",
        "sample_key",
        "sha256",
        "size_bytes",
        "mtime_ns",
        "compressed",
        "record_count",
        "feature_count",
    }
    if revision == 2:
        required_source_columns.add("source_identity")
    missing_sources = sorted(required_source_columns - source_columns)
    if missing_sources:
        raise ValueError(
            "incomplete legacy gbparse index; sources is missing columns: "
            + ", ".join(missing_sources)
        )
    required_xref_columns = {"feature_pk", "namespace", "value", "source_field"}
    missing_xrefs = sorted(required_xref_columns - _table_columns(connection, "xrefs"))
    if missing_xrefs:
        raise ValueError(
            "incomplete legacy gbparse index; xrefs is missing columns: "
            + ", ".join(missing_xrefs)
        )


def _validate_current_schema_shape(connection: sqlite3.Connection) -> None:
    """Validate the structural portion of the revision-3 contract."""

    for table, required in _REQUIRED_COLUMNS.items():
        missing = sorted(required - _table_columns(connection, table))
        if missing:
            raise ValueError(
                f"incomplete gbparse index; table {table!r} is missing columns: "
                + ", ".join(missing)
            )
    source_info = {
        str(row[1]): row for row in connection.execute("PRAGMA table_info(sources)")
    }
    if int(source_info["source_identity"][3]) != 1:
        raise ValueError("invalid gbparse index; sources.source_identity must be NOT NULL")
    source_sql = connection.execute(
        "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = 'sources'"
    ).fetchone()[0]
    if "display_path TEXT NOT NULL UNIQUE" in str(source_sql):
        raise ValueError("invalid gbparse index; display_path must not be unique")
    if "ordinal" not in _table_columns(connection, "xrefs"):
        raise ValueError("invalid gbparse index; xrefs.ordinal is missing")
    _validate_current_index_definitions(connection)


def _integrity_check(connection: sqlite3.Connection) -> None:
    foreign = connection.execute("PRAGMA foreign_key_check").fetchall()
    if foreign:
        raise ValueError(f"index foreign-key check failed: {foreign[:3]}")
    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise ValueError(f"index integrity check failed: {integrity}")


def migrate_schema(connection: sqlite3.Connection) -> dict[str, str]:
    """Migrate revision-1 or revision-2 indexes to revision 3 transactionally."""

    configure_existing_writer_connection(connection)
    try:
        user_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
        metadata = read_metadata(connection)
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a readable gbparse index") from exc
    if metadata.get("schema_version") != INDEX_SCHEMA_VERSION:
        raise ValueError(f"unsupported index schema {metadata.get('schema_version', '<missing>')!r}")
    revision = metadata.get("schema_revision")
    if user_version == SQLITE_USER_VERSION and revision == str(INDEX_SCHEMA_REVISION):
        return metadata
    if (user_version, revision) not in {(1, "1"), (2, "2")}:
        raise ValueError("only revision-1 and revision-2 gbparse indexes can be migrated")

    revision_number = int(revision)
    _validate_legacy_contract(connection, revision_number)
    _validate_existing_index_definitions(connection)
    if revision_number == 1:
        source_rows = connection.execute(
            "SELECT source_pk, display_path, resolved_path_at_index, sample_key, sha256, "
            "size_bytes, mtime_ns, compressed, record_count, feature_count "
            "FROM sources ORDER BY source_pk"
        ).fetchall()
    else:
        source_rows = connection.execute(
            "SELECT source_pk, source_identity, display_path, resolved_path_at_index, sample_key, sha256, "
            "size_bytes, mtime_ns, compressed, record_count, feature_count "
            "FROM sources ORDER BY source_pk"
        ).fetchall()
    migrated_sources: list[tuple[object, ...]] = []
    identities: dict[str, int] = {}
    for row in source_rows:
        if revision_number == 1:
            source_pk, display_path, resolved_path, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count = row
            identity = _legacy_source_identity(str(resolved_path))
        else:
            source_pk, stored_identity, display_path, resolved_path, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count = row
            identity = str(stored_identity)
            if not identity:
                raise ValueError("cannot migrate legacy index: source_identity is empty")
        previous = identities.setdefault(identity, int(source_pk))
        if previous != int(source_pk):
            raise ValueError(
                "cannot migrate legacy index: source rows collapse to one canonical identity"
            )
        migrated_sources.append(
            (
                source_pk,
                identity,
                display_path,
                resolved_path,
                sample_key,
                sha256,
                size_bytes,
                mtime_ns,
                compressed,
                record_count,
                feature_count,
            )
        )

    xref_columns = _table_columns(connection, "xrefs")
    if "ordinal" in xref_columns:
        xref_rows = connection.execute(
            "SELECT feature_pk, namespace, ordinal, value, source_field "
            "FROM xrefs ORDER BY feature_pk, namespace, ordinal, rowid"
        ).fetchall()
    else:
        xref_rows = connection.execute(
            "SELECT feature_pk, namespace, value, source_field "
            "FROM xrefs ORDER BY feature_pk, namespace, rowid"
        ).fetchall()
    next_ordinal: dict[tuple[int, str], int] = {}
    migrated_xrefs: list[tuple[object, ...]] = []
    for row in xref_rows:
        if "ordinal" in xref_columns:
            feature_pk, namespace, _legacy_ordinal, value, source_field = row
        else:
            feature_pk, namespace, value, source_field = row
        key = (int(feature_pk), str(namespace))
        ordinal = next_ordinal.get(key, 0)
        next_ordinal[key] = ordinal + 1
        migrated_xrefs.append((feature_pk, namespace, ordinal, value, source_field))

    counts_before = {
        table: int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
        for table in ("sources", "records", "features", "segments", "qualifiers", "xrefs")
    }
    try:
        connection.execute("PRAGMA foreign_keys = OFF")
        connection.execute("BEGIN IMMEDIATE")
        connection.execute("DROP TABLE xrefs")
        connection.execute(_XREFS_V3_SQL)
        connection.executemany(
            "INSERT INTO xrefs(feature_pk, namespace, ordinal, value, source_field) "
            "VALUES (?, ?, ?, ?, ?)",
            migrated_xrefs,
        )
        connection.execute("DROP TABLE sources")
        connection.execute(_SOURCES_V3_SQL)
        connection.executemany(
            "INSERT INTO sources(source_pk, source_identity, display_path, resolved_path_at_index, sample_key, sha256, size_bytes, mtime_ns, compressed, record_count, feature_count) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            migrated_sources,
        )
        for index_name in _EXPLICIT_INDEX_DEFINITIONS:
            connection.execute(f"DROP INDEX IF EXISTS {index_name}")
        # The revision-3 sample-key index is deliberately installed only after
        # deterministic rekeying has made all stored rows collision-safe.
        _rekey_stored_sources(connection)
        for statement in _CURRENT_INDEX_STATEMENTS:
            connection.execute(statement)
        _validate_current_schema_shape(connection)
        _integrity_check(connection)
        counts_after = {
            table: int(connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
            for table in counts_before
        }
        if counts_after != counts_before:
            raise ValueError(
                f"legacy index row counts changed during migration: {counts_before} -> {counts_after}"
            )
        connection.execute("PRAGMA user_version = 3")
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES ('schema_revision', '3')"
        )
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES ('updated_at', ?)",
            (utc_now(),),
        )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.execute("PRAGMA foreign_keys = ON")
    return read_metadata(connection)


def read_metadata(connection: sqlite3.Connection) -> dict[str, str]:
    try:
        return {str(key): str(value) for key, value in connection.execute("SELECT key, value FROM metadata")}
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a compatible gbparse index: metadata table is unavailable") from exc


_REQUIRED_COLUMNS: dict[str, frozenset[str]] = {
    "metadata": frozenset({"key", "value"}),
    "sources": frozenset(
        {
            "source_pk",
            "source_identity",
            "display_path",
            "resolved_path_at_index",
            "sample_key",
            "sha256",
            "size_bytes",
            "mtime_ns",
            "compressed",
            "record_count",
            "feature_count",
        }
    ),
    "records": frozenset(
        {
            "record_pk",
            "source_pk",
            "record_index",
            "record_id",
            "record_name",
            "description",
            "length",
            "topology",
            "molecule_type",
            "organism",
            "strain",
            "gc_percent",
        }
    ),
    "features": frozenset(
        {
            "feature_pk",
            "record_pk",
            "feature_index",
            "type",
            "locus_tag",
            "gene",
            "product",
            "protein_id",
            "start",
            "end",
            "strand",
            "biological_length",
            "partial",
            "pseudo",
        }
    ),
    "segments": frozenset({"feature_pk", "ordinal", "start", "end"}),
    "qualifiers": frozenset({"feature_pk", "key", "ordinal", "value"}),
    "xrefs": frozenset({"feature_pk", "namespace", "ordinal", "value", "source_field"}),
}


def validate_schema(
    connection: sqlite3.Connection,
    *,
    read_only: bool = True,
) -> dict[str, str]:
    try:
        if read_only:
            configure_reader_connection(connection)
        else:
            configure_existing_writer_connection(connection)
        user_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a readable SQLite gbparse index") from exc
    if user_version != SQLITE_USER_VERSION:
        raise ValueError(f"unsupported gbparse index SQLite user_version {user_version}; expected {SQLITE_USER_VERSION}")
    try:
        tables = {
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        missing_tables = sorted(set(_REQUIRED_COLUMNS) - tables)
        if missing_tables:
            raise ValueError(
                "incomplete gbparse index; missing tables: " + ", ".join(missing_tables)
            )
        for table, required in _REQUIRED_COLUMNS.items():
            columns = {
                str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")
            }
            missing_columns = sorted(required - columns)
            if missing_columns:
                raise ValueError(
                    f"incomplete gbparse index; table {table!r} is missing columns: "
                    + ", ".join(missing_columns)
                )
        index_rows = list(connection.execute("PRAGMA index_list('sources')"))
        sample_key_index = next(
            (row for row in index_rows if str(row[1]) == "uq_sources_sample_key_nocase"),
            None,
        )
        if sample_key_index is None:
            raise ValueError(
                "incomplete gbparse index; required sample-key uniqueness index is missing"
            )
        unique = int(sample_key_index[2])
        partial = int(sample_key_index[4]) if len(sample_key_index) > 4 else 0
        if unique != 1 or partial != 0:
            raise ValueError(
                "invalid gbparse index; sample-key uniqueness index must be "
                "unique and non-partial"
            )
        key_columns = [
            (str(row[2]), str(row[4] or "BINARY").upper())
            for row in connection.execute(
                "PRAGMA index_xinfo('uq_sources_sample_key_nocase')"
            )
            if int(row[5]) == 1
        ]
        if key_columns != [("sample_key", "NOCASE")]:
            raise ValueError(
                "invalid gbparse index; sample-key index must be "
                "sample_key COLLATE NOCASE"
            )
        metadata = read_metadata(connection)
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a readable SQLite gbparse index") from exc
    if metadata.get("schema_version") != INDEX_SCHEMA_VERSION:
        raise ValueError(f"unsupported index schema {metadata.get('schema_version', '<missing>')!r}")
    if metadata.get("schema_revision") != str(INDEX_SCHEMA_REVISION):
        raise ValueError("unsupported gbparse index schema revision")
    return metadata


__all__ = [
    "INDEX_SCHEMA_REVISION",
    "INDEX_SCHEMA_VERSION",
    "SQLITE_USER_VERSION",
    "configure_connection",
    "configure_reader_connection",
    "initialize_schema",
    "migrate_schema",
    "read_metadata",
    "utc_now",
    "validate_schema",
]
