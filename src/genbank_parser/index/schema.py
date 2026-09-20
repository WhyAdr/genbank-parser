"""SQLite schema and compatibility checks for ``gbparse.index.v1``."""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

INDEX_SCHEMA_VERSION = "gbparse.index.v1"
INDEX_SCHEMA_REVISION = 2
SQLITE_USER_VERSION = 2

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


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def configure_connection(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = DELETE")
    connection.execute("PRAGMA synchronous = FULL")


def initialize_schema(connection: sqlite3.Connection, metadata: Mapping[str, str]) -> None:
    configure_connection(connection)
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


def migrate_schema(connection: sqlite3.Connection) -> dict[str, str]:
    """Migrate a revision-1 index to revision 2 inside one transaction."""

    configure_connection(connection)
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
    if user_version != 1 or revision != "1":
        raise ValueError("only revision-1 gbparse indexes can be migrated")

    try:
        connection.execute("BEGIN IMMEDIATE")
        source_rows = connection.execute(
            "SELECT source_pk, resolved_path_at_index FROM sources ORDER BY source_pk"
        ).fetchall()
        identities: dict[str, int] = {}
        for source_pk, resolved_path in source_rows:
            identity = _legacy_source_identity(str(resolved_path))
            previous = identities.setdefault(identity, int(source_pk))
            if previous != int(source_pk):
                raise ValueError(
                    "cannot migrate revision-1 index: source rows collapse to one canonical identity"
                )
        connection.execute("ALTER TABLE sources ADD COLUMN source_identity TEXT")
        connection.executemany(
            "UPDATE sources SET source_identity = ? WHERE source_pk = ?",
            [(_legacy_source_identity(str(path)), int(source_pk)) for source_pk, path in source_rows],
        )
        connection.execute("CREATE UNIQUE INDEX uq_sources_source_identity ON sources(source_identity)")

        connection.execute(
            """
            CREATE TABLE xrefs_v2 (
                feature_pk INTEGER NOT NULL REFERENCES features(feature_pk) ON DELETE CASCADE,
                namespace TEXT NOT NULL,
                ordinal INTEGER NOT NULL,
                value TEXT NOT NULL,
                source_field TEXT NOT NULL,
                PRIMARY KEY(feature_pk, namespace, ordinal),
                UNIQUE(feature_pk, namespace, value, source_field)
            )
            """
        )
        connection.execute(
            """
            INSERT INTO xrefs_v2(feature_pk, namespace, ordinal, value, source_field)
            SELECT feature_pk, namespace,
                   ROW_NUMBER() OVER (PARTITION BY feature_pk, namespace ORDER BY rowid) - 1,
                   value, source_field
            FROM xrefs
            ORDER BY feature_pk, namespace, rowid
            """
        )
        connection.execute("DROP TABLE xrefs")
        connection.execute("ALTER TABLE xrefs_v2 RENAME TO xrefs")
        connection.execute("CREATE INDEX idx_xrefs_namespace_value ON xrefs(namespace, value)")
        connection.execute("PRAGMA user_version = 2")
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES ('schema_revision', '2')"
        )
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES ('updated_at', ?)",
            (utc_now(),),
        )
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return read_metadata(connection)


def read_metadata(connection: sqlite3.Connection) -> dict[str, str]:
    try:
        return {str(key): str(value) for key, value in connection.execute("SELECT key, value FROM metadata")}
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a compatible gbparse index: metadata table is unavailable") from exc


def validate_schema(connection: sqlite3.Connection) -> dict[str, str]:
    try:
        configure_connection(connection)
        user_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
    except sqlite3.DatabaseError as exc:
        raise ValueError("not a readable SQLite gbparse index") from exc
    if user_version != SQLITE_USER_VERSION:
        raise ValueError(f"unsupported gbparse index SQLite user_version {user_version}; expected {SQLITE_USER_VERSION}")
    metadata = read_metadata(connection)
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
    "initialize_schema",
    "migrate_schema",
    "read_metadata",
    "utc_now",
    "validate_schema",
]
