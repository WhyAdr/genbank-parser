"""Compact status and integrity reports for persistent cohort indexes."""

from __future__ import annotations

import sqlite3
from pathlib import Path

from ..discovery import fingerprint_file
from .schema import validate_schema


def inspect_index(
    database: str | Path,
    *,
    check: bool = False,
    verify_sources: bool = False,
) -> dict[str, object]:
    """Return a schema-versioned index status without hashing by default."""

    path = Path(database)
    if not path.is_file():
        raise FileNotFoundError(f"index does not exist: {path}")
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    try:
        metadata = validate_schema(connection, read_only=True)
        source_count = int(connection.execute("SELECT COUNT(*) FROM sources").fetchone()[0])
        record_count = int(connection.execute("SELECT COUNT(*) FROM records").fetchone()[0])
        feature_count = int(connection.execute("SELECT COUNT(*) FROM features").fetchone()[0])
        feature_types = {
            str(row[0]): int(row[1])
            for row in connection.execute("SELECT type, COUNT(*) FROM features GROUP BY type ORDER BY type")
        }
        topologies = {
            str(row[0] or "unknown"): int(row[1])
            for row in connection.execute("SELECT topology, COUNT(*) FROM records GROUP BY topology ORDER BY topology")
        }
        duplicate_hashes = [
            {"sha256": str(row[0]), "source_count": int(row[1])}
            for row in connection.execute(
                "SELECT sha256, COUNT(*) FROM sources GROUP BY sha256 HAVING COUNT(*) > 1 ORDER BY sha256"
            )
        ]
        source_checks: list[dict[str, object]] = []
        if verify_sources:
            for row in connection.execute(
                "SELECT display_path, resolved_path_at_index, sha256 FROM sources ORDER BY display_path"
            ):
                source_path = Path(row[1])
                item: dict[str, object] = {
                    "display_path": row[0],
                    "resolved_path": row[1],
                    "exists": source_path.is_file(),
                    "sha256_matches": False,
                }
                if source_path.is_file():
                    item["sha256_matches"] = fingerprint_file(source_path).sha256 == row[2]
                source_checks.append(item)
        report: dict[str, object] = {
            "schema_version": "gbparse.index-status.v1",
            "database": str(path),
            "metadata": metadata,
            "counts": {
                "sources": source_count,
                "records": record_count,
                "features": feature_count,
            },
            "feature_type_counts": feature_types,
            "topology_counts": topologies,
            "database_size_bytes": path.stat().st_size,
            "duplicate_content_hashes": duplicate_hashes,
        }
        if verify_sources:
            report["source_checks"] = source_checks
        if check:
            foreign_key = connection.execute("PRAGMA foreign_key_check").fetchall()
            quick_check = str(connection.execute("PRAGMA quick_check").fetchone()[0])
            report["checks"] = {
                "quick_check": quick_check,
                "foreign_key_check": "ok" if not foreign_key else [list(row) for row in foreign_key],
            }
        return report
    finally:
        connection.close()


__all__ = ["inspect_index"]
