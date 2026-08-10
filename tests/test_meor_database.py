"""Integrity tests for the externalized MEOR knowledge base."""
from __future__ import annotations

from pathlib import Path
import re

import pytest
import yaml

from genbank_parser.meor.database import MeorDatabaseError, load_meor_database


def test_packaged_meor_database_integrity() -> None:
    database = load_meor_database()
    assert database.catalog_version == "1.1"
    assert len(database.categories) == 9
    assert len(database.markers) == 48
    assert len(database.pathways) == 7
    assert len(database.marker_map) == 48
    assert len(database.category_map) == 9
    provenance = database.provenance["markers"]
    assert set(provenance) == set(database.marker_map)
    known_sources = set(database.provenance["sources"])
    for marker in database.markers:
        assert marker.sources
        assert set(marker.sources) <= known_sources
        assert tuple(provenance[marker.id]["sources"]) == marker.sources
        assert provenance[marker.id]["interpretation"]
        for pattern in (
            *marker.gene_patterns,
            *marker.product_patterns,
            *marker.note_patterns,
        ):
            re.compile(pattern, re.IGNORECASE)
    for pathway in database.pathways:
        for step in pathway.steps:
            assert step.any_of
            assert set(step.any_of) <= set(database.marker_map)


def test_loader_rejects_duplicate_marker_ids(tmp_path: Path) -> None:
    source = Path("src/genbank_parser/data/meor/markers.yaml")
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["markers"].append(dict(payload["markers"][0]))
    custom = tmp_path / "markers.yaml"
    custom.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    with pytest.raises(MeorDatabaseError, match="Duplicate marker IDs"):
        load_meor_database(markers_path=custom)


def test_loader_rejects_invalid_marker_regex(tmp_path: Path) -> None:
    source = Path("src/genbank_parser/data/meor/markers.yaml")
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["markers"][0]["genes"] = ["["]
    custom = tmp_path / "markers.yaml"
    custom.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    with pytest.raises(MeorDatabaseError, match="Invalid regex"):
        load_meor_database(markers_path=custom)


def test_custom_catalog_without_version_is_backward_compatible(tmp_path: Path) -> None:
    source = Path("src/genbank_parser/data/meor/markers.yaml")
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload.pop("catalog_version")
    custom = tmp_path / "markers.yaml"
    custom.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    assert load_meor_database(markers_path=custom).catalog_version == "custom"
