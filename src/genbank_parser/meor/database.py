"""Validated loader for the packaged MEOR marker knowledge base."""
from __future__ import annotations

from importlib import resources
from pathlib import Path
import re
from typing import Any

import yaml

from .models import (
    MeorCategory,
    MeorDatabase,
    MeorMarker,
    MeorPathway,
    MeorPathwayStep,
)


class MeorDatabaseError(ValueError):
    """Raised when MEOR package data violates its schema or references."""


def _read_yaml(path: str | Path | None, packaged_name: str) -> dict[str, Any]:
    try:
        if path is None:
            resource = resources.files("genbank_parser").joinpath(
                "data", "meor", packaged_name
            )
            text = resource.read_text(encoding="utf-8")
        else:
            text = Path(path).read_text(encoding="utf-8")
        payload = yaml.safe_load(text)
    except (OSError, yaml.YAMLError) as exc:
        raise MeorDatabaseError(f"Could not load {packaged_name}: {exc}") from exc
    if not isinstance(payload, dict):
        raise MeorDatabaseError(f"{packaged_name} must contain a YAML mapping")
    if payload.get("schema_version") != 1:
        raise MeorDatabaseError(
            f"Unsupported {packaged_name} schema_version: "
            f"{payload.get('schema_version')!r}"
        )
    return payload


def _string_list(value: Any, field: str, owner: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item for item in value
    ):
        raise MeorDatabaseError(f"{owner}.{field} must be a list of strings")
    return tuple(value)


def _unique_ids(items: list[dict[str, Any]], kind: str) -> None:
    ids = [item.get("id") for item in items]
    if any(not isinstance(item_id, str) or not item_id for item_id in ids):
        raise MeorDatabaseError(f"Every {kind} requires a non-empty string id")
    duplicates = sorted({item_id for item_id in ids if ids.count(item_id) > 1})
    if duplicates:
        raise MeorDatabaseError(f"Duplicate {kind} IDs: {', '.join(duplicates)}")


def load_meor_database(
    markers_path: str | Path | None = None,
    pathways_path: str | Path | None = None,
) -> MeorDatabase:
    marker_data = _read_yaml(markers_path, "markers.yaml")
    pathway_data = _read_yaml(pathways_path, "pathways.yaml")
    provenance = _read_yaml(None, "provenance.yaml")

    raw_categories = marker_data.get("categories")
    raw_markers = marker_data.get("markers")
    raw_pathways = pathway_data.get("pathways")
    if not isinstance(raw_categories, list) or not all(
        isinstance(item, dict) for item in raw_categories
    ):
        raise MeorDatabaseError("markers.yaml categories must be a list of mappings")
    if not isinstance(raw_markers, list) or not all(
        isinstance(item, dict) for item in raw_markers
    ):
        raise MeorDatabaseError("markers.yaml markers must be a list of mappings")
    if not isinstance(raw_pathways, list) or not all(
        isinstance(item, dict) for item in raw_pathways
    ):
        raise MeorDatabaseError("pathways.yaml pathways must be a list of mappings")

    _unique_ids(raw_categories, "category")
    _unique_ids(raw_markers, "marker")
    _unique_ids(raw_pathways, "pathway")

    categories = tuple(
        MeorCategory(
            id=item["id"],
            title=str(item.get("title", "")),
            role=str(item.get("role", "")),
        )
        for item in raw_categories
    )
    category_ids = {category.id for category in categories}

    markers: list[MeorMarker] = []
    for item in raw_markers:
        marker_id = item["id"]
        category_id = item.get("category")
        if category_id not in category_ids:
            raise MeorDatabaseError(
                f"Marker {marker_id} references missing category {category_id!r}"
            )
        gene_patterns = _string_list(item.get("genes", []), "genes", marker_id)
        product_patterns = _string_list(
            item.get("products", []), "products", marker_id
        )
        note_patterns = _string_list(item.get("notes", []), "notes", marker_id)
        for field, patterns in (
            ("genes", gene_patterns),
            ("products", product_patterns),
            ("notes", note_patterns),
        ):
            for pattern in patterns:
                try:
                    re.compile(pattern, re.IGNORECASE)
                except re.error as exc:
                    raise MeorDatabaseError(
                        f"Invalid regex in {marker_id}.{field}: {pattern!r}: {exc}"
                    ) from exc
        markers.append(
            MeorMarker(
                id=marker_id,
                name=str(item.get("name", "")),
                category_id=category_id,
                gene_patterns=gene_patterns,
                ecs=_string_list(item.get("ecs", []), "ecs", marker_id),
                kos=_string_list(item.get("kos", []), "kos", marker_id),
                product_patterns=product_patterns,
                note_patterns=note_patterns,
                sources=_string_list(item.get("sources", []), "sources", marker_id),
            )
        )

    marker_ids = {marker.id for marker in markers}
    pathways: list[MeorPathway] = []
    for item in raw_pathways:
        pathway_id = item["id"]
        category_id = item.get("category")
        if category_id not in category_ids:
            raise MeorDatabaseError(
                f"Pathway {pathway_id} references missing category {category_id!r}"
            )
        raw_steps = item.get("steps")
        if not isinstance(raw_steps, list) or not raw_steps:
            raise MeorDatabaseError(f"Pathway {pathway_id} requires at least one step")
        steps: list[MeorPathwayStep] = []
        for index, step in enumerate(raw_steps, 1):
            if not isinstance(step, dict):
                raise MeorDatabaseError(
                    f"Pathway {pathway_id} step {index} must be a mapping"
                )
            any_of = _string_list(step.get("any_of"), "any_of", pathway_id)
            if not any_of:
                raise MeorDatabaseError(
                    f"Pathway {pathway_id} step {index} requires any_of markers"
                )
            unresolved = sorted(set(any_of) - marker_ids)
            if unresolved:
                raise MeorDatabaseError(
                    f"Pathway {pathway_id} has unresolved markers: "
                    f"{', '.join(unresolved)}"
                )
            steps.append(MeorPathwayStep(str(step.get("name", "")), any_of))
        pathways.append(
            MeorPathway(
                id=pathway_id,
                name=str(item.get("name", "")),
                category_id=category_id,
                steps=tuple(steps),
            )
        )

    if markers_path is None and pathways_path is None:
        expected = (len(categories), len(markers), len(pathways))
        if expected != (9, 48, 7):
            raise MeorDatabaseError(
                "Packaged MEOR catalog must contain 9 categories, 48 markers, "
                f"and 7 pathways; found {expected}"
            )
        provenance_markers = provenance.get("markers", {})
        missing = sorted(marker_ids - set(provenance_markers))
        if missing:
            raise MeorDatabaseError(
                f"Markers missing scientific provenance: {', '.join(missing)}"
            )
        known_sources = set(provenance.get("sources", {}))
        for marker in markers:
            marker_provenance = provenance_markers[marker.id]
            provenance_sources = marker_provenance.get("sources", [])
            if tuple(provenance_sources) != marker.sources:
                raise MeorDatabaseError(
                    f"Marker {marker.id} sources disagree between markers.yaml "
                    "and provenance.yaml"
                )
            unresolved_sources = sorted(set(marker.sources) - known_sources)
            if unresolved_sources:
                raise MeorDatabaseError(
                    f"Marker {marker.id} has unresolved provenance sources: "
                    f"{', '.join(unresolved_sources)}"
                )
            if not marker_provenance.get("interpretation"):
                raise MeorDatabaseError(
                    f"Marker {marker.id} lacks a provenance interpretation"
                )

    return MeorDatabase(
        categories=categories,
        markers=tuple(markers),
        pathways=tuple(pathways),
        provenance=provenance,
    )
