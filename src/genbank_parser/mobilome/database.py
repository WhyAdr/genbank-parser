"""Fail-closed loader for the versioned mobilome evidence catalog."""

from __future__ import annotations

import hashlib
import json
import re
from importlib import resources
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from .models import (
    CatalogResource,
    DisabledAggregateRule,
    ExternalHandoff,
    Facet,
    InferenceComponent,
    InferenceRule,
    MarkerMatcher,
    MobilomeDatabase,
    MobilomeDatabaseError,
    MobilomeMarker,
    ProvenanceSource,
)

_DATABASE_RESOURCE_NAMES = ("markers.yaml", "provenance.yaml", "inference.yaml")
_REPORT_SCHEMA_NAME = "report.schema.json"
_VALID_MATCH_FIELDS = frozenset(
    {
        "feature_type",
        "gene",
        "product",
        "note",
        "db_xref",
        "inference",
        "function",
        "mobile_element_type",
        "regulatory_class",
        "rpt_type",
    }
)
_FIELD_STRENGTH_CEILINGS = {
    "feature_type": 3,
    "regulatory_class": 3,
    "rpt_type": 3,
    "db_xref": 3,
    "gene": 2,
    "product": 2,
    "inference": 2,
    "function": 2,
    "mobile_element_type": 2,
    "note": 1,
}
_VALID_MODES = frozenset({"exact", "casefold_exact", "regex"})
_VALID_STATUSES = frozenset({"supported", "tentative", "insufficient", "conflicting"})
_VALID_SOURCE_TYPES = frozenset(
    {"primary_article", "official_documentation", "curated_database", "local_rule"}
)
_CLASSIFICATION_RULE_IDS = frozenset(
    {
        "source_plasmid_qualifier",
        "source_chromosome_qualifier",
        "record_annotation_plasmid",
        "record_annotation_chromosome",
        "bounded_description_plasmid_token",
        "bounded_description_chromosome_token",
        "bounded_name_plasmid_token",
        "bounded_name_chromosome_token",
    }
)


class _UniqueKeyLoader(yaml.SafeLoader):
    """Safe YAML loader that rejects duplicate mapping keys."""


def _construct_unique_mapping(
    loader: yaml.SafeLoader, node: yaml.Node, deep: bool = False
) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise yaml.constructor.ConstructorError(
                "while constructing a mapping",
                node.start_mark,
                f"found duplicate key {key!r}",
                key_node.start_mark,
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping
)


def _expect_mapping(value: Any, owner: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise MobilomeDatabaseError(f"{owner} must be a mapping")
    return value


def _expect_keys(value: dict[str, Any], owner: str, allowed: set[str]) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise MobilomeDatabaseError(f"{owner} has unknown fields: {', '.join(unknown)}")


def _string(value: Any, owner: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise MobilomeDatabaseError(f"{owner} must be a non-empty string")
    return value.strip()


def _string_list(
    value: Any, owner: str, *, allow_empty: bool = False
) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item.strip() for item in value
    ):
        raise MobilomeDatabaseError(f"{owner} must be a list of non-empty strings")
    result = tuple(item.strip() for item in value)
    if not allow_empty and not result:
        raise MobilomeDatabaseError(f"{owner} must not be empty")
    if len(set(result)) != len(result):
        raise MobilomeDatabaseError(f"{owner} must not contain duplicate values")
    return result


def _bool(value: Any, owner: str) -> bool:
    if not isinstance(value, bool):
        raise MobilomeDatabaseError(f"{owner} must be boolean")
    return value


def _positive_int(value: Any, owner: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise MobilomeDatabaseError(f"{owner} must be a positive integer")
    return value


def _unique_ids(items: list[dict[str, Any]], owner: str) -> None:
    ids = [_string(item.get("id"), f"{owner}.id") for item in items]
    duplicates = sorted({item_id for item_id in ids if ids.count(item_id) > 1})
    if duplicates:
        raise MobilomeDatabaseError(f"Duplicate {owner} IDs: {', '.join(duplicates)}")


def _major_version(value: str, owner: str) -> str:
    match = re.fullmatch(r"(\d+)(?:\.\d+){0,2}(?:[-+][A-Za-z0-9.-]+)?", value)
    if match is None:
        raise MobilomeDatabaseError(f"{owner} must be a semantic-style version")
    return match.group(1)


def _resource_bytes(
    database_dir: Path | None,
) -> tuple[dict[str, bytes], tuple[Path, ...]]:
    if database_dir is None:
        root = resources.files("genbank_parser").joinpath("data", "mobilome")
        try:
            payload = {
                name: root.joinpath(name).read_bytes()
                for name in (*_DATABASE_RESOURCE_NAMES, _REPORT_SCHEMA_NAME)
            }
        except OSError as exc:
            raise MobilomeDatabaseError(
                f"Could not load packaged mobilome resources: {exc}"
            ) from exc
        return payload, ()

    try:
        resolved_dir = database_dir.expanduser().resolve(strict=True)
    except OSError as exc:
        raise MobilomeDatabaseError(f"Could not resolve --database-dir: {exc}") from exc
    if not resolved_dir.is_dir():
        raise MobilomeDatabaseError("--database-dir must name a directory")

    names = {path.name for path in resolved_dir.iterdir() if path.is_file()}
    expected = set(_DATABASE_RESOURCE_NAMES)
    if names != expected:
        raise MobilomeDatabaseError(
            "Custom database directory must contain exactly "
            f"{', '.join(_DATABASE_RESOURCE_NAMES)}"
        )
    paths = tuple(
        (resolved_dir / name).resolve(strict=True) for name in _DATABASE_RESOURCE_NAMES
    )
    try:
        payload = {path.name: path.read_bytes() for path in paths}
        schema = (
            resources.files("genbank_parser")
            .joinpath("data", "mobilome", _REPORT_SCHEMA_NAME)
            .read_bytes()
        )
    except OSError as exc:
        raise MobilomeDatabaseError(
            f"Could not read custom mobilome resources: {exc}"
        ) from exc
    payload[_REPORT_SCHEMA_NAME] = schema
    return payload, paths


def _load_yaml(raw: bytes, name: str) -> dict[str, Any]:
    try:
        loaded = yaml.load(raw.decode("utf-8"), Loader=_UniqueKeyLoader)
    except (UnicodeDecodeError, yaml.YAMLError) as exc:
        raise MobilomeDatabaseError(f"Could not load {name}: {exc}") from exc
    payload = _expect_mapping(loaded, name)
    if payload.get("schema_version") != 1:
        raise MobilomeDatabaseError(
            f"Unsupported {name} schema_version: {payload.get('schema_version')!r}"
        )
    return payload


def _load_schema(raw: bytes) -> None:
    try:
        schema = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise MobilomeDatabaseError(
            f"Could not load {_REPORT_SCHEMA_NAME}: {exc}"
        ) from exc
    if (
        not isinstance(schema, dict)
        or schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema"
    ):
        raise MobilomeDatabaseError(
            "report.schema.json must declare JSON Schema draft 2020-12"
        )
    if (
        schema.get("properties", {}).get("schema_version", {}).get("const")
        != "gbparse.mobilome.v1"
    ):
        raise MobilomeDatabaseError(
            "report.schema.json must require gbparse.mobilome.v1"
        )


def _parse_facets(markers_data: dict[str, Any]) -> tuple[Facet, ...]:
    _expect_keys(
        markers_data,
        "markers.yaml",
        {"schema_version", "catalog_version", "facets", "markers"},
    )
    raw_facets = markers_data.get("facets")
    if not isinstance(raw_facets, list) or not all(
        isinstance(item, dict) for item in raw_facets
    ):
        raise MobilomeDatabaseError("markers.yaml.facets must be a list of mappings")
    _unique_ids(raw_facets, "facet")
    facets: list[Facet] = []
    for raw in raw_facets:
        _expect_keys(
            raw, f"facet {_string(raw.get('id'), 'facet.id')}", {"id", "title"}
        )
        facets.append(
            Facet(
                id=_string(raw["id"], "facet.id"),
                title=_string(raw.get("title"), "facet.title"),
            )
        )
    return tuple(facets)


def _parse_markers(
    markers_data: dict[str, Any], facet_ids: set[str]
) -> tuple[MobilomeMarker, ...]:
    raw_markers = markers_data.get("markers")
    if not isinstance(raw_markers, list) or not all(
        isinstance(item, dict) for item in raw_markers
    ):
        raise MobilomeDatabaseError("markers.yaml.markers must be a list of mappings")
    _unique_ids(raw_markers, "marker")
    markers: list[MobilomeMarker] = []
    allowed = {
        "id",
        "label",
        "facets",
        "feature_types",
        "matchers",
        "requires_non_pseudo",
        "requires_complete",
        "sources",
        "interpretation",
        "limitations",
    }
    for raw_marker in raw_markers:
        marker_id = _string(raw_marker.get("id"), "marker.id")
        _expect_keys(raw_marker, f"marker {marker_id}", allowed)
        facets = _string_list(raw_marker.get("facets"), f"marker {marker_id}.facets")
        missing_facets = sorted(set(facets) - facet_ids)
        if missing_facets:
            raise MobilomeDatabaseError(
                f"Marker {marker_id} references unknown facets: {', '.join(missing_facets)}"
            )
        feature_types = _string_list(
            raw_marker.get("feature_types"), f"marker {marker_id}.feature_types"
        )
        if len({item.casefold() for item in feature_types}) != len(feature_types):
            raise MobilomeDatabaseError(f"Marker {marker_id} duplicates a feature type")
        requires_non_pseudo = _bool(
            raw_marker.get("requires_non_pseudo"),
            f"marker {marker_id}.requires_non_pseudo",
        )
        requires_complete = _bool(
            raw_marker.get("requires_complete"), f"marker {marker_id}.requires_complete"
        )
        if "cds" in {
            item.casefold() for item in feature_types
        } and "pseudogene" not in {item.casefold() for item in feature_types}:
            raise MobilomeDatabaseError(
                f"Marker {marker_id} scans CDS but must also retain pseudogene observations"
            )
        raw_matchers = raw_marker.get("matchers")
        if (
            not isinstance(raw_matchers, list)
            or not raw_matchers
            or not all(isinstance(item, dict) for item in raw_matchers)
        ):
            raise MobilomeDatabaseError(
                f"Marker {marker_id}.matchers must be a non-empty list"
            )
        matchers: list[MarkerMatcher] = []
        for index, raw_matcher in enumerate(raw_matchers):
            owner = f"marker {marker_id}.matchers[{index}]"
            _expect_keys(
                raw_matcher,
                owner,
                {"field", "mode", "patterns", "negative_patterns", "strength"},
            )
            field = _string(raw_matcher.get("field"), f"{owner}.field")
            if field not in _VALID_MATCH_FIELDS:
                raise MobilomeDatabaseError(f"{owner}.field is not supported: {field}")
            mode = _string(raw_matcher.get("mode"), f"{owner}.mode")
            if mode not in _VALID_MODES:
                raise MobilomeDatabaseError(f"{owner}.mode is not supported: {mode}")
            patterns = _string_list(raw_matcher.get("patterns"), f"{owner}.patterns")
            negative_patterns = _string_list(
                raw_matcher.get("negative_patterns"),
                f"{owner}.negative_patterns",
                allow_empty=True,
            )
            strength = _positive_int(raw_matcher.get("strength"), f"{owner}.strength")
            if strength > 3:
                raise MobilomeDatabaseError(f"{owner}.strength must be between 1 and 3")
            if strength > _FIELD_STRENGTH_CEILINGS[field]:
                raise MobilomeDatabaseError(
                    f"{owner}.strength exceeds the {field} evidence ceiling"
                )
            if mode == "regex":
                for pattern in patterns:
                    try:
                        re.compile(pattern)
                    except re.error as exc:
                        raise MobilomeDatabaseError(
                            f"Invalid regex in {owner}: {pattern!r}: {exc}"
                        ) from exc
            for pattern in negative_patterns:
                try:
                    re.compile(pattern)
                except re.error as exc:
                    raise MobilomeDatabaseError(
                        f"Invalid negative regex in {owner}: {pattern!r}: {exc}"
                    ) from exc
            matchers.append(
                MarkerMatcher(
                    field=field,
                    mode=mode,  # type: ignore[arg-type]
                    patterns=patterns,
                    negative_patterns=negative_patterns,
                    strength=strength,
                )
            )
        markers.append(
            MobilomeMarker(
                id=marker_id,
                label=_string(raw_marker.get("label"), f"marker {marker_id}.label"),
                facets=tuple(sorted(facets)),
                feature_types=feature_types,
                matchers=tuple(matchers),
                requires_non_pseudo=requires_non_pseudo,
                requires_complete=requires_complete,
                sources=tuple(
                    sorted(
                        _string_list(
                            raw_marker.get("sources"), f"marker {marker_id}.sources"
                        )
                    )
                ),
                interpretation=_string(
                    raw_marker.get("interpretation"),
                    f"marker {marker_id}.interpretation",
                ),
                limitations=_string_list(
                    raw_marker.get("limitations"), f"marker {marker_id}.limitations"
                ),
            )
        )
    return tuple(markers)


def _parse_provenance(provenance_data: dict[str, Any]) -> tuple[ProvenanceSource, ...]:
    _expect_keys(
        provenance_data,
        "provenance.yaml",
        {"schema_version", "provenance_version", "sources"},
    )
    raw_sources = provenance_data.get("sources")
    if not isinstance(raw_sources, list) or not all(
        isinstance(item, dict) for item in raw_sources
    ):
        raise MobilomeDatabaseError(
            "provenance.yaml.sources must be a list of mappings"
        )
    _unique_ids(raw_sources, "provenance source")
    sources: list[ProvenanceSource] = []
    for raw_source in raw_sources:
        source_id = _string(raw_source.get("id"), "provenance source.id")
        source_type = _string(raw_source.get("type"), f"source {source_id}.type")
        if source_type not in _VALID_SOURCE_TYPES:
            raise MobilomeDatabaseError(
                f"Source {source_id} has unknown type {source_type!r}"
            )
        common = {"id", "type", "title", "url", "supports", "limitations"}
        if source_type == "primary_article":
            allowed = common | {"year", "doi", "pmid", "pmcid"}
            _expect_keys(raw_source, f"source {source_id}", allowed)
            _string(raw_source.get("title"), f"source {source_id}.title")
            _positive_int(raw_source.get("year"), f"source {source_id}.year")
            if not any(raw_source.get(key) for key in ("doi", "pmid", "pmcid", "url")):
                raise MobilomeDatabaseError(
                    f"Source {source_id} needs a stable identifier"
                )
            _string_list(raw_source.get("supports"), f"source {source_id}.supports")
            _string_list(
                raw_source.get("limitations"), f"source {source_id}.limitations"
            )
        elif source_type == "official_documentation":
            allowed = common | {"accessed"}
            _expect_keys(raw_source, f"source {source_id}", allowed)
            _string(raw_source.get("title"), f"source {source_id}.title")
            _string(raw_source.get("url"), f"source {source_id}.url")
            _string(raw_source.get("accessed"), f"source {source_id}.accessed")
            _string_list(raw_source.get("supports"), f"source {source_id}.supports")
            _string_list(
                raw_source.get("limitations"), f"source {source_id}.limitations"
            )
        elif source_type == "curated_database":
            allowed = common | {
                "name",
                "release",
                "database_version",
                "taxonomic_scope",
                "thresholds",
            }
            _expect_keys(raw_source, f"source {source_id}", allowed)
            _string(raw_source.get("name"), f"source {source_id}.name")
            if not raw_source.get("release") and not raw_source.get("database_version"):
                raise MobilomeDatabaseError(
                    f"Source {source_id} needs release or database_version"
                )
            _string(raw_source.get("url"), f"source {source_id}.url")
            _string(
                raw_source.get("taxonomic_scope"), f"source {source_id}.taxonomic_scope"
            )
            _string_list(raw_source.get("supports"), f"source {source_id}.supports")
            _string_list(
                raw_source.get("limitations"), f"source {source_id}.limitations"
            )
        else:
            allowed = {
                "id",
                "type",
                "rule_ids",
                "rationale",
                "adjudicator",
                "date",
                "negative_controls",
                "limitations",
            }
            _expect_keys(raw_source, f"source {source_id}", allowed)
            _string_list(raw_source.get("rule_ids"), f"source {source_id}.rule_ids")
            _string(raw_source.get("rationale"), f"source {source_id}.rationale")
            _string(raw_source.get("adjudicator"), f"source {source_id}.adjudicator")
            _string(raw_source.get("date"), f"source {source_id}.date")
            _string_list(
                raw_source.get("negative_controls"),
                f"source {source_id}.negative_controls",
            )
            _string_list(
                raw_source.get("limitations"), f"source {source_id}.limitations"
            )
        sources.append(
            ProvenanceSource(
                id=source_id, source_type=source_type, payload=dict(raw_source)
            )
        )
    return tuple(sources)


def _parse_components(
    inference_data: dict[str, Any], facet_ids: set[str]
) -> tuple[InferenceComponent, ...]:
    raw_components = inference_data.get("components")
    if not isinstance(raw_components, dict) or not raw_components:
        raise MobilomeDatabaseError(
            "inference.yaml.components must be a non-empty mapping"
        )
    components: list[InferenceComponent] = []
    for component_id in sorted(raw_components):
        raw = _expect_mapping(raw_components[component_id], f"component {component_id}")
        _expect_keys(
            raw,
            f"component {component_id}",
            {"operator", "facets", "min_distinct_features"},
        )
        operator = _string(raw.get("operator"), f"component {component_id}.operator")
        if operator not in {"any", "all"}:
            raise MobilomeDatabaseError(
                f"Component {component_id} has unknown operator"
            )
        facets = _string_list(raw.get("facets"), f"component {component_id}.facets")
        missing = sorted(set(facets) - facet_ids)
        if missing:
            raise MobilomeDatabaseError(
                f"Component {component_id} references unknown facets: {', '.join(missing)}"
            )
        raw_minimums = raw.get("min_distinct_features", {})
        if not isinstance(raw_minimums, dict):
            raise MobilomeDatabaseError(
                f"component {component_id}.min_distinct_features must be a mapping"
            )
        unknown_minimums = sorted(set(raw_minimums) - set(facets))
        if unknown_minimums:
            raise MobilomeDatabaseError(
                f"Component {component_id} has minima for unknown facets: {', '.join(unknown_minimums)}"
            )
        minimums = tuple(
            sorted(
                (facet, _positive_int(count, f"component {component_id}.{facet}"))
                for facet, count in raw_minimums.items()
            )
        )
        components.append(
            InferenceComponent(
                id=_string(component_id, "component id"),
                operator=operator,  # type: ignore[arg-type]
                facets=tuple(sorted(facets)),
                min_distinct_features=minimums,
            )
        )
    return tuple(components)


def _parse_inference_rules(
    inference_data: dict[str, Any],
    component_ids: set[str],
    marker_ids: set[str],
) -> tuple[InferenceRule, ...]:
    raw_rules = inference_data.get("rules")
    if not isinstance(raw_rules, list) or not all(
        isinstance(item, dict) for item in raw_rules
    ):
        raise MobilomeDatabaseError("inference.yaml.rules must be a list of mappings")
    _unique_ids(raw_rules, "inference rule")
    rules: list[InferenceRule] = []
    allowed = {
        "id",
        "kind",
        "required_components",
        "required_marker_ids",
        "same_record",
        "distinct_features",
        "max_circular_gap_bp",
        "status",
        "wording",
        "limitations",
        "sources",
    }
    for raw_rule in raw_rules:
        rule_id = _string(raw_rule.get("id"), "inference rule.id")
        _expect_keys(raw_rule, f"inference rule {rule_id}", allowed)
        raw_components = raw_rule.get("required_components", {})
        if not isinstance(raw_components, dict):
            raise MobilomeDatabaseError(
                f"Rule {rule_id}.required_components must be a mapping"
            )
        parsed_components: list[tuple[str, tuple[str, ...]]] = []
        for role, component_list in raw_components.items():
            role_name = _string(role, f"rule {rule_id}.role")
            component_ids_for_role = _string_list(
                component_list, f"rule {rule_id}.{role_name}"
            )
            unknown_components = sorted(set(component_ids_for_role) - component_ids)
            if unknown_components:
                raise MobilomeDatabaseError(
                    f"Rule {rule_id} references unknown components: {', '.join(unknown_components)}"
                )
            parsed_components.append((role_name, component_ids_for_role))
        required_marker_ids = _string_list(
            raw_rule.get("required_marker_ids", []),
            f"rule {rule_id}.required_marker_ids",
            allow_empty=True,
        )
        unknown_markers = sorted(set(required_marker_ids) - marker_ids)
        if unknown_markers:
            raise MobilomeDatabaseError(
                f"Rule {rule_id} references unknown markers: {', '.join(unknown_markers)}"
            )
        if not parsed_components and not required_marker_ids:
            raise MobilomeDatabaseError(f"Rule {rule_id} requires evidence")
        status = _string(raw_rule.get("status"), f"rule {rule_id}.status")
        if status not in _VALID_STATUSES:
            raise MobilomeDatabaseError(
                f"Rule {rule_id} has unsupported status {status!r}"
            )
        raw_gap = raw_rule.get("max_circular_gap_bp")
        max_gap = (
            None
            if raw_gap is None
            else _positive_int(raw_gap, f"rule {rule_id}.max_circular_gap_bp")
        )
        rules.append(
            InferenceRule(
                id=rule_id,
                kind=_string(raw_rule.get("kind"), f"rule {rule_id}.kind"),
                status=status,  # type: ignore[arg-type]
                wording=_string(raw_rule.get("wording"), f"rule {rule_id}.wording"),
                sources=tuple(
                    sorted(
                        _string_list(raw_rule.get("sources"), f"rule {rule_id}.sources")
                    )
                ),
                limitations=_string_list(
                    raw_rule.get("limitations"), f"rule {rule_id}.limitations"
                ),
                required_components=tuple(sorted(parsed_components)),
                required_marker_ids=tuple(sorted(required_marker_ids)),
                same_record=_bool(
                    raw_rule.get("same_record", False), f"rule {rule_id}.same_record"
                ),
                distinct_features=_bool(
                    raw_rule.get("distinct_features", False),
                    f"rule {rule_id}.distinct_features",
                ),
                max_circular_gap_bp=max_gap,
            )
        )
    return tuple(rules)


def _parse_disabled_rules(
    inference_data: dict[str, Any], facet_ids: set[str], rule_ids: set[str]
) -> tuple[DisabledAggregateRule, ...]:
    raw_rules = inference_data.get("disabled_aggregate_rules")
    if not isinstance(raw_rules, list) or not all(
        isinstance(item, dict) for item in raw_rules
    ):
        raise MobilomeDatabaseError(
            "inference.yaml.disabled_aggregate_rules must be a list of mappings"
        )
    _unique_ids(raw_rules, "disabled aggregate rule")
    parsed: list[DisabledAggregateRule] = []
    for raw in raw_rules:
        rule_id = _string(raw.get("id"), "disabled aggregate rule.id")
        _expect_keys(
            raw,
            f"disabled aggregate rule {rule_id}",
            {
                "id",
                "reason",
                "evidence_retained_as",
                "enablement_requirements",
                "sources",
            },
        )
        if rule_id in rule_ids:
            raise MobilomeDatabaseError(
                f"Disabled aggregate rule {rule_id} is also executable"
            )
        retained = _string_list(
            raw.get("evidence_retained_as"),
            f"disabled rule {rule_id}.evidence_retained_as",
        )
        unknown = sorted(
            set(retained) - facet_ids - {"pxo_numbered_product_annotation"}
        )
        if unknown:
            raise MobilomeDatabaseError(
                f"Disabled aggregate rule {rule_id} has unknown retained evidence: {', '.join(unknown)}"
            )
        parsed.append(
            DisabledAggregateRule(
                rule_id=rule_id,
                reason=_string(raw.get("reason"), f"disabled rule {rule_id}.reason"),
                evidence_retained_as=tuple(sorted(retained)),
                enablement_requirements=_string_list(
                    raw.get("enablement_requirements"),
                    f"disabled rule {rule_id}.enablement_requirements",
                ),
                source_ids=tuple(
                    sorted(
                        _string_list(
                            raw.get("sources"), f"disabled rule {rule_id}.sources"
                        )
                    )
                ),
            )
        )
    return tuple(parsed)


def _parse_handoffs(inference_data: dict[str, Any]) -> tuple[ExternalHandoff, ...]:
    raw_handoffs = inference_data.get("handoffs")
    if not isinstance(raw_handoffs, list) or not all(
        isinstance(item, dict) for item in raw_handoffs
    ):
        raise MobilomeDatabaseError(
            "inference.yaml.handoffs must be a list of mappings"
        )
    _unique_ids(raw_handoffs, "handoff")
    handoffs: list[ExternalHandoff] = []
    for raw in raw_handoffs:
        handoff_id = _string(raw.get("id"), "handoff.id")
        _expect_keys(
            raw,
            f"handoff {handoff_id}",
            {
                "id",
                "tool",
                "purpose",
                "reason",
                "required_input",
                "required_provenance_fields",
                "limitations",
            },
        )
        handoffs.append(
            ExternalHandoff(
                handoff_id=handoff_id,
                tool=_string(raw.get("tool"), f"handoff {handoff_id}.tool"),
                purpose=_string(raw.get("purpose"), f"handoff {handoff_id}.purpose"),
                status="not_run",
                reason=_string(raw.get("reason"), f"handoff {handoff_id}.reason"),
                required_input=_string_list(
                    raw.get("required_input"), f"handoff {handoff_id}.required_input"
                ),
                required_provenance_fields=_string_list(
                    raw.get("required_provenance_fields"),
                    f"handoff {handoff_id}.required_provenance_fields",
                ),
                limitations=_string_list(
                    raw.get("limitations"), f"handoff {handoff_id}.limitations"
                ),
            )
        )
    return tuple(handoffs)


def _validate_source_references(
    markers: tuple[MobilomeMarker, ...],
    rules: tuple[InferenceRule, ...],
    disabled: tuple[DisabledAggregateRule, ...],
    sources: tuple[ProvenanceSource, ...],
) -> None:
    known_sources = {source.id for source in sources}
    for owner, source_ids in [
        *((f"marker {marker.id}", marker.sources) for marker in markers),
        *((f"inference rule {rule.id}", rule.sources) for rule in rules),
        *(
            (f"disabled aggregate rule {rule.rule_id}", rule.source_ids)
            for rule in disabled
        ),
    ]:
        unresolved = sorted(set(source_ids) - known_sources)
        if unresolved:
            raise MobilomeDatabaseError(
                f"{owner} has unresolved sources: {', '.join(unresolved)}"
            )
    known_rules = (
        _CLASSIFICATION_RULE_IDS
        | {rule.id for rule in rules}
        | {rule.rule_id for rule in disabled}
    )
    for source in sources:
        if source.source_type != "local_rule":
            continue
        rule_ids = _string_list(
            source.payload.get("rule_ids", []), f"source {source.id}.rule_ids"
        )
        unresolved = sorted(set(rule_ids) - known_rules)
        if unresolved:
            raise MobilomeDatabaseError(
                f"Local source {source.id} references unknown rules: {', '.join(unresolved)}"
            )


def load_mobilome_database(database_dir: str | Path | None = None) -> MobilomeDatabase:
    """Load one complete packaged or custom mobilome evidence database.

    The three YAML files are intentionally selected as a set.  A custom marker
    resource therefore cannot silently inherit packaged provenance or wording.
    """

    custom_dir = Path(database_dir) if database_dir is not None else None
    raw_resources, source_paths = _resource_bytes(custom_dir)
    markers_data = _load_yaml(raw_resources["markers.yaml"], "markers.yaml")
    provenance_data = _load_yaml(raw_resources["provenance.yaml"], "provenance.yaml")
    inference_data = _load_yaml(raw_resources["inference.yaml"], "inference.yaml")
    _load_schema(raw_resources[_REPORT_SCHEMA_NAME])
    _expect_keys(
        inference_data,
        "inference.yaml",
        {
            "schema_version",
            "inference_version",
            "components",
            "rules",
            "disabled_aggregate_rules",
            "handoffs",
        },
    )

    catalog_version = _string(
        markers_data.get("catalog_version"), "markers.yaml.catalog_version"
    )
    provenance_version = _string(
        provenance_data.get("provenance_version"), "provenance.yaml.provenance_version"
    )
    inference_version = _string(
        inference_data.get("inference_version"), "inference.yaml.inference_version"
    )
    majors = {
        _major_version(catalog_version, "markers.yaml.catalog_version"),
        _major_version(provenance_version, "provenance.yaml.provenance_version"),
        _major_version(inference_version, "inference.yaml.inference_version"),
    }
    if len(majors) != 1:
        raise MobilomeDatabaseError(
            "Catalog, provenance, and inference major versions must agree"
        )

    facets = _parse_facets(markers_data)
    facet_ids = {facet.id for facet in facets}
    markers = _parse_markers(markers_data, facet_ids)
    sources = _parse_provenance(provenance_data)
    components = _parse_components(inference_data, facet_ids)
    rules = _parse_inference_rules(
        inference_data,
        {component.id for component in components},
        {marker.id for marker in markers},
    )
    disabled = _parse_disabled_rules(
        inference_data, facet_ids, {rule.id for rule in rules}
    )
    handoffs = _parse_handoffs(inference_data)
    _validate_source_references(markers, rules, disabled, sources)

    catalog_resources = tuple(
        CatalogResource(
            name=name, sha256=hashlib.sha256(raw_resources[name]).hexdigest()
        )
        for name in (*_DATABASE_RESOURCE_NAMES, _REPORT_SCHEMA_NAME)
    )
    return MobilomeDatabase(
        facets=tuple(sorted(facets, key=lambda item: item.id)),
        markers=tuple(sorted(markers, key=lambda item: item.id)),
        provenance_sources=tuple(sorted(sources, key=lambda item: item.id)),
        components=tuple(sorted(components, key=lambda item: item.id)),
        inference_rules=tuple(sorted(rules, key=lambda item: item.id)),
        disabled_aggregate_rules=tuple(sorted(disabled, key=lambda item: item.rule_id)),
        handoffs=tuple(sorted(handoffs, key=lambda item: item.handoff_id)),
        catalog_version=catalog_version,
        inference_version=inference_version,
        provenance_version=provenance_version,
        resources=catalog_resources,
        source_paths=source_paths,
        database_source="custom" if custom_dir is not None else "packaged",
    )


__all__ = ["MobilomeDatabaseError", "load_mobilome_database"]
