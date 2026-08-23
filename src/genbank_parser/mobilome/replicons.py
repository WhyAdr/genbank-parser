"""Complete record inventory and conservative declared-replicon classification."""

from __future__ import annotations

import re
from collections import defaultdict

from ..model import GenBankDocument, GenBankRecord
from .models import (
    RecordEvidence,
    RepliconClass,
    RepliconClassification,
    RepliconInventory,
    SourceQualifier,
)

_RECORD_TOKEN_PATTERNS = {
    "plasmid": re.compile(r"(?<![A-Za-z0-9_-])plasmid(?![A-Za-z0-9_-])", re.IGNORECASE),
    "chromosome": re.compile(
        r"(?<![A-Za-z0-9_-])chromosome(?![A-Za-z0-9_-])", re.IGNORECASE
    ),
}


def _source_qualifiers(record: GenBankRecord) -> tuple[SourceQualifier, ...]:
    values_by_key: dict[str, list[str]] = defaultdict(list)
    for feature in record.features:
        if feature.type.casefold() != "source":
            continue
        for key, values in feature.qualifiers.items():
            values_by_key[key].extend(str(value) for value in values)
    return tuple(
        SourceQualifier(key=key, values=tuple(sorted(set(values))))
        for key, values in sorted(values_by_key.items())
    )


def collect_record_classification_evidence(
    record: GenBankRecord,
) -> tuple[RecordEvidence, ...]:
    """Collect declarations without using topology, length, or marker content."""

    evidence: list[RecordEvidence] = []
    for feature in record.features:
        if feature.type.casefold() != "source":
            continue
        for key, label, rule_id in (
            ("plasmid", "plasmid", "source_plasmid_qualifier"),
            ("chromosome", "chromosome", "source_chromosome_qualifier"),
        ):
            if key not in feature.qualifiers:
                continue
            values = feature.qualifiers[key] or [""]
            for value in values:
                text = str(value)
                evidence.append(
                    RecordEvidence(
                        rule_id=rule_id,
                        source="source_qualifier",
                        field=key,
                        value=text,
                        matched_text=text,
                        pattern="qualifier-present",
                        source_feature_index=feature.feature_index,
                        interpreted_as=label,  # type: ignore[arg-type]
                        strength=3,
                        source_ids=("insdc-feature-table",),
                    )
                )
    for key, label, rule_id in (
        ("plasmid", "plasmid", "record_annotation_plasmid"),
        ("chromosome", "chromosome", "record_annotation_chromosome"),
    ):
        if key not in record.annotations:
            continue
        raw_values = record.annotations[key]
        values = raw_values if isinstance(raw_values, list) else [raw_values]
        for value in values or [""]:
            text = str(value)
            evidence.append(
                RecordEvidence(
                    rule_id=rule_id,
                    source="record_annotation",
                    field=key,
                    value=text,
                    matched_text=text,
                    pattern="annotation-key-present",
                    source_feature_index=None,
                    interpreted_as=label,  # type: ignore[arg-type]
                    strength=3,
                    source_ids=("insdc-feature-table",),
                )
            )
    for source, field, value in (
        ("record_name", "name", record.name),
        ("record_description", "description", record.description),
    ):
        for label, pattern in _RECORD_TOKEN_PATTERNS.items():
            match = pattern.search(value or "")
            if match is None:
                continue
            rule_id = f"bounded_{field}_{label}_token"
            evidence.append(
                RecordEvidence(
                    rule_id=rule_id,
                    source=source,  # type: ignore[arg-type]
                    field=field,
                    value=value or "",
                    matched_text=match.group(0),
                    pattern=pattern.pattern,
                    source_feature_index=None,
                    interpreted_as=label,  # type: ignore[arg-type]
                    strength=1,
                    source_ids=("local-record-classification-policy",),
                )
            )
    return tuple(
        sorted(
            evidence,
            key=lambda item: (
                item.interpreted_as,
                item.source,
                item.field,
                item.value,
                item.source_feature_index or 0,
                item.rule_id,
            ),
        )
    )


def classify_replicon(record: GenBankRecord) -> RepliconClassification:
    """Describe the input's record declarations; do not predict a replicon."""

    evidence = collect_record_classification_evidence(record)
    plasmid = tuple(item for item in evidence if item.interpreted_as == "plasmid")
    chromosome = tuple(item for item in evidence if item.interpreted_as == "chromosome")
    if plasmid and chromosome:
        return RepliconClassification(
            label="unknown",
            status="conflicting",
            supporting_evidence=(),
            conflicting_evidence=plasmid + chromosome,
            limitations=("Conflicting record-level declarations were retained.",),
        )
    if plasmid or chromosome:
        label: RepliconClass = "plasmid" if plasmid else "chromosome"
        reasons = plasmid or chromosome
        supported = any(
            item.source in {"source_qualifier", "record_annotation"} for item in reasons
        )
        return RepliconClassification(
            label=label,
            status="supported" if supported else "tentative",
            supporting_evidence=reasons,
            conflicting_evidence=(),
            limitations=()
            if supported
            else (
                "Classification is based only on a bounded name or description token.",
            ),
        )
    return RepliconClassification(
        label="unknown",
        status="insufficient",
        supporting_evidence=(),
        conflicting_evidence=(),
        limitations=(
            "No chromosome or plasmid declaration was observed in record metadata.",
        ),
    )


def _topology(value: str | None) -> str:
    if value is None:
        return "unknown"
    folded = value.casefold()
    return folded if folded in {"circular", "linear"} else "unknown"


def inventory_replicons(document: GenBankDocument) -> tuple[RepliconInventory, ...]:
    """Inventory every parsed record, including zero-feature and unresolved records."""

    inventory: list[RepliconInventory] = []
    for record_index, record in enumerate(document.records, 1):
        features = tuple(record.features)
        pseudo_features = tuple(feature for feature in features if feature.is_pseudo)
        cds = tuple(feature for feature in features if feature.type.casefold() == "cds")
        inventory.append(
            RepliconInventory(
                record_index=record_index,
                record_id=record.id,
                name=record.name,
                description=record.description,
                length=record.length,
                topology=_topology(record.topology),  # type: ignore[arg-type]
                source_qualifiers=_source_qualifiers(record),
                classification=classify_replicon(record),
                feature_count=len(features),
                cds_count=len(cds),
                pseudo_feature_count=len(pseudo_features),
                pseudo_cds_count=sum(
                    feature.type.casefold() == "cds" for feature in pseudo_features
                ),
                pseudogene_feature_count=sum(
                    feature.type.casefold() == "pseudogene" for feature in features
                ),
            )
        )
    return tuple(sorted(inventory, key=lambda item: item.record_index))


__all__ = [
    "classify_replicon",
    "collect_record_classification_evidence",
    "inventory_replicons",
]
