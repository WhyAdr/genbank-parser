"""Deterministic JSON, TSV, text, and safe output writing for mobilome reports."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

from .models import (
    CatalogResource,
    DisabledAggregateRule,
    EvidenceReason,
    ExternalHandoff,
    FeatureRef,
    MobilomeHit,
    MobilomeHypothesis,
    MobilomeOutputError,
    MobilomeParameterError,
    MobilomeReport,
    RecordEvidence,
    RepliconAssessment,
    RepliconInventory,
)

TSV_COLUMNS = (
    "schema_version",
    "row_type",
    "record_index",
    "record_id",
    "record_name",
    "record_description",
    "record_length",
    "topology",
    "replicon_class",
    "classification_status",
    "feature_count",
    "cds_count",
    "pseudo_feature_count",
    "pseudo_cds_count",
    "pseudogene_feature_count",
    "feature_index",
    "feature_type",
    "locus_tag",
    "gene",
    "product",
    "start",
    "end",
    "strand",
    "segments_json",
    "wraps_origin",
    "partial",
    "pseudo",
    "hit_id",
    "marker_id",
    "facets_json",
    "reason_id",
    "evidence_field",
    "match_mode",
    "matched_value",
    "matched_text",
    "pattern",
    "strength",
    "eligible",
    "source_ids_json",
    "hypothesis_id",
    "hypothesis_rule_id",
    "hypothesis_kind",
    "hypothesis_status",
    "participants_json",
    "supporting_hit_ids_json",
    "missing_components_json",
    "conflicting_hit_ids_json",
    "summary",
    "limitations_json",
    "disabled_rule_id",
    "disabled_rule_reason",
    "disabled_evidence_retained_as_json",
    "disabled_enablement_requirements_json",
    "handoff_id",
    "handoff_tool",
    "handoff_status",
    "handoff_reason",
    "handoff_provenance_json",
)


def _list(value: Iterable[object]) -> list[object]:
    return list(value)


def _resource_dict(resource: CatalogResource) -> dict[str, object]:
    return {"name": resource.name, "sha256": resource.sha256}


def _record_evidence_dict(evidence: RecordEvidence) -> dict[str, object]:
    return {
        "rule_id": evidence.rule_id,
        "source": evidence.source,
        "field": evidence.field,
        "value": evidence.value,
        "matched_text": evidence.matched_text,
        "pattern": evidence.pattern,
        "source_feature_index": evidence.source_feature_index,
        "interpreted_as": evidence.interpreted_as,
        "strength": evidence.strength,
        "source_ids": _list(evidence.source_ids),
    }


def _inventory_dict(inventory: RepliconInventory) -> dict[str, object]:
    classification = inventory.classification
    return {
        "record_index": inventory.record_index,
        "record_id": inventory.record_id,
        "name": inventory.name,
        "description": inventory.description,
        "length": inventory.length,
        "topology": inventory.topology,
        "source_qualifiers": [
            {"key": qualifier.key, "values": _list(qualifier.values)}
            for qualifier in inventory.source_qualifiers
        ],
        "classification": {
            "label": classification.label,
            "status": classification.status,
            "supporting_evidence": [
                _record_evidence_dict(item)
                for item in classification.supporting_evidence
            ],
            "conflicting_evidence": [
                _record_evidence_dict(item)
                for item in classification.conflicting_evidence
            ],
            "limitations": _list(classification.limitations),
        },
        "feature_count": inventory.feature_count,
        "cds_count": inventory.cds_count,
        "pseudo_feature_count": inventory.pseudo_feature_count,
        "pseudo_cds_count": inventory.pseudo_cds_count,
        "pseudogene_feature_count": inventory.pseudogene_feature_count,
    }


def _feature_dict(feature: FeatureRef) -> dict[str, object]:
    return {
        "record_index": feature.record_index,
        "record_id": feature.record_id,
        "feature_index": feature.feature_index,
        "feature_type": feature.feature_type,
        "locus_tag": feature.locus_tag,
        "gene": feature.gene,
        "product": feature.product,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand,
        "segments": [list(segment) for segment in feature.segments],
        "location_operator": feature.location_operator,
        "biological_length": feature.biological_length,
        "wraps_origin": feature.wraps_origin,
        "partial": feature.is_partial,
        "pseudo": feature.is_pseudo,
    }


def _reason_dict(reason: EvidenceReason) -> dict[str, object]:
    return {
        "reason_id": reason.reason_id,
        "field": reason.field,
        "match_mode": reason.match_mode,
        "matched_value": reason.matched_value,
        "matched_text": reason.matched_text,
        "pattern": reason.pattern,
        "strength": reason.strength,
        "eligible": reason.eligible,
        "source_ids": _list(reason.source_ids),
    }


def _hit_dict(hit: MobilomeHit) -> dict[str, object]:
    return {
        "hit_id": hit.hit_id,
        "marker_id": hit.marker_id,
        "marker_label": hit.marker_label,
        "facets": _list(hit.facets),
        "feature": _feature_dict(hit.feature),
        "reasons": [_reason_dict(reason) for reason in hit.reasons],
        "max_strength": hit.max_strength,
        "eligible_for_inference": hit.eligible_for_inference,
    }


def _hypothesis_dict(hypothesis: MobilomeHypothesis) -> dict[str, object]:
    return {
        "hypothesis_id": hypothesis.hypothesis_id,
        "rule_id": hypothesis.rule_id,
        "kind": hypothesis.kind,
        "status": hypothesis.status,
        "summary": hypothesis.summary,
        "participants": [
            {
                "role": participant.role,
                "record_index": participant.record_index,
                "record_id": participant.record_id,
            }
            for participant in hypothesis.participants
        ],
        "supporting_hit_ids": _list(hypothesis.supporting_hit_ids),
        "missing_components": _list(hypothesis.missing_components),
        "conflicting_hit_ids": _list(hypothesis.conflicting_hit_ids),
        "limitations": _list(hypothesis.limitations),
        "source_ids": _list(hypothesis.source_ids),
    }


def _disabled_rule_dict(rule: DisabledAggregateRule) -> dict[str, object]:
    return {
        "rule_id": rule.rule_id,
        "reason": rule.reason,
        "evidence_retained_as": _list(rule.evidence_retained_as),
        "enablement_requirements": _list(rule.enablement_requirements),
        "source_ids": _list(rule.source_ids),
    }


def _handoff_dict(handoff: ExternalHandoff) -> dict[str, object]:
    return {
        "handoff_id": handoff.handoff_id,
        "tool": handoff.tool,
        "purpose": handoff.purpose,
        "status": handoff.status,
        "reason": handoff.reason,
        "required_input": _list(handoff.required_input),
        "required_provenance_fields": _list(handoff.required_provenance_fields),
        "limitations": _list(handoff.limitations),
    }


def _scanned_replicons(report: MobilomeReport) -> tuple[RepliconAssessment, ...]:
    return report.scanned_replicons or report.replicons


def _summary_dict(report: MobilomeReport) -> dict[str, int]:
    scanned = _scanned_replicons(report)
    hits = tuple(hit for assessment in scanned for hit in assessment.hits)
    hypotheses = (
        tuple(
            hypothesis for assessment in scanned for hypothesis in assessment.hypotheses
        )
        + report.cross_record_hypotheses
    )
    reasons = tuple(reason for hit in hits for reason in hit.reasons)
    return {
        "records_scanned": len(report.inventory),
        "records_reported": len(report.replicons),
        "raw_hits": len(hits),
        "eligible_hits": sum(hit.eligible_for_inference for hit in hits),
        "raw_reasons": len(reasons),
        "eligible_reasons": sum(reason.eligible for reason in reasons),
        "hypotheses": len(hypotheses),
    }


def validate_mobilome_report_semantics(report: MobilomeReport) -> None:
    """Check constraints that are clearer in Python than generic JSON Schema."""

    if len(report.source_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in report.source_sha256
    ):
        raise MobilomeOutputError(
            "Mobilome report source_sha256 must be a lowercase SHA-256"
        )
    inventory_indices = [inventory.record_index for inventory in report.inventory]
    if len(inventory_indices) != len(set(inventory_indices)):
        raise MobilomeOutputError(
            "Mobilome report inventory has duplicate record_index values"
        )
    inventory_by_index = {
        inventory.record_index: inventory for inventory in report.inventory
    }
    detailed_indices = [
        assessment.inventory.record_index for assessment in report.replicons
    ]
    if len(detailed_indices) != len(set(detailed_indices)):
        raise MobilomeOutputError(
            "Mobilome report replicons have duplicate record_index values"
        )
    if not set(detailed_indices).issubset(inventory_by_index):
        raise MobilomeOutputError(
            "Mobilome report detail references unknown inventory records"
        )

    all_assessments = _scanned_replicons(report)
    hit_ids: set[str] = set()
    for assessment in all_assessments:
        if assessment.inventory.record_index not in inventory_by_index:
            raise MobilomeOutputError(
                "Scanned assessment references unknown inventory record"
            )
        for hit in assessment.hits:
            if hit.feature.record_index != assessment.inventory.record_index:
                raise MobilomeOutputError(
                    "Mobilome hit belongs to a different inventory record"
                )
            if hit.hit_id in hit_ids:
                raise MobilomeOutputError("Mobilome report has duplicate hit IDs")
            hit_ids.add(hit.hit_id)
            reason_ids = [reason.reason_id for reason in hit.reasons]
            if not reason_ids or len(reason_ids) != len(set(reason_ids)):
                raise MobilomeOutputError(
                    "Mobilome hit reasons must have unique non-empty IDs"
                )
            if hit.feature.strand not in (-1, 0, 1, None):
                raise MobilomeOutputError(
                    "Mobilome feature strand must be -1, 0, 1, or null"
                )
            if not hit.feature.segments and any(
                value is not None
                for value in (
                    hit.feature.start,
                    hit.feature.end,
                    hit.feature.biological_length,
                )
            ):
                raise MobilomeOutputError(
                    "Unlocatable feature references must use null coordinates"
                )
            for start, end in hit.feature.segments:
                if start < 1 or end < start:
                    raise MobilomeOutputError(
                        "Feature segments must be one-based inclusive"
                    )
    hypotheses = (
        tuple(
            hypothesis
            for assessment in all_assessments
            for hypothesis in assessment.hypotheses
        )
        + report.cross_record_hypotheses
    )
    hypothesis_ids = [hypothesis.hypothesis_id for hypothesis in hypotheses]
    if len(hypothesis_ids) != len(set(hypothesis_ids)):
        raise MobilomeOutputError("Mobilome report has duplicate hypothesis IDs")
    for hypothesis in hypotheses:
        unresolved = set(hypothesis.supporting_hit_ids) - hit_ids
        if unresolved:
            raise MobilomeOutputError(
                f"Hypothesis {hypothesis.hypothesis_id} references unknown hits: {', '.join(sorted(unresolved))}"
            )
        roles = [participant.role for participant in hypothesis.participants]
        if len(roles) != len(set(roles)):
            raise MobilomeOutputError("Hypothesis participant roles must be unique")
        if {"target", "helper"}.issubset(roles):
            if len(hypothesis.participants) < 2:
                raise MobilomeOutputError(
                    "Cross-record hypothesis needs target and helper participants"
                )
            target = next(
                item for item in hypothesis.participants if item.role == "target"
            )
            helper = next(
                item for item in hypothesis.participants if item.role == "helper"
            )
            if target.record_index == helper.record_index:
                raise MobilomeOutputError(
                    "Cross-record target and helper must be distinct records"
                )
        elif roles != ["subject"]:
            raise MobilomeOutputError(
                "Record-local hypotheses require one subject participant"
            )
    for rule in report.disabled_aggregate_rules:
        if not rule.enablement_requirements:
            raise MobilomeOutputError(
                "Disabled aggregate rules need enablement requirements"
            )


def mobilome_report_to_dict(report: MobilomeReport) -> dict[str, object]:
    """Convert a validated report to the schema-versioned canonical mapping."""

    validate_mobilome_report_semantics(report)
    return {
        "schema_version": "gbparse.mobilome.v1",
        "tool": "gbparse",
        "analysis": "mobilome",
        "tool_version": report.tool_version,
        "catalog_version": report.catalog_version,
        "inference_version": report.inference_version,
        "input": {"file": report.source_file, "sha256": report.source_sha256},
        "parameters": {
            "include": report.parameters.include,
            "min_evidence": report.parameters.min_evidence,
            "database_source": report.parameters.database_source,
        },
        "provenance": {
            "resources": [_resource_dict(resource) for resource in report.resources],
            "annotation_pipelines": [
                {
                    "record_indices": _list(item.record_indices),
                    "name": item.name,
                    "version": item.version,
                    "database_version": item.database_version,
                    "evidence_fields": _list(item.evidence_fields),
                }
                for item in report.annotation_provenance
            ],
        },
        "summary": _summary_dict(report),
        "inventory": [_inventory_dict(inventory) for inventory in report.inventory],
        "replicons": [
            {
                "inventory": _inventory_dict(assessment.inventory),
                "hits": [_hit_dict(hit) for hit in assessment.hits],
                "hypotheses": [
                    _hypothesis_dict(hypothesis) for hypothesis in assessment.hypotheses
                ],
            }
            for assessment in report.replicons
        ],
        "cross_record_hypotheses": [
            _hypothesis_dict(hypothesis)
            for hypothesis in report.cross_record_hypotheses
        ],
        "disabled_aggregate_rules": [
            _disabled_rule_dict(rule) for rule in report.disabled_aggregate_rules
        ],
        "handoffs": [_handoff_dict(handoff) for handoff in report.handoffs],
        "limitations": _list(report.limitations),
    }


def render_json(report: MobilomeReport) -> str:
    return (
        json.dumps(mobilome_report_to_dict(report), ensure_ascii=False, indent=2) + "\n"
    )


def _compact_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _blank_row(row_type: str) -> dict[str, object]:
    row: dict[str, object] = {column: "" for column in TSV_COLUMNS}
    row["schema_version"] = "gbparse.mobilome.v1"
    row["row_type"] = row_type
    return row


def _inventory_row(
    inventory: RepliconInventory, row_type: str = "inventory"
) -> dict[str, object]:
    row = _blank_row(row_type)
    row.update(
        {
            "record_index": inventory.record_index,
            "record_id": inventory.record_id,
            "record_name": inventory.name,
            "record_description": inventory.description,
            "record_length": inventory.length,
            "topology": inventory.topology,
            "replicon_class": inventory.classification.label,
            "classification_status": inventory.classification.status,
            "feature_count": inventory.feature_count,
            "cds_count": inventory.cds_count,
            "pseudo_feature_count": inventory.pseudo_feature_count,
            "pseudo_cds_count": inventory.pseudo_cds_count,
            "pseudogene_feature_count": inventory.pseudogene_feature_count,
        }
    )
    return row


def _feature_row_values(feature: FeatureRef) -> dict[str, object]:
    return {
        "record_index": feature.record_index,
        "record_id": feature.record_id,
        "feature_index": feature.feature_index,
        "feature_type": feature.feature_type,
        "locus_tag": feature.locus_tag or "",
        "gene": feature.gene or "",
        "product": feature.product or "",
        "start": feature.start if feature.start is not None else "",
        "end": feature.end if feature.end is not None else "",
        "strand": feature.strand if feature.strand is not None else "",
        "segments_json": _compact_json([list(item) for item in feature.segments]),
        "wraps_origin": str(feature.wraps_origin).lower(),
        "partial": str(feature.is_partial).lower(),
        "pseudo": str(feature.is_pseudo).lower(),
    }


def _reason_row(
    inventory: RepliconInventory, hit: MobilomeHit, reason: EvidenceReason
) -> dict[str, object]:
    row = _inventory_row(inventory, "evidence")
    row.update(_feature_row_values(hit.feature))
    row.update(
        {
            "hit_id": hit.hit_id,
            "marker_id": hit.marker_id,
            "facets_json": _compact_json(hit.facets),
            "reason_id": reason.reason_id,
            "evidence_field": reason.field,
            "match_mode": reason.match_mode,
            "matched_value": reason.matched_value,
            "matched_text": reason.matched_text,
            "pattern": reason.pattern,
            "strength": reason.strength,
            "eligible": str(reason.eligible).lower(),
            "source_ids_json": _compact_json(reason.source_ids),
        }
    )
    return row


def _hypothesis_row(
    hypothesis: MobilomeHypothesis,
    *,
    row_type: str,
    inventory: RepliconInventory | None = None,
) -> dict[str, object]:
    row = (
        _inventory_row(inventory, row_type)
        if inventory is not None
        else _blank_row(row_type)
    )
    row.update(
        {
            "hypothesis_id": hypothesis.hypothesis_id,
            "hypothesis_rule_id": hypothesis.rule_id,
            "hypothesis_kind": hypothesis.kind,
            "hypothesis_status": hypothesis.status,
            "participants_json": _compact_json(
                [
                    {
                        "role": participant.role,
                        "record_index": participant.record_index,
                        "record_id": participant.record_id,
                    }
                    for participant in hypothesis.participants
                ]
            ),
            "supporting_hit_ids_json": _compact_json(hypothesis.supporting_hit_ids),
            "missing_components_json": _compact_json(hypothesis.missing_components),
            "conflicting_hit_ids_json": _compact_json(hypothesis.conflicting_hit_ids),
            "summary": hypothesis.summary,
            "limitations_json": _compact_json(hypothesis.limitations),
            "source_ids_json": _compact_json(hypothesis.source_ids),
        }
    )
    return row


def render_tsv(report: MobilomeReport) -> str:
    validate_mobilome_report_semantics(report)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for inventory in report.inventory:
        writer.writerow(_inventory_row(inventory))
    for assessment in report.replicons:
        for hit in assessment.hits:
            for reason in hit.reasons:
                writer.writerow(_reason_row(assessment.inventory, hit, reason))
        for hypothesis in assessment.hypotheses:
            writer.writerow(
                _hypothesis_row(
                    hypothesis, row_type="hypothesis", inventory=assessment.inventory
                )
            )
    for hypothesis in report.cross_record_hypotheses:
        writer.writerow(_hypothesis_row(hypothesis, row_type="cross_record_hypothesis"))
    for rule in report.disabled_aggregate_rules:
        row = _blank_row("disabled_aggregate_rule")
        row.update(
            {
                "disabled_rule_id": rule.rule_id,
                "disabled_rule_reason": rule.reason,
                "disabled_evidence_retained_as_json": _compact_json(
                    rule.evidence_retained_as
                ),
                "disabled_enablement_requirements_json": _compact_json(
                    rule.enablement_requirements
                ),
                "source_ids_json": _compact_json(rule.source_ids),
            }
        )
        writer.writerow(row)
    for handoff in report.handoffs:
        row = _blank_row("handoff")
        row.update(
            {
                "handoff_id": handoff.handoff_id,
                "handoff_tool": handoff.tool,
                "handoff_status": handoff.status,
                "handoff_reason": handoff.reason,
                "handoff_provenance_json": _compact_json(
                    {
                        "required_input": handoff.required_input,
                        "required_provenance_fields": handoff.required_provenance_fields,
                        "limitations": handoff.limitations,
                    }
                ),
            }
        )
        writer.writerow(row)
    return output.getvalue()


def render_text(report: MobilomeReport) -> str:
    """Render a deterministic, cautious human summary without a score verdict."""

    validate_mobilome_report_semantics(report)
    output = io.StringIO(newline="\n")
    print("Mobilome evidence report", file=output)
    print("Input and provenance", file=output)
    print(f"  Input: {report.source_file}", file=output)
    print(f"  SHA-256: {report.source_sha256}", file=output)
    print(
        f"  Catalog/inference: {report.catalog_version} / {report.inference_version}",
        file=output,
    )
    print("Record inventory", file=output)
    for inventory in report.inventory:
        print(
            f"  [{inventory.record_index}] {inventory.record_id} | "
            f"{inventory.classification.label} ({inventory.classification.status}) | "
            f"{inventory.topology} | {inventory.length} bp | {inventory.cds_count} CDS",
            file=output,
        )
        for limitation in inventory.classification.limitations:
            print(f"    Limitation: {limitation}", file=output)
    print("Replicon assessments", file=output)
    for assessment in report.replicons:
        inventory = assessment.inventory
        print(f"  [{inventory.record_index}] {inventory.record_id}", file=output)
        eligible = [hit for hit in assessment.hits if hit.eligible_for_inference]
        below = [hit for hit in assessment.hits if not hit.eligible_for_inference]
        print("    Observed eligible evidence", file=output)
        if eligible:
            for hit in eligible:
                print(
                    f"      {hit.marker_label} ({hit.marker_id}) on feature "
                    f"{hit.feature.feature_index}",
                    file=output,
                )
        else:
            print("      None at the configured threshold.", file=output)
        print("    Tentative hypotheses", file=output)
        tentative = [
            item for item in assessment.hypotheses if item.status == "tentative"
        ]
        if tentative:
            for item in tentative:
                print(f"      {item.summary} [{item.rule_id}]", file=output)
        else:
            print("      None.", file=output)
        insufficient = [
            item for item in assessment.hypotheses if item.status != "tentative"
        ]
        if insufficient:
            print("    Conflicts and missing components", file=output)
            for item in insufficient:
                missing = ", ".join(item.missing_components) or "none configured"
                print(f"      {item.summary}; missing: {missing}", file=output)
        if below:
            print("    Below-threshold observations", file=output)
            for hit in below:
                print(f"      {hit.marker_label} ({hit.marker_id})", file=output)
    print("Cross-record hypotheses", file=output)
    if report.cross_record_hypotheses:
        for item in report.cross_record_hypotheses:
            print(f"  {item.summary} [{item.rule_id}]", file=output)
            for limitation in item.limitations:
                print(f"    Limitation: {limitation}", file=output)
    else:
        print("  None.", file=output)
    print("Disabled aggregate rules and retained evidence", file=output)
    for rule in report.disabled_aggregate_rules:
        print(f"  {rule.rule_id}: {rule.reason}", file=output)
    print("External handoffs not run", file=output)
    for handoff in report.handoffs:
        print(f"  {handoff.tool}: {handoff.purpose} ({handoff.status})", file=output)
    print("Interpretation limits", file=output)
    for limitation in report.limitations:
        print(f"  - {limitation}", file=output)
    return output.getvalue()


def serialize_mobilome_report(report: MobilomeReport, format_type: str) -> str:
    """Serialize a report after validating the closed output-format set."""

    _validate_format_type(format_type)
    if format_type == "json":
        return render_json(report)
    if format_type == "tsv":
        return render_tsv(report)
    return render_text(report)


def _validate_format_type(format_type: str) -> None:
    if format_type not in {"text", "json", "tsv"}:
        raise MobilomeParameterError("format_type must be 'text', 'json', or 'tsv'")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise MobilomeOutputError(
            f"Could not hash protected input {path}: {exc}"
        ) from exc
    return digest.hexdigest()


def _same_path(left: Path, right: Path) -> bool:
    left_resolved = left.expanduser().resolve(strict=False)
    right_resolved = right.expanduser().resolve(strict=False)
    if left_resolved == right_resolved:
        return True
    try:
        return left.exists() and right.exists() and os.path.samefile(left, right)
    except OSError:
        return False


def _assert_not_protected(output_path: Path, protected_paths: Sequence[Path]) -> None:
    for protected in protected_paths:
        if _same_path(output_path, protected):
            raise MobilomeOutputError(
                f"Output path must not overwrite protected input or database resource: {protected}"
            )


def _atomic_write_text(path: Path, text: str) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise MobilomeOutputError(f"Could not create output directory: {exc}") from exc
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except OSError as exc:
        raise MobilomeOutputError(
            f"Could not write mobilome report safely: {exc}"
        ) from exc
    finally:
        if temporary_name is not None:
            try:
                if os.path.exists(temporary_name):
                    os.unlink(temporary_name)
            except OSError:
                pass


def write_mobilome_report(
    report: MobilomeReport,
    *,
    source_path: str | Path,
    database_paths: Sequence[str | Path],
    format_type: str,
    output_path: str | Path | None = None,
) -> str:
    """Render or atomically write a report while protecting analyzed resources."""

    _validate_format_type(format_type)
    source = Path(source_path)
    if _sha256(source) != report.source_sha256:
        raise MobilomeOutputError(
            "Input changed after analysis; refusing to write report"
        )
    rendered = serialize_mobilome_report(report, format_type)
    if output_path is None:
        return rendered
    output = Path(output_path)
    protected = (source, *(Path(path) for path in database_paths))
    _assert_not_protected(output, protected)
    _atomic_write_text(output, rendered)
    return ""


__all__ = [
    "TSV_COLUMNS",
    "mobilome_report_to_dict",
    "render_json",
    "render_text",
    "render_tsv",
    "serialize_mobilome_report",
    "validate_mobilome_report_semantics",
    "write_mobilome_report",
]
