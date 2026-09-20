"""Canonical feature projections and deterministic tabular serializers."""

from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass

from .io import extract_xref_sources, extract_xrefs
from .model import GenBankFeature, GenBankRecord

FEATURE_SCHEMA_VERSION = "gbparse.feature.v1"

# These are the long-standing ``extract`` columns.  The unified exporter
# reuses them so downstream annotation tables do not need a migration merely
# because the command name changed.
ANNOTATION_COLUMNS = (
    "contig",
    "start",
    "end",
    "strand",
    "type",
    "locus_tag",
    "gene",
    "product",
    "ec_number",
    "cog",
    "kegg_ko",
    "pfam",
    "rfam",
    "db_xref",
    "protein_id",
)


def _json_value(value: object) -> object:
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    return value


@dataclass(frozen=True)
class FeatureRow:
    """A stable, source-aware projection of one physical feature."""

    source: str
    record: str
    record_index: int
    feature_index: int
    type: str
    locus_tag: str | None
    gene: str | None
    product: str | None
    protein_id: str | None
    start: int
    end: int
    strand: int | None
    strand_symbol: str
    length: int
    record_length: int
    topology: str | None
    partial: bool
    pseudo: bool
    segments: tuple[tuple[int, int], ...]
    qualifiers: dict[str, tuple[str, ...]]
    xrefs: dict[str, tuple[str, ...]]
    xref_sources: dict[str, tuple[tuple[str, str], ...]]
    sample_key: str | None = None

    @classmethod
    def from_feature(
        cls,
        feature: GenBankFeature,
        *,
        source: str = "",
    ) -> FeatureRow:
        segments = tuple(feature.join_segments) if feature.is_compound else (
            (feature.start, feature.end),
        )
        xrefs = extract_xrefs(feature)
        xref_sources = extract_xref_sources(feature)
        source_map = {
            category: tuple(
                (item["value"], item["source_field"])
                for item in items
            )
            for category, items in xref_sources.items()
        }
        return cls(
            source=source,
            record=feature.record_id,
            record_index=feature.record_index,
            feature_index=feature.feature_index,
            type=feature.type,
            locus_tag=feature.locus_tag or None,
            gene=feature.gene or None,
            product=feature.product or None,
            protein_id=feature.protein_id or None,
            start=feature.start,
            end=feature.end,
            strand=feature.strand,
            strand_symbol=feature.strand_symbol,
            length=feature.length,
            record_length=feature.record_length,
            topology=feature.topology,
            partial=feature.is_partial,
            pseudo=feature.is_pseudo,
            segments=segments,
            qualifiers={
                key: tuple(str(value) for value in values)
                for key, values in sorted(feature.qualifiers.items())
            },
            xrefs={key: tuple(values) for key, values in sorted(xrefs.items())},
            xref_sources={key: value for key, value in sorted(source_map.items())},
        )

    def value(self, field: str) -> object:
        """Return a queryable field without exposing arbitrary attributes."""

        if field.startswith("qualifier(") and field.endswith(")"):
            key = field[len("qualifier(") : -1].strip().strip("\"'")
            return self.qualifier(key)

        aliases: dict[str, object] = {
            "record": self.record,
            "record_id": self.record,
            "record_index": self.record_index,
            "feature_index": self.feature_index,
            "type": self.type,
            "locus_tag": self.locus_tag,
            "gene": self.gene,
            "product": self.product,
            "protein_id": self.protein_id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand_symbol,
            "strand_value": self.strand,
            "length": self.length,
            "record_length": self.record_length,
            "topology": self.topology,
            "pseudo": self.pseudo,
            "partial": self.partial,
            "ko": self.xrefs.get("kegg_kos", ()),
            "kegg_ko": self.xrefs.get("kegg_kos", ()),
            "ec": self.xrefs.get("ec_numbers", ()),
            "ec_number": self.xrefs.get("ec_numbers", ()),
            "cog": self.xrefs.get("cog_ids", ()),
            "pfam": self.xrefs.get("pfam", ()),
            "rfam": self.xrefs.get("rfam", ()),
            "go": self.xrefs.get("go_terms", ()),
            "go_terms": self.xrefs.get("go_terms", ()),
            "db_xref": self.qualifiers.get("db_xref", ()),
            "source": self.source,
            "sample": self.sample_key,
            "sample_key": self.sample_key,
            "segments": self.segments,
            "qualifiers": self.qualifiers,
            "xrefs": self.xrefs,
            "xref_sources": self.xref_sources,
        }
        return aliases.get(field, None)

    def qualifier(self, key: str) -> tuple[str, ...]:
        return self.qualifiers.get(key, ())

    def to_dict(self, *, selected: tuple[str, ...] | None = None) -> dict[str, object]:
        """Return the schema-complete JSON representation."""

        payload: dict[str, object] = {
            "schema_version": FEATURE_SCHEMA_VERSION,
            "source": self.source,
            "record": self.record,
            "record_index": self.record_index,
            "feature_index": self.feature_index,
            "type": self.type,
            "locus_tag": self.locus_tag,
            "gene": self.gene,
            "product": self.product,
            "protein_id": self.protein_id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand_symbol,
            "strand_value": self.strand,
            "length": self.length,
            "record_length": self.record_length,
            "topology": self.topology,
            "partial": self.partial,
            "pseudo": self.pseudo,
            "segments": [
                {"start": start, "end": end} for start, end in self.segments
            ],
            "qualifiers": {
                key: list(values) for key, values in self.qualifiers.items()
            },
            "xrefs": {key: list(values) for key, values in self.xrefs.items()},
            "xref_sources": {
                key: [
                    {"value": value, "source_field": source_field}
                    for value, source_field in values
                ]
                for key, values in self.xref_sources.items()
            },
        }
        if selected is None:
            return payload
        selected_payload: dict[str, object] = {
            "schema_version": FEATURE_SCHEMA_VERSION,
        }
        for field in selected:
            selected_payload[field] = (
                self.value(field)
                if field not in {"segments", "qualifiers", "xrefs", "xref_sources"}
                else payload[field]
            )
        return selected_payload

    def annotation_dict(self) -> dict[str, str]:
        return {
            "contig": self.record,
            "start": str(self.start),
            "end": str(self.end),
            "strand": self.strand_symbol,
            "type": self.type,
            "locus_tag": self.locus_tag or "",
            "gene": self.gene or "",
            "product": self.product or "",
            "ec_number": ";".join(self.xrefs.get("ec_numbers", ())),
            "cog": ";".join(self.xrefs.get("cog_ids", ())),
            "kegg_ko": ";".join(self.xrefs.get("kegg_kos", ())),
            "pfam": ";".join(self.xrefs.get("pfam", ())),
            "rfam": ";".join(self.xrefs.get("rfam", ())),
            "db_xref": ";".join(self.xrefs.get("db_xrefs", ())),
            "protein_id": self.protein_id or "",
        }


def feature_rows(document: object) -> tuple[FeatureRow, ...]:
    """Build rows in canonical record/feature order."""

    source_value = getattr(document, "source_label", "") or getattr(document, "path", None)
    source = "" if source_value is None else str(source_value)
    rows: list[FeatureRow] = []
    for record in getattr(document, "records", ()):
        for feature in record.features:
            rows.append(FeatureRow.from_feature(feature, source=source))
    return tuple(rows)


def render_annotations_tsv(rows: tuple[FeatureRow, ...] | list[FeatureRow]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=ANNOTATION_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(row.annotation_dict() for row in rows)
    return output.getvalue()


def render_feature_jsonl(rows: tuple[FeatureRow, ...] | list[FeatureRow]) -> str:
    lines = [
        json.dumps(row.to_dict(), ensure_ascii=False, allow_nan=False, sort_keys=True)
        for row in rows
    ]
    return "\n".join(lines) + ("\n" if lines else "")


def render_rows(
    rows: tuple[FeatureRow, ...] | list[FeatureRow],
    *,
    format_type: str,
    selected: tuple[str, ...] | None = None,
) -> str:
    """Render selected query rows in a deterministic machine/human format."""

    if format_type == "json":
        return (
            json.dumps(
                [row.to_dict(selected=selected) for row in rows],
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )
    if format_type == "jsonl":
        return "\n".join(
            json.dumps(
                row.to_dict(selected=selected),
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
            )
            for row in rows
        ) + ("\n" if rows else "")
    fields = selected or (
        "record",
        "feature_index",
        "type",
        "locus_tag",
        "gene",
        "product",
        "start",
        "end",
        "strand",
    )
    output = io.StringIO(newline="")
    delimiter = "\t" if format_type == "tsv" else ","
    writer = csv.writer(output, delimiter=delimiter, lineterminator="\n")
    writer.writerow(fields)
    for row in rows:
        values: list[object] = []
        for field in fields:
            value = row.value(field)
            if isinstance(value, (tuple, list)):
                value = ";".join("" if item is None else str(item) for item in value)
            values.append("" if value is None else value)
        writer.writerow(values)
    if format_type == "text":
        lines = [
            f"{row.record}:{row.start}..{row.end}({row.strand_symbol}) "
            f"{row.type} {row.locus_tag or '-'} {row.gene or '-'} "
            f"{row.product or '-'}"
            for row in rows
        ]
        return "\n".join(lines) + ("\n" if lines else "")
    return output.getvalue()


def record_rows(record: GenBankRecord, *, source: str = "") -> tuple[FeatureRow, ...]:
    return tuple(FeatureRow.from_feature(feature, source=source) for feature in record.features)


__all__ = [
    "ANNOTATION_COLUMNS",
    "FEATURE_SCHEMA_VERSION",
    "FeatureRow",
    "feature_rows",
    "record_rows",
    "render_annotations_tsv",
    "render_feature_jsonl",
    "render_rows",
]
