"""Lossless record inventory, selection, and serialization."""

from __future__ import annotations

import copy
import csv
import hashlib
import io
import json
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from Bio import SeqIO

from .model import GenBankDocument, GenBankRecord

RECORD_SCHEMA_VERSION = "gbparse.record.v1"
MAX_RECORD_REGEX_LENGTH = 1_024


@dataclass(frozen=True)
class RecordRow:
    """Schema-versioned projection of one source GenBank record."""

    source: str
    record_index: int
    id: str
    name: str
    description: str
    length: int
    topology: str | None
    molecule_type: str | None
    division: str | None
    date: str | None
    accessions: tuple[str, ...]
    organism: str | None
    taxonomy: tuple[str, ...]
    strain: str | None
    isolate: str | None
    gc_percent: float
    feature_count: int
    cds_count: int
    rrna_count: int
    trna_count: int
    tmrna_count: int
    ncrna_count: int
    pseudogene_count: int
    has_sequence: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": RECORD_SCHEMA_VERSION,
            "source": self.source,
            "record_index": self.record_index,
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "length": self.length,
            "topology": self.topology,
            "molecule_type": self.molecule_type,
            "division": self.division,
            "date": self.date,
            "accessions": list(self.accessions),
            "organism": self.organism,
            "taxonomy": list(self.taxonomy),
            "strain": self.strain,
            "isolate": self.isolate,
            "gc_percent": self.gc_percent,
            "feature_count": self.feature_count,
            "cds_count": self.cds_count,
            "rrna_count": self.rrna_count,
            "trna_count": self.trna_count,
            "tmrna_count": self.tmrna_count,
            "ncrna_count": self.ncrna_count,
            "pseudogene_count": self.pseudogene_count,
            "has_sequence": self.has_sequence,
        }


@dataclass(frozen=True)
class RecordSelector:
    """Record-level predicates used by list, filter, and split actions."""

    record_ids: tuple[str, ...] = ()
    record_regex: str | None = None
    min_length: int | None = None
    max_length: int | None = None
    topology: str | None = None
    molecule_type: str | None = None
    has_features: tuple[str, ...] = ()
    has_loci: tuple[str, ...] = ()
    invert: bool = False

    @property
    def has_predicates(self) -> bool:
        return bool(
            self.record_ids
            or self.record_regex
            or self.min_length is not None
            or self.max_length is not None
            or self.topology
            or self.molecule_type
            or self.has_features
            or self.has_loci
        )


def _first_qualifier(record: GenBankRecord, key: str) -> str | None:
    for feature in record.features:
        if feature.type.casefold() == "source":
            value = feature.get_qual(key)
            if value:
                return value
    annotation_value = record.annotations.get(key)
    if isinstance(annotation_value, (list, tuple)):
        return str(annotation_value[0]) if annotation_value else None
    return str(annotation_value) if annotation_value else None


def _annotation_values(record: GenBankRecord, *keys: str) -> tuple[str, ...]:
    values: list[str] = []
    for key in keys:
        value = record.annotations.get(key)
        if isinstance(value, (list, tuple)):
            values.extend(str(item) for item in value)
        elif value:
            values.append(str(value))
    return tuple(dict.fromkeys(values))


def record_rows(document: GenBankDocument) -> tuple[RecordRow, ...]:
    source = document.source_label or str(document.path or "")
    rows: list[RecordRow] = []
    for index, record in enumerate(document.records, 1):
        record_index = record.record_index or index
        counts: dict[str, int] = {}
        for feature in record.features:
            key = feature.type.casefold()
            counts[key] = counts.get(key, 0) + 1
        rows.append(
            RecordRow(
                source=source,
                record_index=record_index,
                id=record.id,
                name=record.name,
                description=record.description,
                length=record.length,
                topology=record.topology,
                molecule_type=record.molecule_type,
                division=record.division,
                date=record.date,
                accessions=_annotation_values(record, "accessions", "accession"),
                organism=_first_qualifier(record, "organism"),
                taxonomy=_annotation_values(record, "taxonomy"),
                strain=_first_qualifier(record, "strain"),
                isolate=_first_qualifier(record, "isolate"),
                gc_percent=record.gc_content,
                feature_count=len(record.features),
                cds_count=counts.get("cds", 0),
                rrna_count=counts.get("rrna", 0),
                trna_count=counts.get("trna", 0),
                tmrna_count=counts.get("tmrna", 0),
                ncrna_count=counts.get("ncrna", 0),
                pseudogene_count=sum(feature.is_pseudo for feature in record.features),
                has_sequence=bool(len(record.seq)),
            )
        )
    return tuple(rows)


def _record_matches(record: GenBankRecord, selector: RecordSelector) -> bool:
    if selector.record_ids:
        candidates = {record.id, record.name}
        if not any(value in candidates for value in selector.record_ids):
            return False
    if selector.record_regex and not re.search(selector.record_regex, f"{record.id}\n{record.name}"):
        return False
    if selector.min_length is not None and record.length < selector.min_length:
        return False
    if selector.max_length is not None and record.length > selector.max_length:
        return False
    if selector.topology and (record.topology or "unknown").casefold() != selector.topology.casefold():
        return False
    if selector.molecule_type and (record.molecule_type or "").casefold() != selector.molecule_type.casefold():
        return False
    feature_types = {feature.type.casefold() for feature in record.features}
    if any(value.casefold() not in feature_types for value in selector.has_features):
        return False
    loci = {feature.locus_tag for feature in record.features if feature.locus_tag}
    return not any(value not in loci for value in selector.has_loci)


def select_records(
    document: GenBankDocument,
    selector: RecordSelector | None = None,
) -> tuple[GenBankRecord, ...]:
    """Select records in source order; ``invert`` applies to the final match."""

    active = selector or RecordSelector()
    return tuple(
        record
        for record in document.records
        if (_record_matches(record, active) != active.invert)
    )


def select_exact_records(document: GenBankDocument, selectors: Sequence[str]) -> tuple[GenBankRecord, ...]:
    """Resolve exact IDs first, then exact names, rejecting ambiguity."""

    selected: list[GenBankRecord] = []
    for selector in selectors:
        id_matches = [record for record in document.records if record.id == selector]
        matches = id_matches or [record for record in document.records if record.name == selector]
        if not matches:
            raise ValueError(f"record selector not found: {selector}")
        if len(matches) > 1:
            indices = [str(document.records.index(record) + 1) for record in matches]
            raise ValueError(f"record selector {selector!r} is ambiguous at indices {', '.join(indices)}")
        if matches[0] not in selected:
            selected.append(matches[0])
    order = {id(record): index for index, record in enumerate(document.records)}
    return tuple(sorted(selected, key=lambda record: order[id(record)]))


def render_record_rows(rows: Iterable[RecordRow], format_type: str = "text") -> str:
    rows = tuple(rows)
    if format_type == "json":
        return json.dumps([row.to_dict() for row in rows], ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"
    if format_type == "jsonl":
        return "\n".join(json.dumps(row.to_dict(), ensure_ascii=False, allow_nan=False, sort_keys=True) for row in rows) + ("\n" if rows else "")
    fields = tuple(key for key in RecordRow.__dataclass_fields__ if key != "schema_version")
    if format_type in {"tsv", "csv"}:
        output = io.StringIO(newline="")
        writer = csv.writer(output, delimiter="\t" if format_type == "tsv" else ",", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            payload = row.to_dict()
            writer.writerow(";".join(str(value) for value in payload[field]) if isinstance(payload[field], list) else ("" if payload[field] is None else payload[field]) for field in fields)
        return output.getvalue()
    lines = [
        f"{row.record_index}: {row.id} ({row.length:,} bp, {row.topology or 'unknown'}) "
        f"features={row.feature_count} CDS={row.cds_count} GC={row.gc_percent:.2f}%"
        for row in rows
    ]
    return "\n".join(lines) + ("\n" if lines else "")


def _serialize_records(records: Sequence[GenBankRecord], format_type: str) -> bytes:
    raw_records = []
    for record in records:
        if record.raw_record is None:
            raise ValueError(f"record {record.id!r} has no retained source SeqRecord")
        raw_records.append(copy.deepcopy(record.raw_record))
    if format_type == "genbank":
        handle = io.StringIO()
        SeqIO.write(raw_records, handle, "genbank")
        return handle.getvalue().encode("utf-8")
    if format_type == "fasta":
        handle = io.StringIO()
        SeqIO.write(raw_records, handle, "fasta")
        return handle.getvalue().encode("utf-8")
    raise ValueError(f"unsupported record output format: {format_type}")


def render_selected_records(records: Sequence[GenBankRecord], format_type: str) -> bytes:
    return _serialize_records(records, format_type)


def _safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return cleaned or "record"


def split_record_artifacts(
    records: Sequence[GenBankRecord],
    format_type: str,
    *,
    source: str = "",
) -> dict[str, str | bytes]:
    extension = ".gbff" if format_type == "genbank" else ".fna"
    used: set[str] = set()
    files: dict[str, str | bytes] = {}
    manifest_rows = [[
        "schema_version",
        "source",
        "record_index",
        "record_id",
        "record_name",
        "length",
        "filename",
        "sha256",
    ]]
    for index, record in enumerate(records, 1):
        record_index = record.record_index or index
        base = _safe_filename(record.id)
        filename = f"{base}{extension}"
        if filename.casefold() in used:
            filename = f"{base}-r{index}{extension}"
        while filename.casefold() in used:
            filename = f"{base}-r{index}-{len(used)}{extension}"
        used.add(filename.casefold())
        payload = _serialize_records((record,), format_type)
        files[filename] = payload
        digest = hashlib.sha256(payload).hexdigest()
        source_label = source or record.id
        manifest_rows.append([
            "gbparse.records.v1",
            source_label,
            str(record_index),
            record.id,
            record.name,
            str(record.length),
            filename,
            digest,
        ])
    manifest = io.StringIO(newline="")
    csv.writer(manifest, delimiter="\t", lineterminator="\n").writerows(manifest_rows)
    files["records.manifest.tsv"] = manifest.getvalue()
    return files


__all__ = [
    "MAX_RECORD_REGEX_LENGTH",
    "RECORD_SCHEMA_VERSION",
    "RecordRow",
    "RecordSelector",
    "record_rows",
    "render_record_rows",
    "render_selected_records",
    "select_exact_records",
    "select_records",
    "split_record_artifacts",
]
