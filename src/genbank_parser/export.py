"""Unified interoperable exporters built on the canonical parser model."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from . import __version__
from .gff import convert_to_gff3, render_gff3_document
from .io import read_genbank
from .model import GenBankDocument, GenBankFeature
from .serializers import (
    feature_rows,
    render_annotations_tsv,
    render_feature_jsonl,
)

EXPORT_FORMATS = (
    "annotations-tsv",
    "jsonl",
    "faa",
    "ffn",
    "fna",
    "gff3",
    "bed12",
    "ncbi-table",
)


@dataclass(frozen=True)
class ExportArtifact:
    """Rendered single-file data plus non-fatal diagnostics."""

    content: str
    warnings: tuple[str, ...] = ()
    metadata: dict[str, object] | None = None


def _wrap(sequence: str, width: int) -> list[str]:
    return [sequence[index : index + width] for index in range(0, len(sequence), width)]


def _fasta_id(feature: GenBankFeature) -> str:
    tag = feature.locus_tag or f"feature_{feature.feature_index}"
    return f"{feature.record_id}:{tag}"


def render_faa(document: GenBankDocument, *, strict: bool = False) -> ExportArtifact:
    lines: list[str] = []
    missing = 0
    for feature in document.all_features:
        if feature.type != "CDS":
            continue
        translation = "".join(feature.translation.split())
        if not translation:
            missing += 1
            continue
        header = _fasta_id(feature)
        details = ["type=CDS", f"feature_index={feature.feature_index}"]
        if feature.gene:
            details.append(f"gene={feature.gene}")
        lines.append(f">{header} {' '.join(details)}")
        lines.extend(_wrap(translation, 60))
    if strict and missing:
        raise ValueError(f"{missing} CDS feature(s) have no /translation qualifier")
    if missing:
        warnings = (f"skipped {missing} CDS feature(s) without /translation",)
    else:
        warnings = ()
    content = "\n".join(lines) + ("\n" if lines else "")
    return ExportArtifact(content=content, warnings=warnings, metadata={"skipped": missing})


def render_ffn(document: GenBankDocument, *, strict: bool = False) -> ExportArtifact:
    lines: list[str] = []
    missing = 0
    for record in document.records:
        for feature in record.cds_features:
            if len(record.seq) == 0:
                missing += 1
                continue
            try:
                sequence = str(feature.extract(record.seq)).upper()
            except (AttributeError, TypeError, ValueError):
                missing += 1
                continue
            if not sequence:
                missing += 1
                continue
            lines.append(
                f">{_fasta_id(feature)} type=CDS feature_index={feature.feature_index}"
            )
            lines.extend(_wrap(sequence, 70))
    if strict and missing:
        raise ValueError(f"{missing} CDS feature(s) have no extractable nucleotide sequence")
    warnings = (f"skipped {missing} CDS feature(s) without sequence",) if missing else ()
    content = "\n".join(lines) + ("\n" if lines else "")
    return ExportArtifact(content=content, warnings=warnings, metadata={"skipped": missing})


def render_fna(document: GenBankDocument) -> ExportArtifact:
    lines: list[str] = []
    missing = 0
    for record in document.records:
        sequence = str(record.seq).upper()
        if not sequence:
            missing += 1
            continue
        lines.append(f">{record.id} len={len(sequence)}")
        lines.extend(_wrap(sequence, 70))
    warnings = (f"skipped {missing} record(s) without nucleotide sequence",) if missing else ()
    return ExportArtifact(
        content="\n".join(lines) + ("\n" if lines else ""),
        warnings=warnings,
        metadata={"skipped": missing},
    )


def _bed_name(feature: GenBankFeature, suffix: str = "") -> str:
    base = feature.locus_tag or f"{feature.type}_{feature.feature_index}"
    safe = re.sub(r"[^A-Za-z0-9_.:-]+", "_", f"{feature.record_id}:{base}")
    return safe + suffix


def _bed_line(
    record_id: str,
    start: int,
    end: int,
    name: str,
    strand: str,
    blocks: list[tuple[int, int]],
) -> str:
    if start < 1 or end < start:
        raise ValueError(f"invalid BED source interval {record_id}:{start}..{end}")
    chrom_start = start - 1
    chrom_end = end
    ordered = sorted(blocks)
    block_sizes = [block_end - block_start + 1 for block_start, block_end in ordered]
    block_starts = [block_start - start for block_start, _block_end in ordered]
    if any(size <= 0 for size in block_sizes) or any(offset < 0 for offset in block_starts):
        raise ValueError(f"invalid BED blocks for {record_id}:{start}..{end}")
    bed_strand = strand if strand in {"+", "-"} else "."
    return "\t".join(
        (
            record_id,
            str(chrom_start),
            str(chrom_end),
            name,
            "0",
            bed_strand,
            str(chrom_start),
            str(chrom_end),
            "0",
            str(len(ordered)),
            ",".join(str(size) for size in block_sizes) + ",",
            ",".join(str(offset) for offset in block_starts) + ",",
        )
    )


def render_bed12(document: GenBankDocument) -> ExportArtifact:
    lines: list[str] = []
    warnings: list[str] = []
    for record in document.records:
        for feature in record.features:
            segments = list(feature.join_segments) if feature.is_compound else [(feature.start, feature.end)]
            segments = [(start, end) for start, end in segments if end >= start]
            if not segments:
                continue
            if record.length > 0 and any(
                start < 1 or end > record.length for start, end in segments
            ):
                raise ValueError(
                    f"feature coordinates exceed record bounds: {record.id}:{feature.feature_index}"
                )
            origin_spanning = bool(
                record.topology == "circular"
                and feature.is_compound
                and len(segments) > 1
                and feature.genomic_span >= record.length > 0
            )
            if origin_spanning:
                warnings.append(
                    f"split origin-spanning feature {record.id}:{feature.feature_index} into BED records"
                )
                for part_index, (start, end) in enumerate(sorted(segments), 1):
                    lines.append(
                        _bed_line(
                            record.id,
                            start,
                            end,
                            _bed_name(feature, f"|origin_spanning|part{part_index}"),
                            feature.strand_symbol,
                            [(start, end)],
                        )
                    )
                continue
            start = min(start for start, _end in segments)
            end = max(end for _start, end in segments)
            lines.append(
                _bed_line(
                    record.id,
                    start,
                    end,
                    _bed_name(feature),
                    feature.strand_symbol,
                    segments,
                )
            )
    return ExportArtifact(
        content="\n".join(lines) + ("\n" if lines else ""),
        warnings=tuple(warnings),
        metadata={"records": len(lines)},
    )


_NCBI_SUPPORTED_QUALIFIERS = frozenset(
    {
        "locus_tag",
        "gene",
        "product",
        "protein_id",
        "EC_number",
        "db_xref",
        "note",
        "inference",
        "codon_start",
        "transl_table",
        "translation",
        "pseudo",
        "pseudogene",
        "exception",
        "function",
        "gene_synonym",
        "old_locus_tag",
        "mol_type",
        "organism",
        "strain",
        "isolate",
        "chromosome",
        "plasmid",
        "genome",
        "country",
        "collection_date",
        "host",
        "isolation_source",
        "geo_loc_name",
        "lat_lon",
        "segment",
        "serotype",
    }
)


def _ncbi_location(feature: GenBankFeature) -> list[tuple[int, int]]:
    segments = list(feature.join_segments) if feature.is_compound else [(feature.start, feature.end)]
    if feature.strand == -1:
        return [(end, start) for start, end in reversed(segments)]
    return segments


def render_ncbi_table(
    document: GenBankDocument,
    *,
    source_path: Path | None = None,
) -> dict[str, str | bytes]:
    """Return the candidate table2asn input files for one document."""

    fasta_lines: list[str] = []
    for record in document.records:
        sequence = str(record.seq).upper()
        fasta_lines.append(f">{record.id}")
        fasta_lines.extend(_wrap(sequence, 70))
    table_lines: list[str] = []
    unsupported: list[dict[str, object]] = []
    emitted = 0
    skipped = 0
    qualifier_count = 0
    for record in document.records:
        table_lines.append(f">Feature {record.id}")
        for feature in record.features:
            try:
                locations = _ncbi_location(feature)
            except (TypeError, ValueError):
                skipped += 1
                continue
            if not locations:
                skipped += 1
                continue
            for start, end in locations:
                table_lines.append(f"{start}\t{end}\t{feature.type}")
            emitted += 1
            for key in sorted(feature.qualifiers):
                values = feature.qualifiers[key]
                if key not in _NCBI_SUPPORTED_QUALIFIERS:
                    for value in values:
                        unsupported.append(
                            {
                                "record": record.id,
                                "feature_index": feature.feature_index,
                                "qualifier": key,
                                "value": value,
                            }
                        )
                    continue
                for value in values:
                    qualifier_count += 1
                    table_lines.append(f"\t\t\t/{key}={value}")
    source_hash: str | None = None
    if source_path is not None and source_path.exists() and source_path.is_file():
        digest = hashlib.sha256()
        with source_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        source_hash = digest.hexdigest()
    report = {
        "schema_version": "gbparse.ncbi-export.v1",
        "source": document.source_label,
        "source_sha256": source_hash,
        "gbparse_version": __version__,
        "record_ids": [record.id for record in document.records],
        "record_count": len(document.records),
        "feature_count": document.total_features,
        "qualifier_count": qualifier_count,
        "emitted_features": emitted,
        "skipped_features": skipped,
        "unsupported_qualifiers": unsupported,
        "warnings": [
            "These files are candidate table2asn inputs; run table2asn validation separately.",
            *([f"{len(unsupported)} unsupported qualifier value(s) were reported and omitted"] if unsupported else []),
        ],
    }
    if source_path is not None:
        source_name = source_path.name
        if source_name.casefold().endswith(".gz"):
            source_name = Path(source_name).stem
        prefix = Path(source_name).stem
    else:
        prefix = "genbank-export"
    return {
        f"{prefix}.fsa": "\n".join(fasta_lines) + ("\n" if fasta_lines else ""),
        f"{prefix}.tbl": "\n".join(table_lines) + ("\n" if table_lines else ""),
        f"{prefix}.export.json": json.dumps(
            report, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2
        )
        + "\n",
    }


def render_export(
    source: str | Path | TextIO,
    format_type: str,
    *,
    strict: bool = False,
    include_fasta: bool = False,
) -> ExportArtifact | dict[str, str | bytes]:
    """Render one unified export format from a canonical document."""

    if format_type not in EXPORT_FORMATS:
        raise ValueError(f"unsupported export format: {format_type}")
    document = read_genbank(source)
    rows = feature_rows(document)
    if format_type == "annotations-tsv":
        return ExportArtifact(render_annotations_tsv(rows))
    if format_type == "jsonl":
        return ExportArtifact(render_feature_jsonl(rows))
    if format_type == "faa":
        return render_faa(document, strict=strict)
    if format_type == "ffn":
        return render_ffn(document, strict=strict)
    if format_type == "fna":
        return render_fna(document)
    if format_type == "bed12":
        return render_bed12(document)
    if format_type == "gff3":
        if document.path is None:
            return ExportArtifact(
                render_gff3_document(document, include_fasta=include_fasta)
            )
        return ExportArtifact(convert_to_gff3(document.path, include_fasta=include_fasta))
    return render_ncbi_table(document, source_path=document.path)


__all__ = [
    "EXPORT_FORMATS",
    "ExportArtifact",
    "render_bed12",
    "render_export",
    "render_faa",
    "render_ffn",
    "render_fna",
    "render_ncbi_table",
]
