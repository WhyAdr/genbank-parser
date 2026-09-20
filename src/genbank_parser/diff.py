"""Identity-aware comparison of CDS annotations between two GenBank files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .io import extract_xrefs, read_genbank
from .model import GenBankFeature

DIFF_TSV_COLUMNS = (
    "record",
    "old_locus",
    "new_locus",
    "old_start",
    "old_end",
    "new_start",
    "new_end",
    "old_strand",
    "new_strand",
    "match_tier",
    "classifications",
    "old_gene",
    "new_gene",
    "old_product",
    "new_product",
    "old_kos",
    "new_kos",
    "old_ecs",
    "new_ecs",
)


@dataclass(frozen=True)
class _PafMapping:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int


def _load_paf(path: str | Path) -> tuple[dict[str, _PafMapping], tuple[str, ...]]:
    mappings: dict[str, list[_PafMapping]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 columns")
            try:
                mapping = _PafMapping(
                    query=fields[0],
                    query_length=int(fields[1]),
                    query_start=int(fields[2]),
                    query_end=int(fields[3]),
                    strand=fields[4],
                    target=fields[5],
                    target_length=int(fields[6]),
                    target_start=int(fields[7]),
                    target_end=int(fields[8]),
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid PAF coordinates on line {line_number}") from exc
            if mapping.strand not in {"+", "-"} or not (
                0 <= mapping.query_start < mapping.query_end <= mapping.query_length
                and 0 <= mapping.target_start < mapping.target_end <= mapping.target_length
            ):
                raise ValueError(f"invalid PAF mapping on line {line_number}")
            mappings.setdefault(mapping.query, []).append(mapping)
    unique: dict[str, _PafMapping] = {}
    ambiguous: list[str] = []
    for query, candidates in sorted(mappings.items()):
        if len(candidates) == 1:
            unique[query] = candidates[0]
        else:
            ambiguous.append(query)
    return unique, tuple(ambiguous)


def _map_feature_interval(
    feature: GenBankFeature,
    mapping: _PafMapping | None,
) -> tuple[str, int, int, int | None] | None:
    """Map a wholly contained one-based feature interval through one PAF row."""

    if mapping is None:
        return None
    start0 = feature.start - 1
    end0 = feature.end
    if start0 < mapping.query_start or end0 > mapping.query_end:
        return None
    if mapping.strand == "+":
        mapped_start0 = mapping.target_start + start0 - mapping.query_start
        mapped_end0 = mapping.target_start + end0 - mapping.query_start
    else:
        mapped_start0 = mapping.target_end - (end0 - mapping.query_start)
        mapped_end0 = mapping.target_end - (start0 - mapping.query_start)
    mapped_strand = feature.strand
    if mapped_strand in (1, -1) and mapping.strand == "-":
        mapped_strand = -mapped_strand
    return mapping.target, mapped_start0 + 1, mapped_end0, mapped_strand


def _mapped_coordinates_payload(
    mapped: tuple[str, int, int, int | None] | None,
) -> dict[str, object] | None:
    if mapped is None:
        return None
    record, start, end, strand = mapped
    strand_symbol = "." if strand is None else "+" if strand == 1 else "-" if strand == -1 else "?"
    return {
        "record": record,
        "start": start,
        "end": end,
        "strand": strand_symbol,
        "strand_value": strand,
    }


def _detail_row(detail: dict[str, Any]) -> dict[str, object]:
    old = detail.get("old") or {}
    new = detail.get("new") or {}
    old_xrefs = old.get("xrefs", {})
    new_xrefs = new.get("xrefs", {})
    return {
        "record": detail.get("record", ""),
        "old_locus": detail.get("old_locus", ""),
        "new_locus": detail.get("new_locus", ""),
        "old_start": old.get("start", ""),
        "old_end": old.get("end", ""),
        "new_start": new.get("start", ""),
        "new_end": new.get("end", ""),
        "old_strand": old.get("strand", ""),
        "new_strand": new.get("strand", ""),
        "match_tier": detail.get("match", ""),
        "classifications": ";".join(detail.get("classification", [])),
        "old_gene": old.get("gene", detail.get("old_gene", "")),
        "new_gene": new.get("gene", detail.get("new_gene", "")),
        "old_product": old.get("product", detail.get("old_product", "")),
        "new_product": new.get("product", detail.get("new_product", "")),
        "old_kos": ";".join(old_xrefs.get("kegg_kos", detail.get("old_kos", []))),
        "new_kos": ";".join(new_xrefs.get("kegg_kos", detail.get("new_kos", []))),
        "old_ecs": ";".join(old_xrefs.get("ec_numbers", detail.get("old_ecs", []))),
        "new_ecs": ";".join(new_xrefs.get("ec_numbers", detail.get("new_ecs", []))),
    }


def render_diff(
    result: dict[str, Any],
    format_type: str = "text",
    *,
    max_display: int = 25,
    summary_only: bool = False,
) -> str:
    details = [] if summary_only else list(result.get("changes", []))
    if format_type == "json":
        payload = dict(result)
        payload["summary_only"] = summary_only
        payload["changes"] = details
        payload["result_complete"] = True
        return json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"
    if format_type == "jsonl":
        if summary_only:
            metadata = {
                key: value
                for key, value in result.items()
                if key not in {"changes", "details"}
            }
            metadata.update(
                {
                    "record_type": "summary",
                    "summary_only": True,
                    "result_complete": True,
                }
            )
            return json.dumps(metadata, ensure_ascii=False, allow_nan=False, sort_keys=True) + "\n"
        return "".join(json.dumps(detail, ensure_ascii=False, allow_nan=False, sort_keys=True) + "\n" for detail in details)
    if format_type == "tsv":
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=DIFF_TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(_detail_row(detail) for detail in details)
        return output.getvalue()
    lines = [
        "=" * 70,
        "  GENOME ANNOTATION DIFF REPORT",
        "=" * 70,
        f"  Reference (Old): {result['old_file']}",
        f"  Comparison (New): {result['new_file']}",
        f"  CDS Count      : {result['old_cds_count']} (old) -> {result['new_cds_count']} (new)",
        f"  Matched        : {result['matched_cds']:,}",
        f"  Added / Removed: {result['added_cds']:,} / {result['removed_cds']:,}",
        f"  Boundary shifts: {result['boundary_shifted_cds']:,}",
        f"  Product changes: {result['product_name_changes']:,}",
        f"  Gene changes   : {result['gene_name_changes']:,}",
        f"  KO/EC changes  : {result['xref_changes']:,}",
        f"  Change records : {result['change_count']:,}",
    ]
    if summary_only:
        lines.append("  Detail records omitted by --summary-only")
    else:
        lines.append("")
        for detail in details[:max_display]:
            lines.append(
                f"  [{','.join(detail.get('classification', []))}] "
                f"{detail.get('record', '')} {detail.get('old_locus', '')} -> {detail.get('new_locus', '')}"
            )
        if len(details) > max_display:
            lines.append(f"  ... displayed {max_display} of {len(details)} change records")
    lines.append("=" * 70)
    return "\n".join(lines) + "\n"


def _coordinate_key(feature: GenBankFeature) -> tuple[str, int, int, int | None]:
    return feature.record_id, feature.start, feature.end, feature.strand


def _overlap_fraction(old: GenBankFeature, new: GenBankFeature) -> float:
    if old.record_id != new.record_id or old.strand != new.strand:
        return 0.0
    overlap = max(0, min(old.end, new.end) - max(old.start, new.start) + 1)
    if overlap == 0:
        return 0.0
    old_span = old.end - old.start + 1
    new_span = new.end - new.start + 1
    return min(overlap / old_span, overlap / new_span)


def _translation_hash(feature: GenBankFeature) -> str | None:
    translation = feature.translation.strip()
    if not translation:
        return None
    return hashlib.sha256(translation.encode("utf-8")).hexdigest()


def _pair_features(
    old_features: list[GenBankFeature],
    new_features: list[GenBankFeature],
    coordinate_mappings: dict[str, _PafMapping] | None = None,
) -> tuple[list[tuple[GenBankFeature, GenBankFeature, str]], set[int], set[int]]:
    """Match CDSs by IDs, overlap, then translated sequence identity."""
    unmatched_old = set(range(len(old_features)))
    unmatched_new = set(range(len(new_features)))
    matches: list[tuple[GenBankFeature, GenBankFeature, str]] = []

    def match_by_key(key_fn: Any, tier: str) -> None:
        for old_index in sorted(unmatched_old):
            key = key_fn(old_features[old_index])
            if not key:
                continue
            candidates = [
                index for index in unmatched_new if key_fn(new_features[index]) == key
            ]
            if not candidates:
                continue
            new_index = max(
                candidates,
                key=lambda index: _overlap_fraction(
                    old_features[old_index], new_features[index]
                ),
            )
            unmatched_old.remove(old_index)
            unmatched_new.remove(new_index)
            matches.append((old_features[old_index], new_features[new_index], tier))

    match_by_key(lambda feature: feature.locus_tag, "locus_tag")
    match_by_key(lambda feature: feature.protein_id, "protein_id")

    for old_index in sorted(unmatched_old):
        old_feature = old_features[old_index]
        mapped = _map_feature_interval(
            old_feature,
            (coordinate_mappings or {}).get(old_feature.record_id),
        )
        candidates: list[tuple[int, float]] = []
        for index in unmatched_new:
            new_feature = new_features[index]
            if mapped is not None:
                target_record, mapped_start, mapped_end, mapped_strand = mapped
                if new_feature.record_id != target_record:
                    continue
                if mapped_strand is not None and new_feature.strand != mapped_strand:
                    continue
                overlap = max(0, min(mapped_end, new_feature.end) - max(mapped_start, new_feature.start) + 1)
                if overlap:
                    old_span = mapped_end - mapped_start + 1
                    new_span = new_feature.end - new_feature.start + 1
                    candidates.append((index, min(overlap / old_span, overlap / new_span)))
            else:
                candidates.append((index, _overlap_fraction(old_feature, new_feature)))
        candidates = [(index, score) for index, score in candidates if score >= 0.8]
        if not candidates:
            continue
        new_index, _ = max(candidates, key=lambda item: item[1])
        unmatched_old.remove(old_index)
        unmatched_new.remove(new_index)
        matches.append(
            (old_features[old_index], new_features[new_index], "coordinate_overlap")
        )

    old_by_hash: dict[str, list[int]] = {}
    new_by_hash: dict[str, list[int]] = {}
    for index in unmatched_old:
        digest = _translation_hash(old_features[index])
        if digest:
            old_by_hash.setdefault(digest, []).append(index)
    for index in unmatched_new:
        digest = _translation_hash(new_features[index])
        if digest:
            new_by_hash.setdefault(digest, []).append(index)
    for digest in sorted(set(old_by_hash) & set(new_by_hash)):
        for old_index, new_index in zip(
            sorted(old_by_hash[digest]), sorted(new_by_hash[digest])
        ):
            if old_index in unmatched_old and new_index in unmatched_new:
                unmatched_old.remove(old_index)
                unmatched_new.remove(new_index)
                matches.append(
                    (
                        old_features[old_index],
                        new_features[new_index],
                        "translation_hash",
                    )
                )

    return matches, unmatched_old, unmatched_new


def _classify_pair(
    old: GenBankFeature,
    new: GenBankFeature,
    coordinate_mappings: dict[str, _PafMapping] | None = None,
) -> tuple[list[str], dict[str, Any]]:
    classifications: list[str] = []
    mapped_old = _map_feature_interval(old, (coordinate_mappings or {}).get(old.record_id))
    coordinates_equivalent = _coordinate_key(old) == _coordinate_key(new)
    if mapped_old is not None:
        target_record, mapped_start, mapped_end, mapped_strand = mapped_old
        coordinates_equivalent = (
            target_record == new.record_id
            and mapped_start == new.start
            and mapped_end == new.end
            and mapped_strand == new.strand
        )
    if not coordinates_equivalent:
        classifications.append("boundary_shifted")
    product_changed = old.product.strip() != new.product.strip()
    gene_changed = old.gene.strip() != new.gene.strip()
    old_xrefs = extract_xrefs(old)
    new_xrefs = extract_xrefs(new)
    xref_changes = {
        key: {"old": old_xrefs[key], "new": new_xrefs[key]}
        for key in ("kegg_kos", "ec_numbers")
        if set(old_xrefs[key]) != set(new_xrefs[key])
    }
    if product_changed:
        classifications.append("product_changed")
    if gene_changed:
        classifications.append("gene_changed")
    if xref_changes:
        classifications.append("xref_changed")
    if not classifications:
        classifications.append("unchanged")

    detail: dict[str, Any] = {
        "record": old.record_id,
        "old_locus": old.locus_tag or "-",
        "new_locus": new.locus_tag or "-",
        "old_coordinates": f"{old.start}..{old.end}({old.strand_symbol})",
        "new_coordinates": f"{new.start}..{new.end}({new.strand_symbol})",
        "mapped_old": _mapped_coordinates_payload(mapped_old),
        "match": "",
        "classification": classifications,
        "old": {
            "locus_tag": old.locus_tag or None,
            "gene": old.gene or None,
            "product": old.product or None,
            "start": old.start,
            "end": old.end,
            "strand": old.strand_symbol,
            "xrefs": old_xrefs,
        },
        "new": {
            "locus_tag": new.locus_tag or None,
            "gene": new.gene or None,
            "product": new.product or None,
            "start": new.start,
            "end": new.end,
            "strand": new.strand_symbol,
            "xrefs": new_xrefs,
        },
    }
    if product_changed:
        detail.update(
            {"old_product": old.product or "-", "new_product": new.product or "-"}
        )
    if gene_changed:
        detail.update({"old_gene": old.gene or "-", "new_gene": new.gene or "-"})
    if xref_changes:
        detail["xref_changes"] = xref_changes
    return classifications, detail


def diff_annotations(
    old_filepath: str | Path,
    new_filepath: str | Path,
    format_type: str = "text",
    output_path: str | Path | None = None,
    *,
    summary_only: bool = False,
    max_display: int = 25,
    coordinate_map: str | Path | None = None,
) -> dict[str, Any]:
    old_doc = read_genbank(old_filepath)
    new_doc = read_genbank(new_filepath)
    coordinate_mappings: dict[str, _PafMapping] = {}
    ambiguous_mappings: tuple[str, ...] = ()
    if coordinate_map is not None:
        coordinate_mappings, ambiguous_mappings = _load_paf(coordinate_map)
    old_cds = [feature for feature in old_doc.all_features if feature.type == "CDS"]
    new_cds = [feature for feature in new_doc.all_features if feature.type == "CDS"]

    matches, removed_indices, added_indices = _pair_features(
        old_cds,
        new_cds,
        coordinate_mappings=coordinate_mappings,
    )
    details: list[dict[str, Any]] = []
    classification_counts: dict[str, int] = {}
    for old, new, tier in sorted(
        matches, key=lambda pair: (pair[0].record_id, pair[0].start, pair[0].end)
    ):
        classifications, detail = _classify_pair(old, new, coordinate_mappings)
        detail["match"] = tier
        details.append(detail)
        for classification in classifications:
            classification_counts[classification] = (
                classification_counts.get(classification, 0) + 1
            )

    for index in sorted(removed_indices):
        feature = old_cds[index]
        xrefs = extract_xrefs(feature)
        mapped = _map_feature_interval(
            feature,
            (coordinate_mappings or {}).get(feature.record_id),
        )
        details.append(
            {
                "record": feature.record_id,
                "old_locus": feature.locus_tag or "-",
                "old_coordinates": f"{feature.start}..{feature.end}({feature.strand_symbol})",
                "mapped_old": _mapped_coordinates_payload(mapped),
                "classification": ["removed"],
                "match": "",
                "old": {
                    "locus_tag": feature.locus_tag or None,
                    "gene": feature.gene or None,
                    "product": feature.product or None,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand_symbol,
                    "xrefs": xrefs,
                },
                "new": {},
            }
        )
    for index in sorted(added_indices):
        feature = new_cds[index]
        xrefs = extract_xrefs(feature)
        details.append(
            {
                "record": feature.record_id,
                "new_locus": feature.locus_tag or "-",
                "new_coordinates": f"{feature.start}..{feature.end}({feature.strand_symbol})",
                "classification": ["added"],
                "match": "",
                "old": {},
                "new": {
                    "locus_tag": feature.locus_tag or None,
                    "gene": feature.gene or None,
                    "product": feature.product or None,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand_symbol,
                    "xrefs": xrefs,
                },
                "mapped_old": None,
            }
        )
    classification_counts["removed"] = len(removed_indices)
    classification_counts["added"] = len(added_indices)

    exact_shared = {_coordinate_key(feature) for feature in old_cds} & {
        _coordinate_key(feature) for feature in new_cds
    }
    product_changes = [
        detail
        for detail in details
        if "product_changed" in detail.get("classification", [])
    ]
    gene_changes = [
        detail
        for detail in details
        if "gene_changed" in detail.get("classification", [])
    ]
    xref_changes = [
        detail
        for detail in details
        if "xref_changed" in detail.get("classification", [])
    ]
    boundary_shifts = [
        detail
        for detail in details
        if "boundary_shifted" in detail.get("classification", [])
    ]

    result: dict[str, Any] = {
        "old_file": str(old_filepath),
        "new_file": str(new_filepath),
        "old_cds_count": len(old_cds),
        "new_cds_count": len(new_cds),
        "shared_coordinate_cds": len(exact_shared),
        "matched_cds": len(matches),
        "added_cds": len(added_indices),
        "removed_cds": len(removed_indices),
        "product_name_changes": len(product_changes),
        "gene_name_changes": len(gene_changes),
        "xref_changes": len(xref_changes),
        "boundary_shifted_cds": len(boundary_shifts),
        "classification_counts": classification_counts,
        "changes": details,
        "change_count": len(details),
        "result_complete": True,
        "summary_only": summary_only,
        "coordinate_map": str(coordinate_map) if coordinate_map is not None else None,
        "coordinate_map_summary": {
            "mapped_query_records": sorted(coordinate_mappings),
            "ambiguous_query_records": list(ambiguous_mappings),
            "unmapped_query_records": sorted(
                {
                    feature.record_id
                    for feature in old_cds
                    if feature.record_id not in coordinate_mappings
                    and feature.record_id not in ambiguous_mappings
                }
            ),
            "feature_mapping_failures": sorted(
                f"{feature.record_id}:{feature.feature_index}"
                for feature in old_cds
                if feature.record_id in coordinate_mappings
                and _map_feature_interval(
                    feature, coordinate_mappings[feature.record_id]
                )
                is None
            ),
            "mapping_count": len(coordinate_mappings),
        },
        "details": {
            "product_changes": product_changes,
            "gene_changes": gene_changes,
            "xref_changes": xref_changes,
            "boundary_shifts": boundary_shifts,
            "added": [
                detail
                for detail in details
                if "added" in detail.get("classification", [])
            ],
            "removed": [
                detail
                for detail in details
                if "removed" in detail.get("classification", [])
            ],
        },
    }

    if summary_only:
        result["changes"] = []
        result["details"] = {}

    if format_type in {"json", "jsonl", "tsv", "text"}:
        rendered = render_diff(
            result,
            format_type,
            max_display=max_display,
            summary_only=summary_only,
        )
        if output_path:
            Path(output_path).write_text(rendered, encoding="utf-8", newline="")
        else:
            print(rendered)
        return result
    raise ValueError(f"unsupported diff format: {format_type}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare two annotation versions using identity and overlap matching."
    )
    parser.add_argument("old_file", help="Reference GenBank file")
    parser.add_argument("new_file", help="Updated GenBank file")
    parser.add_argument("--format", choices=["text", "json", "jsonl", "tsv"], default="text")
    parser.add_argument("--output", help="Output file path")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--max-display", type=int, default=25)
    parser.add_argument("--coordinate-map", help="Optional caller-supplied PAF coordinate map")
    args = parser.parse_args()
    diff_annotations(
        args.old_file,
        args.new_file,
        format_type=args.format,
        output_path=args.output,
        summary_only=args.summary_only,
        max_display=args.max_display,
        coordinate_map=args.coordinate_map,
    )


if __name__ == "__main__":
    main()
