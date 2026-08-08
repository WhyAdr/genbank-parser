"""Identity-aware comparison of CDS annotations between two GenBank files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .io import extract_xrefs, read_genbank
from .model import GenBankFeature


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
        candidates = [
            (index, _overlap_fraction(old_features[old_index], new_features[index]))
            for index in unmatched_new
        ]
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
    old: GenBankFeature, new: GenBankFeature
) -> tuple[list[str], dict[str, Any]]:
    classifications: list[str] = []
    if _coordinate_key(old) != _coordinate_key(new):
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
        "match": "",
        "classification": classifications,
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
) -> dict[str, Any]:
    old_doc = read_genbank(old_filepath)
    new_doc = read_genbank(new_filepath)
    old_cds = [feature for feature in old_doc.all_features if feature.type == "CDS"]
    new_cds = [feature for feature in new_doc.all_features if feature.type == "CDS"]

    matches, removed_indices, added_indices = _pair_features(old_cds, new_cds)
    details: list[dict[str, Any]] = []
    classification_counts: dict[str, int] = {}
    for old, new, tier in sorted(
        matches, key=lambda pair: (pair[0].record_id, pair[0].start, pair[0].end)
    ):
        classifications, detail = _classify_pair(old, new)
        detail["match"] = tier
        details.append(detail)
        for classification in classifications:
            classification_counts[classification] = (
                classification_counts.get(classification, 0) + 1
            )

    for index in sorted(removed_indices):
        feature = old_cds[index]
        details.append(
            {
                "record": feature.record_id,
                "old_locus": feature.locus_tag or "-",
                "old_coordinates": f"{feature.start}..{feature.end}({feature.strand_symbol})",
                "classification": ["removed"],
            }
        )
    for index in sorted(added_indices):
        feature = new_cds[index]
        details.append(
            {
                "record": feature.record_id,
                "new_locus": feature.locus_tag or "-",
                "new_coordinates": f"{feature.start}..{feature.end}({feature.strand_symbol})",
                "classification": ["added"],
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
        "changes": details[:200],
        "details": {
            "product_changes": product_changes[:50],
            "gene_changes": gene_changes[:50],
            "xref_changes": xref_changes[:50],
            "boundary_shifts": boundary_shifts[:50],
            "added": [
                detail
                for detail in details
                if "added" in detail.get("classification", [])
            ][:50],
            "removed": [
                detail
                for detail in details
                if "removed" in detail.get("classification", [])
            ][:50],
        },
    }

    if format_type == "json":
        rendered = json.dumps(result, indent=2)
        if output_path:
            Path(output_path).write_text(rendered, encoding="utf-8")
        else:
            print(rendered)
        return result

    print("=" * 70)
    print("  GENOME ANNOTATION DIFF REPORT")
    print("=" * 70)
    print(f"  Reference (Old): {old_filepath}")
    print(f"  Comparison (New): {new_filepath}")
    print(f"  CDS Count      : {len(old_cds)} (old) -> {len(new_cds)} (new)")
    print(f"  Matched        : {len(matches):,}")
    print(f"  Added / Removed: {len(added_indices):,} / {len(removed_indices):,}")
    print(f"  Boundary shifts: {len(boundary_shifts):,}")
    print(f"  Product changes: {len(product_changes):,}")
    print(f"  Gene changes   : {len(gene_changes):,}")
    print(f"  KO/EC changes  : {len(xref_changes):,}")
    print("=" * 70)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare two annotation versions using identity and overlap matching."
    )
    parser.add_argument("old_file", help="Reference GenBank file")
    parser.add_argument("new_file", help="Updated GenBank file")
    parser.add_argument("--format", choices=["text", "json"], default="text")
    parser.add_argument("--output", help="Output file path")
    args = parser.parse_args()
    diff_annotations(
        args.old_file, args.new_file, format_type=args.format, output_path=args.output
    )


if __name__ == "__main__":
    main()
