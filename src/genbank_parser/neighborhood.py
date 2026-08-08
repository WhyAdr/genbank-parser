"""Genomic neighborhood viewer with circular topology support."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .io import read_genbank
from .model import GenBankFeature


def extract_neighborhood(
    filepath: str | Path,
    locus_tag: str,
    window: int = 5,
) -> list[GenBankFeature]:
    if window < 0:
        raise ValueError("window must be non-negative")

    doc = read_genbank(filepath)
    target_match = doc.find_locus(locus_tag)

    if target_match is None:
        # Fallback search by gene name, still preferring a CDS.
        for rec in doc.records:
            gene_matches = [
                f
                for f in rec.features
                if f.gene and f.gene.casefold() == locus_tag.casefold()
            ]
            for f in gene_matches:
                if f.type.casefold() == "cds":
                    target_match = (rec, f)
                    break
            if target_match:
                break
            for f in gene_matches:
                if f.gene and f.gene.casefold() == locus_tag.casefold():
                    target_match = (rec, f)
                    break
            if target_match:
                break

    if target_match is None:
        print(
            f"ERROR: Locus tag or gene '{locus_tag}' not found in {filepath}",
            file=sys.stderr,
        )
        sys.exit(1)

    rec, target_feat = target_match
    cdss = sorted(rec.cds_features, key=lambda f: f.start)
    if not cdss:
        raise ValueError(
            f"Resolved target {locus_tag!r} is on {rec.id!r}, which has no CDS features"
        )

    try:
        target_idx = [f.feature_index for f in cdss].index(target_feat.feature_index)
    except ValueError:
        raise ValueError(
            f"Resolved target {locus_tag!r} is a {target_feat.type!r}, not present in CDS ordering"
        ) from None

    n_feats = len(cdss)
    is_circular = rec.topology == "circular"

    if is_circular:
        # Circular wrap-around indexing
        indices = [
            (target_idx + offset) % n_feats for offset in range(-window, window + 1)
        ]
        neighbor_features = [cdss[idx] for idx in indices]
    else:
        start_idx = max(0, target_idx - window)
        end_idx = min(n_feats, target_idx + window + 1)
        neighbor_features = cdss[start_idx:end_idx]

    print("=" * 90)
    print(f"  GENOMIC NEIGHBORHOOD AROUND {locus_tag}")
    print("=" * 90)
    print(
        f"  Contig       : {rec.id} (topology: {rec.topology or 'linear'}, length: {rec.length:,} bp)"
    )
    print(f"  Window size  : +/- {window} CDSs")
    print()
    print(
        f"{'Strand':6s}  {'Locus Tag':18s}  {'Gene':8s}  {'Start':>10s}  {'End':>10s}  {'Product'}"
    )
    print("-" * 90)

    for f in neighbor_features:
        is_target = f.feature_index == target_feat.feature_index
        pointer = ">> " if is_target else "   "
        tag = f.locus_tag or "-"
        gene = f.gene or "-"
        prod = (f.product or "-")[:40]
        strand_sym = f.strand_symbol
        print(
            f"{pointer}[{strand_sym}]  {tag:18s}  {gene:8s}  {f.start:>10,d}  {f.end:>10,d}  {prod}"
        )

    print("=" * 90)
    return neighbor_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="View genomic neighborhood (+/- N genes) around a target locus."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("locus_tag", help="Target locus tag or gene name")
    parser.add_argument(
        "window",
        nargs="?",
        type=int,
        default=5,
        help="Number of flanking genes (default: 5)",
    )
    args = parser.parse_args()

    extract_neighborhood(args.input, args.locus_tag, window=args.window)


if __name__ == "__main__":
    main()
