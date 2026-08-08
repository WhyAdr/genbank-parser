"""Identify operon candidates: consecutive co-directional genes within gap/overlap bounds."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from .io import read_genbank
from .model import GenBankFeature


def operon_candidates(
    filepath: str | Path,
    max_gap: int = 150,
    min_gap: int = -50,
) -> list[tuple[GenBankFeature, GenBankFeature, int]]:
    doc = read_genbank(filepath)
    all_pairs: list[tuple[GenBankFeature, GenBankFeature, int]] = []
    all_clusters: list[list[GenBankFeature]] = []

    print(f"-- Operon candidates (same-strand, gap between {min_gap} and {max_gap} bp) --")
    print(f"{'Locus A':18s}  {'Gene':6s}  ->  {'Locus B':18s}  {'Gene':6s}  {'Gap':>5s}  Strand  Products")
    print("-" * 110)

    for rec in doc.records:
        cdss = sorted(rec.cds_features, key=lambda f: f.start)
        rec_pairs: list[tuple[GenBankFeature, GenBankFeature, int]] = []

        for i in range(len(cdss) - 1):
            a, b = cdss[i], cdss[i + 1]
            if a.strand == b.strand and a.strand in (1, -1):
                gap = b.start - a.end - 1
                if min_gap <= gap <= max_gap:
                    rec_pairs.append((a, b, gap))
                    all_pairs.append((a, b, gap))

        for a, b, gap in rec_pairs:
            ta = a.locus_tag or '?'
            tb = b.locus_tag or '?'
            ga = a.gene or '-'
            gb = b.gene or '-'
            pa = (a.product or '-')[:30]
            pb = (b.product or '-')[:30]
            print(f"  {ta:18s}  {ga:6s}  ->  {tb:18s}  {gb:6s}  {gap:>4d}bp  {a.strand_symbol:>3s}    {pa} | {pb}")

        # Tight clusters (>= 3 consecutive co-directional genes)
        if rec_pairs:
            current_cluster: list[GenBankFeature] = [rec_pairs[0][0], rec_pairs[0][1]]
            for a, b, _ in rec_pairs[1:]:
                if current_cluster[-1].feature_index == a.feature_index:
                    current_cluster.append(b)
                else:
                    if len(current_cluster) >= 3:
                        all_clusters.append(current_cluster)
                    current_cluster = [a, b]
            if len(current_cluster) >= 3:
                all_clusters.append(current_cluster)

    print(f"\nTotal candidate pairs: {len(all_pairs)}")

    if all_clusters:
        print("\n-- Tight clusters (>=3 consecutive co-directional genes) --")
        for cl in all_clusters:
            tags = [f.locus_tag or '?' for f in cl]
            genes = [f.gene or '-' for f in cl]
            span = f"{cl[0].start:,}..{cl[-1].end:,}  (contig: {cl[0].record_id})"
            print(f"  [{cl[0].strand_symbol}] {' -> '.join(genes)}  ({span})")
            print(f"       Tags: {', '.join(tags)}")
    else:
        print("\n-- Tight clusters: none found --")

    return all_pairs


def main() -> None:
    parser = argparse.ArgumentParser(description="Identify operon candidates: consecutive co-directional genes within gap/overlap bounds.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('max_gap', nargs='?', type=int, default=150, help="Maximum intergenic gap in bp (default: 150)")
    parser.add_argument('--min-gap', type=int, default=-50, help="Minimum gap (allows overlaps, default: -50)")
    args = parser.parse_args()

    operon_candidates(args.input, max_gap=args.max_gap, min_gap=args.min_gap)


if __name__ == '__main__':
    main()
