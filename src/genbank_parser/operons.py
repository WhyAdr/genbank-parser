"""Identify operon candidates: consecutive co-directional genes within gap/overlap bounds."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from .io import read_genbank
from .model import GenBankFeature


@dataclass(frozen=True)
class OperonPair:
    """One adjacent same-strand CDS pair within the requested gap bounds."""

    first: GenBankFeature
    second: GenBankFeature
    gap: int


@dataclass(frozen=True)
class OperonCluster:
    """A connected chain of at least three candidate CDSs."""

    record_id: str
    strand: int
    features: tuple[GenBankFeature, ...]
    gaps: tuple[int, ...]


@dataclass(frozen=True)
class OperonResult:
    """Structured operon-pair and cluster analysis."""

    pairs: tuple[OperonPair, ...]
    clusters: tuple[OperonCluster, ...]


def find_operon_pairs(
    features: list[GenBankFeature],
    max_gap: int = 150,
    min_gap: int = -50,
    *,
    circular: bool = False,
    record_length: int | None = None,
) -> list[tuple[GenBankFeature, GenBankFeature, int]]:
    """Return adjacent, same-strand CDS pairs within the requested gap."""
    if max_gap < min_gap:
        raise ValueError("max_gap must be greater than or equal to min_gap")
    if circular and (record_length is None or record_length <= 0):
        raise ValueError("record_length must be positive for circular pairing")

    cdss = sorted((f for f in features if f.type == "CDS"), key=lambda f: f.start)
    pairs: list[tuple[GenBankFeature, GenBankFeature, int]] = []
    adjacent = list(pairwise(cdss))
    if circular and len(cdss) > 1:
        adjacent.append((cdss[-1], cdss[0]))

    for first, second in adjacent:
        if first.strand == second.strand and first.strand in (1, -1):
            if circular and second.feature_index == cdss[0].feature_index:
                assert record_length is not None
                gap = second.start + record_length - first.end - 1
            else:
                gap = second.start - first.end - 1
            if min_gap <= gap <= max_gap:
                pairs.append((first, second, gap))
    return pairs


def build_operon_result(
    features: list[GenBankFeature],
    max_gap: int = 150,
    min_gap: int = -50,
    *,
    circular: bool = False,
    record_length: int | None = None,
) -> OperonResult:
    """Return structured candidate pairs and connected three-CDS clusters."""
    raw_pairs = find_operon_pairs(
        features,
        max_gap=max_gap,
        min_gap=min_gap,
        circular=circular,
        record_length=record_length,
    )
    pairs = tuple(OperonPair(first, second, gap) for first, second, gap in raw_pairs)
    if not pairs:
        return OperonResult((), ())

    feature_by_index = {
        feature.feature_index: feature
        for pair in pairs
        for feature in (pair.first, pair.second)
    }
    outgoing = {pair.first.feature_index: pair for pair in pairs}
    incoming = {pair.second.feature_index for pair in pairs}
    undirected: dict[int, set[int]] = {}
    for pair in pairs:
        undirected.setdefault(pair.first.feature_index, set()).add(
            pair.second.feature_index
        )
        undirected.setdefault(pair.second.feature_index, set()).add(
            pair.first.feature_index
        )

    clusters: list[OperonCluster] = []
    visited: set[int] = set()
    for seed in sorted(undirected, key=lambda index: feature_by_index[index].start):
        if seed in visited:
            continue
        component: set[int] = set()
        stack = [seed]
        while stack:
            index = stack.pop()
            if index in component:
                continue
            component.add(index)
            stack.extend(undirected.get(index, ()))
        visited.update(component)
        if len(component) < 3:
            continue

        starts = [index for index in component if index not in incoming]
        current = min(
            starts or list(component),
            key=lambda index: feature_by_index[index].start,
        )
        ordered_indices: list[int] = []
        gaps: list[int] = []
        while current not in ordered_indices and current in component:
            ordered_indices.append(current)
            pair = outgoing.get(current)
            if pair is None or pair.second.feature_index not in component:
                break
            gaps.append(pair.gap)
            current = pair.second.feature_index
        for index in sorted(
            component - set(ordered_indices),
            key=lambda item: feature_by_index[item].start,
        ):
            ordered_indices.append(index)

        ordered_features = tuple(feature_by_index[index] for index in ordered_indices)
        clusters.append(
            OperonCluster(
                record_id=ordered_features[0].record_id,
                strand=ordered_features[0].strand or 0,
                features=ordered_features,
                gaps=tuple(gaps[: len(ordered_features) - 1]),
            )
        )
    return OperonResult(pairs, tuple(clusters))


def operon_candidates(
    filepath: str | Path,
    max_gap: int = 150,
    min_gap: int = -50,
) -> list[tuple[GenBankFeature, GenBankFeature, int]]:
    doc = read_genbank(filepath)
    all_pairs: list[tuple[GenBankFeature, GenBankFeature, int]] = []
    all_clusters: list[list[GenBankFeature]] = []

    print(
        f"-- Operon candidates (same-strand, gap between {min_gap} and {max_gap} bp) --"
    )
    print(
        f"{'Locus A':18s}  {'Gene':6s}  ->  {'Locus B':18s}  {'Gene':6s}  {'Gap':>5s}  Strand  Products"
    )
    print("-" * 110)

    for rec in doc.records:
        rec_result = build_operon_result(
            rec.features,
            max_gap=max_gap,
            min_gap=min_gap,
            circular=rec.topology == "circular",
            record_length=rec.length,
        )
        rec_pairs = [
            (pair.first, pair.second, pair.gap) for pair in rec_result.pairs
        ]
        all_pairs.extend(rec_pairs)

        for a, b, gap in rec_pairs:
            ta = a.locus_tag or "?"
            tb = b.locus_tag or "?"
            ga = a.gene or "-"
            gb = b.gene or "-"
            pa = (a.product or "-")[:30]
            pb = (b.product or "-")[:30]
            print(
                f"  {ta:18s}  {ga:6s}  ->  {tb:18s}  {gb:6s}  {gap:>4d}bp  {a.strand_symbol:>3s}    {pa} | {pb}"
            )

        all_clusters.extend([list(cluster.features) for cluster in rec_result.clusters])

    print(f"\nTotal candidate pairs: {len(all_pairs)}")

    if all_clusters:
        print("\n-- Tight clusters (>=3 consecutive co-directional genes) --")
        for cl in all_clusters:
            tags = [f.locus_tag or "?" for f in cl]
            genes = [f.gene or "-" for f in cl]
            span = f"{cl[0].start:,}..{cl[-1].end:,}  (contig: {cl[0].record_id})"
            print(f"  [{cl[0].strand_symbol}] {' -> '.join(genes)}  ({span})")
            print(f"       Tags: {', '.join(tags)}")
    else:
        print("\n-- Tight clusters: none found --")

    return all_pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Identify operon candidates: consecutive co-directional genes within gap/overlap bounds."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument(
        "max_gap",
        nargs="?",
        type=int,
        default=150,
        help="Maximum intergenic gap in bp (default: 150)",
    )
    parser.add_argument(
        "--min-gap",
        type=int,
        default=-50,
        help="Minimum gap (allows overlaps, default: -50)",
    )
    args = parser.parse_args()

    operon_candidates(args.input, max_gap=args.max_gap, min_gap=args.min_gap)


if __name__ == "__main__":
    main()
