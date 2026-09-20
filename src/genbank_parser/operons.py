"""Identify operon candidates: consecutive co-directional genes within gap/overlap bounds."""

from __future__ import annotations

import argparse
import csv
import io
import json
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from .io import read_genbank
from .model import GenBankFeature

SCHEMA_VERSION = "gbparse.operons.v1"
TSV_COLUMNS = (
    "row_type",
    "record",
    "candidate_id",
    "strand",
    "first_feature_index",
    "second_feature_index",
    "feature_index",
    "locus_tag",
    "gene",
    "start",
    "end",
    "gap",
    "candidate",
)


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


def _feature_dict(feature: GenBankFeature) -> dict[str, object]:
    return {
        "record": feature.record_id,
        "record_index": feature.record_index,
        "feature_index": feature.feature_index,
        "locus_tag": feature.locus_tag or None,
        "gene": feature.gene or None,
        "product": feature.product or None,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand_symbol,
        "strand_value": feature.strand,
        "length": feature.length,
        "segments": [
            {"start": start, "end": end}
            for start, end in (feature.join_segments if feature.is_compound else [(feature.start, feature.end)])
        ],
    }


def build_operon_report(
    document: object,
    *,
    max_gap: int = 150,
    min_gap: int = -50,
    min_genes: int = 3,
) -> dict[str, object]:
    """Build a complete record-level candidate-operon report."""

    if min_genes < 2:
        raise ValueError("min_genes must be at least 2")
    records: list[dict[str, object]] = []
    total_pairs = 0
    total_clusters = 0
    for record in getattr(document, "records", ()):
        result = build_operon_result(
            record.features,
            max_gap=max_gap,
            min_gap=min_gap,
            circular=record.topology == "circular",
            record_length=record.length,
            min_genes=min_genes,
        )
        clusters = [
            cluster
            for cluster in result.clusters
            if len(cluster.features) >= min_genes
        ]
        pair_payload = [
            {
                "first": _feature_dict(pair.first),
                "second": _feature_dict(pair.second),
                "gap": pair.gap,
                "candidate": True,
            }
            for pair in result.pairs
        ]
        cluster_payload = [
            {
                "candidate_id": f"{record.id}:candidate_operon_{index}",
                "record": record.id,
                "strand": cluster.strand,
                "strand_symbol": "+" if cluster.strand == 1 else "-",
                "gaps": list(cluster.gaps),
                "features": [_feature_dict(feature) for feature in cluster.features],
                "candidate": True,
            }
            for index, cluster in enumerate(clusters, 1)
        ]
        total_pairs += len(pair_payload)
        total_clusters += len(cluster_payload)
        records.append(
            {
                "record": record.id,
                "record_index": record.features[0].record_index if record.features else None,
                "length": record.length,
                "topology": record.topology or "linear",
                "pairs": pair_payload,
                "clusters": cluster_payload,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "parameters": {
            "max_gap": max_gap,
            "min_gap": min_gap,
            "min_genes": min_genes,
        },
        "record_count": len(records),
        "pair_count": total_pairs,
        "cluster_count": total_clusters,
        "records": records,
        "interpretation": "Proximity-based operon candidates; not evidence of transcription or co-expression.",
    }


def render_operons(report: dict[str, object], format_type: str = "text") -> str:
    """Render a candidate-operon report without performing I/O."""

    if format_type == "json":
        return json.dumps(report, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2) + "\n"
    if format_type == "tsv":
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in report["records"]:  # type: ignore[index]
            for pair in record["pairs"]:  # type: ignore[index]
                first = pair["first"]
                second = pair["second"]
                writer.writerow(
                    {
                        "row_type": "pair",
                        "record": record["record"],
                        "candidate_id": "",
                        "strand": first["strand"],
                        "first_feature_index": first["feature_index"],
                        "second_feature_index": second["feature_index"],
                        "feature_index": "",
                        "locus_tag": f"{first['locus_tag'] or ''};{second['locus_tag'] or ''}",
                        "gene": f"{first['gene'] or ''};{second['gene'] or ''}",
                        "start": min(first["start"], second["start"]),
                        "end": max(first["end"], second["end"]),
                        "gap": pair["gap"],
                        "candidate": "true",
                    }
                )
            for cluster in record["clusters"]:  # type: ignore[index]
                for feature in cluster["features"]:
                    writer.writerow(
                        {
                            "row_type": "cluster_member",
                            "record": record["record"],
                            "candidate_id": cluster["candidate_id"],
                            "strand": cluster["strand_symbol"],
                            "first_feature_index": "",
                            "second_feature_index": "",
                            "feature_index": feature["feature_index"],
                            "locus_tag": feature["locus_tag"] or "",
                            "gene": feature["gene"] or "",
                            "start": feature["start"],
                            "end": feature["end"],
                            "gap": "",
                            "candidate": "true",
                        }
                    )
        return output.getvalue()
    if format_type != "text":
        raise ValueError("operon format must be text, tsv, json, or gff3")
    lines = [
        "=" * 80,
        "  OPERON PROXIMITY CANDIDATES",
        "=" * 80,
        f"  Candidate pairs   : {report['pair_count']}",
        f"  Candidate clusters: {report['cluster_count']}",
        "  Interpretation    : proximity candidates, not evidence of transcription or co-expression.",
        "",
    ]
    for record in report["records"]:  # type: ignore[index]
        lines.append(f"-- Record {record['record']} --")
        for pair in record["pairs"]:  # type: ignore[index]
            first = pair["first"]
            second = pair["second"]
            lines.append(
                f"  candidate pair {first['locus_tag'] or '?'} -> {second['locus_tag'] or '?'} "
                f"gap={pair['gap']} strand={first['strand']}"
            )
        for cluster in record["clusters"]:  # type: ignore[index]
            tags = ", ".join(feature["locus_tag"] or "?" for feature in cluster["features"])
            lines.append(f"  candidate cluster {cluster['candidate_id']}: {tags}")
    lines.append("=" * 80)
    return "\n".join(lines) + "\n"


def render_operons_gff3(report: dict[str, object]) -> str:
    """Render candidate clusters as GFF3 parent/child annotations."""

    lines = ["##gff-version 3"]
    for record in report["records"]:  # type: ignore[index]
        lines.append(f"##sequence-region {record['record']} 1 {record['length']}")
        clusters = record["clusters"]  # type: ignore[index]
        if not clusters:
            continue
        for cluster in clusters:
            features = cluster["features"]
            start = min(feature["start"] for feature in features)
            end = max(feature["end"] for feature in features)
            candidate_id = str(cluster["candidate_id"])
            strand = str(cluster["strand_symbol"])
            lines.append(
                "\t".join(
                    (
                        str(record["record"]),
                        "gbparse",
                        "operon_candidate",
                        str(start),
                        str(end),
                        ".",
                        strand,
                        ".",
                        f"ID={candidate_id};candidate=true",
                    )
                )
            )
            for feature in features:
                child_id = f"{candidate_id}:feature:{feature['feature_index']}"
                attrs = f"ID={child_id};Parent={candidate_id};locus_tag={feature['locus_tag'] or ''};candidate=true"
                lines.append(
                    "\t".join(
                        (
                            str(record["record"]),
                            "gbparse",
                            "CDS",
                            str(feature["start"]),
                            str(feature["end"]),
                            ".",
                            str(feature["strand"]),
                            ".",
                            attrs,
                        )
                    )
                )
    return "\n".join(lines) + "\n"


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

    for genomic_first, genomic_second in adjacent:
        if genomic_first.strand == genomic_second.strand and genomic_first.strand in (1, -1):
            is_origin_edge = circular and genomic_second.feature_index == cdss[0].feature_index
            if is_origin_edge:
                assert record_length is not None
                gap = genomic_second.start + record_length - genomic_first.end - 1
            else:
                # The two intervals are sorted by genomic coordinate.  The
                # physical gap is therefore the same on either strand; only
                # the reported pair order changes for reverse-strand CDSs.
                gap = genomic_second.start - genomic_first.end - 1
            if min_gap <= gap <= max_gap:
                if genomic_first.strand == -1:
                    # A reverse-strand candidate is reported in biological
                    # 5'->3' order.  The gap is symmetric for the two
                    # adjacent intervals, including the circular edge.
                    pairs.append((genomic_second, genomic_first, gap))
                else:
                    pairs.append((genomic_first, genomic_second, gap))
    return pairs


def build_operon_result(
    features: list[GenBankFeature],
    max_gap: int = 150,
    min_gap: int = -50,
    *,
    circular: bool = False,
    record_length: int | None = None,
    min_genes: int = 3,
) -> OperonResult:
    """Return structured candidate pairs and connected CDS clusters."""
    if min_genes < 2:
        raise ValueError("min_genes must be at least 2")
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
        if len(component) < min_genes:
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
        ordered_indices.extend(
            sorted(
                component - set(ordered_indices),
                key=lambda item: feature_by_index[item].start,
            )
        )

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
