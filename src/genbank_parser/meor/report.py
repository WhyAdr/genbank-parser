"""Pure serializers for MEOR reports."""
from __future__ import annotations

import collections
import csv
import io
import json

from .models import MeorReport


TSV_COLUMNS = (
    "contig",
    "locus_tag",
    "gene",
    "start",
    "end",
    "strand",
    "category_key",
    "marker_id",
    "marker_name",
    "weight",
    "reason",
    "product",
    "feature_index",
    "evidence_type",
    "evidence_value",
)


def render_json(report: MeorReport) -> str:
    return json.dumps(report.to_dict(), indent=2) + "\n"


def render_tsv(report: MeorReport) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_COLUMNS)
    for hit in report.hits:
        payload = hit.to_dict()
        writer.writerow([payload[column] for column in TSV_COLUMNS])
    return output.getvalue()


def render_text(report: MeorReport) -> str:
    out = io.StringIO()
    print("=" * 80, file=out)
    print(" MEOR & BIOSURFACTANT DISCOVERY REPORT", file=out)
    print(f" Target File: {report.source_file}", file=out)
    print(f" Total Features Parsed: {report.total_features:,}", file=out)
    print(f" Total MEOR & Biosurfactant Hits: {report.total_hits}", file=out)
    print(f" Co-localized BGC / Operon Candidates: {report.total_clusters}", file=out)
    print("=" * 80, file=out)

    print("\n[1] CATEGORY SUMMARY & HIT COUNTS", file=out)
    print("-" * 94, file=out)
    print(
        f"{'Category':<52} | {'High (W=3)':<10} | {'Med (W=2)':<9} | "
        f"{'Low (W=1)':<9} | {'Total':<6}",
        file=out,
    )
    print("-" * 94, file=out)
    counts = collections.defaultdict(lambda: {1: 0, 2: 0, 3: 0})
    for hit in report.hits:
        counts[hit.category_id][hit.weight] += 1
    for category in report.categories:
        category_counts = counts[category.id]
        total = sum(category_counts.values())
        print(
            f"{category.title[:50]:<52} | {category_counts[3]:<10} | "
            f"{category_counts[2]:<9} | {category_counts[1]:<9} | {total:<6}",
            file=out,
        )
    print("-" * 94, file=out)

    print("\n[2] MEOR PATHWAY COMPLETENESS PROFILE", file=out)
    print("-" * 80, file=out)
    for pathway in report.pathways:
        bar_length = int(pathway.completeness_pct // 10)
        bar = "=" * bar_length + "-" * (10 - bar_length)
        print(
            f" {pathway.pathway:<45} [{bar}] {pathway.completeness_pct:>5.1f}% "
            f"({pathway.steps_present}/{pathway.total_steps} steps)",
            file=out,
        )
        for step in pathway.steps:
            mark = "+" if step["status"] == "PRESENT" else "-"
            print(f"    [{mark}] {step['step']}", file=out)
    print("-" * 80, file=out)

    max_gap = report.parameters["max_gap"]
    print(
        f"\n[3] CO-LOCALIZED MEOR OPERONS & BGC CLUSTERS (Gap <= {max_gap} bp)",
        file=out,
    )
    print("-" * 80, file=out)
    if not report.clusters:
        print(" No co-localized multi-gene MEOR clusters detected.", file=out)
    for cluster in report.clusters:
        print(
            f" Cluster ID : {cluster.cluster_id} "
            f"({cluster.contig}:{cluster.start:,}-{cluster.end:,}, "
            f"span: {cluster.span_bp:,} bp)",
            file=out,
        )
        print(f" Genes ({cluster.gene_count}): {', '.join(cluster.genes)}", file=out)
        print(f" Locus Tags : {', '.join(cluster.locus_tags)}", file=out)
        print(" Members:", file=out)
        for member in cluster.members:
            gene = f"({member.gene})" if member.gene else ""
            product = (
                member.product[:40] + ".."
                if len(member.product) > 40
                else member.product
            )
            print(
                f"   * {member.locus_tag:<14} {member.strand} "
                f"{member.start:>8d}-{member.end:<8d} "
                f"{member.marker_name[:32]:<32} {gene} [{product}]",
                file=out,
            )
        print(file=out)
    print("-" * 80, file=out)

    window_size = report.parameters["window_size"]
    print(
        f"\n[4] GENOMIC DENSITY KARYOGRAMS ({window_size / 1000:g} kb windows)",
        file=out,
    )
    print("-" * 80, file=out)
    for contig, data in report.karyograms.items():
        print(
            f" Contig: {contig} ({data['span_bp']:,} bp, {data['total_hits']} hits)",
            file=out,
        )
        print(
            f" Map: [{data['ascii_map']}] "
            f"(# = hit present in {window_size / 1000:g}kb window)",
            file=out,
        )
    print("=" * 80, file=out)
    return out.getvalue()


def serialize_report(report: MeorReport, format_type: str) -> str:
    if format_type == "json":
        return render_json(report)
    if format_type == "tsv":
        return render_tsv(report)
    if format_type == "text":
        return render_text(report)
    raise ValueError("format_type must be 'text', 'json', or 'tsv'")
