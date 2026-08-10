"""Physical-gene, known-strand spatial clustering for MEOR hits."""
from __future__ import annotations

from .models import MeorCluster, MeorHit


def _is_candidate_cluster(members: list[MeorHit]) -> bool:
    return len({member.feature_index for member in members}) >= 2


def _intervening_gap(left: MeorHit, right: MeorHit) -> int:
    return max(0, right.start - left.end - 1)


def cluster_meor_hits(
    hits: list[MeorHit] | tuple[MeorHit, ...], *, max_gap: int = 200
) -> list[MeorCluster]:
    if max_gap < 0:
        raise ValueError("max_gap must be non-negative")
    if not hits:
        return []

    sorted_hits = sorted(
        hits,
        key=lambda hit: (hit.contig, hit.start, hit.end, hit.feature_index),
    )
    raw_clusters: list[list[MeorHit]] = []
    current = [sorted_hits[0]]
    for hit in sorted_hits[1:]:
        previous = current[-1]
        if (
            hit.contig == previous.contig
            and hit.strand in {"+", "-"}
            and hit.strand == previous.strand
            and _intervening_gap(previous, hit) <= max_gap
        ):
            current.append(hit)
        else:
            if _is_candidate_cluster(current):
                raw_clusters.append(current)
            current = [hit]
    if _is_candidate_cluster(current):
        raw_clusters.append(current)

    clusters: list[MeorCluster] = []
    for index, members in enumerate(raw_clusters, 1):
        feature_indices = {member.feature_index for member in members}
        categories = tuple(sorted({member.category_id for member in members}))
        locus_tags = tuple(dict.fromkeys(member.locus_tag for member in members))
        genes = tuple(
            dict.fromkeys(member.gene or member.marker_id for member in members)
        )
        clusters.append(
            MeorCluster(
                cluster_id=f"MEOR_Cluster_{index:02d}",
                contig=members[0].contig,
                start=members[0].start,
                end=members[-1].end,
                span_bp=members[-1].end - members[0].start + 1,
                gene_count=len(feature_indices),
                hit_count=len(members),
                categories=categories,
                locus_tags=locus_tags,
                genes=genes,
                members=tuple(members),
            )
        )
    return clusters
