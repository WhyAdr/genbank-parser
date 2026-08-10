"""Spatial clustering semantics for MEOR marker hits."""
from __future__ import annotations

from genbank_parser.meor.models import MeorHit
from genbank_parser.meor.spatial import cluster_meor_hits


def hit(
    feature_index: int,
    start: int,
    end: int,
    *,
    contig: str = "c1",
    strand: str = "+",
    locus: str = "NO_LOCUS",
    marker: str = "alkB",
) -> MeorHit:
    return MeorHit(
        feature_index=feature_index,
        contig=contig,
        locus_tag=locus,
        gene="",
        product="",
        start=start,
        end=end,
        strand=strand,
        category_id="cat2_medium_alkane",
        category_title="Medium alkanes",
        marker_id=marker,
        marker_name=marker,
        weight=3,
        evidence_type="gene",
        evidence_value=marker,
    )


def test_gap_boundary_overlap_and_physical_feature_deduplication() -> None:
    exact = cluster_meor_hits([hit(1, 1, 100), hit(2, 301, 400)], max_gap=200)
    assert len(exact) == 1
    assert cluster_meor_hits([hit(1, 1, 100), hit(2, 302, 400)], max_gap=200) == []

    overlapping = cluster_meor_hits([hit(1, 100, 300), hit(2, 200, 400)])
    assert len(overlapping) == 1

    duplicate_marker_hits = cluster_meor_hits(
        [hit(1, 1, 100, marker="alkB"), hit(1, 1, 100, marker="CYP153")]
    )
    assert duplicate_marker_hits == []


def test_contig_and_strand_partitioning() -> None:
    assert cluster_meor_hits([hit(1, 1, 100), hit(2, 150, 250, strand="-")]) == []
    assert cluster_meor_hits([hit(1, 1, 100), hit(2, 150, 250, contig="c2")]) == []
    assert cluster_meor_hits([hit(1, 1, 100, strand="?"), hit(2, 150, 250, strand="?")]) == []
    assert cluster_meor_hits([hit(1, 1, 100, strand="."), hit(2, 150, 250, strand=".")]) == []
