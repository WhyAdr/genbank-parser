"""Genome-level pathway completeness tests."""
from __future__ import annotations

from genbank_parser.meor.database import load_meor_database
from genbank_parser.meor.models import MeorHit
from genbank_parser.meor.pathways import evaluate_meor_pathways


def marker_hit(index: int, marker: str) -> MeorHit:
    return MeorHit(
        feature_index=index,
        contig="c1",
        locus_tag=f"L{index}",
        gene=marker,
        product="",
        start=index * 100,
        end=index * 100 + 50,
        strand="+",
        category_id="cat2_medium_alkane",
        category_title="Medium alkanes",
        marker_id=marker,
        marker_name=marker,
        weight=3,
        evidence_type="gene",
        evidence_value=marker,
    )


def test_pathway_zero_partial_complete_any_of_and_duplicates() -> None:
    database = load_meor_database()
    empty = evaluate_meor_pathways([], database)[0]
    assert (empty.steps_present, empty.completeness_pct) == (0, 0.0)

    partial = evaluate_meor_pathways([marker_hit(1, "CYP153")], database)[0]
    assert partial.steps_present == 1
    assert partial.completeness_pct == 16.7

    ids = ["alkB", "rubAB", "alkJ", "alkH", "alkK", "alkL"]
    complete = evaluate_meor_pathways(
        [marker_hit(index, marker) for index, marker in enumerate(ids, 1)], database
    )[0]
    assert (complete.steps_present, complete.completeness_pct) == (6, 100.0)

    duplicate = evaluate_meor_pathways(
        [marker_hit(1, "alkB"), marker_hit(2, "alkB")], database
    )[0]
    assert duplicate.steps_present == 1
