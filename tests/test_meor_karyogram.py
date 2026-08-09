"""Record-length-aware karyogram tests."""
from __future__ import annotations

from Bio.Seq import Seq

from genbank_parser.meor.karyogram import generate_meor_karyograms
from genbank_parser.meor.models import MeorHit
from genbank_parser.model import GenBankRecord


def karyogram_hit(index: int, contig: str, start: int) -> MeorHit:
    return MeorHit(
        feature_index=index,
        contig=contig,
        locus_tag=str(index),
        gene="alkB",
        product="",
        start=start,
        end=start + 10,
        strand="+",
        category_id="cat2_medium_alkane",
        category_title="Medium alkanes",
        marker_id="alkB",
        marker_name="alkB",
        weight=3,
        evidence_type="gene",
        evidence_value="alkB",
    )


def test_karyogram_uses_record_length_boundaries_and_empty_contigs() -> None:
    records = [
        GenBankRecord("c1", "c1", "", Seq("N" * 120_000), 120_000),
        GenBankRecord("c2", "c2", "", Seq("N" * 60_000), 60_000),
    ]
    result = generate_meor_karyograms(
        records,
        [karyogram_hit(1, "c1", 50_000), karyogram_hit(2, "c1", 50_001)],
        window_size=50_000,
    )
    assert result["c1"] == {
        "span_bp": 120_000,
        "window_size": 50_000,
        "total_hits": 2,
        "ascii_map": "##-",
    }
    assert result["c2"]["ascii_map"] == "--"
    assert result["c2"]["total_hits"] == 0
