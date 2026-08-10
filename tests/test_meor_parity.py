"""Golden parity against genbank-meor v1.1.0 commit 42d0910."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

from genbank_parser.meor import analyze_meor
from genbank_parser.meor.report import render_tsv


LEGACY_HIT_KEYS = (
    "contig",
    "locus_tag",
    "gene",
    "product",
    "start",
    "end",
    "strand",
    "category_key",
    "category_title",
    "marker_id",
    "marker_name",
    "weight",
    "reason",
)


def legacy_hit(hit: dict[str, object]) -> dict[str, object]:
    return {key: hit[key] for key in LEGACY_HIT_KEYS}


def legacy_cluster(cluster: dict[str, object]) -> dict[str, object]:
    keys = (
        "cluster_id",
        "contig",
        "start",
        "end",
        "span_bp",
        "gene_count",
        "hit_count",
        "categories",
        "locus_tags",
        "genes",
    )
    result = {key: cluster[key] for key in keys}
    result["members"] = [legacy_hit(hit) for hit in cluster["members"]]
    return result


def legacy_pathway(pathway: dict[str, object]) -> dict[str, object]:
    keys = (
        "pathway",
        "completeness_pct",
        "steps_present",
        "total_steps",
        "steps",
    )
    return {key: pathway[key] for key in keys}


def test_integrated_engine_matches_legacy_biological_golden() -> None:
    fixture = Path("tests/fixtures/meor_parity.gb")
    golden = json.loads(
        Path("tests/golden/meor_v1_1_0.json").read_text(encoding="utf-8")
    )
    actual = analyze_meor(fixture).to_dict()

    assert actual["total_features"] == golden["total_features"]
    assert actual["total_hits"] == golden["total_hits"]
    assert [legacy_hit(hit) for hit in actual["hits"]] == golden["hits"]
    assert [legacy_pathway(pathway) for pathway in actual["pathways"]] == golden["pathways"]

    # The v1.1.0 golden preserves the legacy off-by-one cluster behavior.
    # Patch 6 uses inclusive GenBank coordinates: 1201 - 1000 - 1 == 200,
    # so legacy clusters 2 and 3 correctly merge at max_gap=200.
    assert golden["total_clusters"] == 3
    assert actual["total_clusters"] == 2
    assert legacy_cluster(actual["clusters"][0]) == golden["clusters"][0]
    assert actual["clusters"][1]["start"] == golden["clusters"][1]["start"]
    assert actual["clusters"][1]["end"] == golden["clusters"][2]["end"]
    assert [
        legacy_hit(member) for member in actual["clusters"][1]["members"]
    ] == golden["clusters"][1]["members"] + golden["clusters"][2]["members"]

    # Record length remains the authority; this fixture's source feature spans it.
    assert actual["karyograms"]["MEOR_CONTIG_1"]["span_bp"] == 2400
    assert golden["karyograms"]["MEOR_CONTIG_1"]["span_bp"] == 2400
    for contig in golden["karyograms"]:
        assert actual["karyograms"][contig]["total_hits"] == golden["karyograms"][contig]["total_hits"]
        assert actual["karyograms"][contig]["ascii_map"] == golden["karyograms"][contig]["ascii_map"]


def test_tsv_preserves_legacy_columns_and_rows() -> None:
    report = analyze_meor("tests/fixtures/meor_parity.gb")
    actual_rows = list(csv.reader(io.StringIO(render_tsv(report)), delimiter="\t"))
    golden_rows = list(
        csv.reader(
            io.StringIO(
                Path("tests/golden/meor_v1_1_0.tsv").read_text(encoding="utf-8")
            ),
            delimiter="\t",
        )
    )
    assert actual_rows[0][:12] == golden_rows[0]
    assert [row[:12] for row in actual_rows[1:]] == golden_rows[1:]
