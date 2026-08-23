"""Field-aware, multi-label evidence-scanner regressions."""

from pathlib import Path

from genbank_parser import read_genbank
from genbank_parser.mobilome import load_mobilome_database
from genbank_parser.mobilome.scanner import scan_mobilome_features


def test_scanner_retains_all_reasons_and_pseudogene_observations() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    hits = scan_mobilome_features(document.all_features, load_mobilome_database())
    by_marker = {hit.marker_id: hit for hit in hits}

    assert by_marker["toxn"].feature.gene == "toxN"
    assert by_marker["toxi"].feature.feature_type == "ncRNA"
    assert by_marker["crispr_array"].reasons[0].field == "note"
    assert {reason.field for reason in by_marker["crispr_array"].reasons} == {
        "note",
        "rpt_type",
    }
    pseudo_reps = [
        hit
        for hit in hits
        if hit.marker_id == "replication_initiation_candidate" and hit.feature.is_pseudo
    ]
    assert len(pseudo_reps) == 3
    assert not any(hit.eligible_for_inference for hit in pseudo_reps)


def test_scanner_preserves_compound_origin_and_multifacet_candidates() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    hits = scan_mobilome_features(document.all_features, load_mobilome_database())
    cas = next(hit for hit in hits if hit.marker_id == "crispr_cas")
    zot = next(hit for hit in hits if hit.marker_id == "zot_like_candidate")

    assert cas.feature.segments == ((571, 600), (1, 30))
    assert cas.feature.wraps_origin
    assert zot.facets == ("phage_regulation", "vf_candidate")
    assert zot.reasons[0].strength == 1


def test_threshold_changes_eligibility_without_erasing_raw_hits() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    database = load_mobilome_database()
    low = scan_mobilome_features(document.all_features, database, min_evidence=1)
    high = scan_mobilome_features(document.all_features, database, min_evidence=3)

    assert [hit.hit_id for hit in low] == [hit.hit_id for hit in high]
    assert sum(hit.eligible_for_inference for hit in low) > sum(
        hit.eligible_for_inference for hit in high
    )
