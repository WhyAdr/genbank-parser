"""Field-aware, multi-label evidence-scanner regressions."""

from dataclasses import replace
from pathlib import Path

from genbank_parser import read_genbank
from genbank_parser.mobilome import load_mobilome_database
from genbank_parser.mobilome.models import MarkerMatcher, MobilomeMarker
from genbank_parser.mobilome.scanner import _match_value, scan_mobilome_features
from genbank_parser.model import GenBankFeature


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


def test_annotation_candidates_preserve_provenance_without_aggregate_calls() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    database = load_mobilome_database()
    hits = scan_mobilome_features(document.all_features, database)
    by_marker = {hit.marker_id: hit for hit in hits}

    zot = by_marker["zot_like_candidate"]
    assert zot.facets == ("phage_regulation", "vf_candidate")
    assert all(reason.source_ids for reason in zot.reasons)
    assert not {item.rule_id for item in database.disabled_aggregate_rules} & set(
        database.inference_rule_map
    )


def test_matching_modes_and_negative_patterns_are_field_local() -> None:
    assert (
        _match_value("RepA", mode="exact", pattern="RepA", negative_patterns=())
        == "RepA"
    )
    assert (
        _match_value("repa", mode="exact", pattern="RepA", negative_patterns=()) is None
    )
    assert (
        _match_value(
            "repa", mode="casefold_exact", pattern="RepA", negative_patterns=()
        )
        == "repa"
    )
    assert (
        _match_value(
            "replication initiator",
            mode="regex",
            pattern=r"\binitiator\b",
            negative_patterns=(),
        )
        == "initiator"
    )
    assert (
        _match_value(
            "transcriptional repressor",
            mode="regex",
            pattern=r"\brepressor\b",
            negative_patterns=(r"\btranscriptional repressor\b",),
        )
        is None
    )


def test_multi_value_qualifiers_remain_traceable_without_inference_inflation() -> None:
    database = load_mobilome_database()
    feature = GenBankFeature(
        record_id="synthetic",
        record_index=1,
        feature_index=1,
        type="CDS",
        location=None,
        qualifiers={
            "gene": ["toxN", "toxN"],
            "product": ["ToxN-family toxin", "unrelated"],
        },
        record_length=100,
    )

    hits = scan_mobilome_features([feature], database)
    toxn = next(hit for hit in hits if hit.marker_id == "toxn")

    assert {reason.field for reason in toxn.reasons} == {"gene", "product"}
    assert len({reason.reason_id for reason in toxn.reasons}) == len(toxn.reasons)
    assert toxn.feature.start is None
    assert toxn.feature.end is None


def test_scanner_keeps_dbxref_product_and_field_scoped_reasons_distinct() -> None:
    database = load_mobilome_database()
    marker = MobilomeMarker(
        id="scanner_contract_marker",
        label="Scanner contract marker",
        facets=("amr_candidate",),
        feature_types=("CDS",),
        matchers=(
            MarkerMatcher(
                field="db_xref",
                mode="regex",
                patterns=(r"AMRFinder",),
                negative_patterns=(),
                strength=2,
            ),
            MarkerMatcher(
                field="product",
                mode="regex",
                patterns=(r"efflux pump", r"multidrug"),
                negative_patterns=(),
                strength=1,
            ),
            MarkerMatcher(
                field="gene",
                mode="casefold_exact",
                patterns=("shared",),
                negative_patterns=(),
                strength=2,
            ),
        ),
        requires_non_pseudo=False,
        requires_complete=False,
        sources=("scanner-test",),
        interpretation="Test-only marker.",
        limitations=("Test-only evidence.",),
    )
    database = replace(database, markers=database.markers + (marker,))
    feature = GenBankFeature(
        record_id="synthetic",
        record_index=1,
        feature_index=1,
        type="CDS",
        location=None,
        qualifiers={
            "gene": ["shared"],
            "product": ["shared efflux pump multidrug"],
            "db_xref": ["AMRFinder:shared"],
        },
        record_length=100,
    )

    hit = next(
        item
        for item in scan_mobilome_features([feature], database)
        if item.marker_id == "scanner_contract_marker"
    )

    assert {reason.field for reason in hit.reasons} == {"db_xref", "product", "gene"}
    assert len([reason for reason in hit.reasons if reason.field == "product"]) == 2
    assert len({reason.reason_id for reason in hit.reasons}) == len(hit.reasons)
    assert hit.max_strength == 2
    assert hit.eligible_for_inference
