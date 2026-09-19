"""C4/C5 spatial-cluster calibration tests for toxin-antitoxin inference."""

from dataclasses import replace
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from genbank_parser.mobilome import analyze_mobilome, load_mobilome_database
from genbank_parser.mobilome.inference import infer_replicon_hypotheses
from genbank_parser.mobilome.replicons import inventory_replicons
from genbank_parser.mobilome.scanner import scan_mobilome_features
from genbank_parser.model import GenBankDocument, GenBankFeature, GenBankRecord

FORBIDDEN = (
    "obligate co-transfer",
    "confirmed conjugative",
    "satellite plasmid",
    "theta replication",
    "rolling-circle replication",
    "confirmed phagemid",
    "post-segregational killing is active",
)


def _synthetic_features(
    specs: tuple[tuple[str, int, int, dict[str, list[str]]], ...],
    *,
    topology: str = "linear",
    length: int = 10_000,
    strands: tuple[int | None, ...] | None = None,
    return_features: bool = False,
):
    features = [
        GenBankFeature(
            record_id="synthetic.1",
            record_index=1,
            feature_index=index,
            type=feature_type,
            location=FeatureLocation(
                start - 1,
                end,
                strand=1 if strands is None else strands[index - 1],
            ),
            qualifiers=qualifiers,
            record_length=length,
            topology=topology,
        )
        for index, (feature_type, start, end, qualifiers) in enumerate(specs, 1)
    ]
    record = GenBankRecord(
        id="synthetic.1",
        name="synthetic",
        description="Synthetic compact-module regression record.",
        seq=Seq("A" * length),
        length=length,
        topology=topology,
        features=features,
    )
    document = GenBankDocument(path=Path("synthetic.gb"), records=[record])
    inventory = inventory_replicons(document)[0]
    database = load_mobilome_database()
    hits = scan_mobilome_features(features, database)
    if return_features:
        return inventory, hits, database, tuple(features)
    return inventory, hits, database


def test_tandem_toxin_antitoxin_array_collapses_into_one_aggregate_cluster() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_ta_cluster.gb"))
    replicon = report.replicons[0]
    clusters = [
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "type_iii_toxin_antitoxin_pair_candidate"
    ]

    # Nine nearby cross-product pairs (3 toxN x 3 toxI within the gap) collapse to one
    # aggregate cluster hypothesis instead of the previous combinatorial list.
    assert len(clusters) == 1
    cluster = clusters[0]
    assert len(cluster.supporting_hit_ids) == 6
    assert (
        cluster.hypothesis_id
        == "h:type_iii_toxin_antitoxin_pair_candidate:1:2:3:4:5:6:7"
    )
    assert cluster.status == "tentative"
    assert any("observed cluster span" in item for item in cluster.limitations)
    assert any("Aggregate cluster" in item for item in cluster.limitations)
    assert "2000" in "".join(cluster.limitations)


def test_isolated_toxin_hit_without_cognate_emits_no_hypothesis() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_ta_cluster.gb"))
    replicon = report.replicons[0]

    toxn_hits = [hit for hit in replicon.hits if hit.marker_id == "toxn"]
    assert len(toxn_hits) == 4
    clusters = [
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "type_iii_toxin_antitoxin_pair_candidate"
    ]
    # The isolated toxN (no toxI within the configured gap) is retained as a
    # hit but never joined to the aggregate cluster.
    isolated = next(hit for hit in toxn_hits if hit.feature.feature_index == 8)
    assert isolated.hit_id not in clusters[0].supporting_hit_ids
    assert clusters[0].hypothesis_id.endswith("2:3:4:5:6:7")


def test_separated_pairs_keep_one_hypothesis_each() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_ta_cluster.gb"))
    replicon = report.replicons[0]
    pairs = [
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "hepn_mnt_annotation_pair_candidate"
    ]

    assert sorted(pair.hypothesis_id for pair in pairs) == [
        "h:hepn_mnt_annotation_pair_candidate:1:11:12",
        "h:hepn_mnt_annotation_pair_candidate:1:9:10",
    ]
    for pair in pairs:
        assert len(pair.supporting_hit_ids) == 2
        assert any("observed gap:" in item for item in pair.limitations)
        assert not any("Aggregate cluster" in item for item in pair.limitations)


def test_tripartite_module_rules_fire_with_required_markers() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_tripartite.gb"))
    replicon = report.replicons[0]
    summaries = {hypothesis.rule_id: hypothesis for hypothesis in replicon.hypotheses}

    expected = {
        "omega_epsilon_zeta_module_candidate": (
            3,
            "h:omega_epsilon_zeta_module_candidate:1:2:3:4",
        ),
        "mqsrac_module_candidate": (3, "h:mqsrac_module_candidate:1:7:8:9"),
        "higba_tac_module_candidate": (3, "h:higba_tac_module_candidate:1:16:17:18"),
        "retron_core_annotation_candidate": (
            2,
            "h:retron_core_annotation_candidate:1:10:11",
        ),
        "type_iii_tenpin_annotation_pair_candidate": (
            2,
            "h:type_iii_tenpin_annotation_pair_candidate:1:12:13",
        ),
        "type_iii_cptin_annotation_pair_candidate": (
            2,
            "h:type_iii_cptin_annotation_pair_candidate:1:14:15",
        ),
    }
    for rule_id, (support, hypothesis_id) in expected.items():
        hypothesis = summaries[rule_id]
        assert len(hypothesis.supporting_hit_ids) == support, rule_id
        assert hypothesis.hypothesis_id == hypothesis_id, rule_id
        assert hypothesis.status == "tentative", rule_id
    assert "Omega-epsilon-zeta annotation-module candidate" in (
        summaries["omega_epsilon_zeta_module_candidate"].summary
    )
    assert "MqsRAC" in summaries["mqsrac_module_candidate"].summary
    assert "TAC" in summaries["higba_tac_module_candidate"].summary
    assert "TenpIN" in summaries["type_iii_tenpin_annotation_pair_candidate"].summary
    assert "CptIN" in summaries["type_iii_cptin_annotation_pair_candidate"].summary
    assert "Retron RT-msr/msd" in summaries["retron_core_annotation_candidate"].summary


def test_incomplete_tripartite_cluster_emits_no_module_hypothesis() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_tripartite.gb"))
    replicon = report.replicons[0]

    # TRI_0004/TRI_0005 (mqsR + mqsA without mqsC, beyond the configured gap
    # from the complete trio) are retained as hits but produce no module
    # hypothesis: the rule requires all three markers.
    mqs_hits = {
        hit.marker_id
        for hit in replicon.hits
        if hit.marker_id in {"mqsr", "mqsa", "mqsc"}
    }
    assert mqs_hits == {"mqsr", "mqsa", "mqsc"}
    mqsrac = [
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "mqsrac_module_candidate"
    ]
    assert len(mqsrac) == 1
    assert len(mqsrac[0].supporting_hit_ids) == 3
    assert mqsrac[0].hypothesis_id == "h:mqsrac_module_candidate:1:7:8:9"


def test_tripartite_adversarial_negatives_never_match_markers() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_tripartite.gb"))
    replicon = report.replicons[0]
    adversarial_locus_tags = {"TRI_0018", "TRI_0019", "TRI_0020"}

    matched = {
        hit.feature.locus_tag
        for hit in replicon.hits
        if hit.feature.locus_tag in adversarial_locus_tags
    }
    assert matched == set()
    summaries = " ".join(
        hypothesis.summary for hypothesis in replicon.hypotheses
    ).casefold()
    assert all(forbidden not in summaries for forbidden in FORBIDDEN)


def test_retron_rule_documents_the_system_specific_effector_limit() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_tripartite.gb"))
    replicon = report.replicons[0]
    retron = next(
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "retron_core_annotation_candidate"
    )

    assert retron.kind == "retron"
    assert any(
        "not a toxin-antitoxin call" in limitation for limitation in retron.limitations
    )
    assert retron.source_ids == ("bobonis-2022-retron",)


def test_tripartite_collision_requires_three_physical_features() -> None:
    inventory, hits, database = _synthetic_features(
        (
            (
                "CDS",
                100,
                200,
                {"product": ["MqsR MqsA fusion protein"]},
            ),
            ("CDS", 300, 400, {"gene": ["mqsC"]}),
        )
    )

    hypotheses = infer_replicon_hypotheses(inventory, hits, database)

    assert not any(item.rule_id == "mqsrac_module_candidate" for item in hypotheses)


def test_pair_rule_also_requires_distinct_physical_features() -> None:
    inventory, hits, database = _synthetic_features(
        (
            (
                "CDS",
                100,
                200,
                {"gene": ["toxN"], "product": ["ToxN ToxI fusion annotation"]},
            ),
        )
    )

    hypotheses = infer_replicon_hypotheses(inventory, hits, database)

    assert not any(
        item.rule_id == "type_iii_toxin_antitoxin_pair_candidate" for item in hypotheses
    )


def test_tripartite_bridge_over_global_extent_is_rejected() -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsR"]}),
            ("CDS", 2601, 2700, {"gene": ["mqsA"]}),
            ("CDS", 5101, 5200, {"gene": ["mqsC"]}),
        )
    )

    hypotheses = infer_replicon_hypotheses(inventory, hits, database)

    assert not any(item.rule_id == "mqsrac_module_candidate" for item in hypotheses)


def test_compact_submodule_survives_an_unrelated_single_linkage_bridge() -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsR"]}),
            ("CDS", 2601, 2700, {"gene": ["mqsA"]}),
            ("CDS", 5101, 5200, {"gene": ["mqsC"]}),
            ("CDS", 6000, 6100, {"gene": ["mqsR"]}),
            ("CDS", 6200, 6300, {"gene": ["mqsA"]}),
            ("CDS", 6400, 6500, {"gene": ["mqsC"]}),
        )
    )

    hypotheses = [
        item
        for item in infer_replicon_hypotheses(inventory, hits, database)
        if item.rule_id == "mqsrac_module_candidate"
    ]

    assert len(hypotheses) == 1
    assert {
        hit_id.split(":f")[1].split(":")[0]
        for hit_id in hypotheses[0].supporting_hit_ids
    } == {
        "4",
        "5",
        "6",
    }


def test_ambiguous_marker_hits_use_a_deterministic_full_matching() -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 100, 200, {"product": ["MqsR MqsA fusion protein"]}),
            ("CDS", 300, 400, {"product": ["MqsA MqsC fusion protein"]}),
            ("CDS", 500, 600, {"gene": ["mqsC"]}),
        )
    )

    hypotheses = [
        item
        for item in infer_replicon_hypotheses(inventory, hits, database)
        if item.rule_id == "mqsrac_module_candidate"
    ]

    assert len(hypotheses) == 1
    assert hypotheses[0].supporting_hit_ids == (
        "r1:f1:mmqsr",
        "r1:f2:mmqsa",
        "r1:f3:mmqsc",
    )


def test_two_complete_compact_modules_remain_two_selected_modules() -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsR"]}),
            ("CDS", 300, 400, {"gene": ["mqsA"]}),
            ("CDS", 500, 600, {"gene": ["mqsC"]}),
            ("CDS", 1000, 1100, {"gene": ["mqsR"]}),
            ("CDS", 1200, 1300, {"gene": ["mqsA"]}),
            ("CDS", 1400, 1500, {"gene": ["mqsC"]}),
        )
    )

    hypotheses = [
        item
        for item in infer_replicon_hypotheses(inventory, hits, database)
        if item.rule_id == "mqsrac_module_candidate"
    ]

    assert len(hypotheses) == 2
    assert [item.supporting_hit_ids for item in hypotheses] == [
        ("r1:f1:mmqsr", "r1:f2:mmqsa", "r1:f3:mmqsc"),
        ("r1:f4:mmqsr", "r1:f5:mmqsa", "r1:f6:mmqsc"),
    ]


def test_compact_module_span_accepts_exact_boundary_and_rejects_boundary_plus_one() -> (
    None
):
    exact_inventory, exact_hits, exact_database = _synthetic_features(
        (
            ("CDS", 100, 100, {"gene": ["mqsR"]}),
            ("CDS", 1100, 1100, {"gene": ["mqsA"]}),
            ("CDS", 3099, 3099, {"gene": ["mqsC"]}),
        )
    )
    over_inventory, over_hits, over_database = _synthetic_features(
        (
            ("CDS", 100, 100, {"gene": ["mqsR"]}),
            ("CDS", 1100, 1100, {"gene": ["mqsA"]}),
            ("CDS", 3100, 3100, {"gene": ["mqsC"]}),
        )
    )

    exact = infer_replicon_hypotheses(exact_inventory, exact_hits, exact_database)
    over = infer_replicon_hypotheses(over_inventory, over_hits, over_database)

    assert any(item.rule_id == "mqsrac_module_candidate" for item in exact)
    assert not any(item.rule_id == "mqsrac_module_candidate" for item in over)


def test_origin_straddling_compact_module_uses_circular_covering_arc() -> None:
    inventory, hits, database = _synthetic_features(
        (
            ("CDS", 9900, 9950, {"gene": ["mqsR"]}),
            ("CDS", 1, 50, {"gene": ["mqsA"]}),
            ("CDS", 100, 150, {"gene": ["mqsC"]}),
        ),
        topology="circular",
        length=10_000,
    )

    hypotheses = infer_replicon_hypotheses(inventory, hits, database)

    assert (
        len([item for item in hypotheses if item.rule_id == "mqsrac_module_candidate"])
        == 1
    )


def test_structural_policies_are_applied_to_the_selected_assignment() -> None:
    inventory, hits, database, features = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsR"]}),
            ("CDS", 300, 400, {"gene": ["mqsA"]}),
            ("CDS", 500, 600, {"gene": ["mqsC"]}),
        ),
        return_features=True,
    )
    base_rule = database.inference_rule_map["mqsrac_module_candidate"]
    structural_rule = replace(
        base_rule,
        strand_policy="same",
        allowed_orders=(("mqsr", "mqsa", "mqsc"),),
        max_intervening_features=0,
        intervening_feature_types=("CDS",),
    )
    database = replace(
        database,
        inference_rules=tuple(
            structural_rule if rule.id == structural_rule.id else rule
            for rule in database.inference_rules
        ),
    )

    hypotheses = infer_replicon_hypotheses(
        inventory,
        hits,
        database,
        features=features,
    )
    module = next(
        item for item in hypotheses if item.rule_id == "mqsrac_module_candidate"
    )

    assert any("strand policy (same) passed" in item for item in module.limitations)
    assert any("biological-order policy passed" in item for item in module.limitations)
    assert any(
        "intervening-feature policy passed" in item for item in module.limitations
    )


def test_unknown_strand_fails_a_non_any_structural_policy() -> None:
    inventory, hits, database, features = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsR"]}),
            ("CDS", 300, 400, {"gene": ["mqsA"]}),
            ("CDS", 500, 600, {"gene": ["mqsC"]}),
        ),
        strands=(1, None, 1),
        return_features=True,
    )
    base_rule = database.inference_rule_map["mqsrac_module_candidate"]
    structural_rule = replace(base_rule, strand_policy="same")
    database = replace(
        database,
        inference_rules=tuple(
            structural_rule if rule.id == structural_rule.id else rule
            for rule in database.inference_rules
        ),
    )

    hypotheses = infer_replicon_hypotheses(
        inventory,
        hits,
        database,
        features=features,
    )

    assert not any(item.rule_id == "mqsrac_module_candidate" for item in hypotheses)


def test_reverse_strand_order_and_interveners_are_normalized() -> None:
    inventory, hits, database, features = _synthetic_features(
        (
            ("CDS", 100, 200, {"gene": ["mqsC"]}),
            ("CDS", 250, 260, {"gene": ["unrelated"]}),
            ("CDS", 300, 400, {"gene": ["mqsA"]}),
            ("CDS", 500, 600, {"gene": ["mqsR"]}),
        ),
        strands=(-1, -1, -1, -1),
        return_features=True,
    )
    base_rule = database.inference_rule_map["mqsrac_module_candidate"]
    structural_rule = replace(
        base_rule,
        strand_policy="same",
        allowed_orders=(("mqsr", "mqsa", "mqsc"),),
        max_intervening_features=1,
        intervening_feature_types=("CDS",),
    )
    database = replace(
        database,
        inference_rules=tuple(
            structural_rule if rule.id == structural_rule.id else rule
            for rule in database.inference_rules
        ),
    )

    hypotheses = infer_replicon_hypotheses(
        inventory,
        hits,
        database,
        features=features,
    )
    module = next(
        item for item in hypotheses if item.rule_id == "mqsrac_module_candidate"
    )

    assert any("biological-order policy passed" in item for item in module.limitations)
    assert any("maximum observed: 1" in item for item in module.limitations)
