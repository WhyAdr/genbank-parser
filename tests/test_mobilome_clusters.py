"""C4/C5 spatial-cluster calibration tests for toxin-antitoxin inference."""

from pathlib import Path

from genbank_parser.mobilome import analyze_mobilome

FORBIDDEN = (
    "obligate co-transfer",
    "confirmed conjugative",
    "satellite plasmid",
    "theta replication",
    "rolling-circle replication",
    "confirmed phagemid",
    "post-segregational killing is active",
)


def test_tandem_toxin_antitoxin_array_collapses_into_one_aggregate_cluster() -> None:
    report = analyze_mobilome(Path("tests/fixtures/mobilome_ta_cluster.gb"))
    replicon = report.replicons[0]
    clusters = [
        hypothesis
        for hypothesis in replicon.hypotheses
        if hypothesis.rule_id == "type_iii_toxin_antitoxin_pair_candidate"
    ]

    # Six cross-product pairs (3 toxN x 2 toxI within the gap) collapse to one
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
        "retron_rt_msdna_module_candidate": (
            2,
            "h:retron_rt_msdna_module_candidate:1:10:11",
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
    assert "Retron RT-msDNA" in summaries["retron_rt_msdna_module_candidate"].summary


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
        if hypothesis.rule_id == "retron_rt_msdna_module_candidate"
    )

    assert any("effector toxin" in limitation for limitation in retron.limitations)
    assert retron.source_ids == ("bobonis-2022-retron",)
