"""Tests for fail-closed mobilome evidence-resource loading."""

import json
from collections.abc import Callable
from copy import deepcopy
from importlib import resources
from pathlib import Path

import pytest
import yaml

from genbank_parser.mobilome import load_mobilome_database
from genbank_parser.mobilome.database import MobilomeDatabaseError


def _copy_custom_database(destination: Path) -> None:
    source = resources.files("genbank_parser").joinpath("data", "mobilome")
    for name in ("markers.yaml", "provenance.yaml", "inference.yaml"):
        (destination / name).write_bytes(source.joinpath(name).read_bytes())


def test_packaged_database_is_versioned_hashed_and_schema_valid() -> None:
    database = load_mobilome_database()
    schema_path = resources.files("genbank_parser").joinpath(
        "data", "mobilome", "report.schema.json"
    )
    schema = json.loads(schema_path.read_text(encoding="utf-8"))

    assert database.catalog_version == "1.4.0"
    assert database.inference_version == "1.4.0"
    assert database.provenance_version == "1.4.0"
    assert database.database_source == "packaged"
    assert database.source_paths == ()
    assert [resource.name for resource in database.resources] == [
        "markers.yaml",
        "provenance.yaml",
        "inference.yaml",
        "report.schema.json",
    ]
    assert all(len(resource.sha256) == 64 for resource in database.resources)
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"


def test_provenance_sources_use_verified_citations() -> None:
    database = load_mobilome_database()
    by_id = {source.id: source for source in database.provenance_sources}

    assert "anjum-2014-pxo1" not in by_id
    akhtar = by_id["akhtar-2012-pxo1"]
    assert akhtar.payload["year"] == 2012
    assert akhtar.payload["title"].startswith("Two independent replicons")
    assert (
        "MOB-suite: software tools"
        in by_id["robertson-2018-mob-suite"].payload["title"]
    )
    assert (
        by_id["yao-2020-hepn-mnt"]
        .payload["title"]
        .startswith("Novel polyadenylylation")
    )
    assert (
        by_id["blower-2012-toxin"]
        .payload["title"]
        .endswith("encoded in chromosomal and plasmid genomes")
    )
    for anchor in (
        "garcillan-barcia-2009-relaxases",
        "christie-2025-t4ss",
        "makarova-2020-crispr",
        "siguier-2006-isfinder",
        "bouet-2019-partition",
        "schwengers-2021-bakta",
        "garcillan-barcia-2025-extended-mobility",
        "johnson-2015-ice",
        "mazel-2006-integrons",
        "qiu-2022-ta-classification",
        "yao-2026-mnt-hepn",
        "liu-2022-vfdb",
        "partridge-2018-mge-amr",
        "arndt-2016-phaster",
        "camargo-2024-genomad",
        "xie-2017-isescan",
        "johansson-2021-mobileelementfinder",
        "anand-2008-repx",
        "tinsley-2006-repx",
        "makarova-2025-crispr-classification",
    ):
        assert anchor in by_id
    assert by_id["schwengers-2021-bakta"].payload["pmid"] == "34739369"
    assert by_id["anand-2008-repx"].payload["pmid"] == "18179418"
    assert by_id["tinsley-2006-repx"].payload["pmid"] == "16585744"
    assert by_id["christie-2025-t4ss"].payload["pmid"] == "41474020"
    assert by_id["makarova-2025-crispr-classification"].payload["pmid"] == "41198952"
    assert (
        by_id["garcillan-barcia-2025-extended-mobility"].payload["pmid"] == "40694848"
    )
    assert by_id["mazel-2006-integrons"].payload["doi"] == "10.1038/nrmicro1462"
    assert (
        by_id["johansson-2021-mobileelementfinder"]
        .payload["title"]
        .endswith("web tool: MobileElementFinder")
    )


def test_amr_and_vf_citations_are_anchored_to_their_own_markers() -> None:
    database = load_mobilome_database()
    markers = database.marker_map

    assert markers["generic_amr_candidate"].sources == (
        "feldgarden-2021-amrfinderplus",
        "local-annotation-candidate-policy",
        "partridge-2018-mge-amr",
    )
    assert markers["generic_vf_candidate"].sources == (
        "chen-2005-vfdb",
        "liu-2022-vfdb",
        "local-annotation-candidate-policy",
    )
    # The AMRFinderPlus and VFDB citations must never support phage markers.
    assert markers["phage_packaging_candidate"].sources == (
        "local-annotation-candidate-policy",
    )
    assert markers["phage_structure_candidate"].sources == (
        "local-annotation-candidate-policy",
    )


def test_ta_and_mobility_rules_carry_their_configured_sources() -> None:
    database = load_mobilome_database()
    rules = database.inference_rule_map

    assert rules["hepn_mnt_annotation_pair_candidate"].sources == (
        "qiu-2022-ta-classification",
        "yao-2020-hepn-mnt",
        "yao-2026-mnt-hepn",
    )
    assert rules["type_iii_toxin_antitoxin_pair_candidate"].sources == (
        "blower-2012-toxin",
        "qiu-2022-ta-classification",
    )
    assert rules["mobilizable_core_annotation_candidate"].sources == (
        "garcillan-barcia-2025-extended-mobility",
        "smillie-2010-mobility",
    )
    assert any(
        "Type VII" in limitation
        for limitation in rules["hepn_mnt_annotation_pair_candidate"].limitations
    )


def test_complete_custom_database_is_loaded_as_one_resource_set(tmp_path: Path) -> None:
    _copy_custom_database(tmp_path)

    database = load_mobilome_database(tmp_path)

    assert database.database_source == "custom"
    assert tuple(path.name for path in database.source_paths) == (
        "markers.yaml",
        "provenance.yaml",
        "inference.yaml",
    )


def test_database_rejects_incomplete_or_invalid_custom_resources(
    tmp_path: Path,
) -> None:
    _copy_custom_database(tmp_path)
    (tmp_path / "extra.txt").write_text("unexpected", encoding="utf-8")
    with pytest.raises(MobilomeDatabaseError, match="exactly"):
        load_mobilome_database(tmp_path)

    (tmp_path / "extra.txt").unlink()
    markers_path = tmp_path / "markers.yaml"
    payload = yaml.safe_load(markers_path.read_text(encoding="utf-8"))
    payload["markers"][0]["matchers"][0]["strength"] = 4
    markers_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    with pytest.raises(MobilomeDatabaseError, match="strength"):
        load_mobilome_database(tmp_path)


def test_aggregate_candidate_claims_are_explicitly_disabled() -> None:
    database = load_mobilome_database()
    disabled = {item.rule_id for item in database.disabled_aggregate_rules}
    executable = set(database.inference_rule_map)

    assert disabled == {
        "phage_module_bearing_plasmid_candidate",
        "pxo_like_annotation_pattern_candidate",
    }
    assert not disabled & executable


def _mutated_database(
    tmp_path: Path,
    resource_name: str,
    mutate: Callable[[dict], None],
) -> None:
    _copy_custom_database(tmp_path)
    path = tmp_path / resource_name
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    mutate(payload)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


@pytest.mark.parametrize(
    ("resource_name", "mutate", "message"),
    [
        (
            "markers.yaml",
            lambda p: p["facets"].append(deepcopy(p["facets"][0])),
            "Duplicate facet IDs",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"].append(deepcopy(p["markers"][0])),
            "Duplicate marker IDs",
        ),
        (
            "provenance.yaml",
            lambda p: p["sources"].append(deepcopy(p["sources"][0])),
            "Duplicate provenance source IDs",
        ),
        (
            "inference.yaml",
            lambda p: p["rules"].append(deepcopy(p["rules"][0])),
            "Duplicate inference rule IDs",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0]["facets"].append("missing_facet"),
            "unknown facets",
        ),
        (
            "inference.yaml",
            lambda p: next(
                rule for rule in p["rules"] if rule.get("required_marker_ids")
            )["required_marker_ids"].append("missing_marker"),
            "unknown markers",
        ),
        (
            "inference.yaml",
            lambda p: p["rules"][0]["required_components"]["subject_record"].append(
                "missing_component"
            ),
            "unknown components",
        ),
        (
            "provenance.yaml",
            lambda p: p["sources"][-1]["rule_ids"].append("missing_rule"),
            "unknown rules",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0]["matchers"][0].update(
                {"mode": "regex", "patterns": ["["]}
            ),
            "Invalid regex",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0]["matchers"][0].update({"strength": 0}),
            "positive integer",
        ),
        (
            "markers.yaml",
            lambda p: next(
                marker
                for marker in p["markers"]
                if marker["id"] == "crispr_repeat_array"
            )["matchers"][1].update({"strength": 2}),
            "evidence ceiling",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0]["matchers"][0].update({"patterns": []}),
            "must not be empty",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0].update({"matchers": []}),
            "non-empty list",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0].update({"limitations": []}),
            "must not be empty",
        ),
        (
            "markers.yaml",
            lambda p: p["markers"][0].update({"id": "replication_origin"}),
            "shadows an existing facet ID",
        ),
        (
            "provenance.yaml",
            lambda p: p["sources"][-1].pop("rationale"),
            "rationale",
        ),
        (
            "provenance.yaml",
            lambda p: p.update({"provenance_version": "2.0.0"}),
            "major versions",
        ),
    ],
)
def test_database_rejects_fail_closed_matrix(
    tmp_path: Path,
    resource_name: str,
    mutate: Callable[[dict], None],
    message: str,
) -> None:
    _mutated_database(tmp_path, resource_name, mutate)

    with pytest.raises(MobilomeDatabaseError, match=message):
        load_mobilome_database(tmp_path)


def test_database_rejects_duplicate_component_mapping_key(tmp_path: Path) -> None:
    _copy_custom_database(tmp_path)
    path = tmp_path / "inference.yaml"
    raw = path.read_text(encoding="utf-8")
    raw = raw.replace("  mobilizable_core:\n", "  replication_candidate:\n", 1)
    path.write_text(raw, encoding="utf-8")

    with pytest.raises(MobilomeDatabaseError, match="duplicate key"):
        load_mobilome_database(tmp_path)


def test_database_rejects_custom_directory_missing_one_resource(tmp_path: Path) -> None:
    _copy_custom_database(tmp_path)
    (tmp_path / "inference.yaml").unlink()

    with pytest.raises(MobilomeDatabaseError, match="exactly"):
        load_mobilome_database(tmp_path)
