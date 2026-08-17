"""Core genomic-neighborhood selection regressions."""

from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pytest

from genbank_parser.model import GenBankFeature, GenBankRecord
from genbank_parser.neighborhood import (
    SCHEMA_VERSION,
    build_neighborhood,
    serialize_neighborhood,
    write_neighborhood,
)
from genbank_parser.discover import load_ruleset, match_feature_rules
from genbank_parser.io import read_genbank
from genbank_parser.spatial import TargetNotFoundError, select_cds_window


def _feature(index: int, start: int, end: int, *, strand: int = 1) -> GenBankFeature:
    return GenBankFeature(
        record_id="CIRC",
        record_index=0,
        feature_index=index,
        type="CDS",
        location=FeatureLocation(start - 1, end, strand=strand),
        qualifiers={"locus_tag": [f"CDS_{index}"]},
        record_length=1000,
        topology="circular",
    )


def _record(*, circular: bool = True) -> GenBankRecord:
    features = [
        _feature(1, 20, 80),
        _feature(2, 200, 260),
        _feature(3, 500, 560),
        _feature(4, 900, 960),
    ]
    return GenBankRecord(
        id="CIRC",
        name="CIRC",
        description="",
        seq=Seq("A" * 1000),
        length=1000,
        topology="circular" if circular else "linear",
        features=features,
    )


def test_circular_window_crosses_origin_once() -> None:
    record = _record()
    result = select_cds_window(record, record.features[0], 1)
    assert [feature.locus_tag for feature in result.features] == [
        "CDS_4",
        "CDS_1",
        "CDS_2",
    ]
    assert result.wraps_origin
    assert sum(feature.feature_index == 1 for feature in result.features) == 1


def test_oversized_circular_window_has_no_duplicates() -> None:
    record = _record()
    result = select_cds_window(record, record.features[1], 50)
    indices = [feature.feature_index for feature in result.features]
    assert len(indices) == len(set(indices)) == 4
    assert indices == [1, 2, 3, 4]


def test_linear_boundary_and_zero_windows() -> None:
    record = _record(circular=False)
    assert [f.feature_index for f in select_cds_window(record, record.features[0], 2).features] == [1, 2, 3]
    assert [f.feature_index for f in select_cds_window(record, record.features[-1], 2).features] == [2, 3, 4]
    assert [f.feature_index for f in select_cds_window(record, record.features[2], 0).features] == [3]


def _write_circular_genbank(path) -> None:
    path.write_text(
        """LOCUS       CIRC                     1000 bp    DNA     circular BCT 17-AUG-2026
DEFINITION  Synthetic circular neighborhood.
ACCESSION   CIRC
FEATURES             Location/Qualifiers
     CDS             20..80
                     /locus_tag="CDS_1"
                     /gene="dup"
                     /product="alpha protein"
     CDS             200..260
                     /locus_tag="CDS_2"
                     /gene="dup"
                     /product="beta protein"
     CDS             complement(500..560)
                     /locus_tag="CDS_3"
                     /gene="third"
                     /product="gamma protein"
     CDS             900..960
                     /locus_tag="CDS_4"
                     /gene="fourth"
                     /product="delta protein"
ORIGIN
        1 """
        + "a" * 1000
        + "\n//\n",
        encoding="utf-8",
    )


def test_build_neighborhood_is_pure_and_resolves_duplicate_gene_deterministically(
    tmp_path, capsys
) -> None:
    source = tmp_path / "circular.gb"
    _write_circular_genbank(source)
    result = build_neighborhood(source, "dup", 1)
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""
    assert result.target.feature.locus_tag == "CDS_1"
    assert [item.feature.locus_tag for item in result.features] == [
        "CDS_4",
        "CDS_1",
        "CDS_2",
    ]
    assert result.wraps_origin
    assert [item.local_start for item in result.features] == sorted(
        item.local_start for item in result.features
    )


def test_neighborhood_serializers_are_schema_versioned_and_stable(tmp_path) -> None:
    source = tmp_path / "circular.gb"
    _write_circular_genbank(source)
    result = build_neighborhood(source, "CDS_1", 1)
    json_text = serialize_neighborhood(result, "json")
    tsv_text = serialize_neighborhood(result, "tsv")
    assert f'"schema_version": "{SCHEMA_VERSION}"' in json_text
    assert tsv_text.splitlines()[0].startswith("record\torder\tis_target")
    assert len(tsv_text.splitlines()) == 4

    output = tmp_path / "neighborhood.json"
    assert write_neighborhood(result, format_type="json", output_path=output) == json_text
    assert output.read_text(encoding="utf-8") == json_text
    with pytest.raises(ValueError, match="cannot overwrite"):
        write_neighborhood(result, output_path=source)


def test_missing_target_raises_domain_error(tmp_path) -> None:
    source = tmp_path / "circular.gb"
    _write_circular_genbank(source)
    with pytest.raises(TargetNotFoundError, match="was not found"):
        build_neighborhood(source, "absent", 1)


def test_canonical_rules_context_features_and_operon_links(
    context_neighborhood_gbff,
) -> None:
    document = read_genbank(context_neighborhood_gbff)
    matches = match_feature_rules(
        document.records[0].cds_features[0], load_ruleset("mobilome")
    )
    assert [(match.rule_id, match.weight) for match in matches] == [
        ("transposase", 3)
    ]

    result = build_neighborhood(
        context_neighborhood_gbff,
        "CTX_001",
        1,
        include_feature_types=("CDS", "tRNA"),
        ruleset="mobilome",
        show_operons=True,
        operon_gap=60,
    )
    assert [item.feature.type for item in result.features] == [
        "CDS",
        "CDS",
        "tRNA",
        "CDS",
    ]
    assert result.target.rule_matches[0].rule_id == "transposase"
    assert len(result.operon_links) == 2
    payload = result.to_dict()
    assert payload["ruleset"] == "mobilome"
    assert payload["features"][1]["rule_matches"][0]["rule"] == "transposase"
