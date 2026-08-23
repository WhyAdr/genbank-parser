"""Adjudicated typed-parser facts used by the native mobilome subsystem."""

from pathlib import Path

from genbank_parser import read_genbank

REPORT_SCHEMA_VERSION = "gbparse.mobilome.v1"
CATALOG_VERSION = "1.0.0"
INFERENCE_VERSION = "1.0.0"
FORBIDDEN_CLAIMS = (
    "obligate co-mobilisation",
    "obligate co-mobilization",
    "obligate co-transfer",
    "guaranteed transfer",
    "confirmed conjugative",
    "diagnostic pXO",
    "proves pXO ancestry",
    "confirmed phagemid",
    "confirmed resistance",
    "confirmed virulence",
    "active communication",
    "satellite plasmid",
    "evolutionary transition",
    "theta replication",
    "rolling-circle replication",
)


def test_mobilome_inventory_fixture_preserves_record_identity() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_inventory.gb"))

    assert len(document.records) == 11
    assert len(document.all_features) == 14
    duplicate_ids = [
        record.id for record in document.records if record.id == "MOB_DUP.1"
    ]
    assert duplicate_ids == ["MOB_DUP.1", "MOB_DUP.1"]
    assert document.records[8].topology is None
    assert document.records[8].features == []


def test_mobilome_evidence_fixture_uses_canonical_pseudo_and_location_semantics() -> (
    None
):
    document = read_genbank(Path("tests/fixtures/mobilome_evidence.gb"))
    record = document.records[0]
    pseudo = {
        feature.feature_index: feature
        for feature in record.features
        if feature.is_pseudo
    }

    assert {feature.type.casefold() for feature in pseudo.values()} == {
        "cds",
        "pseudogene",
    }
    assert len(pseudo) == 3
    wrapped = next(feature for feature in record.features if feature.gene == "cas9")
    assert wrapped.join_segments == [(571, 600), (1, 30)]
    assert wrapped.length == 60
    assert wrapped.is_compound
    assert next(
        feature for feature in record.features if feature.gene == "repFR"
    ).is_partial
