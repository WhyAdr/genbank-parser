"""Record inventory and declared-replicon classification regressions."""

from pathlib import Path

from genbank_parser import read_genbank
from genbank_parser.mobilome.replicons import (
    classify_replicon,
    inventory_replicons,
)
from genbank_parser.model import GenBankRecord


def test_inventory_uses_record_declarations_not_topology_or_marker_content() -> None:
    inventory = inventory_replicons(
        read_genbank(Path("tests/fixtures/mobilome_inventory.gb"))
    )
    by_id = {
        item.record_id: item for item in inventory if item.record_id != "MOB_DUP.1"
    }

    assert by_id["MOB_CHR.1"].classification.label == "chromosome"
    assert by_id["MOB_CHR.1"].classification.status == "supported"
    assert by_id["MOB_PLASMID.1"].classification.label == "plasmid"
    assert by_id["MOB_TENT.1"].classification.status == "tentative"
    assert by_id["MOB_UNKNOWN.1"].classification.label == "unknown"
    assert by_id["MOB_LINEAR.1"].classification.label == "plasmid"
    assert by_id["MOB_LINEAR.1"].topology == "linear"
    assert by_id["MOB_CONFLICT.1"].classification.status == "conflicting"
    assert len(by_id["MOB_CONFLICT.1"].classification.conflicting_evidence) == 2
    assert by_id["MOB_EMPTY.1"].topology == "unknown"
    assert by_id["MOB_EMPTY.1"].feature_count == 0
    assert (
        by_id["MOB_NEG_PLASMID.1"].classification
        == by_id["MOB_NEG_CHROMOSOME.1"].classification
    )
    assert by_id["MOB_NEG_PLASMID.1"].classification.label == "unknown"
    assert by_id["MOB_NEG_PLASMID.1"].classification.status == "insufficient"


def test_structured_negative_controls_do_not_classify_product_phrases() -> None:
    inventory = inventory_replicons(
        read_genbank(Path("tests/fixtures/mobilome_inventory.gb"))
    )
    negatives = {
        item.record_id: item
        for item in inventory
        if item.record_id in {"MOB_NEG_PLASMID.1", "MOB_NEG_CHROMOSOME.1"}
    }

    assert set(negatives) == {"MOB_NEG_PLASMID.1", "MOB_NEG_CHROMOSOME.1"}
    assert all(item.classification.label == "unknown" for item in negatives.values())
    assert all(
        item.classification.status == "insufficient" for item in negatives.values()
    )


def test_inventory_preserves_duplicate_record_ids_by_record_index() -> None:
    inventory = inventory_replicons(
        read_genbank(Path("tests/fixtures/mobilome_inventory.gb"))
    )
    duplicates = [item for item in inventory if item.record_id == "MOB_DUP.1"]

    assert [item.record_index for item in duplicates] == [10, 11]
    assert all(item.classification.label == "plasmid" for item in duplicates)


def test_defensive_record_annotation_classification_path_is_supported() -> None:
    record = GenBankRecord(
        id="synthetic",
        name="synthetic",
        description="",
        seq="A" * 20,
        length=20,
        annotations={"plasmid": "pDEFENSIVE"},
    )

    classification = classify_replicon(record)

    assert classification.label == "plasmid"
    assert classification.status == "supported"
