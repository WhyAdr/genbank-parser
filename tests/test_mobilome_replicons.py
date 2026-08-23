"""Record inventory and declared-replicon classification regressions."""

from pathlib import Path

from genbank_parser import read_genbank
from genbank_parser.mobilome.replicons import inventory_replicons


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


def test_inventory_preserves_duplicate_record_ids_by_record_index() -> None:
    inventory = inventory_replicons(
        read_genbank(Path("tests/fixtures/mobilome_inventory.gb"))
    )
    duplicates = [item for item in inventory if item.record_id == "MOB_DUP.1"]

    assert [item.record_index for item in duplicates] == [8, 9]
    assert all(item.classification.label == "plasmid" for item in duplicates)
