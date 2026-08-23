"""Cautious record-local inference before spatial and cross-record rules."""

from pathlib import Path

from genbank_parser import read_genbank
from genbank_parser.mobilome.database import load_mobilome_database
from genbank_parser.mobilome.inference import infer_replicon_hypotheses
from genbank_parser.mobilome.replicons import inventory_replicons
from genbank_parser.mobilome.scanner import scan_mobilome_features


def test_empty_record_gets_only_a_calibrated_replication_assessment() -> None:
    document = read_genbank(Path("tests/fixtures/mobilome_inference.gb"))
    database = load_mobilome_database()
    inventory = inventory_replicons(document)
    hits = scan_mobilome_features(document.all_features, database)

    empty = inventory[2]
    hypotheses = infer_replicon_hypotheses(empty, hits, database)

    assert [item.summary for item in hypotheses] == [
        "Replication evidence insufficient; mechanism unresolved"
    ]
    assert "satellite" not in hypotheses[0].summary.casefold()
