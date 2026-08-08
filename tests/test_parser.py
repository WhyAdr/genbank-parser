"""Test parser IO and cross-reference extraction."""
from pathlib import Path

from genbank_parser.io import extract_xrefs, parse_features, read_genbank


def test_parse_features_flat_list(simple_cds_gbff: Path) -> None:
    features = parse_features(simple_cds_gbff)
    assert len(features) >= 4  # source, gene, CDS, gene, CDS
    types = [f.type for f in features]
    assert 'source' in types
    assert 'CDS' in types


def test_extract_xrefs_categories(simple_cds_gbff: Path) -> None:
    doc = read_genbank(simple_cds_gbff)
    fwd_cds = doc.records[0].cds_features[0]

    xrefs = extract_xrefs(fwd_cds)
    assert '1.1.1.1' in xrefs['ec_numbers']
    assert 'COG0596' in xrefs['cog_ids']
    assert 'K00844' in xrefs['kegg_kos']
    assert 'PF00067' in xrefs['pfam']
    assert 'GO:0008152' in xrefs['go_terms']


def test_partial_features(special_cds_gbff: Path) -> None:
    doc = read_genbank(special_cds_gbff)
    part_cds = doc.records[0].find_locus("PARTIAL_001")
    assert part_cds is not None
    assert part_cds.is_partial is True
    assert part_cds.is_partial_start is True
    assert part_cds.is_partial_end is True


def test_pseudogene_detection(special_cds_gbff: Path) -> None:
    doc = read_genbank(special_cds_gbff)
    pseudo_cds = doc.records[0].find_locus("PSEUDO_001")
    assert pseudo_cds is not None
    assert pseudo_cds.is_pseudo is True
