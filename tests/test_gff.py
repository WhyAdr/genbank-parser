"""Test GFF3 conversion, phase calculation, and compound location semantics."""
from pathlib import Path

from genbank_parser.gff import convert_to_gff3


def test_ordinary_feature_not_treated_as_compound(simple_cds_gbff: Path) -> None:
    gff_out = convert_to_gff3(simple_cds_gbff)
    lines = [ln for ln in gff_out.splitlines() if ln and not ln.startswith('#')]

    # Should contain exactly the annotated features without artificial parent/mRNA sub-splits
    feature_types = [ln.split('\t')[2] for ln in lines]
    assert 'mRNA' not in feature_types
    assert 'gene' in feature_types
    assert 'CDS' in feature_types

    # ##sequence-region should cover full 100 bp
    assert "##sequence-region CONTIG_1 1 100" in gff_out


def test_compound_joined_gff3_export_and_phases(compound_joined_gbff: Path) -> None:
    gff_out = convert_to_gff3(compound_joined_gbff)
    lines = [ln for ln in gff_out.splitlines() if ln and not ln.startswith('#')]

    # The ordinary feature (ORDINARY_001) should NOT have mRNA parent
    ordinary_lines = [ln for ln in lines if 'ORDINARY_001' in ln]
    assert len(ordinary_lines) == 1
    assert ordinary_lines[0].split('\t')[2] == 'CDS'

    # The forward join (JOIN_FWD_001) should have gene, mRNA, and 2 CDS segments
    fwd_join_lines = [ln for ln in lines if 'JOIN_FWD_001' in ln]
    types_fwd = [ln.split('\t')[2] for ln in fwd_join_lines]
    assert types_fwd == ['gene', 'mRNA', 'CDS', 'CDS']

    # CDS segments should have phase
    cds_segs = [ln for ln in fwd_join_lines if ln.split('\t')[2] == 'CDS']
    assert len(cds_segs) == 2
    assert cds_segs[0].split('\t')[7] == '0'
    assert cds_segs[1].split('\t')[7] == '0'

    # Reverse join (JOIN_REV_001)
    rev_join_lines = [ln for ln in lines if 'JOIN_REV_001' in ln]
    assert len(rev_join_lines) == 4
    rev_cds_segs = [ln for ln in rev_join_lines if ln.split('\t')[2] == 'CDS']
    assert len(rev_cds_segs) == 2
    assert rev_cds_segs[0].split('\t')[6] == '-'
    assert rev_cds_segs[1].split('\t')[6] == '-'
