"""Test sequence extraction and codon usage calculation."""
from pathlib import Path

from genbank_parser.codon import analyze_codon_usage
from genbank_parser.io import read_genbank
from genbank_parser.sequence import parse_sequences


def test_sequence_extraction_compound_vs_genomic_slice(compound_joined_gbff: Path) -> None:
    doc = read_genbank(compound_joined_gbff)
    rec = doc.records[0]

    fwd_join = rec.find_locus("JOIN_FWD_001")
    assert fwd_join is not None

    # Biological extraction via Biopython location extract
    extracted_seq = str(fwd_join.extract(rec.seq)).upper()
    assert len(extracted_seq) == 30  # 15 + 15 bp
    # Segments: 10..24 (atgaaagtattatgg) + 34..48 (gcaggcctgattact)
    assert extracted_seq == "ATGAAAGTATTATGGGCAGGCCTGATTACT"

    # Genomic slice without join awareness would have included intervening 9 bp (39 bp total)
    raw_slice = str(rec.seq[fwd_join.start - 1:fwd_join.end]).upper()
    assert len(raw_slice) == 39
    assert "CCCCCCCCC" in raw_slice
    assert "CCCCCCCCC" not in extracted_seq


def test_codon_usage_calculation(simple_cds_gbff: Path) -> None:
    stats = analyze_codon_usage(simple_cds_gbff, min_len_aa=5)
    assert stats['total_codons'] >= 20
    assert 'rscu' in stats
    assert 'gc3s' in stats
    # Check that known codons exist in codon_counts
    assert 'ATG' in stats['codon_counts']
