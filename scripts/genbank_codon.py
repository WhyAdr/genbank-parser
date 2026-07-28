#!/usr/bin/env python3
"""Codon usage analysis and RSCU (Relative Synonymous Codon Usage) calculation.

Requires a full .gbff file with ORIGIN blocks. Uses /translation qualifiers
for amino acid mapping and the .ffn sequences for codon counting.

Outputs:
  - Per-codon frequency table
  - RSCU values per synonymous codon family
  - High-frequency codon summary (useful for heterologous expression planning)

Usage:
    python genbank_codon.py INPUT.gbff [--min-len 100] [--output codon_usage.tsv]
"""
import sys, re, argparse, collections, csv
from genbank_parser import parse_features, get_qual

# Standard genetic code (NCBI transl_table 11, most bacteria)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C',
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Group codons by amino acid
AA_SYNONYMS = collections.defaultdict(list)
for codon, aa in CODON_TABLE.items():
    if aa != '*':
        AA_SYNONYMS[aa].append(codon)


def _extract_sequences(filepath):
    """Extract genome sequences keyed by contig name from ORIGIN blocks."""
    seqs = {}
    current_locus = None
    in_origin = False
    seq_parts = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('LOCUS'):
                if current_locus and seq_parts:
                    seqs[current_locus] = ''.join(seq_parts).upper()
                current_locus = line.split()[1]
                seq_parts = []
                in_origin = False
            elif line.startswith('ORIGIN'):
                in_origin = True
            elif line.startswith('//'):
                if current_locus and seq_parts:
                    seqs[current_locus] = ''.join(seq_parts).upper()
                seq_parts = []
                in_origin = False
            elif in_origin:
                seq_parts.append(re.sub(r'[\d\s]', '', line))
    return seqs


def codon_usage(filepath, min_len=100, output_tsv=None):
    features = parse_features(filepath)
    genome   = _extract_sequences(filepath)

    if not genome:
        print("No ORIGIN sequences found. A full .gbff file is required.")
        sys.exit(1)

    cdss = [f for f in features if f['type'] == 'CDS']
    codon_counts = collections.Counter()
    genes_counted = 0

    for f in cdss:
        contig_seq = genome.get(f['contig'])
        if not contig_seq:
            continue
        # Extract coding sequence (0-based slicing; GenBank coords are 1-based)
        s, e = f['start'] - 1, f['end']
        cds_seq = contig_seq[s:e]
        if f['strand'] == '-':
            cds_seq = cds_seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))
        if len(cds_seq) < min_len or len(cds_seq) % 3 != 0:
            continue
        for i in range(0, len(cds_seq) - 3, 3):   # exclude stop codon
            codon = cds_seq[i:i + 3]
            if codon in CODON_TABLE:
                codon_counts[codon] += 1
        genes_counted += 1

    if not codon_counts:
        print("No valid CDS sequences extracted -- check sequence availability.")
        sys.exit(1)

    total_codons = sum(codon_counts.values())

    # RSCU: observed / expected (= count / (total_synonyms_for_aa * uniform_expected))
    rows = []
    for aa in sorted(AA_SYNONYMS):
        synonyms    = AA_SYNONYMS[aa]
        aa_total    = sum(codon_counts[c] for c in synonyms)
        n_syn       = len(synonyms)
        expected_ea = aa_total / n_syn if aa_total > 0 else 0
        for codon in sorted(synonyms):
            count = codon_counts[codon]
            rscu  = (count / expected_ea) if expected_ea > 0 else 0.0
            freq_per_k = 1000 * count / total_codons if total_codons > 0 else 0.0
            rows.append({
                'codon':       codon,
                'amino_acid':  aa,
                'count':       count,
                'freq_per_k':  round(freq_per_k, 3),
                'rscu':        round(rscu, 4),
                'n_synonyms':  n_syn,
            })

    # Console summary
    print(f"Codon Usage Analysis")
    print(f"  Genes counted : {genes_counted}")
    print(f"  Total codons  : {total_codons:,}")
    print()
    print(f"  {'Codon':6s}  {'AA':4s}  {'Count':>8s}  {'Freq/k':>8s}  {'RSCU':>8s}")
    print("  " + "-" * 44)
    for r in rows:
        flag = "  <-- preferred" if r['rscu'] >= 1.5 else ""
        print(f"  {r['codon']:6s}  {r['amino_acid']:4s}  {r['count']:>8,d}  "
              f"{r['freq_per_k']:>8.3f}  {r['rscu']:>8.4f}{flag}")

    # Most preferred codon per amino acid
    print()
    print("-- Preferred Codon per Amino Acid (highest RSCU) --")
    for aa in sorted(AA_SYNONYMS):
        aa_rows = [r for r in rows if r['amino_acid'] == aa]
        best    = max(aa_rows, key=lambda r: r['rscu'])
        print(f"  {aa} :  {best['codon']}  (RSCU={best['rscu']:.3f})")

    if output_tsv:
        headers = ['codon', 'amino_acid', 'count', 'freq_per_k', 'rscu', 'n_synonyms']
        with open(output_tsv, 'w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nFull table written to {output_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute codon usage frequency and RSCU from a full GenBank file."
    )
    parser.add_argument('input',     help="Input GenBank file with ORIGIN blocks")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Minimum CDS length in bp (default: 100)")
    parser.add_argument('--output',  help="Optional TSV output path")
    args = parser.parse_args()
    codon_usage(args.input, args.min_len, args.output)


if __name__ == '__main__':
    main()
