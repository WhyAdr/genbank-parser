#!/usr/bin/env python3
"""Extract nucleotide sequences from the ORIGIN block of a GenBank file.
Outputs genome FASTA (.fna) and optionally CDS nucleotide slices (.ffn).
Reports GC content and coding density."""
import sys, re
from genbank_parser import parse_features, get_qual


COMPLEMENT = str.maketrans('ATCGatcg', 'TAGCtagc')


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_sequences(filepath):
    """Parse ORIGIN blocks. Returns dict: contig_name -> sequence (uppercase)."""
    sequences = {}
    current_locus = '_unknown_'
    current_seq_lines = []
    in_origin = False

    locus_re = re.compile(r'^LOCUS\s+(\S+)')

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')

            lm = locus_re.match(line)
            if lm:
                current_locus = lm.group(1)
                continue

            if line.startswith('ORIGIN'):
                in_origin = True
                current_seq_lines = []
                continue

            if line.startswith('//'):
                if in_origin and current_seq_lines:
                    seq = ''.join(current_seq_lines).upper()
                    sequences[current_locus] = seq
                in_origin = False
                current_seq_lines = []
                continue

            if in_origin:
                cleaned = re.sub(r'[\s\d]', '', line)
                if cleaned:
                    current_seq_lines.append(cleaned)

    if in_origin and current_seq_lines:
        seq = ''.join(current_seq_lines).upper()
        sequences[current_locus] = seq

    return sequences


def extract_sequences(filepath, out_fna=None, out_ffn=None):
    sequences = parse_sequences(filepath)

    if not sequences:
        print("No ORIGIN sequences found in file.")
        print("(This file may be a feature-table-only snippet.)")
        return

    base = filepath.rsplit('.', 1)[0]
    if out_fna is None:
        out_fna = base + '.fna'
    if out_ffn is None:
        out_ffn = base + '.ffn'

    total_bp = 0
    total_gc = 0
    with open(out_fna, 'w') as fh:
        for contig, seq in sequences.items():
            fh.write(f">{contig} len={len(seq)}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i+70] + '\n')
            total_bp += len(seq)
            total_gc += seq.count('G') + seq.count('C')

    gc_pct = 100 * total_gc / total_bp if total_bp > 0 else 0

    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']
    cds_written = 0
    cds_total_bp = 0

    with open(out_ffn, 'w') as fh:
        for f in cdss:
            contig = f['contig']
            if contig not in sequences:
                continue
            seq = sequences[contig]
            start = f['start'] - 1
            end = f['end']
            if start < 0 or end > len(seq):
                continue

            nt_seq = seq[start:end]
            if f['strand'] == '-':
                nt_seq = reverse_complement(nt_seq)

            tag = get_qual(f, 'locus_tag', 'unknown')
            gene = get_qual(f, 'gene')
            product = get_qual(f, 'product')
            header = f">{tag}"
            if gene:
                header += f" gene={gene}"
            header += f" product={product} [{contig}:{f['start']}..{f['end']}({f['strand']})]"

            fh.write(header + '\n')
            for i in range(0, len(nt_seq), 70):
                fh.write(nt_seq[i:i+70] + '\n')
            cds_written += 1
            cds_total_bp += len(nt_seq)

    coding_density = 100 * cds_total_bp / total_bp if total_bp > 0 else 0

    print("=" * 70)
    print("  SEQUENCE EXTRACTION REPORT")
    print("=" * 70)
    print(f"  File           : {filepath}")
    print(f"  Contigs        : {len(sequences)}")
    print(f"  Total length   : {total_bp:,} bp")
    print(f"  GC content     : {gc_pct:.1f}%")
    print()
    print(f"  CDS extracted  : {cds_written}")
    print(f"  Coding bp      : {cds_total_bp:,}")
    print(f"  Coding density : {coding_density:.1f}%")
    print()
    print(f"  Genome FASTA   : {out_fna}")
    print(f"  CDS nucleotide : {out_ffn}")
    print("=" * 70)


if __name__ == '__main__':
    filepath = sys.argv[1] if len(sys.argv) > 1 else input("File path: ")
    out_fna = sys.argv[2] if len(sys.argv) > 2 else None
    out_ffn = sys.argv[3] if len(sys.argv) > 3 else None
    extract_sequences(filepath, out_fna, out_ffn)
