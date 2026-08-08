"""Codon usage analysis and RSCU (Relative Synonymous Codon Usage) calculation."""
from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path
import sys
from typing import Any

from Bio.Data import CodonTable
from Bio.Seq import Seq

from .io import read_genbank
from .model import GenBankDocument, GenBankFeature

CODON_TABLE_11 = {
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

AA_SYNONYMS = collections.defaultdict(list)
for codon, aa in CODON_TABLE_11.items():
    if aa != '*':
        AA_SYNONYMS[aa].append(codon)


def analyze_codon_usage(
    filepath: str | Path,
    min_len_aa: int = 100,
    output_path: str | Path | None = None,
    include_pseudo: bool = False,
) -> dict[str, Any]:
    doc = read_genbank(filepath)

    has_seqs = any(len(rec.seq) > 0 for rec in doc.records)
    if not has_seqs:
        print("ERROR: No ORIGIN sequences found. Codon usage requires nucleotide sequences.", file=sys.stderr)
        sys.exit(1)

    codon_counts: dict[str, int] = collections.defaultdict(int)
    aa_counts: dict[str, int] = collections.defaultdict(int)
    total_codons = 0
    cds_evaluated = 0

    pos1_gc = 0
    pos2_gc = 0
    pos3_gc = 0
    pos3_syn_gc = 0
    pos3_syn_total = 0

    for rec in doc.records:
        rec_seq = rec.seq
        for f in rec.cds_features:
            if not include_pseudo and f.is_pseudo:
                continue

            extracted_nt = str(f.extract(rec_seq)).upper()
            offset = f.codon_start - 1
            coding_nt = extracted_nt[offset:]

            # Require minimum length and modulo 3
            if len(coding_nt) < min_len_aa * 3 or len(coding_nt) % 3 != 0:
                continue

            cds_evaluated += 1
            for i in range(0, len(coding_nt) - 2, 3):
                codon = coding_nt[i:i + 3]
                if len(codon) == 3 and all(c in 'ACGT' for c in codon):
                    aa = CODON_TABLE_11.get(codon, '?')
                    codon_counts[codon] += 1
                    total_codons += 1

                    if aa != '*':
                        aa_counts[aa] += 1

                    # Positional GC
                    if codon[0] in 'GC':
                        pos1_gc += 1
                    if codon[1] in 'GC':
                        pos2_gc += 1
                    if codon[2] in 'GC':
                        pos3_gc += 1

                    # Synonymous GC3s
                    if aa not in ('M', 'W', '*', '?') and len(AA_SYNONYMS[aa]) > 1:
                        pos3_syn_total += 1
                        if codon[2] in 'GC':
                            pos3_syn_gc += 1

    # Calculate RSCU
    # RSCU = (Count of codon * number of synonymous codons) / sum of counts for that amino acid
    rscu: dict[str, float] = {}
    for aa, synonyms in AA_SYNONYMS.items():
        family_total = sum(codon_counts[c] for c in synonyms)
        n_syn = len(synonyms)
        for c in synonyms:
            if family_total > 0:
                rscu[c] = (codon_counts[c] * n_syn) / family_total
            else:
                rscu[c] = 0.0

    # Summary positional stats
    gc1 = (100.0 * pos1_gc / total_codons) if total_codons > 0 else 0.0
    gc2 = (100.0 * pos2_gc / total_codons) if total_codons > 0 else 0.0
    gc3 = (100.0 * pos3_gc / total_codons) if total_codons > 0 else 0.0
    gc3s = (100.0 * pos3_syn_gc / pos3_syn_total) if pos3_syn_total > 0 else 0.0

    print("=" * 70)
    print("  CODON USAGE & RSCU ANALYSIS")
    print("=" * 70)
    print(f"  File            : {filepath}")
    print(f"  CDSs evaluated  : {cds_evaluated} (min length: {min_len_aa} aa)")
    print(f"  Total codons    : {total_codons:,}")
    print(f"  GC1 / GC2 / GC3 : {gc1:.1f}% / {gc2:.1f}% / {gc3:.1f}%")
    print(f"  GC3s (synonymous): {gc3s:.1f}%")
    print()
    print(f"  {'AA':2s}  {'Codon':5s}  {'Count':>7s}  {'per 1k':>7s}  {'RSCU':>6s}  Visual")
    print("  " + "-" * 55)

    rows: list[dict[str, Any]] = []
    for aa in sorted(AA_SYNONYMS.keys()):
        for codon in sorted(AA_SYNONYMS[aa]):
            cnt = codon_counts[codon]
            per_k = (cnt / total_codons * 1000) if total_codons > 0 else 0.0
            r_val = rscu.get(codon, 0.0)
            bar = '#' * int(r_val * 5)
            flag = ' *' if r_val >= 1.5 else ''
            print(f"  {aa:2s}  {codon:5s}  {cnt:>7,d}  {per_k:>7.1f}  {r_val:>6.2f}  {bar}{flag}")
            rows.append({
                'AminoAcid': aa,
                'Codon': codon,
                'Count': cnt,
                'PerThousand': f"{per_k:.2f}",
                'RSCU': f"{r_val:.3f}",
            })

    if output_path:
        out_p = Path(output_path)
        with out_p.open('w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=['AminoAcid', 'Codon', 'Count', 'PerThousand', 'RSCU'], delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nWrote codon usage table to {output_path}")

    return {
        'total_codons': total_codons,
        'cds_evaluated': cds_evaluated,
        'gc1': gc1,
        'gc2': gc2,
        'gc3': gc3,
        'gc3s': gc3s,
        'codon_counts': dict(codon_counts),
        'rscu': rscu,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate codon usage bias, RSCU, and GC position statistics.")
    parser.add_argument('input', help="Input GenBank file with ORIGIN sequences")
    parser.add_argument('--min-len', type=int, default=100, help="Minimum CDS length in amino acids (default: 100)")
    parser.add_argument('--output', help="Output TSV file path")
    parser.add_argument('--include-pseudo', action='store_true', help="Include pseudogenes in codon counts")
    args = parser.parse_args()

    analyze_codon_usage(args.input, min_len_aa=args.min_len, output_path=args.output, include_pseudo=args.include_pseudo)


if __name__ == '__main__':
    main()
