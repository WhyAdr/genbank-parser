"""Extract nucleotide sequences from GenBank records.
Outputs genome FASTA (.fna) and CDS nucleotide sequences (.ffn) via biological extraction.
Reports GC content and coding density."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from Bio.Seq import Seq

from .io import read_genbank


def reverse_complement(seq_str: str) -> str:
    """Return reverse complement of a nucleotide sequence string."""
    return str(Seq(seq_str).reverse_complement())


def parse_sequences(filepath: str | Path) -> dict[str, str]:
    """Return dict mapping contig/record ID to full nucleotide sequence string (uppercase)."""
    doc = read_genbank(filepath)
    return {rec.id: str(rec.seq).upper() for rec in doc.records if len(rec.seq) > 0}


def extract_sequences(
    filepath: str | Path,
    out_fna: str | Path | None = None,
    out_ffn: str | Path | None = None,
) -> None:
    doc = read_genbank(filepath)

    has_seqs = any(len(rec.seq) > 0 for rec in doc.records)
    if not has_seqs:
        print("No nucleotide sequences found in file.")
        print("(This file may be a feature-table-only snippet without ORIGIN blocks.)")
        return

    base = os.path.splitext(str(filepath))[0]
    if out_fna is None:
        out_fna = base + ".fna"
    if out_ffn is None:
        out_ffn = base + ".ffn"

    total_acgt = 0
    total_gc = 0
    total_bp = 0

    with open(out_fna, "w", encoding="utf-8") as fh:
        for rec in doc.records:
            seq_str = str(rec.seq).upper()
            fh.write(f">{rec.id} len={len(seq_str)}\n")
            fh.writelines(
                seq_str[i : i + 70] + "\n" for i in range(0, len(seq_str), 70)
            )
            total_bp += len(seq_str)
            total_gc += seq_str.count("G") + seq_str.count("C")
            total_acgt += (
                seq_str.count("A")
                + seq_str.count("C")
                + seq_str.count("G")
                + seq_str.count("T")
            )

    gc_pct = (100.0 * total_gc / total_acgt) if total_acgt > 0 else 0.0

    cds_written = 0
    cds_total_bp = 0

    with open(out_ffn, "w", encoding="utf-8") as fh:
        for rec in doc.records:
            rec_seq = rec.seq
            for f in rec.cds_features:
                # Biological extraction respecting joins and strand
                nt_seq = str(f.extract(rec_seq)).upper()
                if not nt_seq:
                    continue

                tag = f.locus_tag or "unknown"
                gene = f.gene
                product = f.product
                header = f">{tag}"
                if gene:
                    header += f" gene={gene}"
                header += f" product={product} [{rec.id}:{f.start}..{f.end}({f.strand_symbol})]"

                fh.write(header + "\n")
                fh.writelines(
                    nt_seq[i : i + 70] + "\n" for i in range(0, len(nt_seq), 70)
                )
                cds_written += 1
                cds_total_bp += len(nt_seq)

    cds_feature_bp_sum = sum(record.cds_feature_bp_sum for record in doc.records)
    nonredundant_coding_bp = sum(
        record.nonredundant_coding_bp for record in doc.records
    )
    coding_density = (
        (100.0 * nonredundant_coding_bp / total_bp) if total_bp > 0 else 0.0
    )

    print("=" * 70)
    print("  SEQUENCE EXTRACTION REPORT")
    print("=" * 70)
    print(f"  File           : {filepath}")
    print(f"  Contigs        : {len(doc.records)}")
    print(f"  Total length   : {total_bp:,} bp")
    print(f"  GC content     : {gc_pct:.1f}%")
    print()
    print(f"  CDS extracted  : {cds_written}")
    print(f"  CDS bp sum     : {cds_feature_bp_sum:,}")
    print(f"  Nonredundant bp: {nonredundant_coding_bp:,}")
    print(f"  Coding density : {coding_density:.1f}%")
    print()
    print(f"  Genome FASTA   : {out_fna}")
    print(f"  CDS nucleotide : {out_ffn}")
    print("=" * 70)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract genome FASTA and per-CDS nucleotide sequences."
    )
    parser.add_argument("input", help="Input GenBank file with ORIGIN sequences")
    parser.add_argument("--fna", help="Genome FASTA output path (default: <input>.fna)")
    parser.add_argument(
        "--ffn", help="CDS nucleotide FASTA output path (default: <input>.ffn)"
    )
    args = parser.parse_args()
    extract_sequences(args.input, args.fna, args.ffn)


if __name__ == "__main__":
    main()
