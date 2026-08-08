"""Extract genomic sub-regions with optional coordinate rebasing for comparative tools."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .io import read_genbank
from .model import GenBankDocument, GenBankFeature


def extract_region(
    filepath: str | Path,
    locus_tag: str | None = None,
    record_id: str | None = None,
    start: int | None = None,
    end: int | None = None,
    flank_genes: int = 0,
    flank_bp: int = 0,
    rebase: bool = False,
    output_path: str | Path | None = None,
) -> SeqRecord:
    doc = read_genbank(filepath)

    target_rec = None
    target_start = start
    target_end = end

    if locus_tag:
        match = doc.find_locus(locus_tag)
        if match is None:
            print(f"ERROR: Locus tag '{locus_tag}' not found in {filepath}", file=sys.stderr)
            sys.exit(1)
        rec, feat = match
        target_rec = rec
        cdss = sorted(rec.cds_features, key=lambda f: f.start)

        if flank_genes > 0:
            try:
                idx = [f.feature_index for f in cdss].index(feat.feature_index)
                st_idx = max(0, idx - flank_genes)
                en_idx = min(len(cdss) - 1, idx + flank_genes)
                target_start = cdss[st_idx].start
                target_end = cdss[en_idx].end
            except ValueError:
                target_start = feat.start
                target_end = feat.end
        else:
            target_start = feat.start
            target_end = feat.end

        if flank_bp > 0:
            target_start = max(1, target_start - flank_bp)
            target_end = min(rec.length if rec.length > 0 else target_end + flank_bp, target_end + flank_bp)

    elif record_id:
        target_rec = doc.get_record(record_id)
        if target_rec is None:
            print(f"ERROR: Record '{record_id}' not found in {filepath}", file=sys.stderr)
            sys.exit(1)
        target_start = start or 1
        target_end = end or target_rec.length
    else:
        target_rec = doc.records[0]
        target_start = start or 1
        target_end = end or target_rec.length

    target_start = max(1, target_start)
    target_end = max(target_start, target_end)

    # Slice sub-sequence
    full_seq = target_rec.seq
    has_seq = len(full_seq) > 0
    sub_seq = full_seq[target_start - 1:target_end] if has_seq else Seq("")

    # Filter features inside region
    region_features: list[SeqFeature] = []
    rebase_offset = (target_start - 1) if rebase else 0

    for f in target_rec.features:
        # Check feature overlap with [target_start, target_end]
        if f.end < target_start or f.start > target_end:
            continue

        if rebase:
            # Rebase coordinates relative to 1
            new_st = max(0, f.start - target_start)
            new_en = min(target_end - target_start + 1, f.end - target_start + 1)
            new_loc = FeatureLocation(new_st, new_en, strand=f.strand)
        else:
            new_st = max(target_start - 1, f.start - 1)
            new_en = min(target_end, f.end)
            new_loc = FeatureLocation(new_st, new_en, strand=f.strand)

        bio_feat = SeqFeature(
            location=new_loc,
            type=f.type,
            qualifiers=dict(f.qualifiers),
        )
        region_features.append(bio_feat)

    sub_rec_id = f"{target_rec.id}_region_{target_start}_{target_end}"
    sub_record = SeqRecord(
        seq=sub_seq,
        id=sub_rec_id,
        name=sub_rec_id[:16],
        description=f"Extracted region {target_rec.id}:{target_start}..{target_end} (rebase={rebase})",
        features=region_features,
    )
    sub_record.annotations['molecule_type'] = 'DNA'

    print("=" * 70)
    print("  GENOMIC REGION EXTRACTION")
    print("=" * 70)
    print(f"  Parent Record   : {target_rec.id}")
    print(f"  Region Span     : {target_start:,} .. {target_end:,} ({target_end - target_start + 1:,} bp)")
    print(f"  Rebased coords  : {rebase}")
    print(f"  Features sliced : {len(region_features)}")
    print("=" * 70)

    if output_path:
        out_p = Path(output_path)
        out_fmt = 'genbank'
        if out_p.suffix.lower() in ('.fasta', '.fna', '.fa'):
            out_fmt = 'fasta'
        with out_p.open('w', encoding='utf-8') as fh:
            SeqIO.write(sub_record, fh, out_fmt)
        print(f"Wrote region file to {output_path}")

    return sub_record


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract genomic regions with optional flanking windows and coordinate rebasing.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--locus', help="Target locus tag or gene name")
    parser.add_argument('--record', help="Target record/contig ID")
    parser.add_argument('--start', type=int, help="Start coordinate (1-based)")
    parser.add_argument('--end', type=int, help="End coordinate (1-based)")
    parser.add_argument('--flank-genes', type=int, default=0, help="Number of flanking genes around target locus")
    parser.add_argument('--flank-bp', type=int, default=0, help="Number of flanking bp around target locus")
    parser.add_argument('--rebase', action='store_true', help="Rebase coordinates to start at 1")
    parser.add_argument('--output', help="Output file path (.gbk, .gbff, or .fna)")
    args = parser.parse_args()

    extract_region(
        args.input,
        locus_tag=args.locus,
        record_id=args.record,
        start=args.start,
        end=args.end,
        flank_genes=args.flank_genes,
        flank_bp=args.flank_bp,
        rebase=args.rebase,
        output_path=args.output,
    )


if __name__ == '__main__':
    main()
