"""Phylogenomic marker gene extractor."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys
from typing import Any

from .io import get_qual, read_genbank
from .model import GenBankFeature

RIBOSOMAL = {
    'rpsB': '30S ribosomal protein S2',
    'rpsC': '30S ribosomal protein S3',
    'rpsD': '30S ribosomal protein S4',
    'rpsE': '30S ribosomal protein S5',
    'rpsG': '30S ribosomal protein S7',
    'rpsH': '30S ribosomal protein S8',
    'rpsI': '30S ribosomal protein S9',
    'rpsJ': '30S ribosomal protein S10',
    'rpsK': '30S ribosomal protein S11',
    'rpsL': '30S ribosomal protein S12',
    'rpsM': '30S ribosomal protein S13',
    'rplA': '50S ribosomal protein L1',
    'rplB': '50S ribosomal protein L2',
    'rplC': '50S ribosomal protein L3',
    'rplD': '50S ribosomal protein L4',
    'rplE': '50S ribosomal protein L5',
    'rplF': '50S ribosomal protein L6',
    'rplK': '50S ribosomal protein L11',
    'rplL': '50S ribosomal protein L12',
    'rplN': '50S ribosomal protein L14',
    'rplP': '50S ribosomal protein L16',
    'rplR': '50S ribosomal protein L18',
    'rplV': '50S ribosomal protein L22',
    'rplW': '50S ribosomal protein L23',
}

HOUSEKEEPING = {
    'recA': 'DNA recombinase A',
    'gyrB': 'DNA gyrase subunit B',
    'rpoB': 'RNA polymerase subunit beta',
    'dnaK': 'DnaK chaperone',
    'groEL': 'GroEL chaperonin (Hsp60)',
    'tsf': 'elongation factor Ts',
    'tuf': 'elongation factor Tu',
    'fusA': 'elongation factor G',
}


def _match_marker(f: GenBankFeature, marker_name: str, desc: str) -> bool:
    gene = f.gene
    prod = f.product

    # Priority 1: Exact /gene match
    if gene and gene.casefold() == marker_name.casefold():
        return True

    # Priority 2: Exact /product match
    if prod and prod.casefold() == desc.casefold():
        return True

    # Priority 3: Regex match ONLY if /gene is empty
    if not gene and prod:
        pattern = re.compile(re.escape(desc) + r'(?![\w\-])', re.IGNORECASE)
        if pattern.search(prod):
            return True

    return False


def extract_phylogenomic_markers(
    filepath: str | Path,
    marker_set: str = 'all',
    output_dir: str | Path | None = None,
) -> dict[str, GenBankFeature]:
    doc = read_genbank(filepath)
    cdss = [f for f in doc.all_features if f.type == 'CDS']

    targets: dict[str, str] = {}
    if marker_set in ('core', 'all'):
        targets.update(RIBOSOMAL)
    if marker_set in ('housekeeping', 'all'):
        targets.update(HOUSEKEEPING)

    found_markers: dict[str, GenBankFeature] = {}

    for marker_name, desc in targets.items():
        for f in cdss:
            if _match_marker(f, marker_name, desc):
                found_markers[marker_name] = f
                break

    print("=" * 70)
    print("  PHYLOGENOMIC MARKER EXTRACTION")
    print("=" * 70)
    print(f"  File             : {filepath}")
    print(f"  Marker panel     : {marker_set} ({len(targets)} candidate genes)")
    print(f"  Markers found    : {len(found_markers)} / {len(targets)} ({100.0 * len(found_markers) / len(targets):.1f}%)")
    print()

    print(f"{'Marker':10s}  {'Status':8s}  {'Locus Tag':18s}  {'Length':>8s}  {'Product'}")
    print("-" * 70)
    for m, desc in sorted(targets.items()):
        if m in found_markers:
            f = found_markers[m]
            tlen = f"{len(f.translation)} aa" if f.translation else f"{f.length} bp"
            tag = f.locus_tag or '-'
            prod = (f.product or desc)[:30]
            print(f"  {m:10s}  FOUND     {tag:18s}  {tlen:>8s}  {prod}")
        else:
            print(f"  {m:10s}  MISSING   {'-':18s}  {'-':>8s}  {desc[:30]}")
    print("=" * 70)

    if output_dir and found_markers:
        out_d = Path(output_dir)
        out_d.mkdir(parents=True, exist_ok=True)
        for m, f in found_markers.items():
            if f.translation:
                fa_path = out_d / f"{m}.faa"
                fa_path.write_text(f">{m} {f.locus_tag} {f.product}\n{f.translation}\n", encoding='utf-8')
        print(f"\nWrote individual marker FASTA files to {output_dir}")

    return found_markers


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract phylogenomic core & housekeeping marker genes.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--markers', choices=['core', 'housekeeping', 'all'], default='all', help="Marker panel to query")
    parser.add_argument('--output-dir', help="Directory to write individual marker protein FASTA files")
    args = parser.parse_args()

    extract_phylogenomic_markers(args.input, marker_set=args.markers, output_dir=args.output_dir)


if __name__ == '__main__':
    main()
