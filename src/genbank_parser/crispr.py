"""CRISPR array and Cas gene cluster detector."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from .io import get_qual, read_genbank
from .model import GenBankFeature

CAS_KEYWORDS = [
    "cas1", "cas2", "cas4", "cas3", "cas5", "cas6", "cas7", "cas8",
    "cas9", "cpf1", "cas12", "cas10", "csm", "cmr", "cas13", "c2c2",
    "crispr-associated", "crispr associated",
]

REPEAT_KEYWORDS = [
    "crispr", "direct repeat", "palindromic repeat",
]


def _is_cas(f: GenBankFeature) -> bool:
    gene = (f.gene or '').lower()
    product = (f.product or '').lower()
    text = f"{gene} {product}"
    return any(kw in text for kw in CAS_KEYWORDS)


def _is_crispr_array(f: GenBankFeature) -> bool:
    if f.type != 'repeat_region':
        return False
    product = (f.product or '').lower()
    notes = ' '.join(f.qualifiers.get('note', [])).lower()
    text = f"{product} {notes}"
    return any(kw in text for kw in REPEAT_KEYWORDS)


def detect_crispr(filepath: str | Path, window: int = 15000) -> dict[str, Any]:
    doc = read_genbank(filepath)
    all_features = doc.all_features

    arrays = [f for f in all_features if _is_crispr_array(f)]
    cas_cdss = [f for f in all_features if f.type == 'CDS' and _is_cas(f)]

    print("=" * 70)
    print("  CRISPR ARRAY / Cas GENE REPORT")
    print("=" * 70)
    print(f"  File          : {filepath}")
    print(f"  CRISPR arrays : {len(arrays)}")
    print(f"  Cas CDSs      : {len(cas_cdss)}")
    print()

    if not arrays and not cas_cdss:
        print("  No CRISPR-associated features detected.")
        print("=" * 70)
        return {'arrays': [], 'cas_genes': [], 'colocalized': []}

    if arrays:
        print("-- CRISPR Arrays --")
        for a in arrays:
            note = '; '.join(a.qualifiers.get('note', [])) or '-'
            print(f"  {a.record_id:20s}  {a.start:>10,d}..{a.end:>10,d}  {a.strand_symbol}  note: {note[:60]}")
        print()

    if cas_cdss:
        print("-- Cas Genes --")
        for c in cas_cdss:
            tag = c.locus_tag or '-'
            gene = c.gene or '-'
            prod = (c.product or '-')[:40]
            print(f"  {c.record_id:20s}  {c.start:>10,d}..{c.end:>10,d}  {c.strand_symbol}  {tag:16s}  {gene:8s}  {prod}")
        print()

    # Spatial co-localization
    colocalized: list[dict[str, Any]] = []
    for a in arrays:
        nearby_cas = [
            c for c in cas_cdss
            if c.record_id == a.record_id and abs(c.start - a.end) <= window
        ]
        if nearby_cas:
            colocalized.append({'array': a, 'cas_genes': nearby_cas})

    if colocalized:
        print(f"-- Spatial Co-localization (within {window:,} bp) --")
        for pair in colocalized:
            arr = pair['array']
            print(f"  Array at {arr.record_id}:{arr.start}..{arr.end} linked with {len(pair['cas_genes'])} Cas gene(s):")
            for c in pair['cas_genes']:
                dist = min(abs(c.start - arr.end), abs(arr.start - c.end))
                print(f"    -> {c.locus_tag or '-'} ({c.gene or '-'}) {c.product} [dist: {dist:,} bp]")
        print()

    print("=" * 70)
    return {
        'arrays': [a.to_dict() for a in arrays],
        'cas_genes': [c.to_dict() for c in cas_cdss],
        'colocalized': len(colocalized),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Detect CRISPR repeat arrays and Cas gene clusters.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--window', type=int, default=15000, help="Window size in bp for array-Cas linking (default: 15000)")
    args = parser.parse_args()

    detect_crispr(args.input, window=args.window)


if __name__ == '__main__':
    main()
