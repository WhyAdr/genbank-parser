#!/usr/bin/env python3
"""CRISPR array and Cas gene detector.

Scans a GenBank file for:
  1. repeat_region features with CRISPR-associated annotations
  2. Known Cas protein-coding genes by keyword/gene name
  3. Spatial co-localization of CRISPR arrays with Cas gene clusters

Usage:
    python genbank_crispr.py INPUT.gbff [--window 15000]
"""
import sys, re, argparse, collections
from genbank_parser import parse_features, get_qual

# Cas gene keywords (covers Type I-VI systems)
CAS_KEYWORDS = [
    # Universal / signature
    "cas1", "cas2", "cas4",
    # Type I
    "cas3", "cas5", "cas6", "cas7", "cas8",
    # Type II
    "cas9", "cpf1", "cas12",
    # Type III
    "cas10", "csm", "cmr",
    # Type V/VI
    "cas13", "c2c2",
    # Generic
    "crispr-associated", "crispr associated",
]

REPEAT_KEYWORDS = [
    "crispr", "direct repeat", "palindromic repeat",
]


def _is_cas(f):
    gene    = get_qual(f, 'gene').lower()
    product = get_qual(f, 'product').lower()
    text    = gene + " " + product
    return any(kw in text for kw in CAS_KEYWORDS)


def _is_crispr_array(f):
    if f['type'] != 'repeat_region':
        return False
    product = get_qual(f, 'product').lower()
    note    = ' '.join(f['qualifiers'].get('note', [])).lower()
    text    = product + " " + note
    return any(kw in text for kw in REPEAT_KEYWORDS)


def detect_crispr(filepath, window=15000):
    features = parse_features(filepath)

    arrays   = [f for f in features if _is_crispr_array(f)]
    cas_cdss = [f for f in features if f['type'] == 'CDS' and _is_cas(f)]

    print("=" * 65)
    print("  CRISPR ARRAY / Cas GENE REPORT")
    print("=" * 65)
    print(f"  File          : {filepath}")
    print(f"  CRISPR arrays : {len(arrays)}")
    print(f"  Cas CDS       : {len(cas_cdss)}")
    print()

    if not arrays and not cas_cdss:
        print("  No CRISPR-associated features detected.")
        return

    if arrays:
        print("-- CRISPR Arrays --")
        for a in arrays:
            note = '; '.join(a['qualifiers'].get('note', [])) or '-'
            print(f"  {a['contig']:20s}  {a['start']:>10,d}..{a['end']:>10,d}  "
                  f"{a['strand']}  note: {note[:60]}")
        print()

    if cas_cdss:
        print("-- Cas Genes --")
        for c in cas_cdss:
            print(f"  {get_qual(c, 'locus_tag'):18s}  {get_qual(c, 'gene'):10s}  "
                  f"{c['contig']:20s}  {c['start']:>10,d}..{c['end']:>10,d}  "
                  f"{get_qual(c, 'product')[:50]}")
        print()

    # Spatial co-localization: for each CRISPR array, find Cas genes within window
    if arrays and cas_cdss:
        print(f"-- Co-localization (Cas genes within {window // 1000} kb of array) --")
        found_any = False
        for arr in arrays:
            nearby = [
                c for c in cas_cdss
                if c['contig'] == arr['contig']
                and abs(c['start'] - arr['start']) <= window
            ]
            if nearby:
                found_any = True
                arr_len   = arr['end'] - arr['start'] + 1
                print(f"\n  Array: {arr['contig']}:{arr['start']:,}..{arr['end']:,} "
                      f"({arr_len:,} bp)")
                for c in nearby:
                    dist = c['start'] - arr['start']
                    sign = "+" if dist >= 0 else ""
                    print(f"    {sign}{dist // 1000}kb  {get_qual(c, 'locus_tag'):18s}  "
                          f"{get_qual(c, 'gene'):10s}  {get_qual(c, 'product')[:45]}")
        if not found_any:
            print("  No Cas genes found within window of any CRISPR array.")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Detect CRISPR arrays and Cas genes in a GenBank file."
    )
    parser.add_argument('input',    help="Input GenBank file")
    parser.add_argument('--window', type=int, default=15000,
                        help="Co-localization window in bp (default: 15000)")
    args = parser.parse_args()
    detect_crispr(args.input, args.window)


if __name__ == '__main__':
    main()
