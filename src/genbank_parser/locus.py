"""Single-locus deep-dive: display all qualifiers and cross-references for a target locus."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature


def inspect_locus(filepath: str | Path, locus_tag: str) -> GenBankFeature | None:
    doc = read_genbank(filepath)
    match = doc.find_locus(locus_tag)

    if match is None:
        # Check by gene name
        for rec in doc.records:
            for f in rec.features:
                if f.gene and f.gene.casefold() == locus_tag.casefold():
                    match = (rec, f)
                    break
            if match:
                break

    if match is None:
        print(f"ERROR: Locus tag or gene '{locus_tag}' not found in {filepath}", file=sys.stderr)
        sys.exit(1)

    rec, f = match
    xrefs = extract_xrefs(f)

    print("=" * 70)
    print(f"  FEATURE DEEP-DIVE: {locus_tag}")
    print("=" * 70)
    print(f"  Record / Contig   : {rec.id} (length: {rec.length:,} bp)")
    print(f"  Feature type      : {f.type}")
    print(f"  Coordinates       : {f.start:,} .. {f.end:,} ({f.strand_symbol})")
    print(f"  Biological length : {f.length:,} bp ({f.length // 3} aa)")
    print(f"  Genomic span      : {f.genomic_span:,} bp")
    print(f"  Compound / Join   : {f.is_compound}")
    if f.is_compound:
        print(f"  Join segments     : {f.join_segments}")
    print(f"  Partial coords    : start={f.is_partial_start}, end={f.is_partial_end}")
    print(f"  Pseudogene        : {f.is_pseudo}")
    print()
    print("-- Standard Qualifiers --")
    if f.gene:
        print(f"  Gene symbol       : {f.gene}")
    if f.product:
        print(f"  Product           : {f.product}")
    if f.protein_id:
        print(f"  Protein ID        : {f.protein_id}")
    print(f"  Codon start       : {f.codon_start}")
    print(f"  Transl table      : {f.transl_table}")
    print()
    print("-- Cross-References --")
    if xrefs['kegg_kos']:
        print(f"  KEGG KO           : {', '.join(xrefs['kegg_kos'])}")
    if xrefs['ec_numbers']:
        print(f"  EC Number         : {', '.join(xrefs['ec_numbers'])}")
    if xrefs['cog_ids']:
        print(f"  COG IDs           : {', '.join(xrefs['cog_ids'])}")
    if xrefs['pfam']:
        print(f"  Pfam              : {', '.join(xrefs['pfam'])}")
    if xrefs['rfam']:
        print(f"  Rfam              : {', '.join(xrefs['rfam'])}")
    if xrefs['go_terms']:
        print(f"  GO Terms          : {', '.join(xrefs['go_terms'])}")
    if xrefs['db_xrefs']:
        print(f"  Other db_xrefs    : {', '.join(xrefs['db_xrefs'][:10])}")
    print()
    print("-- All Qualifiers --")
    for k, v in sorted(f.qualifiers.items()):
        if k == 'translation':
            trans = v[0]
            print(f"  /{k:16s} : {trans[:40]}... (length: {len(trans)} aa)")
        else:
            for item in v:
                print(f"  /{k:16s} : {item}")
    print("=" * 70)

    return f


def main() -> None:
    parser = argparse.ArgumentParser(description="Deep-dive inspection of a single locus tag or gene.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('locus_tag', help="Target locus tag or gene name")
    args = parser.parse_args()

    inspect_locus(args.input, args.locus_tag)


if __name__ == '__main__':
    main()
