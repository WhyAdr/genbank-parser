"""Annotation diff mode: compare two annotation versions of the same genome."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankDocument, GenBankFeature


def diff_annotations(
    old_filepath: str | Path,
    new_filepath: str | Path,
    format_type: str = 'text',
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    doc_old = read_genbank(old_filepath)
    doc_new = read_genbank(new_filepath)

    cdss_old = [f for f in doc_old.all_features if f.type == 'CDS']
    cdss_new = [f for f in doc_new.all_features if f.type == 'CDS']

    # Map by coordinate interval (record_id, start, end, strand)
    coords_old = {(f.record_id, f.start, f.end, f.strand): f for f in cdss_old}
    coords_new = {(f.record_id, f.start, f.end, f.strand): f for f in cdss_new}

    shared_coords = set(coords_old.keys()) & set(coords_new.keys())
    added_coords = set(coords_new.keys()) - set(coords_old.keys())
    removed_coords = set(coords_old.keys()) - set(coords_new.keys())

    product_changes: list[dict[str, Any]] = []
    gene_changes: list[dict[str, Any]] = []
    xref_changes: list[dict[str, Any]] = []

    for k in sorted(shared_coords):
        fo = coords_old[k]
        fn = coords_new[k]

        # Product changes
        if (fo.product or '').strip() != (fn.product or '').strip():
            product_changes.append({
                'coords': f"{k[0]}:{k[1]}..{k[2]}({fo.strand_symbol})",
                'old_locus': fo.locus_tag or '-',
                'new_locus': fn.locus_tag or '-',
                'old_product': fo.product or '-',
                'new_product': fn.product or '-',
            })

        # Gene changes
        if (fo.gene or '').strip() != (fn.gene or '').strip():
            gene_changes.append({
                'coords': f"{k[0]}:{k[1]}..{k[2]}({fo.strand_symbol})",
                'old_gene': fo.gene or '-',
                'new_gene': fn.gene or '-',
            })

        # KO changes
        ko_old = set(extract_xrefs(fo)['kegg_kos'])
        ko_new = set(extract_xrefs(fn)['kegg_kos'])
        if ko_old != ko_new:
            xref_changes.append({
                'coords': f"{k[0]}:{k[1]}..{k[2]}({fo.strand_symbol})",
                'type': 'KEGG KO',
                'old': list(ko_old),
                'new': list(ko_new),
            })

    result = {
        'old_file': str(old_filepath),
        'new_file': str(new_filepath),
        'old_cds_count': len(cdss_old),
        'new_cds_count': len(cdss_new),
        'shared_coordinate_cds': len(shared_coords),
        'added_cds': len(added_coords),
        'removed_cds': len(removed_coords),
        'product_name_changes': len(product_changes),
        'gene_name_changes': len(gene_changes),
        'xref_changes': len(xref_changes),
        'details': {
            'product_changes': product_changes[:50],
            'gene_changes': gene_changes[:50],
            'xref_changes': xref_changes[:50],
        },
    }

    if format_type == 'json':
        out_json = json.dumps(result, indent=2)
        if output_path:
            Path(output_path).write_text(out_json, encoding='utf-8')
        else:
            print(out_json)
        return result

    print("=" * 70)
    print("  GENOME ANNOTATION DIFF REPORT")
    print("=" * 70)
    print(f"  Reference (Old): {old_filepath}")
    print(f"  Comparison (New): {new_filepath}")
    print()
    print(f"  CDS Count      : {len(cdss_old)} (old) -> {len(cdss_new)} (new) [delta: {len(cdss_new) - len(cdss_old):+d}]")
    print(f"  Shared Coords  : {len(shared_coords):,}")
    print(f"  Added CDSs     : {len(added_coords):,}")
    print(f"  Removed CDSs   : {len(removed_coords):,}")
    print(f"  Product Changes: {len(product_changes):,}")
    print(f"  Gene Changes   : {len(gene_changes):,}")
    print(f"  Xref Changes   : {len(xref_changes):,}")
    print()

    if product_changes:
        print("-- Sample Product Name Changes (first 10) --")
        for pc in product_changes[:10]:
            print(f"  {pc['coords']:25s}  '{pc['old_product'][:25]}' -> '{pc['new_product'][:25]}'")
        print()

    print("=" * 70)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two annotation versions of the same genome.")
    parser.add_argument('old_file', help="Reference GenBank file")
    parser.add_argument('new_file', help="Updated / alternative GenBank file")
    parser.add_argument('--format', choices=['text', 'json'], default='text', help="Output format")
    parser.add_argument('--output', help="Output file path")
    args = parser.parse_args()

    diff_annotations(args.old_file, args.new_file, format_type=args.format, output_path=args.output)


if __name__ == '__main__':
    main()
