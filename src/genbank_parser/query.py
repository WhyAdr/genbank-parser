"""Feature search and query engine supporting gene, product, KO, COG, EC, Pfam and regex filters."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import sys
from typing import Any

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankDocument, GenBankFeature


def search_features(
    filepath: str | Path,
    gene: str | None = None,
    product: str | None = None,
    ko: str | None = None,
    ec: str | None = None,
    cog: str | None = None,
    pfam: str | None = None,
    ftype: str | None = None,
    gene_regex: str | None = None,
    product_regex: str | None = None,
    format_type: str = 'text',
    output_path: str | Path | None = None,
) -> list[dict[str, Any]]:
    doc = read_genbank(filepath)
    results: list[dict[str, Any]] = []

    gene_pat = re.compile(gene_regex, re.IGNORECASE) if gene_regex else None
    prod_pat = re.compile(product_regex, re.IGNORECASE) if product_regex else None

    for f in doc.all_features:
        if ftype and f.type.casefold() != ftype.casefold():
            continue

        f_gene = f.gene or ''
        f_prod = f.product or ''
        xrefs = extract_xrefs(f)

        if gene and gene.casefold() != f_gene.casefold():
            continue
        if product and product.casefold() not in f_prod.casefold():
            continue
        if gene_pat and not gene_pat.search(f_gene):
            continue
        if prod_pat and not prod_pat.search(f_prod):
            continue
        if ko and not any(ko.casefold() == k.casefold() for k in xrefs['kegg_kos']):
            continue
        if ec and not any(ec in e for e in xrefs['ec_numbers']):
            continue
        if cog and not any(cog.casefold() == c.casefold() for c in xrefs['cog_ids']):
            continue
        if pfam and not any(pfam.casefold() == p.casefold() for p in xrefs['pfam']):
            continue

        res = {
            'record': f.record_id,
            'locus_tag': f.locus_tag or '-',
            'type': f.type,
            'gene': f.gene or '-',
            'product': f.product or '-',
            'start': f.start,
            'end': f.end,
            'strand': f.strand_symbol,
            'length': f.length,
            'kegg_ko': ';'.join(xrefs['kegg_kos']),
            'ec_number': ';'.join(xrefs['ec_numbers']),
            'cog': ';'.join(xrefs['cog_ids']),
            'pfam': ';'.join(xrefs['pfam']),
        }
        results.append(res)

    if format_type == 'json':
        out_json = json.dumps(results, indent=2)
        if output_path:
            Path(output_path).write_text(out_json, encoding='utf-8')
        else:
            print(out_json)
        return results

    if format_type in ('tsv', 'csv'):
        delim = '\t' if format_type == 'tsv' else ','
        fields = ['record', 'locus_tag', 'type', 'gene', 'product', 'start', 'end', 'strand', 'length', 'kegg_ko', 'ec_number', 'cog', 'pfam']
        if output_path:
            with Path(output_path).open('w', newline='', encoding='utf-8') as fh:
                writer = csv.DictWriter(fh, fieldnames=fields, delimiter=delim)
                writer.writeheader()
                writer.writerows(results)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter=delim)
            writer.writeheader()
            writer.writerows(results)
        return results

    # Text mode
    print("=" * 80)
    print(f"  FEATURE SEARCH RESULTS: {filepath}")
    print("=" * 80)
    print(f"  Matches found: {len(results)}")
    print()
    if results:
        print(f"{'Record':16s}  {'Strand':6s}  {'Start':>9s}..{'End':<9s}  {'Locus Tag':16s}  {'Gene':8s}  {'Product'}")
        print("-" * 80)
        for r in results:
            coords = f"{r['start']:,}..{r['end']:,}"
            print(f"  {r['record']:16s}  [{r['strand']}]    {r['start']:>9,d}..{r['end']:<9,d}  {r['locus_tag']:16s}  {r['gene']:8s}  {r['product'][:35]}")
        print("=" * 80)
    else:
        print("  (No features matched search criteria)")
        print("=" * 80)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Search GenBank annotations by gene, product, KO, EC, Pfam, or regex.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--gene', help="Exact gene name match (e.g. 'ladA')")
    parser.add_argument('--product', help="Substring product match (e.g. 'monooxygenase')")
    parser.add_argument('--ko', help="KEGG KO match (e.g. 'K20938')")
    parser.add_argument('--ec', help="EC number match (e.g. '1.14.14.1')")
    parser.add_argument('--cog', help="COG match (e.g. 'COG0596')")
    parser.add_argument('--pfam', help="Pfam ID match (e.g. 'PF00067')")
    parser.add_argument('--feature', dest='ftype', help="Feature type (e.g. 'CDS', 'tRNA')")
    parser.add_argument('--gene-regex', help="Regex pattern on /gene")
    parser.add_argument('--product-regex', help="Regex pattern on /product")
    parser.add_argument('--format', choices=['text', 'tsv', 'csv', 'json'], default='text', help="Output format")
    parser.add_argument('--output', help="Output file path")
    args = parser.parse_args()

    search_features(
        args.input,
        gene=args.gene,
        product=args.product,
        ko=args.ko,
        ec=args.ec,
        cog=args.cog,
        pfam=args.pfam,
        ftype=args.ftype,
        gene_regex=args.gene_regex,
        product_regex=args.product_regex,
        format_type=args.format,
        output_path=args.output,
    )


if __name__ == '__main__':
    main()
