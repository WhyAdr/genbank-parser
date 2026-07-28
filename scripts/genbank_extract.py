#!/usr/bin/env python3
"""Extract a tab-delimited annotation summary from a GenBank feature table.

Columns are semantically typed: GO_terms contains only GO: identifiers,
COG only COG identifiers, etc. Previously the GO_terms column was incorrectly
populated from the entire /db_xref list.
"""
import sys, csv, io, argparse
from genbank_parser import parse_features, get_qual, extract_xrefs


def extract_summary(filepath, output_tsv):
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']

    headers = [
        'locus_tag', 'gene', 'contig', 'start', 'end', 'strand',
        'length_bp', 'length_aa', 'product',
        'EC_number', 'COG', 'KEGG', 'GO_terms', 'Pfam', 'RFAM',
        'db_xrefs', 'inference',
    ]

    rows = []
    for f in cdss:
        xr   = extract_xrefs(f)
        trans = get_qual(f, 'translation')

        rows.append({
            'locus_tag':  get_qual(f, 'locus_tag'),
            'gene':       get_qual(f, 'gene'),
            'contig':     f['contig'],
            'start':      f['start'],
            'end':        f['end'],
            'strand':     f['strand'],
            'length_bp':  f['end'] - f['start'] + 1,
            'length_aa':  len(trans),
            'product':    get_qual(f, 'product'),
            'EC_number':  '; '.join(xr['ec_numbers']),
            'COG':        '; '.join(xr['cog_ids']),
            'KEGG':       '; '.join(xr['kegg_kos']),
            'GO_terms':   '; '.join(xr['go_terms']),
            'Pfam':       '; '.join(xr['pfam']),
            'RFAM':       '; '.join(xr['rfam']),
            'db_xrefs':   '; '.join(xr['db_xrefs']),
            'inference':  '; '.join(f['qualifiers'].get('inference', [])),
        })

    with open(output_tsv, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} CDS records to {output_tsv}")

    # Compact preview
    buf = io.StringIO()
    preview_cols = ['locus_tag', 'gene', 'contig', 'product', 'COG', 'EC_number', 'GO_terms']
    writer2 = csv.DictWriter(buf, fieldnames=preview_cols, delimiter='\t')
    writer2.writeheader()
    for r in rows[:15]:
        writer2.writerow({k: r[k] for k in preview_cols})
    print("\n-- Preview (first 15 CDS) --")
    print(buf.getvalue())


def main():
    parser = argparse.ArgumentParser(
        description="Extract tab-delimited annotation summary from a GenBank file."
    )
    parser.add_argument('input',  help="Input GenBank file (.gbff/.gbk/.gb)")
    parser.add_argument('output', nargs='?', help="Output TSV path (default: <input>_summary.tsv)")
    args = parser.parse_args()

    import os
    outfile = args.output or os.path.splitext(args.input)[0] + '_summary.tsv'
    extract_summary(args.input, outfile)


if __name__ == '__main__':
    main()
