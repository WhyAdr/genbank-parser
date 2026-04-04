#!/usr/bin/env python3
"""Extract a tab-delimited annotation summary from a GenBank feature table."""
import sys, csv, io
from genbank_parser import parse_features, get_qual, get_notes


def extract_summary(filepath, output_tsv):
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']

    headers = [
        'locus_tag', 'gene', 'start', 'end', 'strand', 'length_bp', 'length_aa',
        'product', 'EC_number', 'COG', 'KEGG', 'GO_terms', 'RefSeq', 'inference'
    ]

    rows = []
    for f in cdss:
        notes = f['qualifiers'].get('note', [])
        cog = next((n.split(':')[1] for n in notes if n.startswith('COG:COG')), '')
        kegg = next((n.split(':')[1] for n in notes if n.startswith('KEGG:')), '')
        refseq = next((n.split(':')[1] for n in notes if n.startswith('RefSeq:')), '')
        go = ';'.join(f['qualifiers'].get('db_xref', []))
        trans = get_qual(f, 'translation')

        rows.append({
            'locus_tag': get_qual(f, 'locus_tag'),
            'gene': get_qual(f, 'gene'),
            'start': f['start'],
            'end': f['end'],
            'strand': f['strand'],
            'length_bp': f['end'] - f['start'] + 1,
            'length_aa': len(trans),
            'product': get_qual(f, 'product'),
            'EC_number': get_qual(f, 'EC_number'),
            'COG': cog,
            'KEGG': kegg,
            'GO_terms': go,
            'RefSeq': refseq,
            'inference': '; '.join(f['qualifiers'].get('inference', [])),
        })

    with open(output_tsv, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} CDS records to {output_tsv}")

    # Print a compact preview (first 15 rows)
    buf = io.StringIO()
    preview_cols = ['locus_tag', 'gene', 'product', 'COG', 'EC_number']
    writer2 = csv.DictWriter(buf, fieldnames=preview_cols, delimiter='\t')
    writer2.writeheader()
    for r in rows[:15]:
        writer2.writerow({k: r[k] for k in preview_cols})
    print("\n-- Preview (first 15 CDS) --")
    print(buf.getvalue())


if __name__ == '__main__':
    infile = sys.argv[1] if len(sys.argv) > 1 else input("File path: ")
    outfile = sys.argv[2] if len(sys.argv) > 2 else infile.rsplit('.', 1)[0] + '_summary.tsv'
    extract_summary(infile, outfile)
