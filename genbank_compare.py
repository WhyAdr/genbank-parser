#!/usr/bin/env python3
"""
Comparative mode: scan multiple GenBank files for user-specified gene markers.
Outputs a presence/absence matrix as TSV.

Usage:
    python genbank_compare.py <dir_or_glob> <gene_list> [output.tsv]

    gene_list: comma-separated gene names, EC numbers, COG IDs, or KEGG orthologs.
    Examples:
        "ladA,ssuD,ntaA"
        "EC:1.14.14.28,EC:1.14.14.5"
        "COG:COG2141,KEGG:K20938"
"""
import sys, os, glob, csv
from genbank_parser import parse_features, get_qual, get_notes


def scan_genome(filepath, markers):
    """Scan a single GenBank file for presence of each marker.
    Returns dict: marker -> list of matching locus tags.
    """
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']
    results = {m: [] for m in markers}

    for f in cdss:
        tag = get_qual(f, 'locus_tag', '?')
        gene = get_qual(f, 'gene').lower()
        product = get_qual(f, 'product').lower()
        ec = get_qual(f, 'EC_number')
        notes = f['qualifiers'].get('note', [])

        for m in markers:
            ml = m.lower()
            matched = False

            if gene == ml:
                matched = True
            elif m.startswith('EC:') and ec == m[3:]:
                matched = True
            elif m.startswith('COG:') and any(m[4:] in n for n in notes):
                matched = True
            elif m.startswith('KEGG:') and any(m[5:] in n for n in notes):
                matched = True
            elif ml in product:
                matched = True

            if matched:
                results[m].append(tag)

    return results


def compare_genomes(input_path, markers_str, output_tsv):
    markers = [m.strip() for m in markers_str.split(',')]

    if os.path.isdir(input_path):
        files = sorted(
            glob.glob(os.path.join(input_path, '*.gbff')) +
            glob.glob(os.path.join(input_path, '*.gbk')) +
            glob.glob(os.path.join(input_path, '*.gb')) +
            glob.glob(os.path.join(input_path, '*.txt'))
        )
    else:
        files = sorted(glob.glob(input_path))

    if not files:
        print(f"No GenBank files found in '{input_path}'")
        sys.exit(1)

    print(f"Scanning {len(files)} genomes for {len(markers)} markers...")

    rows = []
    for fp in files:
        genome = os.path.basename(fp).rsplit('.', 1)[0]
        results = scan_genome(fp, markers)
        row = {'genome': genome}
        for m in markers:
            hits = results[m]
            row[m] = f"{len(hits)} ({','.join(hits[:3])})" if hits else '0'
        rows.append(row)
        status = '  '.join(f"{m}={'Y' if results[m] else 'N'}" for m in markers)
        print(f"  {genome}: {status}")

    headers = ['genome'] + markers
    with open(output_tsv, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote presence/absence matrix to {output_tsv}")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python genbank_compare.py <dir_or_glob> <gene1,gene2,...> [output.tsv]")
        sys.exit(1)
    input_path = sys.argv[1]
    markers = sys.argv[2]
    output = sys.argv[3] if len(sys.argv) > 3 else 'gene_matrix.tsv'
    compare_genomes(input_path, markers, output)
