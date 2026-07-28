#!/usr/bin/env python3
"""Comparative mode: scan multiple GenBank files for user-specified gene markers.
Outputs a presence/absence matrix as TSV.

Usage:
    python genbank_compare.py <dir_or_glob> <gene_list> [output.tsv]

    gene_list: comma-separated gene names, EC numbers, COG IDs, or KEGG orthologs.
    Examples:
        "ladA,ssuD,ntaA"
        "EC:1.14.14.28,EC:1.14.14.5"
        "COG:COG2141,KEGG:K20938"
"""
import sys, os, glob, csv, re, argparse
from genbank_parser import parse_features, get_qual, get_notes, extract_xrefs


def _word_in(needle, haystack):
    """True if *needle* appears as a whole word in *haystack* (case-insensitive)."""
    return bool(re.search(r'\b' + re.escape(needle) + r'\b', haystack, re.IGNORECASE))


def scan_genome(filepath, markers):
    """Scan a single GenBank file for presence of each marker.
    Returns dict: marker -> list of matching locus tags.
    """
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']
    results = {m: [] for m in markers}

    for f in cdss:
        tag     = get_qual(f, 'locus_tag', '?')
        gene    = get_qual(f, 'gene').lower()
        product = get_qual(f, 'product').lower()
        xr      = extract_xrefs(f)

        for m in markers:
            ml      = m.lower()
            matched = False

            if gene == ml:
                # Exact gene name match
                matched = True
            elif m.startswith('EC:') and any(m[3:] == ec for ec in xr['ec_numbers']):
                matched = True
            elif m.startswith('COG:') and any(m[4:] == cog or m[4:] in cog for cog in xr['cog_ids']):
                matched = True
            elif m.startswith('KEGG:') and any(m[5:] == ko for ko in xr['kegg_kos']):
                matched = True
            elif not m.startswith(('EC:', 'COG:', 'KEGG:')):
                # Word-boundary product search: "alkB" must not match "walkBase"
                matched = _word_in(ml, product)

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
        genome  = os.path.basename(fp).rsplit('.', 1)[0]
        results = scan_genome(fp, markers)
        row     = {'genome': genome}
        for m in markers:
            hits   = results[m]
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


def main():
    parser = argparse.ArgumentParser(
        description="Scan multiple GenBank files for gene markers; output presence/absence matrix."
    )
    parser.add_argument('input',   help="Directory of GenBank files or a glob pattern")
    parser.add_argument('markers', help="Comma-separated markers (gene names, EC:x, COG:x, KEGG:x)")
    parser.add_argument('output',  nargs='?', default='gene_matrix.tsv',
                        help="Output TSV path (default: gene_matrix.tsv)")
    args = parser.parse_args()
    compare_genomes(args.input, args.markers, args.output)


if __name__ == '__main__':
    main()
