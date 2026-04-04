#!/usr/bin/env python3
"""Extract all qualifiers for a single locus tag -- the surgical scalpel."""
import sys
from genbank_parser import parse_features, get_qual


def extract_locus(filepath, target_tag):
    features = parse_features(filepath)
    hits = [f for f in features if f['type'] == 'CDS'
            and target_tag in f['qualifiers'].get('locus_tag', [])]

    if not hits:
        print(f"Locus tag '{target_tag}' not found among CDS features.")
        sys.exit(1)

    for f in hits:
        q = f['qualifiers']
        loc = f"complement({f['start']}..{f['end']})" if f['strand'] == '-' else f"{f['start']}..{f['end']}"
        print(f"-- {target_tag} --")
        print(f"  Location : {loc}")
        print(f"  Type     : {f['type']}")
        print()

        identity_keys = ['gene', 'product', 'protein_id', 'locus_tag']
        function_keys = ['EC_number', 'codon_start', 'transl_table']
        evidence_keys = ['inference']

        print("  -- Identity --")
        for k in identity_keys:
            for v in q.get(k, []):
                print(f"    /{k}=\"{v}\"")

        print("  -- Function --")
        for k in function_keys:
            for v in q.get(k, []):
                print(f"    /{k}=\"{v}\"")

        print("  -- Cross-references --")
        for v in q.get('note', []):
            print(f"    /note=\"{v}\"")
        for v in q.get('db_xref', []):
            print(f"    /db_xref=\"{v}\"")

        print("  -- Evidence --")
        for k in evidence_keys:
            for v in q.get(k, []):
                print(f"    /{k}=\"{v}\"")

        seq = get_qual(f, 'translation')
        if seq:
            print()
            print(f"  -- Protein ({len(seq)} aa) --")
            print(f"    {seq[:60]}{'...' if len(seq) > 60 else ''}")

        covered = set(identity_keys + function_keys + evidence_keys +
                       ['note', 'db_xref', 'translation'])
        remaining = {k: v for k, v in q.items() if k not in covered}
        if remaining:
            print("  -- Other --")
            for k, vals in remaining.items():
                for v in vals:
                    print(f"    /{k}=\"{v}\"")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python genbank_locus.py <file> <locus_tag>")
        sys.exit(1)
    extract_locus(sys.argv[1], sys.argv[2])
