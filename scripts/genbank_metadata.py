#!/usr/bin/env python3
"""Extract top-level GenBank metadata (LOCUS, DEFINITION, SOURCE, ORGANISM)
without parsing features or sequences. Fast scan for large files."""
import sys, re


def extract_metadata(filepath):
    records = []
    current = {}

    locus_re = re.compile(
        r'^LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(\S+)\s+(\w+)\s+(\w+)\s+(.*)'
    )

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')

            lm = locus_re.match(line)
            if lm:
                if current:
                    records.append(current)
                current = {
                    'locus': lm.group(1),
                    'length': int(lm.group(2)),
                    'mol_type': lm.group(3),
                    'topology': lm.group(4),
                    'division': lm.group(5),
                    'date': lm.group(6).strip(),
                    'definition': '',
                    'accession': '',
                    'version': '',
                    'source': '',
                    'organism': '',
                    'strain': '',
                }
                continue

            if not current:
                continue

            if line.startswith('DEFINITION'):
                current['definition'] = line[12:].strip()
            elif line.startswith('ACCESSION'):
                current['accession'] = line[12:].strip()
            elif line.startswith('VERSION'):
                current['version'] = line[12:].strip()
            elif line.startswith('SOURCE'):
                current['source'] = line[12:].strip()
            elif line.startswith('  ORGANISM'):
                current['organism'] = line[12:].strip()

            if '/strain=' in line:
                m = re.search(r'/strain="([^"]+)"', line)
                if m and not current.get('strain'):
                    current['strain'] = m.group(1)

    if current:
        records.append(current)

    print("=" * 70)
    print("  GENBANK METADATA REPORT")
    print("=" * 70)
    print(f"  File    : {filepath}")
    print(f"  Records : {len(records)}")
    print()

    total_len = 0
    for i, rec in enumerate(records, 1):
        total_len += rec['length']
        size = rec['length']
        if size >= 1_000_000:
            size_str = f"{size/1_000_000:.2f} Mb"
        elif size >= 1_000:
            size_str = f"{size/1_000:.1f} kb"
        else:
            size_str = f"{size} bp"

        print(f"  Record {i}: {rec['locus']}")
        print(f"    Length     : {size_str} ({rec['length']:,} bp)")
        print(f"    Topology   : {rec['topology']}")
        print(f"    Mol. type  : {rec['mol_type']}")
        print(f"    Definition : {rec['definition'][:80]}")
        if rec['organism']:
            print(f"    Organism   : {rec['organism']}")
        if rec['strain']:
            print(f"    Strain     : {rec['strain']}")
        print()

    if total_len >= 1_000_000:
        total_str = f"{total_len/1_000_000:.2f} Mb"
    else:
        total_str = f"{total_len:,} bp"
    print(f"  Total: {len(records)} record(s), {total_str}")
    if len(records) > 1:
        print(f"  WARNING: Multi-contig file. Downstream scripts must respect "
              f"contig boundaries.")
    print("=" * 70)

    return records

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract LOCUS header metadata from a GenBank file.")
    parser.add_argument('input', help="Input GenBank file")
    args = parser.parse_args()
    extract_metadata(args.input)

if __name__ == '__main__':
    main()
