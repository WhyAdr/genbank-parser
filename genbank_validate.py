#!/usr/bin/env python3
"""GenBank feature table structural validator."""
import sys, collections
from genbank_parser import parse_features, get_qual


def validate(filepath):
    features = parse_features(filepath)
    if not features:
        print("ERROR: No features parsed. Check file format.")
        sys.exit(1)

    type_counts = collections.Counter(f['type'] for f in features)
    strand_counts = collections.Counter(f['strand'] for f in features)
    all_quals = set()
    for f in features:
        all_quals.update(f['qualifiers'].keys())

    cdss = [f for f in features if f['type'] == 'CDS']

    # Coordinate span
    all_starts = [f['start'] for f in features]
    all_ends   = [f['end']   for f in features]
    span_min, span_max = min(all_starts), max(all_ends)

    # CDS length distribution
    cds_lens = [f['end'] - f['start'] + 1 for f in cdss]

    # Locus tags
    locus_tags = []
    for f in features:
        lt = f['qualifiers'].get('locus_tag', [])
        locus_tags.extend(lt)
    unique_tags = set(locus_tags)
    dup_tags = [t for t, c in collections.Counter(locus_tags).items() if c > 2]

    # Named genes vs hypothetical
    named = [f for f in cdss if f['qualifiers'].get('gene')]
    products = [get_qual(f, 'product') for f in cdss]
    hypothetical = [p for p in products
                    if 'hypothetical' in p.lower()
                    or 'domain-containing' in p.lower()
                    or 'DUF' in p]

    print("=" * 60)
    print("GENBANK FEATURE TABLE -- STRUCTURAL REPORT")
    print("=" * 60)
    print(f"File: {filepath}")
    print(f"Total features parsed: {len(features)}")
    print()
    print("-- Feature type counts --")
    for ft, c in type_counts.most_common():
        print(f"  {ft:20s}  {c}")
    print()
    print("-- Strand distribution --")
    for s, c in strand_counts.items():
        print(f"  {s}  {c}")
    print()
    print("-- Coordinate span --")
    print(f"  Min start : {span_min:,}")
    print(f"  Max end   : {span_max:,}")
    print(f"  Span      : {span_max - span_min + 1:,} bp")
    print()
    print("-- CDS statistics --")
    print(f"  Count        : {len(cdss)}")
    if cds_lens:
        print(f"  Min length   : {min(cds_lens):,} bp")
        print(f"  Max length   : {max(cds_lens):,} bp")
        print(f"  Mean length  : {sum(cds_lens)/len(cds_lens):,.0f} bp")
        print(f"  Median length: {sorted(cds_lens)[len(cds_lens)//2]:,} bp")
    print(f"  Named genes  : {len(named)}")
    print(f"  Hypothetical / DUF: {len(hypothetical)}")
    print()
    print("-- Locus tags --")
    print(f"  Total occurrences : {len(locus_tags)}")
    print(f"  Unique tags       : {len(unique_tags)}")
    if dup_tags:
        print(f"  WARNING: Tags with >2 occurrences (possible issue): {dup_tags[:10]}")
    if unique_tags:
        tags_sorted = sorted(unique_tags)
        print(f"  First tag : {tags_sorted[0]}")
        print(f"  Last tag  : {tags_sorted[-1]}")
    print()
    print("-- Qualifier keys found --")
    print(f"  {', '.join(sorted(all_quals))}")
    print()

    # Validation checks
    issues = []
    for f in cdss:
        if not f['qualifiers'].get('product'):
            issues.append(f"  Line {f['line']}: CDS at {f['start']}..{f['end']} missing /product")
        if not f['qualifiers'].get('translation'):
            issues.append(f"  Line {f['line']}: CDS at {f['start']}..{f['end']} missing /translation")
        tlen = len(get_qual(f, 'translation'))
        nlen = f['end'] - f['start'] + 1
        expected_aa = (nlen // 3) - 1  # minus stop codon
        if tlen > 0 and abs(tlen - expected_aa) > 2:
            issues.append(
                f"  Line {f['line']}: CDS {get_qual(f, 'locus_tag', '?')} -- "
                f"nucleotide length {nlen} predicts ~{expected_aa} aa, got {tlen} aa"
            )
    if issues:
        print(f"-- Validation issues ({len(issues)}) --")
        for iss in issues[:20]:
            print(iss)
    else:
        print("-- Validation: ALL CHECKS PASSED --")
    print("=" * 60)


if __name__ == '__main__':
    validate(sys.argv[1] if len(sys.argv) > 1 else input("File path: "))
