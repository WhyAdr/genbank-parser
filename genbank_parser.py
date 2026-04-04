#!/usr/bin/env python3
"""
Canonical GenBank feature table parser.
Single source of truth -- all analysis scripts import from here.

Handles:
  - Multi-contig/multi-record files (LOCUS splitting)
  - complement() locations
  - Partial features (>start, <end)
  - Multi-line qualifier continuation (19-22 space indent)
  - Graceful skip + stderr warning for join()/order() complex locations

Does NOT handle:
  - join(), order(), complement(join(...)) locations (warns, skips)
"""
import re, collections, sys


def parse_features(filepath):
    """Parse a GenBank-style feature table into structured records.

    Returns a list of dicts, each with keys:
        type (str), start (int), end (int), strand ('+'/'-'),
        contig (str), qualifiers (defaultdict(list)), line (int)

    Each feature includes a 'contig' key from its parent LOCUS record,
    preventing cross-contig clustering in downstream analyses.
    """
    features = []
    current_contig = '_unknown_'
    current = None
    skipped = collections.defaultdict(int)
    in_origin = False

    locus_re = re.compile(r'^LOCUS\s+(\S+)')

    feature_re = re.compile(
        r'^     (gene|CDS|tRNA|rRNA|tmRNA|ncRNA|regulatory|misc_feature'
        r'|repeat_region|source|mRNA|sig_peptide|mat_peptide'
        r'|misc_binding|mobile_element|oriT|gap)'
        r'\s+(complement\()?<?(\d+)\.\.>?(\d+)\)?\s*$'
    )

    complex_re = re.compile(r'^     \w+\s+.*(?:join|order)\(')
    qualifier_re = re.compile(r'^\s{19,22}/(\w+)(?:="?(.*?)"?\s*)?$')
    continuation_re = re.compile(r'^\s{19,22}(?!/)\S')

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\r\n')

            if line.startswith('ORIGIN'):
                in_origin = True
                if current:
                    features.append(current)
                    current = None
                continue

            if line.startswith('//'):
                in_origin = False
                if current:
                    features.append(current)
                    current = None
                continue

            if in_origin:
                continue

            lm = locus_re.match(line)
            if lm:
                if current:
                    features.append(current)
                    current = None
                current_contig = lm.group(1)
                continue

            if complex_re.match(line):
                skipped[current_contig] += 1
                if current:
                    features.append(current)
                    current = None
                continue

            fm = feature_re.match(line)
            if fm:
                if current:
                    features.append(current)
                current = {
                    'type': fm.group(1),
                    'start': int(fm.group(3)),
                    'end': int(fm.group(4)),
                    'strand': '-' if fm.group(2) else '+',
                    'contig': current_contig,
                    'qualifiers': collections.defaultdict(list),
                    'line': lineno,
                }
                continue

            qm = qualifier_re.match(line)
            if qm and current:
                key, val = qm.group(1), qm.group(2) or ''
                current['qualifiers'][key].append(val)
                continue

            if current and continuation_re.match(line) and current['qualifiers']:
                last_key = list(current['qualifiers'].keys())[-1]
                current['qualifiers'][last_key][-1] += line.strip()

    if current:
        features.append(current)

    if skipped:
        total = sum(skipped.values())
        print(f"WARNING: Skipped {total} complex-location features "
              f"(join/order not supported):", file=sys.stderr)
        for contig, count in skipped.items():
            print(f"  {contig}: {count} skipped", file=sys.stderr)

    return features


def get_qual(feature, key, default=''):
    """Convenience: get first qualifier value for a key, or default."""
    vals = feature['qualifiers'].get(key, [])
    return vals[0] if vals else default


def get_notes(feature, prefix):
    """Extract /note values matching a prefix (e.g., 'COG:', 'KEGG:')."""
    return [n for n in feature['qualifiers'].get('note', []) if n.startswith(prefix)]
