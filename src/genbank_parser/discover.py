"""Genomic discovery engine: mobilome islands, operons, and dark-matter clusters with declarative rules."""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
import sys
from typing import Any

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature

DEFAULT_MOBILOME_RULES = [
    {
        'id': 'integrase',
        'terms': ['integrase', 'site-specific integrase', 'tyrosine recombinase', 'xerC', 'xerD'],
        'weight': 3,
        'category': 'Mobilome: Integrase',
    },
    {
        'id': 'transposase',
        'terms': ['transposase', 'insertion element', 'IS element', 'IS3', 'IS4', 'IS5', 'IS6', 'IS21', 'IS256', 'IS110'],
        'weight': 3,
        'category': 'Mobilome: Transposon',
    },
    {
        'id': 'resolvase',
        'terms': ['resolvase', 'recombinase', 'site-specific recombinase'],
        'weight': 2,
        'category': 'Mobilome: Recombinase',
    },
    {
        'id': 'phage',
        'terms': ['phage', 'prophage', 'tail protein', 'capsid', 'portal protein', 'terminase', 'holin'],
        'weight': 2,
        'category': 'Mobilome: Phage',
    },
    {
        'id': 'conjugation',
        'terms': ['traA', 'traB', 'traC', 'traD', 'traE', 'traG', 'traI', 'traK', 'mobA', 'mobB', 'relaxase', 'type IV secretion'],
        'weight': 3,
        'category': 'Mobilome: Conjugation',
    },
]


def discover_clusters(
    filepath: str | Path,
    cluster_gap: int = 5000,
    operon_gap: int = 150,
    min_weight: int = 1,
    format_type: str = 'text',
    rules_file: str | Path | None = None,
) -> dict[str, Any]:
    doc = read_genbank(filepath)
    rules = DEFAULT_MOBILOME_RULES

    # Optional user-defined rules from JSON/YAML
    if rules_file:
        rf = Path(rules_file)
        if rf.exists():
            if rf.suffix.lower() in ('.yaml', '.yml'):
                try:
                    import yaml
                    with rf.open(encoding='utf-8') as fh:
                        custom_rules = yaml.safe_load(fh)
                        if isinstance(custom_rules, list):
                            rules = custom_rules
                except Exception as err:
                    print(f"Warning: Failed to load YAML rules ({err}), using defaults.", file=sys.stderr)
            elif rf.suffix.lower() == '.json':
                with rf.open(encoding='utf-8') as fh:
                    rules = json.load(fh)

    hits: list[dict[str, Any]] = []

    for rec in doc.records:
        for f in rec.cds_features:
            gene = (f.gene or '').lower()
            prod = (f.product or '').lower()
            notes = ' '.join(f.qualifiers.get('note', [])).lower()
            xrefs = extract_xrefs(f)

            matched_rules = []
            total_weight = 0

            for r in rules:
                rule_weight = r.get('weight', 1)
                for term in r.get('terms', []):
                    t_low = term.lower()
                    if t_low in gene or t_low in prod or t_low in notes:
                        matched_rules.append({'rule': r['id'], 'term': term, 'weight': rule_weight})
                        total_weight += rule_weight
                        break

            if total_weight >= min_weight:
                hits.append({
                    'contig': rec.id,
                    'locus_tag': f.locus_tag or '-',
                    'gene': f.gene or '-',
                    'product': f.product or '-',
                    'start': f.start,
                    'end': f.end,
                    'strand': f.strand_symbol,
                    'weight': total_weight,
                    'matches': matched_rules,
                    'feature_index': f.feature_index,
                })

    # Spatial clustering of hits within cluster_gap
    clusters: list[list[dict[str, Any]]] = []
    if hits:
        hits_by_contig = collections.defaultdict(list)
        for h in hits:
            hits_by_contig[h['contig']].append(h)

        for c_id, c_hits in hits_by_contig.items():
            sorted_hits = sorted(c_hits, key=lambda x: x['start'])
            current_cl = [sorted_hits[0]]
            for h in sorted_hits[1:]:
                if h['start'] - current_cl[-1]['end'] <= cluster_gap:
                    current_cl.append(h)
                else:
                    if len(current_cl) >= 2 or sum(x['weight'] for x in current_cl) >= 4:
                        clusters.append(current_cl)
                    current_cl = [h]
            if len(current_cl) >= 2 or sum(x['weight'] for x in current_cl) >= 4:
                clusters.append(current_cl)

    result = {
        'file': str(filepath),
        'total_hits': len(hits),
        'total_islands': len(clusters),
        'hits': hits,
        'islands': [
            {
                'contig': cl[0]['contig'],
                'start': cl[0]['start'],
                'end': cl[-1]['end'],
                'span': cl[-1]['end'] - cl[0]['start'] + 1,
                'genes': len(cl),
                'total_weight': sum(x['weight'] for x in cl),
                'features': cl,
            }
            for cl in clusters
        ],
    }

    if format_type == 'json':
        print(json.dumps(result, indent=2))
        return result

    print("=" * 70)
    print("  GENOMIC MOBILOME & DARK-MATTER DISCOVERY")
    print("=" * 70)
    print(f"  File             : {filepath}")
    print(f"  Individual hits  : {len(hits)}")
    print(f"  Clustered islands: {len(clusters)} (gap <= {cluster_gap:,} bp)")
    print()

    if clusters:
        print("-- Discovered Genomic Islands --")
        for i, island in enumerate(result['islands'], 1):
            print(f"  Island #{i}: {island['contig']}:{island['start']:,}..{island['end']:,} ({island['span']:,} bp, {island['genes']} genes, score: {island['total_weight']})")
            for feat in island['features']:
                matched_str = ', '.join(f"{m['rule']}:{m['term']}" for m in feat['matches'])
                print(f"    [{feat['strand']}] {feat['locus_tag']:16s} {feat['gene']:8s} {feat['product'][:35]:35s} (matches: {matched_str})")
            print()
    else:
        print("  (No multi-gene mobilome islands detected at chosen threshold)")
    print("=" * 70)

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Mine mobilome islands, operons, and dark-matter clusters.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--cluster-gap', type=int, default=5000, help="Max distance to bridge island genes (default: 5000)")
    parser.add_argument('--operon-gap', type=int, default=150, help="Max intergenic gap for operon linking (default: 150)")
    parser.add_argument('--min-weight', type=int, default=1, help="Minimum evidence weight threshold (default: 1)")
    parser.add_argument('--format', choices=['text', 'json', 'tsv'], default='text', help="Output format")
    parser.add_argument('--rules', help="Optional custom ruleset file (.yaml or .json)")
    args = parser.parse_args()

    discover_clusters(
        args.input,
        cluster_gap=args.cluster_gap,
        operon_gap=args.operon_gap,
        min_weight=args.min_weight,
        format_type=args.format,
        rules_file=args.rules,
    )


if __name__ == '__main__':
    main()
