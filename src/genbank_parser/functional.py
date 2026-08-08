"""Unified functional profiling: COG distribution + metabolic pathway completeness."""
from __future__ import annotations

import argparse
import collections
import csv
import io
import json
from pathlib import Path
import sys
from typing import Any

from .io import extract_xrefs, get_qual, read_genbank
from .model import GenBankFeature

COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking',
    'O': 'Post-translational modification',
    'X': 'Mobilome: prophages, transposons',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis',
    'R': 'General function prediction only',
    'S': 'Function unknown',
}

PATHWAYS = {
    "Glycolysis (Embden-Meyerhof)": [
        {"name": "Glucokinase / Hexokinase", "refs": ["K00845", "K00844", "EC:2.7.1.2", "EC:2.7.1.1"]},
        {"name": "Phosphoglucose isomerase", "refs": ["K01810", "EC:5.3.1.9"]},
        {"name": "6-Phosphofructokinase", "refs": ["K00850", "K21071", "EC:2.7.1.11"]},
        {"name": "Fructose-bisphosphate aldolase", "refs": ["K01623", "K01624", "K11645", "EC:4.1.2.13"]},
        {"name": "Triosephosphate isomerase", "refs": ["K01803", "EC:5.3.1.1"]},
        {"name": "Glyceraldehyde-3-phosphate dehydrogenase", "refs": ["K00134", "K10705", "EC:1.2.1.12"]},
        {"name": "Phosphoglycerate kinase", "refs": ["K00927", "EC:2.7.2.3"]},
        {"name": "Phosphoglycerate mutase", "refs": ["K01834", "K01837", "K15633", "EC:5.4.2.11", "EC:5.4.2.12"]},
        {"name": "Enolase", "refs": ["K01689", "EC:4.2.1.11"]},
        {"name": "Pyruvate kinase", "refs": ["K00873", "EC:2.7.1.40"]},
    ],
    "TCA Cycle": [
        {"name": "Citrate synthase", "refs": ["K01647", "K01659", "EC:2.3.3.1", "EC:2.3.3.5", "EC:2.3.3.16"]},
        {"name": "Aconitate hydratase", "refs": ["K01681", "K01682", "EC:4.2.1.3"]},
        {"name": "Isocitrate dehydrogenase", "refs": ["K00031", "K00030", "EC:1.1.1.42", "EC:1.1.1.41"]},
        {"name": "2-Oxoglutarate dehydrogenase", "refs": ["K00164", "K00174", "EC:1.2.4.2", "EC:1.2.7.3"]},
        {"name": "Succinyl-CoA synthetase", "refs": ["K01899", "K01900", "K01902", "EC:6.2.1.4", "EC:6.2.1.5"]},
        {"name": "Succinate dehydrogenase", "refs": ["K00234", "K00235", "K00236", "K00237", "EC:1.3.5.1"]},
        {"name": "Fumarate hydratase", "refs": ["K01676", "K01677", "K01678", "K01679", "EC:4.2.1.2"]},
        {"name": "Malate dehydrogenase", "refs": ["K00024", "K00025", "K00026", "EC:1.1.1.37", "EC:1.1.1.38", "EC:1.1.5.4"]},
    ],
    "Glyoxylate Bypass": [
        {"name": "Isocitrate lyase", "refs": ["K01637", "EC:4.1.3.1"]},
        {"name": "Malate synthase", "refs": ["K01638", "EC:2.3.3.9"]},
    ],
    "Fatty Acid Beta-Oxidation": [
        {"name": "Acyl-CoA synthetase", "refs": ["K01897", "EC:6.2.1.3"]},
        {"name": "Acyl-CoA dehydrogenase", "refs": ["K00249", "K00248", "K00255", "EC:1.3.8.7", "EC:1.3.8.1", "EC:1.3.8.8"]},
        {"name": "Enoyl-CoA hydratase", "refs": ["K01692", "EC:4.2.1.17"]},
        {"name": "3-Hydroxyacyl-CoA dehydrogenase", "refs": ["K00022", "EC:1.1.1.35"]},
        {"name": "Acetyl-CoA C-acyltransferase (thiolase)", "refs": ["K00632", "EC:2.3.1.16"]},
    ],
}


def analyze_functional(
    filepath: str | Path,
    format_type: str = 'tsv',
    pathways_only: bool = False,
    cog_only: bool = False,
) -> dict[str, Any]:
    doc = read_genbank(filepath)
    cdss = [f for f in doc.all_features if f.type == 'CDS']

    # Cross-reference harvesting
    all_kos: set[str] = set()
    all_ecs: set[str] = set()
    all_cogs: list[str] = []

    for f in cdss:
        xr = extract_xrefs(f)
        all_kos.update(xr['kegg_kos'])
        all_ecs.update(f"EC:{ec}" for ec in xr['ec_numbers'])
        all_cogs.extend(xr['cog_ids'])

    all_annotated = all_kos | all_ecs

    # COG breakdown
    cog_counts: dict[str, int] = collections.defaultdict(int)
    for c in all_cogs:
        # Bakta / standard COGs may have category letter or ID
        cog_counts[c] += 1

    # Pathway evaluations
    pathway_results: dict[str, Any] = {}
    for p_name, steps in PATHWAYS.items():
        present_steps = 0
        step_details = []
        for s in steps:
            matched = [r for r in s['refs'] if r in all_annotated]
            if matched:
                present_steps += 1
            step_details.append({
                'name': s['name'],
                'matched': matched,
                'status': 'PRESENT' if matched else 'ABSENT',
            })
        pct = (100.0 * present_steps / len(steps)) if steps else 0.0
        pathway_results[p_name] = {
            'steps_present': present_steps,
            'steps_total': len(steps),
            'completeness_pct': round(pct, 1),
            'steps': step_details,
        }

    report = {
        'file': str(filepath),
        'total_cds': len(cdss),
        'unique_kos': len(all_kos),
        'unique_ecs': len(all_ecs),
        'total_cogs': len(all_cogs),
        'pathways': pathway_results,
    }

    if format_type == 'json':
        print(json.dumps(report, indent=2))
        return report

    # Formatted TSV / Text
    if not pathways_only:
        print("=" * 70)
        print("  FUNCTIONAL ANNOTATION SUMMARY")
        print("=" * 70)
        print(f"  File         : {filepath}")
        print(f"  Total CDSs   : {len(cdss):,}")
        print(f"  Unique KOs   : {len(all_kos):,}")
        print(f"  Unique ECs   : {len(all_ecs):,}")
        print(f"  COG entries  : {len(all_cogs):,}")
        print()

    if not cog_only:
        print("-- Pathway Completeness Profiles --")
        for p_name, p_data in pathway_results.items():
            bar = '#' * int(p_data['completeness_pct'] / 10)
            print(f"  {p_name:35s} {p_data['steps_present']:2d}/{p_data['steps_total']:2d} ({p_data['completeness_pct']:>5.1f}%)  {bar}")
        print("=" * 70)

    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified functional profiling (COG + metabolic pathway completeness).")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--format', choices=['tsv', 'json'], default='tsv', help="Output format (default: tsv)")
    parser.add_argument('--pathways-only', action='store_true', help="Skip COG section")
    parser.add_argument('--cog-only', action='store_true', help="Skip pathway section")
    args = parser.parse_args()

    analyze_functional(args.input, format_type=args.format, pathways_only=args.pathways_only, cog_only=args.cog_only)


if __name__ == '__main__':
    main()
