#!/usr/bin/env python3
"""Unified functional profiling: COG distribution + metabolic pathway completeness.

Replaces genbank_cog.py and adds pathway completeness mapping using
KEGG KOs and EC numbers extracted via extract_xrefs().

Usage:
    python genbank_functional.py INPUT.gbff                    # TSV (default)
    python genbank_functional.py INPUT.gbff --format json      # JSON
    python genbank_functional.py INPUT.gbff --pathways-only    # skip COG section
    python genbank_functional.py INPUT.gbff --cog-only         # skip pathway section
"""
import sys, os, argparse, collections, json, csv, io
from genbank_parser import parse_features, get_qual, extract_xrefs

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
        {"name": "Glucokinase / Hexokinase",
         "refs": ["K00845", "K00844", "EC:2.7.1.2", "EC:2.7.1.1"]},
        {"name": "Phosphoglucose isomerase",
         "refs": ["K01810", "EC:5.3.1.9"]},
        {"name": "6-Phosphofructokinase",
         "refs": ["K00850", "K21071", "EC:2.7.1.11"]},
        {"name": "Fructose-bisphosphate aldolase",
         "refs": ["K01623", "K01624", "K11645", "EC:4.1.2.13"]},
        {"name": "Triosephosphate isomerase",
         "refs": ["K01803", "EC:5.3.1.1"]},
        {"name": "Glyceraldehyde-3-phosphate dehydrogenase",
         "refs": ["K00134", "K10705", "EC:1.2.1.12"]},
        {"name": "Phosphoglycerate kinase",
         "refs": ["K00927", "EC:2.7.2.3"]},
        {"name": "Phosphoglycerate mutase",
         "refs": ["K01834", "K01837", "K15633", "EC:5.4.2.11", "EC:5.4.2.12"]},
        {"name": "Enolase",
         "refs": ["K01689", "EC:4.2.1.11"]},
        {"name": "Pyruvate kinase",
         "refs": ["K00873", "EC:2.7.1.40"]},
    ],
    "TCA Cycle": [
        {"name": "Citrate synthase",
         "refs": ["K01647", "K01659", "EC:2.3.3.1", "EC:2.3.3.5", "EC:2.3.3.16"]},
        {"name": "Aconitate hydratase",
         "refs": ["K01681", "K01682", "EC:4.2.1.3", "EC:4.2.1.99"]},
        {"name": "Isocitrate dehydrogenase",
         "refs": ["K00030", "K00031", "EC:1.1.1.41", "EC:1.1.1.42"]},
        {"name": "2-Oxoglutarate dehydrogenase (E1)",
         "refs": ["K00164", "K00165", "EC:1.2.4.2"]},
        {"name": "Dihydrolipoamide succinyltransferase (E2)",
         "refs": ["K00658", "EC:2.3.1.61"]},
        {"name": "Succinyl-CoA synthetase",
         "refs": ["K01899", "K01900", "K01902", "K01903", "EC:6.2.1.4", "EC:6.2.1.5"]},
        {"name": "Succinate dehydrogenase",
         "refs": ["K00234", "K00235", "K00236", "K00237", "EC:1.3.5.1", "EC:1.3.5.4"]},
        {"name": "Fumarate hydratase",
         "refs": ["K01676", "K01679", "EC:4.2.1.2"]},
        {"name": "Malate dehydrogenase",
         "refs": ["K00024", "K00025", "K00026", "EC:1.1.1.37"]},
    ],
    "Pentose Phosphate Pathway": [
        {"name": "Glucose-6-phosphate dehydrogenase",
         "refs": ["K00036", "EC:1.1.1.49"]},
        {"name": "6-Phosphogluconolactonase",
         "refs": ["K01053", "K07404", "EC:3.1.1.31"]},
        {"name": "6-Phosphogluconate dehydrogenase",
         "refs": ["K00033", "EC:1.1.1.44"]},
        {"name": "Ribulose-5-phosphate 3-epimerase",
         "refs": ["K01783", "EC:5.1.3.1"]},
        {"name": "Ribose-5-phosphate isomerase",
         "refs": ["K01807", "K01808", "EC:5.3.1.6"]},
        {"name": "Transketolase",
         "refs": ["K00615", "EC:2.2.1.1"]},
        {"name": "Transaldolase",
         "refs": ["K00616", "K13810", "EC:2.2.1.2"]},
    ],
    "Entner-Doudoroff Pathway": [
        {"name": "Glucose-6-phosphate dehydrogenase",
         "refs": ["K00036", "EC:1.1.1.49"]},
        {"name": "6-Phosphogluconolactonase",
         "refs": ["K01053", "K07404", "EC:3.1.1.31"]},
        {"name": "Phosphogluconate dehydratase",
         "refs": ["K01690", "EC:4.2.1.12"]},
        {"name": "KDPG aldolase",
         "refs": ["K01625", "EC:4.1.2.14"]},
    ],
    "Glyoxylate Bypass": [
        {"name": "Isocitrate lyase",
         "refs": ["K01637", "EC:4.1.3.1"]},
        {"name": "Malate synthase",
         "refs": ["K01638", "EC:2.3.3.9"]},
    ],
    "Shikimate Pathway": [
        {"name": "DAHP synthase",
         "refs": ["K01626", "K03856", "K13853", "EC:2.5.1.54"]},
        {"name": "3-Dehydroquinate synthase",
         "refs": ["K01735", "K13829", "K13830", "EC:4.2.3.4"]},
        {"name": "3-Dehydroquinate dehydratase",
         "refs": ["K03785", "K03786", "EC:4.2.1.10"]},
        {"name": "Shikimate dehydrogenase",
         "refs": ["K00014", "K13830", "K13832", "EC:1.1.1.25"]},
        {"name": "Shikimate kinase",
         "refs": ["K00891", "K13829", "EC:2.7.1.71"]},
        {"name": "EPSP synthase",
         "refs": ["K00800", "K13830", "EC:2.5.1.19"]},
        {"name": "Chorismate synthase",
         "refs": ["K01736", "EC:4.2.3.5"]},
    ],
    "Peptidoglycan Biosynthesis": [
        {"name": "MurA (UDP-GlcNAc enolpyruvyl transferase)",
         "refs": ["K00790", "EC:2.5.1.7"]},
        {"name": "MurB (UDP-MurNAc reductase)",
         "refs": ["K00075", "EC:1.3.1.98", "EC:1.1.1.158"]},
        {"name": "MurC (L-Ala ligase)",
         "refs": ["K01924", "EC:6.3.2.8"]},
        {"name": "MurD (D-Glu ligase)",
         "refs": ["K01925", "EC:6.3.2.9"]},
        {"name": "MurE (m-DAP/L-Lys ligase)",
         "refs": ["K01928", "EC:6.3.2.13"]},
        {"name": "MurF (D-Ala-D-Ala ligase)",
         "refs": ["K01929", "EC:6.3.2.10"]},
        {"name": "MraY (phospho-MurNAc-pentapeptide transferase)",
         "refs": ["K01000", "EC:2.7.8.13"]},
        {"name": "MurG (GlcNAc transferase -> Lipid II)",
         "refs": ["K02563", "EC:2.4.1.227"]},
    ],
    "Mixed-Acid Fermentation": [
        {"name": "Pyruvate formate-lyase",
         "refs": ["K00656", "EC:2.3.1.54"]},
        {"name": "Lactate dehydrogenase",
         "refs": ["K00016", "K03778", "EC:1.1.1.27", "EC:1.1.1.28"]},
        {"name": "Phosphotransacetylase",
         "refs": ["K00625", "K15024", "EC:2.3.1.8"]},
        {"name": "Acetate kinase",
         "refs": ["K00925", "EC:2.7.2.1"]},
        {"name": "Alcohol dehydrogenase",
         "refs": ["K04072", "K00001", "K13953", "EC:1.1.1.1"]},
        {"name": "Fumarate reductase",
         "refs": ["K00244", "K00245", "K00246", "K00247", "EC:1.3.5.4"]},
    ],
    "Assimilatory Nitrate Reduction": [
        {"name": "Assimilatory nitrate reductase",
         "refs": ["K00367", "K10534", "K10535", "EC:1.7.99.4", "EC:1.7.7.2"]},
        {"name": "Assimilatory nitrite reductase",
         "refs": ["K00362", "K00363", "EC:1.7.1.4", "EC:1.7.1.15"]},
    ],
    "Assimilatory Sulfate Reduction": [
        {"name": "Sulfate adenylyltransferase",
         "refs": ["K00958", "K00955", "EC:2.7.7.4"]},
        {"name": "Adenylylsulfate kinase",
         "refs": ["K00860", "EC:2.7.1.25"]},
        {"name": "PAPS reductase",
         "refs": ["K00390", "EC:1.8.4.8"]},
        {"name": "Sulfite reductase (NADPH)",
         "refs": ["K00392", "K00380", "K00381", "EC:1.8.1.2"]},
    ],
}

def cog_distribution(cdss):
    """Extract COG category counts from parsed CDS features."""
    cog_cats = collections.Counter()
    for f in cdss:
        all_notes = f['qualifiers'].get('note', []) + f['qualifiers'].get('db_xref', [])
        for n in all_notes:
            if n.startswith('COG:') and len(n.split(':')[1]) <= 2 \
               and not n.split(':')[1].startswith('COG'):
                for ch in n.split(':')[1]:
                    cog_cats[ch] += 1
    return cog_cats


def collect_genome_markers(cdss):
    """Single-pass extraction of all KO and EC identifiers from CDS features."""
    ko_set = set()
    ec_set = set()
    for f in cdss:
        xr = extract_xrefs(f)
        ko_set.update(xr['kegg_kos'])
        ec_set.update(xr['ec_numbers'])
    return ko_set, ec_set


def assess_pathways(ko_set, ec_set):
    """Evaluate pathway completeness using OR-semantics per step.

    A step is 'present' if ANY of its refs match either the genome's
    KO set or EC set. Returns a list of pathway result dicts.
    """
    results = []
    for pathway_name, steps in PATHWAYS.items():
        present = []
        missing = []
        for step in steps:
            kos_in_step = {r for r in step['refs'] if r.startswith('K')}
            ecs_in_step = {r.replace('EC:', '') for r in step['refs']
                          if r.startswith('EC:')}
            found = bool(kos_in_step & ko_set) or bool(ecs_in_step & ec_set)
            if found:
                present.append(step['name'])
            else:
                missing.append(step['name'])
        completeness = len(present) / len(steps) if steps else 0.0
        results.append({
            'pathway': pathway_name,
            'total_steps': len(steps),
            'present': len(present),
            'missing_count': len(missing),
            'completeness': round(completeness, 3),
            'present_steps': present,
            'missing_steps': missing,
        })
    return results


def emit_tsv(cog_cats, pathway_results, total_cds, mode):
    """Emit structured TSV to stdout."""
    if mode != 'pathways-only':
        print("## COG Functional Category Distribution")
        print("category\tcount\tdescription")
        for cat, count in cog_cats.most_common():
            desc = COG_CATEGORIES.get(cat, 'Unknown')
            print(f"{cat}\t{count}\t{desc}")
        assigned = sum(cog_cats.values())
        print(f"#\tCDS_with_COG={assigned}\tCDS_without_COG={total_cds - assigned}")
        print()

    if mode != 'cog-only':
        print("## Pathway Completeness")
        print("pathway\ttotal_steps\tpresent\tmissing\tcompleteness\tmissing_steps")
        for r in pathway_results:
            missing_str = '; '.join(r['missing_steps']) if r['missing_steps'] else '-'
            print(f"{r['pathway']}\t{r['total_steps']}\t{r['present']}\t"
                  f"{r['missing_count']}\t{r['completeness']:.1%}\t{missing_str}")


def emit_json(cog_cats, pathway_results, total_cds, mode):
    """Emit structured JSON to stdout."""
    result = {}
    if mode != 'pathways-only':
        assigned = sum(cog_cats.values())
        result['cog'] = {
            'total_cds': total_cds,
            'assigned': assigned,
            'unassigned': total_cds - assigned,
            'categories': {cat: {'count': count,
                                 'description': COG_CATEGORIES.get(cat, 'Unknown')}
                           for cat, count in cog_cats.most_common()},
        }
    if mode != 'cog-only':
        result['pathways'] = pathway_results
    print(json.dumps(result, indent=2, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(
        description="Unified functional profiling: COG distribution "
                    "+ metabolic pathway completeness.")
    parser.add_argument('input', help="Input GenBank file")
    parser.add_argument('--format', choices=['tsv', 'json'], default='tsv')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--pathways-only', action='store_true',
                       help="Skip COG section, report pathways only")
    group.add_argument('--cog-only', action='store_true',
                       help="Skip pathway section, report COG only")
    args = parser.parse_args()

    features = parse_features(args.input)
    cdss = [f for f in features if f['type'] == 'CDS']
    mode = 'pathways-only' if args.pathways_only else \
           'cog-only' if args.cog_only else 'full'

    cog_cats = cog_distribution(cdss)
    ko_set, ec_set = collect_genome_markers(cdss)
    pathway_results = assess_pathways(ko_set, ec_set)

    if args.format == 'json':
        emit_json(cog_cats, pathway_results, len(cdss), mode)
    else:
        emit_tsv(cog_cats, pathway_results, len(cdss), mode)


if __name__ == '__main__':
    main()
