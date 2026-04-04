#!/usr/bin/env python3
"""Summarize COG functional category distribution from parsed features."""
import sys, collections
from genbank_parser import parse_features

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

def cog_distribution(filepath):
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']

    cog_cats = collections.Counter()
    for f in cdss:
        notes = f['qualifiers'].get('note', [])
        for n in notes:
            if n.startswith('COG:') and len(n.split(':')[1]) <= 2 and not n.split(':')[1].startswith('COG'):
                for ch in n.split(':')[1]:
                    cog_cats[ch] += 1

    print("-- COG Functional Category Distribution --")
    print(f"{'Cat':>4}  {'Count':>5}  Description")
    print("-" * 60)
    for cat, count in cog_cats.most_common():
        desc = COG_CATEGORIES.get(cat, 'Unknown')
        print(f"  {cat:>2}  {count:>5}  {desc}")
    print()
    uncategorized = len(cdss) - sum(cog_cats.values())
    print(f"  CDS with COG assignment : {sum(cog_cats.values())}")
    print(f"  CDS without COG         : {uncategorized}")


if __name__ == '__main__':
    cog_distribution(sys.argv[1])
