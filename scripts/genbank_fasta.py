#!/usr/bin/env python3
"""Export all CDS translations as a multi-FASTA file."""
import sys
from genbank_parser import parse_features, get_qual


def export_fasta(filepath, outfasta):
    features = parse_features(filepath)
    cdss = [f for f in features if f['type'] == 'CDS']
    count = 0
    with open(outfasta, 'w') as fh:
        for f in cdss:
            tag = get_qual(f, 'locus_tag', 'unknown')
            gene = get_qual(f, 'gene')
            product = get_qual(f, 'product')
            seq = get_qual(f, 'translation')
            if seq:
                header = f">{tag}"
                if gene:
                    header += f" gene={gene}"
                header += f" product={product}"
                fh.write(header + '\n')
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i+70] + '\n')
                count += 1
    print(f"Exported {count} protein sequences to {outfasta}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Export all CDS as protein FASTA.")
    parser.add_argument('input',  help="Input GenBank file")
    parser.add_argument('output', nargs='?', help="Output FASTA path (default: <input>_proteins.faa)")
    args = parser.parse_args()
    import os
    outfile = args.output or os.path.splitext(args.input)[0] + '_proteins.faa'
    export_fasta(args.input, outfile)

if __name__ == '__main__':
    main()
