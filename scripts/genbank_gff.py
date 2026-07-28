#!/usr/bin/env python3
"""Export a GenBank file to GFF3 format.

GFF3 is the standard interchange format accepted by most genome browsers
(IGV, JBrowse, UCSC) and annotation pipelines (MAKER, Augustus, etc.).

Handles:
  - All feature types parsed by genbank_parser.py
  - join() features: emitted as a parent mRNA/CDS with sub-features
  - Multi-contig files: correct seqid per feature
  - Proper GFF3 attribute quoting and URL-encoding

Usage:
    python genbank_gff.py INPUT.gbff [OUTPUT.gff3]
    python genbank_gff.py INPUT.gbff --include-fasta   # append ##FASTA block
"""
import sys, re, argparse, urllib.parse
from genbank_parser import parse_features, get_qual

GFF3_HEADER = "##gff-version 3"

# GFF3 attribute values that must be URL-encoded per spec
_ENCODE_CHARS = re.compile(r'[;\s=%&,]')


def _encode(s):
    """Percent-encode GFF3 attribute value special characters."""
    return _ENCODE_CHARS.sub(lambda m: urllib.parse.quote(m.group()), s)


def _qual_list(feature, key):
    return feature['qualifiers'].get(key, [])


def _build_attributes(feature, feat_id, parent_id=None):
    """Assemble GFF3 attribute string for a feature."""
    attrs = []

    if feat_id:
        attrs.append(f"ID={_encode(feat_id)}")
    if parent_id:
        attrs.append(f"Parent={_encode(parent_id)}")

    locus_tag = get_qual(feature, 'locus_tag')
    if locus_tag:
        attrs.append(f"locus_tag={_encode(locus_tag)}")

    gene = get_qual(feature, 'gene')
    if gene:
        attrs.append(f"gene={_encode(gene)}")

    product = get_qual(feature, 'product')
    if product:
        attrs.append(f"product={_encode(product)}")

    # EC_number
    for ec in _qual_list(feature, 'EC_number'):
        attrs.append(f"ec_number={_encode(ec)}")

    # db_xref
    for xref in _qual_list(feature, 'db_xref'):
        attrs.append(f"Dbxref={_encode(xref)}")

    # note
    for note in _qual_list(feature, 'note'):
        attrs.append(f"Note={_encode(note)}")

    # inference
    for inf in _qual_list(feature, 'inference'):
        attrs.append(f"inference={_encode(inf)}")

    return ';'.join(attrs)


def _gff3_line(seqid, source, ftype, start, end, strand, attrs, score='.', phase='.'):
    """Produce one GFF3 data line (1-based, inclusive coords)."""
    strand_char = strand if strand in ('+', '-') else '.'
    return '\t'.join([
        seqid, source, ftype,
        str(start), str(end),
        score, strand_char, str(phase),
        attrs,
    ])


def convert_to_gff3(filepath, output_path=None, include_fasta=False):
    features = parse_features(filepath)

    lines = [GFF3_HEADER]

    # Emit ##sequence-region directives per contig
    contig_extents = {}
    for f in features:
        c = f['contig']
        if c not in contig_extents:
            contig_extents[c] = [f['start'], f['end']]
        else:
            contig_extents[c][0] = min(contig_extents[c][0], f['start'])
            contig_extents[c][1] = max(contig_extents[c][1], f['end'])
    for contig, (lo, hi) in contig_extents.items():
        lines.append(f"##sequence-region {contig} {lo} {hi}")

    # Counters for unique IDs
    id_counter = {}

    def _unique_id(prefix):
        id_counter[prefix] = id_counter.get(prefix, 0) + 1
        return f"{prefix}_{id_counter[prefix]:05d}"

    for f in features:
        ftype    = f['type']
        seqid    = f['contig']
        strand   = f['strand']
        source   = 'genbank_gff'
        locus_tag = get_qual(f, 'locus_tag') or _unique_id(ftype)

        # join() features: emit as top-level gene + per-segment sub-features
        if 'join_segments' in f:
            gene_id = f"gene:{locus_tag}"
            gene_attrs = _build_attributes(f, gene_id)
            # Emit spanning gene record
            lines.append(_gff3_line(seqid, source, 'gene',
                                     f['start'], f['end'], strand, gene_attrs))
            # Emit parent mRNA/feature record
            parent_type = 'mRNA' if ftype in ('CDS', 'mRNA') else ftype
            parent_id   = f"{parent_type}:{locus_tag}"
            parent_attrs = _build_attributes(f, parent_id, gene_id)
            lines.append(_gff3_line(seqid, source, parent_type,
                                     f['start'], f['end'], strand, parent_attrs))
            # Emit one GFF3 record per join segment
            for seg_idx, (seg_s, seg_e) in enumerate(f['join_segments'], 1):
                seg_id    = f"{ftype}:{locus_tag}.{seg_idx}"
                seg_attrs = f"ID={seg_id};Parent={parent_id}"
                # Phase: estimate for CDS (0 for first segment, carried forward)
                phase = 0
                lines.append(_gff3_line(seqid, source, ftype,
                                         seg_s, seg_e, strand, seg_attrs,
                                         phase=phase))
        else:
            feat_id   = f"{ftype}:{locus_tag}"
            feat_attrs = _build_attributes(f, feat_id)
            lines.append(_gff3_line(seqid, source, ftype,
                                     f['start'], f['end'], strand, feat_attrs))

    # Optional ##FASTA block (genome sequences)
    if include_fasta:
        genome_seqs = _read_origin(filepath)
        if genome_seqs:
            lines.append('##FASTA')
            for contig, seq in genome_seqs.items():
                lines.append(f">{contig}")
                for i in range(0, len(seq), 60):
                    lines.append(seq[i:i + 60])

    output = '\n'.join(lines) + '\n'

    if output_path:
        with open(output_path, 'w', encoding='utf-8') as fh:
            fh.write(output)
        feat_count = sum(1 for ln in lines if ln and not ln.startswith('#'))
        print(f"Wrote {feat_count} GFF3 records to {output_path}")
    else:
        sys.stdout.write(output)


def _read_origin(filepath):
    """Read ORIGIN blocks and return dict[contig -> sequence]."""
    seqs = {}
    current = None
    parts   = []
    in_orig = False
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('LOCUS'):
                if current and parts:
                    seqs[current] = ''.join(parts).upper()
                current = line.split()[1]
                parts   = []
                in_orig = False
            elif line.startswith('ORIGIN'):
                in_orig = True
            elif line.startswith('//'):
                if current and parts:
                    seqs[current] = ''.join(parts).upper()
                parts   = []
                in_orig = False
            elif in_orig:
                parts.append(re.sub(r'[\d\s]', '', line))
    return seqs


def main():
    parser = argparse.ArgumentParser(
        description="Convert a GenBank file to GFF3 format."
    )
    parser.add_argument('input',           help="Input GenBank file (.gbff/.gbk/.gb)")
    parser.add_argument('output', nargs='?', help="Output GFF3 path (default: stdout)")
    parser.add_argument('--include-fasta', action='store_true',
                        help="Append ##FASTA block from ORIGIN sequences")
    args = parser.parse_args()
    convert_to_gff3(args.input, args.output, args.include_fasta)


if __name__ == '__main__':
    main()
