"""Export a GenBank file to standard GFF3 format."""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
import urllib.parse
from typing import Any

from .io import get_qual, read_genbank
from .model import GenBankDocument, GenBankFeature

GFF3_HEADER = "##gff-version 3"
_ENCODE_CHARS = re.compile(r'[;\s=%&,]')


def _encode(s: str) -> str:
    """Percent-encode GFF3 attribute value special characters."""
    return _ENCODE_CHARS.sub(lambda m: urllib.parse.quote(m.group()), str(s))


def _qual_list(feature: GenBankFeature | dict[str, Any], key: str) -> list[str]:
    if hasattr(feature, 'get_quals'):
        return feature.get_quals(key)
    return feature.get('qualifiers', {}).get(key, [])


def _build_attributes(feature: GenBankFeature | dict[str, Any], feat_id: str, parent_id: str | None = None) -> str:
    """Assemble GFF3 attribute string for a feature."""
    attrs: list[str] = []

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

    for ec in _qual_list(feature, 'EC_number'):
        attrs.append(f"ec_number={_encode(ec)}")

    for xref in _qual_list(feature, 'db_xref'):
        attrs.append(f"Dbxref={_encode(xref)}")

    for note in _qual_list(feature, 'note'):
        attrs.append(f"Note={_encode(note)}")

    for inf in _qual_list(feature, 'inference'):
        attrs.append(f"inference={_encode(inf)}")

    return ';'.join(attrs)


def _gff3_line(
    seqid: str,
    source: str,
    ftype: str,
    start: int,
    end: int,
    strand: str,
    attrs: str,
    score: str = '.',
    phase: str | int = '.',
) -> str:
    """Produce one GFF3 data line (1-based, inclusive coordinates)."""
    strand_char = strand if strand in ('+', '-') else '.'
    return '\t'.join([
        seqid,
        source,
        ftype,
        str(start),
        str(end),
        score,
        strand_char,
        str(phase),
        attrs,
    ])


def convert_to_gff3(
    filepath: str | Path,
    output_path: str | Path | None = None,
    include_fasta: bool = False,
) -> str:
    """Convert GenBank file to GFF3 format preserving compound locations and sequence regions."""
    doc = read_genbank(filepath)
    lines: list[str] = [GFF3_HEADER]

    # Emit ##sequence-region directives using true record lengths
    for rec in doc.records:
        rec_len = rec.length if rec.length > 0 else (max([f.end for f in rec.features], default=1))
        lines.append(f"##sequence-region {rec.id} 1 {rec_len}")

    id_counter: dict[str, int] = {}

    def _unique_id(prefix: str) -> str:
        id_counter[prefix] = id_counter.get(prefix, 0) + 1
        return f"{prefix}_{id_counter[prefix]:05d}"

    source = 'genbank_gff'

    for rec in doc.records:
        seqid = rec.id
        for f in rec.features:
            ftype = f.type
            strand = f.strand_symbol
            tag = f.locus_tag
            base_id = tag if tag else _unique_id(ftype)

            # Compound locations (join/order): emit parent feature + child segments
            if f.is_compound and len(f.join_segments) > 0:
                gene_id = f"gene:{base_id}"
                gene_attrs = _build_attributes(f, gene_id)
                lines.append(_gff3_line(seqid, source, 'gene', f.start, f.end, strand, gene_attrs))

                parent_type = 'mRNA' if ftype in ('CDS', 'mRNA') else ftype
                parent_id = f"{parent_type}:{base_id}"
                parent_attrs = _build_attributes(f, parent_id, gene_id)
                lines.append(_gff3_line(seqid, source, parent_type, f.start, f.end, strand, parent_attrs))

                # Order segments 5' -> 3' for biological phase calculation
                segments = list(f.join_segments)
                if strand == '-':
                    # 5' to 3' order on negative strand is descending coordinate order
                    segments_ordered = sorted(segments, key=lambda s: s[0], reverse=True)
                else:
                    segments_ordered = sorted(segments, key=lambda s: s[0])

                current_phase = 0
                for seg_idx, (seg_s, seg_e) in enumerate(segments_ordered, 1):
                    seg_len = seg_e - seg_s + 1
                    seg_id = f"{ftype}:{base_id}.{seg_idx}"
                    seg_attrs = f"ID={_encode(seg_id)};Parent={_encode(parent_id)}"
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            ftype,
                            seg_s,
                            seg_e,
                            strand,
                            seg_attrs,
                            phase=current_phase if ftype == 'CDS' else '.',
                        )
                    )
                    # Next phase: (3 - ((L - P) % 3)) % 3
                    if ftype == 'CDS':
                        current_phase = (3 - ((seg_len - current_phase) % 3)) % 3
            else:
                # Ordinary non-compound feature
                feat_id = f"{ftype}:{base_id}"
                feat_attrs = _build_attributes(f, feat_id)
                lines.append(_gff3_line(seqid, source, ftype, f.start, f.end, strand, feat_attrs, phase=0 if ftype == 'CDS' else '.'))

    # Optional ##FASTA block from record sequences
    if include_fasta:
        has_seq = any(len(r.seq) > 0 for r in doc.records)
        if has_seq:
            lines.append('##FASTA')
            for rec in doc.records:
                if len(rec.seq) > 0:
                    lines.append(f">{rec.id}")
                    seq_str = str(rec.seq)
                    for i in range(0, len(seq_str), 60):
                        lines.append(seq_str[i:i + 60])

    output = '\n'.join(lines) + '\n'

    if output_path:
        out_p = Path(output_path)
        out_p.write_text(output, encoding='utf-8')
        feat_count = sum(1 for ln in lines if ln and not ln.startswith('#'))
        print(f"Wrote {feat_count} GFF3 records to {output_path}")

    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a GenBank file to GFF3 format.")
    parser.add_argument('input', help="Input GenBank file (.gbff/.gbk/.gb)")
    parser.add_argument('output', nargs='?', help="Output GFF3 path (default: stdout)")
    parser.add_argument('--include-fasta', action='store_true', help="Append ##FASTA block from record sequences")
    args = parser.parse_args()

    out = convert_to_gff3(args.input, args.output, args.include_fasta)
    if not args.output:
        sys.stdout.write(out)


if __name__ == '__main__':
    main()
