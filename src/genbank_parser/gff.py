"""Export a GenBank file to standard GFF3 format."""

from __future__ import annotations

import argparse
import re
import sys
import urllib.parse
from pathlib import Path
from typing import Any

from .io import get_qual, read_genbank
from .model import GenBankFeature

GFF3_HEADER = "##gff-version 3"
_ENCODE_CHARS = re.compile(r"[;\s=%&,]")


def _encode(s: str) -> str:
    """Percent-encode GFF3 attribute value special characters."""
    return _ENCODE_CHARS.sub(lambda m: urllib.parse.quote(m.group()), str(s))


def _qual_list(feature: GenBankFeature | dict[str, Any], key: str) -> list[str]:
    if hasattr(feature, "get_quals"):
        return feature.get_quals(key)
    return feature.get("qualifiers", {}).get(key, [])


def _build_attributes(
    feature: GenBankFeature | dict[str, Any],
    feat_id: str,
    parent_id: str | None = None,
) -> str:
    """Assemble GFF3 attributes for a feature."""
    attrs: list[str] = []
    if feat_id:
        attrs.append(f"ID={_encode(feat_id)}")
    if parent_id:
        attrs.append(f"Parent={_encode(parent_id)}")

    for key, output_key in (
        ("locus_tag", "locus_tag"),
        ("gene", "gene"),
        ("product", "product"),
    ):
        value = get_qual(feature, key)
        if value:
            attrs.append(f"{output_key}={_encode(value)}")
    for ec in _qual_list(feature, "EC_number"):
        attrs.append(f"ec_number={_encode(ec)}")
    for xref in _qual_list(feature, "db_xref"):
        attrs.append(f"Dbxref={_encode(xref)}")
    for note in _qual_list(feature, "note"):
        attrs.append(f"Note={_encode(note)}")
    for inf in _qual_list(feature, "inference"):
        attrs.append(f"inference={_encode(inf)}")
    return ";".join(attrs)


def _gff3_line(
    seqid: str,
    source: str,
    ftype: str,
    start: int,
    end: int,
    strand: str,
    attrs: str,
    score: str = ".",
    phase: str | int = ".",
) -> str:
    """Produce one GFF3 data line (1-based, inclusive coordinates)."""
    strand_char = strand if strand in ("+", "-") else "."
    return "\t".join(
        [
            seqid,
            source,
            ftype,
            str(start),
            str(end),
            score,
            strand_char,
            str(phase),
            attrs,
        ]
    )


def convert_to_gff3(
    filepath: str | Path,
    output_path: str | Path | None = None,
    include_fasta: bool = False,
) -> str:
    """Convert GenBank to GFF3, preserving compound segments and safe IDs."""
    doc = read_genbank(filepath)
    lines: list[str] = [GFF3_HEADER]
    for rec in doc.records:
        rec_len = (
            rec.length
            if rec.length > 0
            else max((f.end for f in rec.features), default=1)
        )
        lines.append(f"##sequence-region {rec.id} 1 {rec_len}")

    emitted_ids: set[str] = set()
    id_suffixes: dict[str, int] = {}

    def claim(desired: str) -> str:
        """Claim a globally unique ID, retaining the readable base when possible."""
        if desired not in emitted_ids:
            emitted_ids.add(desired)
            return desired
        index = id_suffixes.get(desired, 1) + 1
        candidate = f"{desired}_{index}"
        while candidate in emitted_ids:
            index += 1
            candidate = f"{desired}_{index}"
        id_suffixes[desired] = index
        emitted_ids.add(candidate)
        return candidate

    # Reserve IDs for real gene features before a CDS can request a synthetic
    # hierarchy.  This prevents a compound CDS from duplicating gene:TAG.
    gene_ids_by_feature: dict[int, str] = {}
    gene_ids_by_tag: dict[tuple[str, str], str] = {}
    for rec in doc.records:
        for feature in rec.features:
            if feature.type == "gene" and feature.locus_tag:
                gene_id = claim(f"gene:{feature.locus_tag}")
                gene_ids_by_feature[feature.feature_index] = gene_id
                gene_ids_by_tag.setdefault((rec.id, feature.locus_tag), gene_id)

    source = "genbank_gff"
    for rec in doc.records:
        seqid = rec.id
        for feature in rec.features:
            ftype = feature.type
            strand = feature.strand_symbol
            tag = feature.locus_tag
            base = tag or claim(ftype)
            real_gene_id = gene_ids_by_tag.get((rec.id, tag)) if tag else None

            if feature.is_compound and feature.join_segments:
                origin_spanning = bool(
                    rec.topology == "circular"
                    and feature.genomic_span >= rec.length > 0
                )

                if ftype == "gene":
                    parent_id = gene_ids_by_feature.get(feature.feature_index) or claim(
                        f"gene:{base}"
                    )
                    parent_attrs = _build_attributes(feature, parent_id)
                    if origin_spanning:
                        parent_attrs += ";origin_spanning=true"
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            "gene",
                            feature.start,
                            feature.end,
                            strand,
                            parent_attrs,
                        )
                    )
                elif real_gene_id:
                    parent_id = claim(f"mRNA:{base}")
                    parent_attrs = _build_attributes(feature, parent_id, real_gene_id)
                    if origin_spanning:
                        parent_attrs += ";origin_spanning=true"
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            "mRNA",
                            feature.start,
                            feature.end,
                            strand,
                            parent_attrs,
                        )
                    )
                else:
                    gene_id = claim(f"gene:{base}")
                    gene_attrs = _build_attributes(feature, gene_id)
                    if origin_spanning:
                        gene_attrs += ";origin_spanning=true"
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            "gene",
                            feature.start,
                            feature.end,
                            strand,
                            gene_attrs,
                        )
                    )
                    parent_id = (
                        claim(f"mRNA:{base}")
                        if ftype == "CDS"
                        else claim(f"{ftype}:{base}")
                    )
                    parent_attrs = _build_attributes(feature, parent_id, gene_id)
                    if origin_spanning:
                        parent_attrs += ";origin_spanning=true"
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            "mRNA" if ftype == "CDS" else ftype,
                            feature.start,
                            feature.end,
                            strand,
                            parent_attrs,
                        )
                    )

                segments = list(feature.join_segments)
                if strand == "-":
                    segments = sorted(segments, key=lambda item: item[0], reverse=True)
                else:
                    segments = sorted(segments, key=lambda item: item[0])

                current_phase = feature.codon_start - 1 if ftype == "CDS" else "."
                for segment_index, (segment_start, segment_end) in enumerate(
                    segments, 1
                ):
                    segment_id = claim(f"{ftype}:{base}.{segment_index}")
                    segment_attrs = _build_attributes(feature, segment_id, parent_id)
                    lines.append(
                        _gff3_line(
                            seqid,
                            source,
                            ftype,
                            segment_start,
                            segment_end,
                            strand,
                            segment_attrs,
                            phase=current_phase if ftype == "CDS" else ".",
                        )
                    )
                    if ftype == "CDS":
                        segment_length = segment_end - segment_start + 1
                        current_phase = (
                            3 - ((segment_length - int(current_phase)) % 3)
                        ) % 3
                continue

            if ftype == "gene" and tag:
                feature_id = gene_ids_by_feature[feature.feature_index]
            else:
                feature_id = claim(f"{ftype}:{base}")
            ordinary_parent_id = (
                real_gene_id if ftype == "CDS" and real_gene_id else None
            )
            attrs = _build_attributes(feature, feature_id, ordinary_parent_id)
            phase = feature.codon_start - 1 if ftype == "CDS" else "."
            lines.append(
                _gff3_line(
                    seqid,
                    source,
                    ftype,
                    feature.start,
                    feature.end,
                    strand,
                    attrs,
                    phase=phase,
                )
            )

    if include_fasta and any(len(record.seq) > 0 for record in doc.records):
        lines.append("##FASTA")
        for rec in doc.records:
            if len(rec.seq) > 0:
                lines.append(f">{rec.id}")
                seq_str = str(rec.seq)
                lines.extend(seq_str[i : i + 60] for i in range(0, len(seq_str), 60))

    output = "\n".join(lines) + "\n"
    if output_path:
        out_p = Path(output_path)
        out_p.write_text(output, encoding="utf-8")
        feat_count = sum(1 for line in lines if line and not line.startswith("#"))
        print(f"Wrote {feat_count} GFF3 records to {output_path}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert GenBank file to standard GFF3 format."
    )
    parser.add_argument("input", help="Input GenBank file (.gbff/.gbk/.gb)")
    parser.add_argument("output", nargs="?", help="Output GFF3 path (default: stdout)")
    parser.add_argument(
        "--include-fasta", action="store_true", help="Append ##FASTA block"
    )
    args = parser.parse_args()
    output = convert_to_gff3(args.input, args.output, args.include_fasta)
    if not args.output:
        sys.stdout.write(output)


if __name__ == "__main__":
    main()
