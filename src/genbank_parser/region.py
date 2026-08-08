"""Extract valid local GenBank sub-regions while preserving feature locations."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .io import read_genbank
from .model import GenBankFeature


def _combine_location(
    parts: list[FeatureLocation], template: object
) -> FeatureLocation | CompoundLocation:
    """Build a location using the source operator and the surviving parts."""
    if len(parts) == 1:
        return parts[0]
    return CompoundLocation(parts, operator=getattr(template, "operator", "join"))


def _clip_and_shift_location(
    location: FeatureLocation | CompoundLocation,
    region_start_1based: int,
    region_end_1based: int,
) -> FeatureLocation | CompoundLocation | None:
    """Intersect a linear location with a region and shift it to local coordinates.

    Coordinates in this helper are 1-based inclusive at the API boundary and
    0-based half-open internally, matching Biopython.  Clipping necessarily
    makes fuzzy boundary positions exact; the original location is retained by
    :func:`extract_region` in a qualifier.
    """
    region_start = region_start_1based - 1
    region_end = region_end_1based
    surviving: list[FeatureLocation] = []
    source_parts = list(getattr(location, "parts", (location,)))
    for part in source_parts:
        start = max(region_start, int(part.start))
        end = min(region_end, int(part.end))
        if end > start:
            surviving.append(
                FeatureLocation(
                    start - region_start, end - region_start, strand=part.strand
                )
            )

    if not surviving:
        return None
    return _combine_location(surviving, location)


def _wrap_sequence(sequence: Seq, start_1based: int, end_1based: int) -> Seq:
    """Extract an unwrapped interval from a circular sequence."""
    if not sequence or end_1based < start_1based:
        return Seq("")
    text = str(sequence)
    length = len(text)
    count = end_1based - start_1based + 1
    start = (start_1based - 1) % length
    repeated = text[start:] + text * ((count + length - 1) // length + 1)
    return Seq(repeated[:count])


def _circular_location(
    location: FeatureLocation | CompoundLocation,
    region_start_1based: int,
    region_end_1based: int,
    record_length: int,
) -> FeatureLocation | CompoundLocation | None:
    """Clip a location against an unwrapped circular window."""
    if record_length <= 0:
        return None
    region_start = region_start_1based - 1
    region_end = region_end_1based
    mapped: list[FeatureLocation] = []

    for part in list(getattr(location, "parts", (location,))):
        part_start = int(part.start)
        part_end = int(part.end)
        lower = math.floor((region_start - part_end) / record_length) - 1
        upper = math.ceil((region_end - part_start) / record_length) + 1
        candidates: list[FeatureLocation] = []
        for copy_index in range(lower, upper + 1):
            copy_start = part_start + copy_index * record_length
            copy_end = part_end + copy_index * record_length
            start = max(region_start, copy_start)
            end = min(region_end, copy_end)
            if end > start:
                candidates.append(
                    FeatureLocation(
                        start - region_start, end - region_start, strand=part.strand
                    )
                )
        if part.strand == -1:
            candidates.sort(key=lambda p: int(p.start), reverse=True)
        else:
            candidates.sort(key=lambda p: int(p.start))
        mapped.extend(candidates)

    if not mapped:
        return None
    return _combine_location(mapped, location)


def _resolve_gene_window(
    rec_features: list[GenBankFeature],
    target: GenBankFeature,
    flank_genes: int,
    circular: bool,
    record_length: int,
) -> tuple[int, int]:
    """Return an unwrapped feature window around a resolved CDS."""
    cdss = sorted((f for f in rec_features if f.type == "CDS"), key=lambda f: f.start)
    if target.type != "CDS" or not cdss:
        raise ValueError("Resolved target is not a CDS in a record with CDS features")
    try:
        index = [f.feature_index for f in cdss].index(target.feature_index)
    except ValueError:
        raise ValueError("Resolved target is not present in CDS ordering") from None

    if flank_genes <= 0:
        return target.start, target.end

    if not circular:
        left = cdss[max(0, index - flank_genes)]
        right = cdss[min(len(cdss) - 1, index + flank_genes)]
        return left.start, right.end

    left = cdss[(index - flank_genes) % len(cdss)]
    right = cdss[(index + flank_genes) % len(cdss)]
    end = right.end
    if end < left.start:
        end += record_length
    return left.start, end


def extract_region(
    filepath: str | Path,
    locus_tag: str | None = None,
    record_id: str | None = None,
    start: int | None = None,
    end: int | None = None,
    flank_genes: int = 0,
    flank_bp: int = 0,
    rebase: bool = False,
    output_path: str | Path | None = None,
) -> SeqRecord:
    """Extract a standalone region with feature coordinates local to its sequence.

    ``rebase`` is retained for CLI compatibility.  Local coordinates are now
    always required for a valid standalone GenBank record; when it is false the
    original parent coordinates are preserved as metadata instead.
    """
    if flank_genes < 0 or flank_bp < 0:
        raise ValueError("flank_genes and flank_bp must be non-negative")

    doc = read_genbank(filepath)
    target_rec = None
    raw_start: int
    raw_end: int

    if locus_tag:
        match = doc.find_locus(locus_tag)
        if match is None:
            # A gene-name lookup follows the same CDS preference as locus tags.
            for rec in doc.records:
                gene_matches = [
                    f
                    for f in rec.features
                    if f.gene and f.gene.casefold() == locus_tag.casefold()
                ]
                feat = next((f for f in gene_matches if f.type == "CDS"), None)
                if feat is not None:
                    match = (rec, feat)
                    break
        if match is None:
            print(
                f"ERROR: Locus tag '{locus_tag}' not found in {filepath}",
                file=sys.stderr,
            )
            raise ValueError(f"Locus tag '{locus_tag}' not found")
        target_rec, feat = match
        circular = target_rec.topology == "circular"
        raw_start, raw_end = _resolve_gene_window(
            target_rec.features, feat, flank_genes, circular, target_rec.length
        )
    elif record_id:
        target_rec = doc.get_record(record_id)
        if target_rec is None:
            raise ValueError(f"Record '{record_id}' not found in {filepath}")
        circular = target_rec.topology == "circular"
        raw_start = start if start is not None else 1
        raw_end = end if end is not None else target_rec.length
        if circular and raw_end < raw_start:
            raw_end += target_rec.length
    else:
        if not doc.records:
            raise ValueError(f"No records found in {filepath}")
        target_rec = doc.records[0]
        circular = target_rec.topology == "circular"
        raw_start = start if start is not None else 1
        raw_end = end if end is not None else target_rec.length
        if circular and raw_end < raw_start:
            raw_end += target_rec.length

    assert target_rec is not None
    if flank_bp:
        raw_start -= flank_bp
        raw_end += flank_bp

    if not circular:
        raw_start = max(1, raw_start)
        raw_end = min(target_rec.length, max(raw_start, raw_end))
    elif target_rec.length <= 0:
        raise ValueError(
            "Cannot extract a circular region from a record without a sequence length"
        )

    if raw_end < raw_start:
        raise ValueError("Region end must not precede region start")

    if circular:
        sub_seq = _wrap_sequence(target_rec.seq, raw_start, raw_end)
    else:
        sub_seq = target_rec.seq[raw_start - 1 : raw_end] if target_rec.seq else Seq("")

    region_features: list[SeqFeature] = []
    for feature in target_rec.features:
        if circular:
            new_location = _circular_location(
                feature.location, raw_start, raw_end, target_rec.length
            )
        else:
            new_location = _clip_and_shift_location(
                feature.location, raw_start, raw_end
            )
        if new_location is None:
            continue

        qualifiers = {key: list(values) for key, values in feature.qualifiers.items()}
        qualifiers.setdefault("original_record", [target_rec.id])
        qualifiers.setdefault("original_location", [str(feature.location)])
        qualifiers.setdefault("original_start", [str(feature.start)])
        qualifiers.setdefault("original_end", [str(feature.end)])
        region_features.append(
            SeqFeature(location=new_location, type=feature.type, qualifiers=qualifiers)
        )

    span = raw_end - raw_start + 1
    sub_rec_id = f"{target_rec.id}_region_{raw_start}_{raw_end}"
    sub_record = SeqRecord(
        seq=sub_seq,
        id=sub_rec_id,
        name=sub_rec_id[:16],
        description=(
            f"Extracted local region from {target_rec.id}:{raw_start}..{raw_end}"
            f" (rebase flag={rebase})"
        ),
        features=region_features,
    )
    sub_record.annotations["molecule_type"] = "DNA"
    sub_record.annotations["topology"] = "linear"
    sub_record.annotations["parent_record"] = target_rec.id
    sub_record.annotations["parent_start"] = raw_start
    sub_record.annotations["parent_end"] = raw_end
    sub_record.annotations["parent_topology"] = target_rec.topology or "linear"

    print("=" * 70)
    print("  GENOMIC REGION EXTRACTION")
    print("=" * 70)
    print(f"  Parent Record   : {target_rec.id}")
    print(f"  Region Span     : {raw_start:,} .. {raw_end:,} ({span:,} bp)")
    print("  Feature coords  : local to extracted sequence")
    print(f"  Rebase flag     : {rebase} (compatibility flag)")
    print(f"  Features sliced : {len(region_features)}")
    print("=" * 70)

    if output_path:
        out_p = Path(output_path)
        out_fmt = (
            "fasta" if out_p.suffix.lower() in (".fasta", ".fna", ".fa") else "genbank"
        )
        with out_p.open("w", encoding="utf-8") as handle:
            SeqIO.write(sub_record, handle, out_fmt)
        print(f"Wrote region file to {output_path}")

    return sub_record


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract genomic regions with local feature coordinates."
    )
    parser.add_argument("input", help="Input GenBank file")
    parser.add_argument("--locus", help="Target locus tag or gene name")
    parser.add_argument("--record", help="Target record/contig ID")
    parser.add_argument("--start", type=int, help="Start coordinate (1-based)")
    parser.add_argument("--end", type=int, help="End coordinate (1-based)")
    parser.add_argument(
        "--flank-genes", type=int, default=0, help="Number of flanking CDSs"
    )
    parser.add_argument("--flank-bp", type=int, default=0, help="Number of flanking bp")
    parser.add_argument(
        "--rebase",
        action="store_true",
        help="Compatibility flag; output coordinates are always local",
    )
    parser.add_argument("--output", help="Output file path (.gbk, .gbff, or .fna)")
    args = parser.parse_args()

    extract_region(
        args.input,
        locus_tag=args.locus,
        record_id=args.record,
        start=args.start,
        end=args.end,
        flank_genes=args.flank_genes,
        flank_bp=args.flank_bp,
        rebase=args.rebase,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
