"""Genetic-code-aware codon usage analysis and RSCU calculation."""

from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path
from typing import Any

from Bio.Data import CodonTable

from .io import read_genbank


def _table_for_feature(feature: Any) -> tuple[int, Any] | None:
    raw = feature.get_qual("transl_table")
    table_id: int | None = 11 if not raw else None
    if raw:
        try:
            table_id = int(raw)
        except ValueError:
            table_id = None
    if table_id is None or table_id not in CodonTable.unambiguous_dna_by_id:
        return None
    return table_id, CodonTable.unambiguous_dna_by_id[table_id]


def _synonymous_codons(table: Any) -> dict[str, list[str]]:
    families: dict[str, list[str]] = collections.defaultdict(list)
    for codon, amino_acid in table.forward_table.items():
        families[amino_acid].append(codon)
    return dict(families)


def _table_metrics(codon_counts: dict[str, int], table: Any) -> dict[str, float]:
    rscu: dict[str, float] = {}
    for synonyms in _synonymous_codons(table).values():
        family_total = sum(codon_counts.get(codon, 0) for codon in synonyms)
        for codon in synonyms:
            rscu[codon] = (
                codon_counts.get(codon, 0) * len(synonyms) / family_total
                if family_total
                else 0.0
            )
    return rscu


def analyze_codon_usage(
    filepath: str | Path,
    min_len_aa: int = 100,
    output_path: str | Path | None = None,
    include_pseudo: bool = False,
) -> dict[str, Any]:
    """Analyze sense codons separately for every translation table encountered."""
    doc = read_genbank(filepath)
    if not any(len(record.seq) > 0 for record in doc.records):
        print(
            "ERROR: No ORIGIN sequences found. Codon usage requires nucleotide sequences.",
            file=sys.stderr,
        )
        raise ValueError("Codon usage requires nucleotide sequences")

    codon_counts: dict[str, int] = collections.defaultdict(int)
    table_counts: dict[int, dict[str, int]] = collections.defaultdict(
        lambda: collections.defaultdict(int)
    )
    table_objects: dict[int, Any] = {}
    total_codons = 0
    sense_codons = 0
    terminal_stops = 0
    internal_stops = 0
    ambiguous_codons = 0
    cds_evaluated = 0
    cds_excluded_partial = 0
    cds_excluded_pseudo = 0
    cds_excluded_nontriplet = 0
    unknown_translation_tables: set[str] = set()
    translation_tables: set[int] = set()

    pos1_gc = pos2_gc = pos3_gc = 0
    pos3_syn_gc = pos3_syn_total = 0

    for record in doc.records:
        for feature in record.cds_features:
            if feature.is_pseudo and not include_pseudo:
                cds_excluded_pseudo += 1
                continue
            if feature.is_partial:
                cds_excluded_partial += 1
                continue

            table_info = _table_for_feature(feature)
            if table_info is None:
                raw_table = feature.get_qual("transl_table") or "11"
                unknown_translation_tables.add(raw_table)
                continue
            table_id, table = table_info
            translation_tables.add(table_id)
            table_objects[table_id] = table

            coding_nt = str(feature.extract(record.seq)).upper()[
                feature.codon_start - 1 :
            ]
            if len(coding_nt) < min_len_aa * 3:
                cds_excluded_nontriplet += 1
                continue
            if len(coding_nt) % 3:
                cds_excluded_nontriplet += 1
                continue

            cds_evaluated += 1
            counts_for_table = table_counts[table_id]
            stop_codons = set(table.stop_codons)
            codons = [coding_nt[i : i + 3] for i in range(0, len(coding_nt), 3)]
            for index, codon in enumerate(codons):
                if not all(base in "ACGT" for base in codon):
                    ambiguous_codons += 1
                    continue
                if codon in stop_codons:
                    if index == len(codons) - 1:
                        terminal_stops += 1
                    else:
                        internal_stops += 1
                    continue

                amino_acid = table.forward_table.get(codon)
                if amino_acid is None:
                    ambiguous_codons += 1
                    continue
                counts_for_table[codon] += 1
                codon_counts[codon] += 1
                total_codons += 1
                sense_codons += 1
                if codon[0] in "GC":
                    pos1_gc += 1
                if codon[1] in "GC":
                    pos2_gc += 1
                if codon[2] in "GC":
                    pos3_gc += 1
                synonyms = _synonymous_codons(table).get(amino_acid, [])
                if len(synonyms) > 1:
                    pos3_syn_total += 1
                    if codon[2] in "GC":
                        pos3_syn_gc += 1

    if unknown_translation_tables:
        print(
            "WARNING: ignored unknown translation table(s): "
            + ", ".join(sorted(unknown_translation_tables)),
            file=sys.stderr,
        )

    rscu_by_table = {
        str(table_id): _table_metrics(
            dict(table_counts[table_id]), table_objects[table_id]
        )
        for table_id in sorted(translation_tables)
    }
    # Keep the original flat fields useful for the common single-table case.
    rscu = next(iter(rscu_by_table.values()), {}) if len(rscu_by_table) == 1 else {}
    gc1 = 100.0 * pos1_gc / sense_codons if sense_codons else 0.0
    gc2 = 100.0 * pos2_gc / sense_codons if sense_codons else 0.0
    gc3 = 100.0 * pos3_gc / sense_codons if sense_codons else 0.0
    gc3s = 100.0 * pos3_syn_gc / pos3_syn_total if pos3_syn_total else 0.0

    print("=" * 70)
    print("  CODON USAGE & RSCU ANALYSIS")
    print("=" * 70)
    print(f"  File            : {filepath}")
    print(f"  CDSs evaluated  : {cds_evaluated} (min length: {min_len_aa} aa)")
    print(
        f"  Translation tbl : {', '.join(map(str, sorted(translation_tables))) or 'none'}"
    )
    print(f"  Sense codons    : {sense_codons:,}")
    print(f"  Terminal stops  : {terminal_stops:,}")
    print(f"  GC1 / GC2 / GC3 : {gc1:.1f}% / {gc2:.1f}% / {gc3:.1f}%")
    print(f"  GC3s (synonymous): {gc3s:.1f}%")
    print()

    rows: list[dict[str, Any]] = []
    for table_id in sorted(translation_tables):
        table = table_objects[table_id]
        table_rscu = rscu_by_table[str(table_id)]
        for amino_acid, synonyms in sorted(_synonymous_codons(table).items()):
            for codon in sorted(synonyms):
                count = table_counts[table_id].get(codon, 0)
                per_k = count / total_codons * 1000 if total_codons else 0.0
                r_value = table_rscu.get(codon, 0.0)
                print(
                    f"  {table_id:>3d}  {amino_acid:2s}  {codon:5s}  {count:>7,d}  {per_k:>7.1f}  {r_value:>6.2f}"
                )
                rows.append(
                    {
                        "TranslationTable": table_id,
                        "AminoAcid": amino_acid,
                        "Codon": codon,
                        "Count": count,
                        "PerThousand": f"{per_k:.2f}",
                        "RSCU": f"{r_value:.3f}",
                    }
                )

    if output_path:
        with Path(output_path).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "TranslationTable",
                    "AminoAcid",
                    "Codon",
                    "Count",
                    "PerThousand",
                    "RSCU",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nWrote codon usage table to {output_path}")

    return {
        "total_codons": total_codons,
        "sense_codons": sense_codons,
        "terminal_stop_codons": terminal_stops,
        "internal_stop_codons": internal_stops,
        "ambiguous_codons": ambiguous_codons,
        "cds_evaluated": cds_evaluated,
        "cds_excluded_partial": cds_excluded_partial,
        "cds_excluded_pseudo": cds_excluded_pseudo,
        "cds_excluded_nontriplet": cds_excluded_nontriplet,
        "translation_tables_encountered": sorted(translation_tables),
        "unknown_translation_tables": sorted(unknown_translation_tables),
        "gc1": gc1,
        "gc2": gc2,
        "gc3": gc3,
        "gc3s": gc3s,
        "codon_counts": dict(codon_counts),
        "codon_counts_by_table": {str(k): dict(v) for k, v in table_counts.items()},
        "rscu": rscu,
        "rscu_by_table": rscu_by_table,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate genetic-code-aware codon usage and RSCU statistics."
    )
    parser.add_argument("input", help="Input GenBank file with ORIGIN sequences")
    parser.add_argument(
        "--min-len", type=int, default=100, help="Minimum CDS length in amino acids"
    )
    parser.add_argument("--output", help="Output TSV file path")
    parser.add_argument(
        "--include-pseudo",
        action="store_true",
        help="Include pseudogenes in codon counts",
    )
    args = parser.parse_args()
    analyze_codon_usage(
        args.input,
        min_len_aa=args.min_len,
        output_path=args.output,
        include_pseudo=args.include_pseudo,
    )


if __name__ == "__main__":
    main()
