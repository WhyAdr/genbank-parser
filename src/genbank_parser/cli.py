"""Unified command-line interface (gbparse) for GenBank feature table parsing and mining."""
from __future__ import annotations

import argparse
import sys
from typing import Sequence

from .bakta import batch_summary
from .codon import analyze_codon_usage
from .compare import compare_genomes
from .crispr import detect_crispr
from .diff import diff_annotations
from .discover import discover_clusters
from .extract import export_annotations_tsv
from .fasta import export_protein_fasta
from .functional import analyze_functional
from .gff import convert_to_gff3
from .locus import inspect_locus
from .metadata import extract_metadata
from .neighborhood import extract_neighborhood
from .phylo import extract_phylogenomic_markers
from .query import search_features
from .region import extract_region
from .sequence import extract_sequences
from .validate import validate


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gbparse",
        description="Unified GenBank feature parser, validation, and genomic evidence engine.",
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True, help="Analysis command to run")

    # 1. validate
    p_val = subparsers.add_parser("validate", help="Validate GenBank structure and biological translation semantics")
    p_val.add_argument("input", help="Input GenBank file")
    p_val.add_argument("--json", action="store_true", help="Output findings as JSON")

    # 2. summary / metadata
    p_meta = subparsers.add_parser("summary", aliases=["metadata"], help="Extract record metadata and summary stats")
    p_meta.add_argument("input", help="Input GenBank file")

    # 3. extract
    p_ext = subparsers.add_parser("extract", help="Export tab-delimited annotation TSV")
    p_ext.add_argument("input", help="Input GenBank file")
    p_ext.add_argument("output", nargs="?", help="Output TSV path (default: stdout)")

    # 4. search
    p_srch = subparsers.add_parser("search", help="Search features by gene, product, KO, EC, Pfam, or regex")
    p_srch.add_argument("input", help="Input GenBank file")
    p_srch.add_argument("--gene", help="Exact gene name")
    p_srch.add_argument("--product", help="Product substring")
    p_srch.add_argument("--ko", help="KEGG KO identifier")
    p_srch.add_argument("--ec", help="EC number")
    p_srch.add_argument("--cog", help="COG identifier")
    p_srch.add_argument("--pfam", help="Pfam domain ID")
    p_srch.add_argument("--feature", dest="ftype", help="Feature type (e.g. CDS, tRNA)")
    p_srch.add_argument("--gene-regex", help="Regex on /gene")
    p_srch.add_argument("--product-regex", help="Regex on /product")
    p_srch.add_argument("--format", choices=["text", "tsv", "csv", "json"], default="text")
    p_srch.add_argument("--output", help="Output file path")

    # 5. locus
    p_loc = subparsers.add_parser("locus", help="Deep-dive single-locus qualifier viewer")
    p_loc.add_argument("input", help="Input GenBank file")
    p_loc.add_argument("locus_tag", help="Target locus tag or gene name")

    # 6. neighborhood
    p_neigh = subparsers.add_parser("neighborhood", help="Extract genomic neighborhood window around target locus")
    p_neigh.add_argument("input", help="Input GenBank file")
    p_neigh.add_argument("locus_tag", help="Target locus tag or gene name")
    p_neigh.add_argument("window", nargs="?", type=int, default=5, help="Window size (+/- N genes, default: 5)")

    # 7. region
    p_reg = subparsers.add_parser("region", help="Extract genomic sub-region with coordinate rebasing")
    p_reg.add_argument("input", help="Input GenBank file")
    p_reg.add_argument("--locus", help="Target locus tag")
    p_reg.add_argument("--record", help="Record ID")
    p_reg.add_argument("--start", type=int, help="Start coordinate")
    p_reg.add_argument("--end", type=int, help="End coordinate")
    p_reg.add_argument("--flank-genes", type=int, default=0, help="Flanking genes")
    p_reg.add_argument("--flank-bp", type=int, default=0, help="Flanking bp")
    p_reg.add_argument("--rebase", action="store_true", help="Rebase coordinates to start at 1")
    p_reg.add_argument("--output", help="Output file path (.gbk/.fna)")

    # 8. fasta
    p_fa = subparsers.add_parser("fasta", help="Export all CDS translations as protein FASTA")
    p_fa.add_argument("input", help="Input GenBank file")
    p_fa.add_argument("output", nargs="?", help="Output FASTA path (default: stdout)")

    # 9. sequence
    p_seq = subparsers.add_parser("sequence", help="Extract genome nucleotide (.fna) and CDS nucleotide (.ffn) sequences")
    p_seq.add_argument("input", help="Input GenBank file")
    p_seq.add_argument("--fna", help="Genome FASTA output path")
    p_seq.add_argument("--ffn", help="CDS nucleotide output path")

    # 10. codon
    p_cod = subparsers.add_parser("codon", help="Calculate codon usage, RSCU, and positional GC metrics")
    p_cod.add_argument("input", help="Input GenBank file")
    p_cod.add_argument("--min-len", type=int, default=100, help="Min CDS length in aa")
    p_cod.add_argument("--output", help="Output TSV path")
    p_cod.add_argument("--include-pseudo", action="store_true", help="Include pseudogenes")

    # 11. functional
    p_func = subparsers.add_parser("functional", help="Profile COG distributions and metabolic pathway completeness")
    p_func.add_argument("input", help="Input GenBank file")
    p_func.add_argument("--format", choices=["tsv", "json"], default="tsv")
    p_func.add_argument("--pathways-only", action="store_true", help="Skip COG section")
    p_func.add_argument("--cog-only", action="store_true", help="Skip pathway section")

    # 12. discover
    p_disc = subparsers.add_parser("discover", help="Mine mobilome islands, operons, and dark-matter clusters")
    p_disc.add_argument("input", help="Input GenBank file")
    p_disc.add_argument("--cluster-gap", type=int, default=5000)
    p_disc.add_argument("--operon-gap", type=int, default=150)
    p_disc.add_argument("--min-weight", type=int, default=1)
    p_disc.add_argument("--format", choices=["text", "json", "tsv"], default="text")
    p_disc.add_argument("--rules", help="Custom YAML/JSON ruleset file")

    # 13. compare
    p_comp = subparsers.add_parser("compare", help="Multi-genome marker presence/absence matrix")
    p_comp.add_argument("genomes", nargs="+", help="GenBank files or directories")
    p_comp.add_argument("--targets", required=True, help="Comma-separated marker names (e.g. 'ladA,ssuD,K20938')")
    p_comp.add_argument("--output", help="Output matrix TSV path")

    # 14. diff
    p_diff = subparsers.add_parser("diff", help="Compare two annotation versions of the same genome")
    p_diff.add_argument("old_file", help="Reference / older GenBank file")
    p_diff.add_argument("new_file", help="Updated / newer GenBank file")
    p_diff.add_argument("--format", choices=["text", "json"], default="text")
    p_diff.add_argument("--output", help="Output file path")

    # 15. phylo
    p_phy = subparsers.add_parser("phylo", help="Extract phylogenomic core and housekeeping marker genes")
    p_phy.add_argument("input", help="Input GenBank file")
    p_phy.add_argument("--markers", choices=["core", "housekeeping", "all"], default="all")
    p_phy.add_argument("--output-dir", help="Output directory for marker FASTAs")

    # 16. crispr
    p_cris = subparsers.add_parser("crispr", help="Detect CRISPR repeat arrays and Cas gene clusters")
    p_cris.add_argument("input", help="Input GenBank file")
    p_cris.add_argument("--window", type=int, default=15000, help="Window for array-Cas linking")

    # 17. gff
    p_gff = subparsers.add_parser("gff", help="Convert GenBank file to standard GFF3 format")
    p_gff.add_argument("input", help="Input GenBank file")
    p_gff.add_argument("output", nargs="?", help="Output GFF3 path")
    p_gff.add_argument("--include-fasta", action="store_true", help="Append ##FASTA block")

    # 18. batch-summary
    p_bat = subparsers.add_parser("batch-summary", help="Parse Bakta summaries across multiple isolates")
    p_bat.add_argument("inputs", nargs="+", help="GenBank files or directories")
    p_bat.add_argument("--csv", default="bakta_summary.csv")
    p_bat.add_argument("--tsv", default="bakta_summary.tsv")
    p_bat.add_argument("--md", default="bakta_summary.md")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)

    cmd = args.subcommand
    if cmd == "validate":
        validate(args.input, json_mode=args.json)
    elif cmd in ("summary", "metadata"):
        extract_metadata(args.input)
    elif cmd == "extract":
        export_annotations_tsv(args.input, args.output)
    elif cmd == "search":
        search_features(
            args.input,
            gene=args.gene,
            product=args.product,
            ko=args.ko,
            ec=args.ec,
            cog=args.cog,
            pfam=args.pfam,
            ftype=args.ftype,
            gene_regex=args.gene_regex,
            product_regex=args.product_regex,
            format_type=args.format,
            output_path=args.output,
        )
    elif cmd == "locus":
        inspect_locus(args.input, args.locus_tag)
    elif cmd == "neighborhood":
        extract_neighborhood(args.input, args.locus_tag, window=args.window)
    elif cmd == "region":
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
    elif cmd == "fasta":
        export_protein_fasta(args.input, args.output)
    elif cmd == "sequence":
        extract_sequences(args.input, out_fna=args.fna, out_ffn=args.ffn)
    elif cmd == "codon":
        analyze_codon_usage(args.input, min_len_aa=args.min_len, output_path=args.output, include_pseudo=args.include_pseudo)
    elif cmd == "functional":
        analyze_functional(args.input, format_type=args.format, pathways_only=args.pathways_only, cog_only=args.cog_only)
    elif cmd == "discover":
        discover_clusters(
            args.input,
            cluster_gap=args.cluster_gap,
            operon_gap=args.operon_gap,
            min_weight=args.min_weight,
            format_type=args.format,
            rules_file=args.rules,
        )
    elif cmd == "compare":
        targets = [t.strip() for t in args.targets.split(",") if t.strip()]
        compare_genomes(args.genomes, targets, output_path=args.output)
    elif cmd == "diff":
        diff_annotations(args.old_file, args.new_file, format_type=args.format, output_path=args.output)
    elif cmd == "phylo":
        extract_phylogenomic_markers(args.input, marker_set=args.markers, output_dir=args.output_dir)
    elif cmd == "crispr":
        detect_crispr(args.input, window=args.window)
    elif cmd == "gff":
        out = convert_to_gff3(args.input, args.output, include_fasta=args.include_fasta)
        if not args.output:
            sys.stdout.write(out)
    elif cmd == "batch-summary":
        batch_summary(args.inputs, csv_out=args.csv, tsv_out=args.tsv, md_out=args.md)

    return 0


if __name__ == "__main__":
    sys.exit(main())
