"""Unified, pipe-safe command-line interface for ``gbparse``."""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .bakta import batch_summary
from .cli_io import (
    InputError,
    OutputError,
    eprint,
    paths_same,
    publish_directory,
    reject_input_output_collision,
    write_text,
)
from .codon import analyze_codon_usage, render_codon_tsv
from .compare import compare_data, render_compare, render_evidence_tsv
from .crispr import build_crispr_report, render_crispr_report
from .diff import diff_annotations
from .discover import discover_clusters
from .export import (
    EXPORT_FORMATS,
    ExportArtifact,
    render_export,
)
from .functional import build_functional_report, render_functional
from .gff import convert_to_gff3
from .io import GenBankInputError, read_genbank
from .locus import build_locus_report, render_locus_report
from .meor import analyze_meor
from .meor.database import load_meor_database
from .meor.report import serialize_report
from .metadata import build_metadata_report, render_metadata
from .mobilome import MobilomeError, analyze_mobilome, load_mobilome_database
from .mobilome.report import serialize_mobilome_report
from .neighborhood import build_neighborhood, serialize_neighborhood
from .operons import build_operon_report, render_operons, render_operons_gff3
from .phylo import (
    build_phylogenomic_report,
    render_phylogenomic_report,
)
from .query import (
    QueryExpressionError,
    parse_select_fields,
    query_features,
    search_features,
)
from .region import extract_region, render_region_record
from .sequence import extract_sequences
from .serializers import feature_rows, render_annotations_tsv, render_rows
from .validate import (
    build_validation_report,
    render_validation_json,
    render_validation_text,
)

EXIT_OK = 0
EXIT_VALIDATION = 1
EXIT_USAGE = 2
EXIT_INPUT = 3
EXIT_OUTPUT = 4


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _add_output_arguments(
    parser: argparse.ArgumentParser,
    *,
    formats: Sequence[str] | None = None,
    default_format: str | None = None,
    output_dir: bool = False,
) -> None:
    if formats is not None:
        parser.add_argument("--format", choices=list(formats), default=default_format or formats[0])
    parser.add_argument("--output", help="Output path (default: stdout)")
    if output_dir:
        parser.add_argument("--output-dir", help="Output directory for multi-file artifacts")
    parser.add_argument("--force", action="store_true", help="Replace an existing output")


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gbparse",
        description="Unified GenBank feature parser, validation, and genomic evidence engine.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="subcommand", required=True, help="Analysis command to run")

    p_val = subparsers.add_parser("validate", help="Validate GenBank structure and translation semantics")
    p_val.add_argument("input")
    _add_output_arguments(p_val, formats=("text", "json"), default_format="text")
    p_val.add_argument("--json", action="store_true", help="Compatibility alias for --format json")
    p_val.add_argument("--fail-on", choices=("never", "error", "warning", "info"), default="never")

    p_meta = subparsers.add_parser("summary", aliases=["metadata"], help="Extract record metadata and summary stats")
    p_meta.add_argument("input")
    _add_output_arguments(p_meta, formats=("text", "tsv", "json"), default_format="text")

    p_ext = subparsers.add_parser("extract", help="Export tab-delimited annotation TSV")
    p_ext.add_argument("input")
    p_ext.add_argument("output", nargs="?", help="Compatibility positional output path")
    p_ext.add_argument("--force", action="store_true")

    p_srch = subparsers.add_parser("search", help="Search features by annotation fields")
    p_srch.add_argument("input")
    p_srch.add_argument("--gene")
    p_srch.add_argument("--product")
    p_srch.add_argument("--ko")
    p_srch.add_argument("--ec")
    p_srch.add_argument("--cog")
    p_srch.add_argument("--pfam")
    p_srch.add_argument("--feature", dest="ftype")
    p_srch.add_argument("--gene-regex")
    p_srch.add_argument("--product-regex")
    p_srch.add_argument("--format", choices=("text", "tsv", "csv", "json"), default="text")
    p_srch.add_argument("--output")
    p_srch.add_argument("--force", action="store_true")

    p_loc = subparsers.add_parser("locus", help="Deep-dive a single locus")
    p_loc.add_argument("input")
    p_loc.add_argument("locus_tag")
    _add_output_arguments(p_loc, formats=("text", "tsv", "json"), default_format="text")

    p_neigh = subparsers.add_parser("neighborhood", help="Extract a genomic neighborhood")
    p_neigh.add_argument("input")
    p_neigh.add_argument("locus_tag")
    p_neigh.add_argument("window", nargs="?", type=_nonnegative_int, default=5)
    p_neigh.add_argument("--format", choices=("text", "tsv", "json"), default="text")
    p_neigh.add_argument("--output")
    p_neigh.add_argument("--force", action="store_true")
    p_neigh.add_argument("--visualize", action="store_true")
    p_neigh.add_argument("--viz-output")
    p_neigh.add_argument("--label-mode", choices=("auto", "gene", "locus_tag", "product", "none"), default="auto")
    p_neigh.add_argument("--include-feature-types", default="CDS")
    p_neigh.add_argument("--color-by", choices=("default", "ruleset"), default="default")
    p_neigh.add_argument("--ruleset", choices=("mobilome", "xenobiotics"), default="mobilome")
    p_neigh.add_argument("--show-operons", action="store_true")
    p_neigh.add_argument("--operon-gap", type=_nonnegative_int, default=150)

    p_reg = subparsers.add_parser("region", help="Extract a genomic sub-region")
    p_reg.add_argument("input")
    p_reg.add_argument("--locus")
    p_reg.add_argument("--record")
    p_reg.add_argument("--start", type=int)
    p_reg.add_argument("--end", type=int)
    p_reg.add_argument("--flank-genes", type=int, default=0)
    p_reg.add_argument("--flank-bp", type=int, default=0)
    p_reg.add_argument("--rebase", action="store_true")
    p_reg.add_argument("--output")
    p_reg.add_argument("--force", action="store_true")

    p_fa = subparsers.add_parser("fasta", help="Export annotated CDS translations")
    p_fa.add_argument("input")
    p_fa.add_argument("output", nargs="?")
    p_fa.add_argument("--force", action="store_true")

    p_seq = subparsers.add_parser("sequence", help="Extract genome and CDS nucleotide FASTA")
    p_seq.add_argument("input")
    p_seq.add_argument("--fna")
    p_seq.add_argument("--ffn")
    p_seq.add_argument("--force", action="store_true")

    p_cod = subparsers.add_parser("codon", help="Calculate codon usage and GC metrics")
    p_cod.add_argument("input")
    p_cod.add_argument("--min-len", type=int, default=100)
    p_cod.add_argument("--output")
    p_cod.add_argument("--include-pseudo", action="store_true")
    p_cod.add_argument("--force", action="store_true")

    p_func = subparsers.add_parser("functional", help="Profile COGs and pathway completeness")
    p_func.add_argument("input")
    _add_output_arguments(p_func, formats=("text", "tsv", "json"), default_format="tsv")
    p_func.add_argument("--pathways-only", action="store_true")
    p_func.add_argument("--cog-only", action="store_true")

    p_disc = subparsers.add_parser("discover", help="Scan annotation-supported islands")
    p_disc.add_argument("input")
    p_disc.add_argument("--cluster-gap", type=_nonnegative_int, default=5000)
    p_disc.add_argument("--operon-gap", type=_nonnegative_int, default=150)
    p_disc.add_argument("--min-weight", type=int, default=1)
    p_disc.add_argument("--format", choices=("text", "json", "tsv"), default="text")
    p_disc.add_argument("--ruleset", choices=("mobilome", "xenobiotics"), default="mobilome")
    p_disc.add_argument("--rules", dest="rules_file")
    p_disc.add_argument("--output")
    p_disc.add_argument("--force", action="store_true")

    p_comp = subparsers.add_parser("compare", help="Compare typed annotation targets across genomes")
    p_comp.add_argument("genomes", nargs="+")
    target_group = p_comp.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--targets")
    target_group.add_argument("--targets-file")
    p_comp.add_argument("--mode", choices=("count", "presence", "status"), default="count")
    p_comp.add_argument("--format", choices=("text", "tsv", "json"), default="text")
    p_comp.add_argument("--output")
    p_comp.add_argument("--evidence-output")
    p_comp.add_argument("--force", action="store_true")

    p_diff = subparsers.add_parser("diff", help="Compare two annotation versions")
    p_diff.add_argument("old_file")
    p_diff.add_argument("new_file")
    p_diff.add_argument("--format", choices=("text", "json", "jsonl", "tsv"), default="text")
    p_diff.add_argument("--summary-only", action="store_true")
    p_diff.add_argument("--max-display", type=_nonnegative_int, default=25)
    p_diff.add_argument("--coordinate-map")
    p_diff.add_argument("--output")
    p_diff.add_argument("--force", action="store_true")

    p_phy = subparsers.add_parser("phylo", help="Extract annotation-based marker candidates")
    p_phy.add_argument("input")
    p_phy.add_argument("--markers", choices=("core", "housekeeping", "all"), default="all")
    p_phy.add_argument("--min-length", type=_positive_int, default=50)
    p_phy.add_argument("--output")
    p_phy.add_argument("--output-dir")
    p_phy.add_argument("--format", choices=("text", "tsv", "json"), default="text")
    p_phy.add_argument("--force", action="store_true")

    p_cris = subparsers.add_parser("crispr", help="Scan CRISPR/Cas annotations")
    p_cris.add_argument("input")
    p_cris.add_argument("--window", type=_nonnegative_int, default=15000)
    _add_output_arguments(p_cris, formats=("text", "tsv", "json"), default_format="text")

    p_gff = subparsers.add_parser("gff", help="Convert GenBank to GFF3")
    p_gff.add_argument("input")
    p_gff.add_argument("output", nargs="?")
    p_gff.add_argument("--include-fasta", action="store_true")
    p_gff.add_argument("--force", action="store_true")

    p_bat = subparsers.add_parser("batch-summary", help="Parse Bakta summaries")
    p_bat.add_argument("inputs", nargs="+")
    p_bat.add_argument("--csv", default="bakta_summary.csv")
    p_bat.add_argument("--tsv", default="bakta_summary.tsv")
    p_bat.add_argument("--md", default="bakta_summary.md")

    p_meor = subparsers.add_parser("meor", help="Scan MEOR and biosurfactant genomic potential")
    p_meor.add_argument("input")
    p_meor.add_argument("--format", choices=("text", "json", "tsv"), default="text")
    p_meor.add_argument("--min-weight", type=int, choices=(1, 2, 3), default=1)
    p_meor.add_argument("--max-gap", type=_nonnegative_int, default=200)
    p_meor.add_argument("--window-size", type=_positive_int, default=50000)
    p_meor.add_argument("--output")
    p_meor.add_argument("--force", action="store_true")
    p_meor.add_argument("--markers")
    p_meor.add_argument("--pathways")

    p_mob = subparsers.add_parser("mobilome", help="Build a mobilome evidence report")
    p_mob.add_argument("input")
    p_mob.add_argument("--format", choices=("text", "json", "tsv"), default="text")
    p_mob.add_argument("--include", choices=("all", "chromosome", "plasmid", "unknown"), default="all")
    p_mob.add_argument("--min-evidence", type=int, choices=(1, 2, 3), default=1)
    p_mob.add_argument("--database-dir")
    p_mob.add_argument("--output")
    p_mob.add_argument("--force", action="store_true")

    p_query = subparsers.add_parser("query", help="Safely query annotated features")
    p_query.add_argument("input")
    p_query.add_argument("--where", required=True)
    p_query.add_argument("--select")
    p_query.add_argument("--format", choices=("text", "tsv", "csv", "json", "jsonl"), default="text")
    p_query.add_argument("--emit", choices=("rows", "faa", "ffn"), default="rows")
    p_query.add_argument("--strict", action="store_true")
    p_query.add_argument("--output")
    p_query.add_argument("--force", action="store_true")

    p_export = subparsers.add_parser("export", help="Export interoperable GenBank representations")
    p_export.add_argument("input")
    p_export.add_argument("--format", choices=list(EXPORT_FORMATS), required=True)
    p_export.add_argument("--output")
    p_export.add_argument("--output-dir")
    p_export.add_argument("--include-fasta", action="store_true")
    p_export.add_argument("--strict", action="store_true")
    p_export.add_argument("--force", action="store_true")

    p_operon = subparsers.add_parser("operons", help="Report circular-aware operon proximity candidates")
    p_operon.add_argument("input")
    p_operon.add_argument("--max-gap", type=int, default=150)
    p_operon.add_argument("--min-gap", type=int, default=-50)
    p_operon.add_argument("--min-genes", type=_positive_int, default=3)
    p_operon.add_argument("--format", choices=("text", "tsv", "json", "gff3"), default="text")
    p_operon.add_argument("--output")
    p_operon.add_argument("--force", action="store_true")

    return parser


def _write_or_stdout(
    rendered: str,
    output: str | Path | None,
    *,
    force: bool,
    inputs: Sequence[str | Path | None],
) -> None:
    write_text(rendered, output, force=force, inputs=inputs)


def _query_fasta(matches: Sequence[object], emit: str, *, strict: bool) -> tuple[str, tuple[str, ...]]:
    lines: list[str] = []
    warnings: list[str] = []
    for match in matches:
        feature = match.feature  # type: ignore[attr-defined]
        record = match.record  # type: ignore[attr-defined]
        identifier = f"{feature.record_id}:{feature.locus_tag or f'feature_{feature.feature_index}'}"
        if emit == "faa":
            sequence = "".join(feature.translation.split())
            if not sequence:
                warnings.append(f"skipped {identifier}: missing /translation")
                continue
        else:
            try:
                sequence = str(feature.extract(record.seq)).upper() if len(record.seq) else ""
            except (AttributeError, TypeError, ValueError):
                sequence = ""
            if not sequence:
                warnings.append(f"skipped {identifier}: no nucleotide sequence")
                continue
        lines.append(f">{identifier} feature_index={feature.feature_index}")
        width = 60 if emit == "faa" else 70
        lines.extend(sequence[index : index + width] for index in range(0, len(sequence), width))
    if strict and warnings:
        raise ValueError("; ".join(warnings))
    return "\n".join(lines) + ("\n" if lines else ""), tuple(warnings)


def _dispatch(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    cmd = args.subcommand

    if cmd == "validate":
        report = build_validation_report(args.input)
        rendered = render_validation_json(report) if (args.json or args.format == "json") else render_validation_text(report)
        _write_or_stdout(rendered, args.output, force=args.force, inputs=(args.input,))
        if args.fail_on != "never":
            levels = {"error": {"ERROR"}, "warning": {"ERROR", "WARNING"}, "info": {"ERROR", "WARNING", "INFO"}}
            if any(finding.severity in levels[args.fail_on] for finding in report.findings):
                return EXIT_VALIDATION
        return EXIT_OK

    if cmd in {"summary", "metadata"}:
        report = build_metadata_report(args.input)
        _write_or_stdout(render_metadata(report, args.format), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "extract":
        output = args.output
        document = read_genbank(args.input)
        rendered = render_annotations_tsv(feature_rows(document))
        _write_or_stdout(rendered, output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "search":
        reject_input_output_collision((args.input,), args.output)
        if args.output and Path(args.output).exists() and not args.force:
            raise OutputError(f"output already exists; use --force: {args.output}")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
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
            )
        _write_or_stdout(buffer.getvalue(), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "locus":
        report = build_locus_report(args.input, args.locus_tag)
        _write_or_stdout(render_locus_report(report, args.format), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "neighborhood":
        if args.output and args.viz_output and Path(args.output).resolve() == Path(args.viz_output).resolve():
            parser.error("Neighborhood data and figure outputs must use different paths")
        feature_types = tuple(item.strip() for item in args.include_feature_types.split(",") if item.strip())
        if not feature_types:
            parser.error("--include-feature-types must contain at least one type")
        result = build_neighborhood(
            args.input,
            args.locus_tag,
            window=args.window,
            include_feature_types=feature_types,
            ruleset=args.ruleset if args.color_by == "ruleset" else None,
            show_operons=args.show_operons,
            operon_gap=args.operon_gap,
        )
        _write_or_stdout(serialize_neighborhood(result, args.format), args.output, force=args.force, inputs=(args.input,))
        if args.visualize or args.viz_output:
            from .visualization.neighborhood import (
                VisualizationDependencyError,
                default_visualization_path,
                render_neighborhood,
            )

            try:
                viz_path = Path(args.viz_output) if args.viz_output else default_visualization_path(result)
                reject_input_output_collision((args.input,), viz_path)
                if viz_path.exists() and not args.force:
                    raise OutputError(f"visualization output already exists; use --force: {viz_path}")
                render_neighborhood(
                    result,
                    args.viz_output,
                    label_mode=args.label_mode,
                    color_mode=args.color_by,
                    force=args.force,
                )
            except (VisualizationDependencyError, ValueError) as exc:
                parser.error(str(exc))
        return EXIT_OK

    if cmd == "region":
        if args.output:
            reject_input_output_collision((args.input,), args.output)
            if Path(args.output).exists() and not args.force:
                raise OutputError(f"output already exists; use --force: {args.output}")
        report_buffer = io.StringIO()
        with contextlib.redirect_stdout(report_buffer):
            region_record = extract_region(
                args.input,
                locus_tag=args.locus,
                record_id=args.record,
                start=args.start,
                end=args.end,
                flank_genes=args.flank_genes,
                flank_bp=args.flank_bp,
                rebase=args.rebase,
            )
        if args.output:
            _write_or_stdout(
                render_region_record(region_record, args.output),
                args.output,
                force=args.force,
                inputs=(args.input,),
            )
            eprint(report_buffer.getvalue().rstrip())
        else:
            sys.stdout.write(report_buffer.getvalue())
        return EXIT_OK

    if cmd == "fasta":
        artifact = render_export(args.input, "faa")
        assert isinstance(artifact, ExportArtifact)
        for warning in artifact.warnings:
            eprint(f"WARNING: {warning}")
        _write_or_stdout(artifact.content, args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "sequence":
        paths = [path for path in (args.fna, args.ffn) if path]
        for path in paths:
            reject_input_output_collision((args.input,), path)
            if Path(path).exists() and not args.force:
                raise OutputError(f"output already exists; use --force: {path}")
        report_buffer = io.StringIO()
        with contextlib.redirect_stdout(report_buffer):
            extract_sequences(
                args.input,
                out_fna=args.fna,
                out_ffn=args.ffn,
                force=args.force,
            )
        eprint(report_buffer.getvalue().rstrip())
        return EXIT_OK

    if cmd == "codon":
        if args.output:
            reject_input_output_collision((args.input,), args.output)
            if Path(args.output).exists() and not args.force:
                raise OutputError(f"output already exists; use --force: {args.output}")
        report_buffer = io.StringIO()
        with contextlib.redirect_stdout(report_buffer):
            report = analyze_codon_usage(
                args.input,
                min_len_aa=args.min_len,
                output_path=None,
                include_pseudo=args.include_pseudo,
            )
        if args.output:
            _write_or_stdout(
                render_codon_tsv(report),
                args.output,
                force=args.force,
                inputs=(args.input,),
            )
            eprint(report_buffer.getvalue().rstrip())
        else:
            sys.stdout.write(report_buffer.getvalue())
        return EXIT_OK

    if cmd == "functional":
        report = build_functional_report(args.input, pathways_only=args.pathways_only, cog_only=args.cog_only)
        _write_or_stdout(render_functional(report, args.format), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "discover":
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            result = discover_clusters(args.input, cluster_gap=args.cluster_gap, operon_gap=args.operon_gap, min_weight=args.min_weight, format_type=args.format, rules_file=args.rules_file, ruleset=args.ruleset)
        _write_or_stdout(buffer.getvalue(), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "compare":
        targets = [item.strip() for item in args.targets.split(",") if item.strip()] if args.targets else None
        if args.output and args.evidence_output and paths_same(args.output, args.evidence_output):
            raise OutputError("matrix and evidence outputs must use different paths")
        matrix, evidence = compare_data(args.genomes, targets=targets, targets_file=args.targets_file, mode=args.mode)
        _write_or_stdout(render_compare(matrix, evidence, format_type=args.format, mode=args.mode), args.output, force=args.force, inputs=args.genomes)
        if args.evidence_output:
            _write_or_stdout(render_evidence_tsv(evidence), args.evidence_output, force=args.force, inputs=args.genomes)
        return EXIT_OK

    if cmd == "diff":
        if args.old_file == "-" and args.new_file == "-":
            raise InputError("diff accepts at most one stdin input")
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            diff_annotations(args.old_file, args.new_file, format_type=args.format, summary_only=args.summary_only, max_display=args.max_display, coordinate_map=args.coordinate_map)
        _write_or_stdout(buffer.getvalue(), args.output, force=args.force, inputs=(args.old_file, args.new_file))
        return EXIT_OK

    if cmd == "phylo":
        report = build_phylogenomic_report(args.input, marker_set=args.markers, min_length=args.min_length)
        if args.output:
            _write_or_stdout(render_phylogenomic_report(report, args.format), args.output, force=args.force, inputs=(args.input,))
        elif not args.output_dir:
            sys.stdout.write(render_phylogenomic_report(report, args.format))
        if args.output_dir:
            files: dict[str, str] = {}
            for marker, marker_data in report["markers"].items():
                if marker_data["hits"]:
                    files[f"{marker}.faa"] = "".join(
                        f">{marker} {hit['locus_tag'] or ''} {hit['record']}\n{hit['translation']}\n"
                        for hit in marker_data["hits"]
                    )
            publish_directory(files, args.output_dir, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "crispr":
        report = build_crispr_report(args.input, window=args.window)
        _write_or_stdout(render_crispr_report(report, args.format), args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "gff":
        if args.output:
            _write_or_stdout(convert_to_gff3(args.input, include_fasta=args.include_fasta), args.output, force=args.force, inputs=(args.input,))
        else:
            sys.stdout.write(convert_to_gff3(args.input, include_fasta=args.include_fasta))
        return EXIT_OK

    if cmd == "batch-summary":
        batch_summary(args.inputs, csv_out=args.csv, tsv_out=args.tsv, md_out=args.md)
        return EXIT_OK

    if cmd == "meor":
        database = load_meor_database(args.markers, args.pathways)
        report = analyze_meor(args.input, min_weight=args.min_weight, max_gap=args.max_gap, window_size=args.window_size, marker_database=database)
        try:
            reject_input_output_collision((args.input, args.markers, args.pathways), args.output)
        except OutputError as exc:
            parser.error(str(exc))
        _write_or_stdout(serialize_report(report, args.format), args.output, force=args.force, inputs=(args.input, args.markers, args.pathways))
        return EXIT_OK

    if cmd == "mobilome":
        if args.input == "-":
            parser.error("mobilome requires a path-backed input for provenance hashing")
        try:
            database = load_mobilome_database(args.database_dir)
        except MobilomeError as exc:
            parser.error(str(exc))
        try:
            report = analyze_mobilome(args.input, include=args.include, min_evidence=args.min_evidence, database=database)
            rendered = serialize_mobilome_report(report, args.format)
        except MobilomeError as exc:
            parser.error(str(exc))
        try:
            reject_input_output_collision((args.input, *database.source_paths), args.output)
        except OutputError as exc:
            parser.error(str(exc))
        _write_or_stdout(rendered, args.output, force=args.force, inputs=(args.input, *database.source_paths))
        return EXIT_OK

    if cmd == "query":
        matches = query_features(args.input, args.where)
        rows = tuple(match.row for match in matches)
        selected = parse_select_fields(args.select)
        if args.emit == "rows":
            rendered = render_rows(rows, format_type=args.format, selected=selected)
        else:
            rendered, warnings = _query_fasta(matches, args.emit, strict=args.strict)
            for warning in warnings:
                eprint(f"WARNING: {warning}")
        _write_or_stdout(rendered, args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "operons":
        document = read_genbank(args.input)
        report = build_operon_report(
            document,
            max_gap=args.max_gap,
            min_gap=args.min_gap,
            min_genes=args.min_genes,
        )
        rendered = render_operons_gff3(report) if args.format == "gff3" else render_operons(report, args.format)
        _write_or_stdout(rendered, args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    if cmd == "export":
        if args.output and args.output_dir:
            parser.error("--output and --output-dir are mutually exclusive")
        if args.format == "ncbi-table":
            if not args.output_dir or args.output:
                parser.error("ncbi-table requires --output-dir and does not support --output")
            artifact = render_export(args.input, args.format, strict=args.strict)
            assert isinstance(artifact, dict)
            publish_directory(artifact, args.output_dir, force=args.force, inputs=(args.input,))
            return EXIT_OK
        artifact = render_export(args.input, args.format, strict=args.strict, include_fasta=args.include_fasta)
        assert isinstance(artifact, ExportArtifact)
        for warning in artifact.warnings:
            eprint(f"WARNING: {warning}")
        if args.output_dir:
            suffix = {"annotations-tsv": ".tsv", "jsonl": ".jsonl", "faa": ".faa", "ffn": ".ffn", "fna": ".fna", "gff3": ".gff3", "bed12": ".bed"}[args.format]
            stem = Path(args.input).name if args.input != "-" else "genbank"
            if stem.casefold().endswith(".gz"):
                stem = Path(stem).stem
            stem = Path(stem).stem
            publish_directory(
                {f"{stem}{suffix}": artifact.content},
                args.output_dir,
                force=args.force,
                inputs=(args.input,),
            )
        else:
            _write_or_stdout(artifact.content, args.output, force=args.force, inputs=(args.input,))
        return EXIT_OK

    raise AssertionError(f"unhandled subcommand: {cmd}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        return _dispatch(args, parser)
    except (GenBankInputError, InputError, OSError, FileNotFoundError, QueryExpressionError, ValueError) as exc:
        eprint(f"ERROR: {exc}")
        return EXIT_INPUT
    except OutputError as exc:
        eprint(f"ERROR: {exc}")
        return EXIT_OUTPUT
    except MobilomeError as exc:
        eprint(f"ERROR: {exc}")
        return EXIT_INPUT


if __name__ == "__main__":
    sys.exit(main())
