# Implementation Plan: Bakta Multi-Isolate Genome Summary Parser

## Refined scope

Create `scripts/parse_bakta_summaries.py`, a command-line tool that produces one comparison table per set of Bakta-annotated GenBank files.

The original proposal assumed a fixed Bakta `.txt`/`.tsv` layout and reimplemented a GenBank parser. Those assumptions are not reliable across Bakta versions and conflict with this repository's canonical Biopython parser. The implementation will therefore require `.gbff`, `.gbk`, or `.gb` input for feature-level metrics. When a same-named Bakta `.txt` exists, its published headline counts, genome size, and GC are preferred; GBFF is authoritative for the per-CDS qualifier metrics (UniRef, RefSeq, COG, and KEGG/KO), and is the complete fallback when no summary exists.

## Outputs and invocation

```powershell
python scripts/parse_bakta_summaries.py -i ./bakta-results
python scripts/parse_bakta_summaries.py -i C14-NMZ.gbff PF_NNT_reoriented.gbff `
  --output-csv summary.csv --output-tsv summary.tsv --output-md summary.md
```

The tool recursively discovers GenBank files under each input directory, writes CSV, TSV, and Markdown outputs, and exits with a clear error if no compatible input is found. Outputs are deterministically ordered by sample name.

## Metric definitions

| Output column | Definition |
|---|---|
| Genome Size (Mbp) | Sibling Bakta `Length`/`Size` when available; otherwise sum of sequence lengths, then `source` feature spans. |
| Total Gene/Feature Count | Sum of the Bakta headline feature categories, including sORFs and pseudogenes, when available. Without a summary, count annotation records except `source` and redundant `gene` records, plus pseudogene-only genes. |
| G + C content (%) | Sibling Bakta `GC` when available; otherwise GC over unambiguous A/C/G/T bases from `ORIGIN`. |
| CDS, rRNA, tRNA / tmRNA | Published Bakta headline counts when available; otherwise counts of the corresponding GenBank feature types. |
| UniRef, RefSeq, COG, KEGG/KO | Number of distinct CDSs carrying a matching identifier in `db_xref` or `note`; COG and KO extraction reuses `extract_xrefs()`. |
| Genes without function prediction | CDS product exactly equal to `hypothetical protein`, case-insensitively. |
| Pseudogenes | Unique `pseudogene` features or `gene`/`CDS` records marked by a `pseudo` or `pseudogene` qualifier. |
| Regulatory ncRNAs | `regulatory` feature records, which are how Bakta encodes leaders, riboswitches, and related elements. |

## Implementation details

1. Reuse `parse_features()`, `get_qual()`, and `extract_xrefs()` from `genbank_parser.py`; reuse `parse_sequences()` from `genbank_sequence.py` for sequence-derived values.
2. Count per-CDS assignments, never raw identifier occurrences, so multiple qualifiers do not inflate a metric.
3. Parse the stable `Key: value` lines found in Bakta text reports (`Length`, `GC`, RNAs, CDSs, sORFs, pseudogenes, origins, and CRISPR arrays). This preserves Bakta's intended distinction between CDSs and sORFs; for example, the supplied C14 report lists 4,626 CDSs while its GBFF has 4,636 CDS records.
4. Do not claim that TSV-only input supports the full report: Bakta TSV schemas are version-dependent and lack the robust feature semantics available in GBFF.
5. Render Markdown safely by escaping table delimiters; write all text outputs in UTF-8 with explicit newlines.
6. Report per-file parsing failures to stderr, continue with other samples, and return a nonzero status only if no records can be created.

## Verification

1. Compile the script with `python -m py_compile`.
2. Run it against the two workspace GBFF fixtures.
3. Confirm all three outputs have two rows and consistent headers. For C14, confirm that its supplied text report is reflected (5.051 Mbp, 59.80% GC, 4,626 CDS); confirm GBFF-derived fallback values against the canonical validator where no sibling summary exists.
4. Exercise the `--help` interface and an invalid input path to verify clear CLI errors.
