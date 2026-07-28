---
name: genbank-feature-parser
description: Parse GenBank feature table files computationally. Validates structure, extracts annotations, and runs diagnostic bioinformatic analyses without loading the entire file into the AI context window. Use this skill whenever the user asks to parse, validate, summarize, compare, or analyze .gbff, .gbk, .gb, or .txt GenBank files; when they want gene neighborhoods, operon candidates, COG/GO/KEGG/EC annotation tables, codon usage, phylogenetic markers, CRISPR detection, or GFF3 export; or when they mention tools like Bakta, Prokka, tbl2asn, or NCBI GenBank in the context of annotation analysis.
---

# GenBank Feature Table Parser

## Purpose

Parse, validate, and analyze GenBank feature table files (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** — never loading the full file content into the AI context window. All logic lives in standalone Python scripts under `scripts/`.

## Epistemic Discipline

1. **Scripts parse; the AI reads summaries.** Never `view_file` on a GenBank file — write a parser and read its output.
2. **One canonical parser, reused everywhere.** All scripts import from `genbank_parser.py` (powered by Biopython `Bio.SeqIO`). Never reconstruct `parse_features()` ad-hoc.
3. **Compact structured output only.** Scripts emit diagnostic tables, not data dumps.

---

## Script Inventory

All scripts reside in the `scripts/` subdirectory. Invoke them in-place; Python resolves the `genbank_parser` import automatically.

| Script | Purpose | Typical Invocation |
|---|---|---|
| `genbank_parser.py` | **Canonical parser + `extract_xrefs()` helper** | (import only) |
| `genbank_validate.py` | Structural validation report | `python scripts/genbank_validate.py INPUT.gbff` |
| `genbank_metadata.py` | LOCUS header / organism summary | `python scripts/genbank_metadata.py INPUT.gbff` |
| `genbank_sequence.py` | ORIGIN → `.fna` + `.ffn` + GC% | `python scripts/genbank_sequence.py INPUT.gbff [--fna out.fna] [--ffn out.ffn]` |
| `genbank_extract.py` | Tab-delimited annotation TSV | `python scripts/genbank_extract.py INPUT.gbff [output.tsv]` |
| `genbank_functional.py` | Unified functional + pathway completeness | `python scripts/genbank_functional.py INPUT.gbff [--format json]` |
| `genbank_fasta.py` | Export all CDS as protein FASTA | `python scripts/genbank_fasta.py INPUT.gbff [proteins.faa]` |
| `genbank_neighborhood.py` | Gene neighborhood viewer (+/- N genes) | `python scripts/genbank_neighborhood.py INPUT.gbff LOCUS_TAG [window]` |
| `genbank_locus.py` | Single-locus deep-dive (all qualifiers) | `python scripts/genbank_locus.py INPUT.gbff LOCUS_TAG` |
| `genbank_operons.py` | Operon candidates (same-strand, gap) | `python scripts/genbank_operons.py INPUT.gbff [max_gap]` |
| `genbank_compare.py` | Multi-genome presence/absence matrix | `python scripts/genbank_compare.py ./genomes/ "ladA,ssuD" [matrix.tsv]` |
| `genbank_discover.py` | Discovery mode — keyword + spatial scan | `python scripts/genbank_discover.py INPUT.gbff [--cluster-gap 5000] [--operon-gap 150] [--min-weight 1] [--format text\|json\|tsv]` |
| `genbank_meor.py` | MEOR & biosurfactant discovery engine | `python scripts/genbank_meor.py INPUT.gbff [--min-weight 1] [--max-gap 200] [--format text|json|tsv]` |
| `genbank_crispr.py` | CRISPR array + Cas gene detection | `python scripts/genbank_crispr.py INPUT.gbff [--window 15000]` |
| `genbank_codon.py` | Codon usage + RSCU table | `python scripts/genbank_codon.py INPUT.gbff [--min-len 100] [--output codon.tsv]` |
| `genbank_phylo.py` | Phylogenetic marker gene extraction | `python scripts/genbank_phylo.py INPUT.gbff [--markers all] [--output-dir ./markers/]` |
| `genbank_gff.py` | Export to GFF3 format | `python scripts/genbank_gff.py INPUT.gbff [OUTPUT.gff3] [--include-fasta]` |

All scripts use `argparse` — run any script with `--help` for full usage.

---

## Supported Formats

Targets GenBank feature tables as produced by **Bakta**, **Prokka**, or NCBI's `tbl2asn`.

```
<feature_key>       <location>
                     /qualifier="value"
```

### Feature keys
The parser is **permissive** — it accepts any feature key at the canonical 5-space indent, including `pseudogene`, `C_region`, `V_segment`, `D_segment`, `misc_RNA`, and any custom key. There is no hardcoded whitelist.

### Common qualifiers
`/locus_tag`, `/gene`, `/product`, `/translation`, `/protein_id`, `/note`, `/inference`, `/db_xref`, `/EC_number`

---

## Parser Capabilities & Limits

### What the parser handles
- Simple locations: `1000..2000`, `complement(1000..2000)`
- Partial features: `<1000..>2000`
- **Single-line `join()` / `order()`**: `join(100..200,300..400)`, `complement(join(...))` — parsed to `start=min`, `end=max`, with a `join_segments` list preserved for downstream use (e.g. GFF3 export)
- Multi-contig files: tracks `LOCUS` lines, injects `feature['contig']` key
- ORIGIN blocks: detected and skipped during feature parsing
- Multi-line qualifiers: robust continuation handling (Biopython)
- Encoding: `utf-8` with `errors='replace'` for Windows compatibility

### What the parser does NOT handle
- **Nested `join()`** — e.g. `complement(join(join(...)))` (warns on stderr, skipped)
- **Multi-line `join()` openers** — location that wraps across lines (warns on stderr, skipped)
- Features spanning the replication origin

For prokaryotic Bakta/Prokka output, these never appear. For eukaryotic data, check stderr warning counts.

### Multi-contig / multi-record files

The parser handles `.gbff` files with multiple `LOCUS` records (draft assemblies, chromosome + plasmids). Each feature gets a `contig` key. **All downstream scripts must respect contig boundaries** — never group genes from different contigs.

### Cross-reference semantics (`extract_xrefs`)

The helper `extract_xrefs(feature)` in `genbank_parser.py` returns a typed dict:

| Key | Prefix matched | Source qualifiers |
|---|---|---|
| `go_terms` | `GO:` | `/db_xref`, `/note` |
| `cog_ids` | `COG` | `/db_xref`, `/note` |
| `kegg_kos` | `K\d{5}` | `/db_xref`, `/note` |
| `pfam` | `PF\d` | `/db_xref`, `/note` |
| `rfam` | `RFAM` | `/db_xref`, `/note` |
| `ec_numbers` | `/EC_number` → `EC:` in `/note` | Priority: `/EC_number` first (INSDC); fallback to `EC:` in `/note` (Bakta style) |
| `db_xrefs` | everything else | `/db_xref` only |

All downstream scripts use `extract_xrefs()` — the old behavior of dumping the entire `/db_xref` list into a `GO_terms` column has been corrected.

### Feature dict structure

```python
{
    'type':          'CDS',
    'start':         1000,
    'end':           2000,
    'strand':        '+',
    'contig':        'contig_1',         # from parent LOCUS
    'qualifiers':    defaultdict(list),  # key -> [values]
    'line':          42,                 # source line number (unique, used as identity)
    'join_segments': [(1000, 1500), (1700, 2000)],  # only if join/order
}
```

---

## Workflow

### Step 0: Locate scripts

```bash
SCRIPTPATH="/path/to/skill/scripts"  # adjust to actual install location
```

### Step 1a: Structural validation
```bash
python $SCRIPTPATH/genbank_validate.py INPUT.gbff
```
Reports: feature type counts, strand distribution, CDS statistics, locus tag integrity.

### Step 1b: Metadata (full `.gbff` only)
```bash
python $SCRIPTPATH/genbank_metadata.py INPUT.gbff
```
Reports: contig count, sizes, topology, organism, strain.

### Step 1c: Sequence extraction (full `.gbff` only)
```bash
python $SCRIPTPATH/genbank_sequence.py INPUT.gbff [--fna genome.fna] [--ffn cds.ffn]
```
Produces: genome FASTA (`.fna`), per-CDS nucleotide slices (`.ffn`), GC% and coding density stats.

For feature-table-only snippets (no ORIGIN), prints a graceful message and exits.

### Step 2: Annotation summary
```bash
python $SCRIPTPATH/genbank_extract.py INPUT.gbff annotations.tsv
```
TSV columns: `locus_tag, gene, contig, start, end, strand, length_bp, length_aa, product, EC_number, COG, KEGG, GO_terms, Pfam, RFAM, db_xrefs, inference`

### Step 3: Diagnostic analyses

| Goal | Script |
|---|---|
| COG and pathway completeness | `genbank_functional.py` |
| Export protein FASTA for BLAST/HMMscan | `genbank_fasta.py` |
| Gene neighborhood for a target locus | `genbank_neighborhood.py` |
| All qualifiers for one locus tag | `genbank_locus.py` |
| Operon candidates (co-directional, gapped) | `genbank_operons.py` |
| CRISPR arrays + Cas gene co-localization | `genbank_crispr.py` |
| Codon usage + RSCU for expression planning | `genbank_codon.py` |
| Phylogenetic marker genes (ribosomal/housekeeping) | `genbank_phylo.py` |
| Export to GFF3 for genome browsers | `genbank_gff.py` |

### Step 4: External bioinformatic hooks

| Analysis | Tool | Command |
|---|---|---|
| BLAST search | `blastp` | `blastp -query proteins.faa -db nr -outfmt 6 -max_target_seqs 5 -evalue 1e-5` |
| Domain annotation | `hmmscan` | `hmmscan --tblout domains.tsv Pfam-A.hmm proteins.faa` |
| Signal peptides | `signalp` | `signalp -fasta proteins.faa -org gram- -format short` |
| Transmembrane | `tmhmm` | `tmhmm --fasta proteins.faa` |
| Antibiotic resistance | `abricate` | `abricate --db resfinder proteins.faa` |

### Step 5: Comparative mode
```bash
python $SCRIPTPATH/genbank_compare.py ./genomes/ "ladA,ssuD,EC:1.14.14.28,KEGG:K20938" matrix.tsv
```
Gene name matching uses **word-boundary regex** — `alkB` will not match `chalkBoard`.

### Step 6: Discovery mode

For uncharacterized isolates, run `genbank_discover.py` as a "metal detector":
1. **Weighted keyword scan** — curated dictionaries with confidence levels (3/2/1)
2. **Spatial clustering** — groups hits into candidate genomic islands (default 5 kb gap)
3. **Dark matter operon detection** — pure (100% hypothetical) and mixed (>=60%) clusters
4. **Mobilome-adjacent island detection** — functional genes within 10 kb of transposase/integrase

```bash
python $SCRIPTPATH/genbank_discover.py INPUT.gbff
python $SCRIPTPATH/genbank_discover.py INPUT.gbff --cluster-gap 10000 --operon-gap 200
python $SCRIPTPATH/genbank_discover.py INPUT.gbff --min-weight 2       # suppress low-confidence hits
python $SCRIPTPATH/genbank_discover.py INPUT.gbff --format json        # structured output
python $SCRIPTPATH/genbank_discover.py INPUT.gbff --format tsv         # tabular output
```
Runs: keyword scan (word-boundary regex) → spatial clustering (per-contig) → dark-matter operon detection → mobilome-adjacent island detection.

Output includes: contig summary, genomic density karyogram, category hits with strand/confidence/multi-category flags, dark matter operons, and mobilome-adjacent functional islands.

---

## Key Principles

1. **Never load the entire file into the AI context window.** Always use scripts to parse and summarize.
2. **Scripts first, raw content never.** Parse computationally, then read the structured output.
3. **Validate first.** Run Step 1a before any downstream analysis.
4. **One canonical parser.** All scripts import from `genbank_parser.py`. Never reconstruct ad-hoc.
5. **Respect contig boundaries.** Every spatial operation (neighborhood, operons, clustering, mobilome proximity) filters by `feature['contig']` before comparing positions. Cross-contig contamination is a hard bug — all scripts enforce this.
6. **Value semantics for cross-references.** Use `extract_xrefs()` — never dump raw `/db_xref` lists into GO or KEGG columns.
7. **Keep reports concise.** Diagnostic tables, not data dumps.
8. **Know the parser's limits.** Nested/multi-line `join()` warns and skips. Verify counts for eukaryotic data.

---

## Quick-Reference Command Sequence

```bash
SCRIPTPATH="path/to/scripts"

# 1. Validate + metadata
python $SCRIPTPATH/genbank_validate.py  INPUT.gbff
python $SCRIPTPATH/genbank_metadata.py  INPUT.gbff

# 2. Sequence extraction (full .gbff only)
python $SCRIPTPATH/genbank_sequence.py  INPUT.gbff

# 3. Annotation table (semantically typed columns)
python $SCRIPTPATH/genbank_extract.py   INPUT.gbff  annotations.tsv

# 4. Diagnostics
python $SCRIPTPATH/genbank_functional.py   INPUT.gbff
python $SCRIPTPATH/genbank_fasta.py        INPUT.gbff  proteins.faa
python $SCRIPTPATH/genbank_neighborhood.py INPUT.gbff  LOCUS_TAG 5
python $SCRIPTPATH/genbank_locus.py        INPUT.gbff  LOCUS_TAG
python $SCRIPTPATH/genbank_operons.py      INPUT.gbff  150

# 5. New analysis scripts
python $SCRIPTPATH/genbank_crispr.py       INPUT.gbff  --window 15000
python $SCRIPTPATH/genbank_codon.py        INPUT.gbff  --output codon.tsv
python $SCRIPTPATH/genbank_phylo.py        INPUT.gbff  --markers all --output-dir ./markers/
python $SCRIPTPATH/genbank_gff.py          INPUT.gbff  OUTPUT.gff3

# 6. Comparative scan
python $SCRIPTPATH/genbank_compare.py      ./genomes/  "ladA,ssuD,alkB"  matrix.tsv

# 7. Discovery mode
python $SCRIPTPATH/genbank_discover.py     INPUT.gbff
python $SCRIPTPATH/genbank_discover.py     INPUT.gbff --format json > report.json
python $SCRIPTPATH/genbank_discover.py     INPUT.gbff --format tsv  > hits.tsv
python $SCRIPTPATH/genbank_discover.py     INPUT.gbff --min-weight 3  # high-confidence only
```

---

## Pressure-Testing Protocol

### 1. Multi-contig contig-boundary test
Use a `.gbff` with multiple `LOCUS` records. Verify `genbank_neighborhood.py` shows only genes from the same contig as the target. Verify `genbank_operons.py` never pairs genes across contigs.

### 2. `join()` location test
Create a CDS with `join(100..200,300..400)`. Verify: no stderr warning, `join_segments` present, `start=100`, `end=400`. Verify GFF3 export emits mRNA + per-segment CDS records.

### 3. Permissive feature key test
Add a `pseudogene` and `V_segment` feature. Verify they appear in `parse_features()` output.

### 4. Cross-reference semantic test
CDS with `/EC_number="1.2.3.4"`, `/db_xref="GO:0001234"`, `/db_xref="UniProtKB:P12345"`, `/note="COG:COG0001"`. Verify `extract_xrefs()` routes each to the correct bucket and `db_xrefs` contains only `UniProtKB:P12345`.

### 5. EC number priority test
CDS with only `/note="EC:1.14.14.28"` (Bakta style, no `/EC_number`). Verify `ec_numbers=['1.14.14.28']`. Then add `/EC_number="1.14.14.28"` and verify it still works (INSDC style takes priority).

### 6. Word-boundary matching test
Gene with product `"chalkBoard reductase"`. Search for `alkB` in `genbank_compare.py`. Verify zero hits. Search for a product `"alkB monooxygenase"`. Verify one hit.

### 7. Full `.gbff` sequence extraction
Run `genbank_sequence.py` on a complete file. Verify `.fna` and `.ffn` produced with correct GC% and coding density.

### 8. Large file performance
File with 5,000+ CDS. Parser completes in < 5 seconds; `genbank_extract.py` TSV written without memory errors.

### 9. Edge cases
- Empty file — validator should error clearly
- Header only, no features — parser returns `[]`
- Feature with no qualifiers — parser still includes it
- Very large file (10k+ CDS) — parser completes in seconds

---

## Changelog (v2)

| Issue | Fix |
|---|---|
| #1 GO_terms was entire `/db_xref` list | `extract_xrefs()` now routes by value prefix; new columns: Pfam, RFAM, db_xrefs |
| #2 Cross-contig neighborhood contamination | `genbank_neighborhood.py` now filters by target contig before sorting |
| #3 Cross-contig dark matter in discover | `spatial_cluster()` groups by contig; `_dark_matter_clusters()` processes each contig independently |
| #4 No join() support | `_join_re` parses single-line join/order; segments preserved in `join_segments` |
| #5 No GFF3 support | New `genbank_gff.py` with join segment expansion and `##FASTA` option |
| #6 Hardcoded feature whitelist | `_feature_re` now uses `\w+` — accepts any feature key |
| #7 No CRISPR detection | New `genbank_crispr.py` (Types I–VI Cas + array co-localization) |
| #8 No codon usage | New `genbank_codon.py` (codon frequencies + RSCU per synonymous family) |
| #9 No phylogenetic markers | New `genbank_phylo.py` (30+ ribosomal + housekeeping markers, concatenated export) |
| #10 `__pycache__` in skill zip | Removed from distribution; exclude on re-package |
| #11 Product matching too loose | `genbank_compare.py` uses `\b` word-boundary regex |
| #12 No argparse | All 15 scripts now use `argparse` with `--help` and type validation |
| #13 `is` identity in operon clusters | Uses `feature['line']` (stable unique integer) for identity comparison |

## Changelog (v3)

| Issue | Fix |
|---|---|
| #14 `mera` keyword matches every `-merase` suffix | Removed `("mera", 3)`; word-boundary regex `\b` now used for all keywords. Added `merT/merC/merD/merP` for full operon coverage |
| #15 `benzoate` triggers menaquinone/siderophore false positives | Replaced with `benzoate 1,2-dioxygenase`, `benzoate catabolism`, `benzoyl-coa`; added `menaquinone`/`o-succinylbenzoate` to exclude-list |
| #16 `domain-containing protein` classified as hypothetical | `is_hypothetical()` refactored: only flags empty, bare "protein", hypothetical, uncharacterized, DUF\d+ |
| #17 Operon gap excludes overlapping gene pairs | Removed `0 <=` floor; overlapping IS200/IS605 TnpA/TnpB pairs now detected |
| #18 Mobilome distance uses start-to-start | New `interval_dist()` computes minimum endpoint-to-endpoint distance |
| #19 Mixed dark matter threshold asymmetric (4 vs 3) | Harmonised to `>= 3` genes for both pure and mixed dark matter |
| #20 No error handling on file parsing | File-existence check + try/except around `parse_features()` |
| #21 Product truncation drops diagnostic info | `_trunc()` helper adds `..` ellipsis; replaces all bare `[:45]` slicing |
| #22 No strand info in output | Strand `>`/`<` symbols on all gene output lines |
| #23 No multi-category flagging | Cross-category `*` tag highlights genes hitting multiple categories |
| #24 No structured output | `--format json` and `--format tsv` modes for programmatic downstream use |
| #25 No `--min-weight` filter | `--min-weight {1,2,3}` suppresses low-confidence hits |
| #26 No contig summary | Per-contig header with span, CDS count, hit density |
| #27 No spatial density overview | ASCII karyogram (50 kb windows) visualises hit clustering across genome |
| #28 Mobilome exclude-list for housekeeping | `rRNA/tRNA/DNA methyltransferase` excluded from Mobilome category |

## Changelog (v4)

| Issue | Fix |
|---|---|
| #29 KEGG `extract_xrefs` bug | Fixed `_kegg_ko_re` regex to correctly parse `KEGG:Kxxxxx` format in `/note` (restored visibility of ~900 KOs) |
| #30 No pathway completeness | Merged `genbank_cog.py` into new `genbank_functional.py` with 10 metabolic pathways evaluated via dual KO/EC semantics |

## Changelog (v5 — Biopython migration)

| Issue | Fix |
|---|---|
| #31 Regex parser fails on multi-line qualifier continuation | Replaced entire `parse_features()` with Biopython `Bio.SeqIO.parse()` |
| #32 Partial feature markers (`<`, `>`) silently discarded | Biopython preserves `BeforePosition`/`AfterPosition`; new `partial_start`/`partial_end` flags |
| #33 `genbank_validate.py` false positives on partial CDS | Skip length validation when `partial_start` or `partial_end` is `True` |
| #34 `genbank_compare.py` COG/KEGG lookup from `/note` only | Now uses `extract_xrefs()` to search both `/db_xref` and `/note` |
| #35 `genbank_functional.py` COG categories from `/note` only | `cog_distribution()` now checks both `/db_xref` and `/note` |
| #36 `genbank_locus.py` CDS-only filter | Now searches all feature types carrying `locus_tag` |
| #37 `reverse_complement` only maps ACGT | Uses `Bio.Seq.Seq.reverse_complement()` for full IUPAC support |
| #38 GC% includes ambiguous bases in denominator | Excludes non-ACGT bases from `total_acgt` denominator |
| #39 Output path splitting fragile | Replaced `rsplit('.', 1)` with `os.path.splitext()` across extract, fasta, sequence |
