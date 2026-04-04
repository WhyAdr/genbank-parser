---
name: genbank-feature-parser
description: Parse GenBank feature table files computationally. Validates structure, extracts annotations, and runs diagnostic bioinformatic analyses without loading the entire file into the context window.
---

# GenBank Feature Table Parser

## Purpose

Parse, validate, and analyze GenBank feature table files (`.gb`, `.gbk`, `.gbff`, `.txt`) **computationally** -- avoiding the need to load the full file content into the AI context window. All logic lives in standalone Python scripts under `scripts/`.

## Epistemic Discipline

1. **Scripts parse; the AI reads summaries.** Never `view_file` on a GenBank file -- write a parser and read its output instead.
2. **One canonical parser, reused everywhere.** All scripts import from `genbank_parser.py`. Never reconstruct `parse_features()` ad-hoc.
3. **Compact structured output only.** Scripts emit diagnostic tables, not data dumps.

---

## Script Inventory

All scripts are located in the `scripts/` subdirectory of this skill. You can run them directly from their location (Python will automatically resolve the `genbank_parser.py` import), or you can copy them to your working directory.

| Script | Purpose | Usage (Run in-place) |
|---|---|---|
| `genbank_parser.py` | **Canonical parser module** | (import only, not run directly) |
| `genbank_validate.py` | Structural validation report | `python scripts/genbank_validate.py INPUT.gbff` |
| `genbank_metadata.py` | LOCUS header / organism summary | `python scripts/genbank_metadata.py INPUT.gbff` |
| `genbank_sequence.py` | ORIGIN block -> `.fna` + `.ffn` + GC% | `python scripts/genbank_sequence.py INPUT.gbff` |
| `genbank_extract.py` | Tab-delimited annotation summary (TSV) | `python scripts/genbank_extract.py INPUT.gbff output.tsv` |
| `genbank_cog.py` | COG functional category distribution | `python scripts/genbank_cog.py INPUT.gbff` |
| `genbank_fasta.py` | Export all CDS as protein FASTA | `python scripts/genbank_fasta.py INPUT.gbff proteins.faa` |
| `genbank_neighborhood.py` | Gene neighborhood viewer (+/- N genes) | `python scripts/genbank_neighborhood.py INPUT.gbff LOCUS_TAG 5` |
| `genbank_locus.py` | Single-locus deep dive (all qualifiers)| `python scripts/genbank_locus.py INPUT.gbff LOCUS_TAG` |
| `genbank_operons.py` | Operon candidates (same-strand, gap) | `python scripts/genbank_operons.py INPUT.gbff 150` |
| `genbank_compare.py` | Multi-genome matrix | `python scripts/genbank_compare.py ./genomes/ "ladA,ssuD" matrix.tsv` |
| `genbank_discover.py` | Discovery mode -- metal detector | `python scripts/genbank_discover.py INPUT.gbff` |

---

## Supported Format

This skill targets GenBank feature table blocks as produced by **Bakta**, **Prokka**, or NCBI's `tbl2asn`:

```
<feature_key>       <location>
                     /qualifier="value"
                     /qualifier="value"
```

### Common feature keys
`gene`, `CDS`, `tRNA`, `rRNA`, `tmRNA`, `ncRNA`, `regulatory`, `misc_feature`, `repeat_region`

### Common qualifiers
`/locus_tag`, `/gene`, `/product`, `/translation`, `/protein_id`, `/note`, `/inference`, `/db_xref`, `/EC_number`

---

## Parser Capabilities & Limits

### What the parser handles
- Simple locations: `1000..2000`, `complement(1000..2000)`
- Partial features: `<1000..>2000`
- Multi-contig files: tracks `LOCUS` lines, injects `feature['contig']` key
- ORIGIN blocks: detected and skipped during feature parsing
- Multi-line qualifiers: robust 19-22 space continuation (handles tool variations)
- Encoding: `utf-8` with `errors='replace'` for Windows compatibility

### What the parser does NOT handle
- `join()`, `order()`, `complement(join(...))` -- **warned on stderr**, skipped
- Features spanning the replication origin

For prokaryotic Bakta/Prokka output, complex locations almost never appear. For eukaryotic/phage data, check the stderr warning output and compare feature counts.

### Multi-contig / multi-record files

The parser handles `.gbff` files with multiple `LOCUS` records (draft assemblies, chromosome + plasmids). Each feature gets a `contig` key. **All downstream scripts must respect contig boundaries** -- never group genes from different contigs.

### Feature dict structure

```python
{
    'type': 'CDS',
    'start': 1000,
    'end': 2000,
    'strand': '+',
    'contig': 'contig_1',       # from parent LOCUS
    'qualifiers': defaultdict(list),  # key -> [values]
    'line': 42,                 # source line number
}
```

---

## Workflow

### Step 0: Identify Script Location

All scripts reside in `.agent/skills/genbank-parser/scripts/`. You can invoke them in-place directly from this path. There is no need to copy `genbank_parser.py` because Python automatically resolves it when you run the other scripts from the same directory.

### Step 1a: Structural validation

Run `genbank_validate.py` to get a structural report: feature type counts, strand distribution, CDS statistics, locus tag integrity, and validation checks.

### Step 1b: Metadata extraction (full `.gbff` only)

Run `genbank_metadata.py` to scan LOCUS headers without parsing features. Reports contig count, sizes, topology, organism, and strain.

### Step 1c: Sequence extraction (full `.gbff` only)

Run `genbank_sequence.py` to extract:
- Genome FASTA (`.fna`) from `ORIGIN` blocks
- Per-CDS nucleotide slices (`.ffn`) with reverse-complement for minus-strand genes
- GC content and coding density statistics

For feature-table-only snippets (no ORIGIN), prints a graceful message and exits.

### Step 2: Annotation summary

Run `genbank_extract.py` to produce a TSV with columns: locus_tag, gene, start, end, strand, length_bp, length_aa, product, EC_number, COG, KEGG, GO_terms, RefSeq, inference.

### Step 3: Diagnostic analyses

Choose based on the user's question:

| Goal | Script |
|---|---|
| COG distribution | `genbank_cog.py` |
| Export protein FASTA for BLAST/HMMscan | `genbank_fasta.py` |
| Gene neighborhood for a target locus | `genbank_neighborhood.py` |
| All qualifiers for one locus tag | `genbank_locus.py` |
| Operon candidates (co-directional, gapped) | `genbank_operons.py` |

### Step 4: External bioinformatic hooks

Once proteins are exported via Step 3, invoke if available:

| Analysis | Tool | Command |
|---|---|---|
| BLAST search | `blastp` | `blastp -query proteins.faa -db nr -outfmt 6 -max_target_seqs 5 -evalue 1e-5` |
| Domain annotation | `hmmscan` | `hmmscan --tblout domains.tsv Pfam-A.hmm proteins.faa` |
| Signal peptides | `signalp` | `signalp -fasta proteins.faa -org gram- -format short` |
| Transmembrane | `tmhmm` | `tmhmm --fasta proteins.faa` |
| Antibiotic resistance | `abricate` | `abricate --db resfinder proteins.faa` |

### Step 5: Comparative mode

Run `genbank_compare.py` to scan multiple genomes for user-specified markers. Supports gene names, EC numbers, COG IDs, and KEGG orthologs:

```bash
python genbank_compare.py ./genomes/ "ladA,ssuD,EC:1.14.14.28,KEGG:K20938" matrix.tsv
```

### Step 6: Discovery mode

For uncharacterized isolates, run `genbank_discover.py` as a "metal detector":
1. **Weighted keyword scan** -- curated dictionaries with confidence levels (3/2/1)
2. **Spatial clustering** -- groups hits into candidate genomic islands (default 5 kb gap)
3. **Dark matter operon detection** -- pure (100% hypothetical) and mixed (>=60%) clusters
4. **Mobilome-adjacent island detection** -- functional genes within 10 kb of transposase/integrase

```bash
python genbank_discover.py INPUT.gbff
python genbank_discover.py INPUT.gbff --cluster-gap 10000 --operon-gap 200
```

---

## Key Principles

1. **Never load the entire file into the AI context window.** Always use scripts to parse and summarize.
2. **Scripts first, raw content never.** Parse computationally, then read the structured output.
3. **Validate first.** Run Step 1a before any downstream analysis.
4. **One canonical parser.** All scripts import from `genbank_parser.py`. Never reconstruct ad-hoc.
5. **Respect contig boundaries.** Always check `a['contig'] == b['contig']` before spatial comparisons.
6. **Keep reports concise.** Diagnostic tables, not data dumps.
7. **Know the parser's limits.** `join()`/`order()` are warned and skipped. Verify counts for eukaryotic data.

---

## Quick-Reference Command Sequence

```bash
# Define script path for convenience
SCRIPTPATH=".agent/skills/genbank-parser/scripts"

# 1a. Validate structure
python $SCRIPTPATH/genbank_validate.py     INPUT.gbff

# 1b. Metadata (contig count, organism, sizes)
python $SCRIPTPATH/genbank_metadata.py     INPUT.gbff

# 1c. Sequence extraction (genome + CDS nucleotides)
python $SCRIPTPATH/genbank_sequence.py     INPUT.gbff

# 2. Annotation table
python $SCRIPTPATH/genbank_extract.py      INPUT.gbff  annotations.tsv

# 3a-e. Diagnostics
python $SCRIPTPATH/genbank_cog.py          INPUT.gbff
python $SCRIPTPATH/genbank_fasta.py        INPUT.gbff  proteins.faa
python $SCRIPTPATH/genbank_neighborhood.py INPUT.gbff  LOCUS_TAG 5
python $SCRIPTPATH/genbank_locus.py        INPUT.gbff  LOCUS_TAG
python $SCRIPTPATH/genbank_operons.py      INPUT.gbff  150

# 5. Comparative scan
python $SCRIPTPATH/genbank_compare.py      ./genomes/  "ladA,ssuD,alkB"  matrix.tsv

# 6. Discovery mode
python $SCRIPTPATH/genbank_discover.py     INPUT.gbff
```

---

## Pressure-Testing Protocol

Before trusting the parser on a new type of input, verify with these edge cases:

### 1. Multi-contig test
Use any `.gbff` with multiple `LOCUS` records. Verify that `genbank_metadata.py` shows multiple records and `genbank_operons.py` does NOT merge genes across contigs.

### 2. Complex location test
Create a file with `join(100..200,300..400)` features. The parser should print a stderr warning with the count of skipped features.

### 3. Whitespace continuation test
Create a CDS with `/translation` continuation lines using 19 or 21 leading spaces. The parser should still concatenate them correctly.

### 4. Full GenBank with ORIGIN
Run `genbank_sequence.py` on a complete `.gbff`. Verify it produces both `.fna` and `.ffn` files with correct GC% and coding density.

### 5. Edge cases
- Empty file -- validator should error clearly
- Header only, no features -- parser returns `[]`
- Feature with no qualifiers -- parser still includes it
- Very large file (10k+ CDS) -- parser completes in seconds
