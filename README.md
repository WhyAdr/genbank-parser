# GenBank Feature Parser & Annotation Engine

A Biopython-powered genome-annotation query engine, validation suite, and CLI toolset (`gbparse`) for parsing, validating, analyzing, and mining GenBank flatfiles (`.gbff`, `.gbk`, `.gb`, `.txt`).

[![CI](https://github.com/WhyAdr/genbank-parser/actions/workflows/ci.yml/badge.svg)](https://github.com/WhyAdr/genbank-parser/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

---

## Key Capabilities & Biological Semantics

- **Standard-Compliant Core**: Built on **Biopython (`Bio.SeqIO`)** with a typed dataclass model (`GenBankDocument`, `GenBankRecord`, `GenBankFeature`).
- **Faithful Biological Locations**: Preserves `FeatureLocation` and `CompoundLocation` (`join`/`order`), calculating biological feature length (`len(location)`) and extracting biological coding sequences via `SeqFeature.extract()`.
- **GFF3 Export**: Correctly computes CDS translation phase ($0, 1, 2$) across multi-exon/joined segments in 5' $\rightarrow$ 3' transcription order, emits complete `##sequence-region` extents, and handles ordinary features without artificial parent splits.
- **QC & Semantic Validation**: Verifies translation integrity against genetic codes (`transl_table`) and `codon_start` offsets, with structured severity findings (`ERROR`, `WARNING`, `INFO`) and pseudogene tolerance.
- **Unified CLI**: Provides `gbparse` subcommands for feature search, valid local sub-region extraction, annotation diffing, genetic-code-aware codon usage, annotation-based candidate phylogenetic markers, CRISPR/Cas annotation scanning, and declarative discovery.
- **MEOR Evidence Engine**: Scans 48 curated markers across 9 hydrocarbon-degradation, biosurfactant, and bio-emulsifier categories; evaluates 7 genome-level pathway models; and reports same-contig, known-strand candidate clusters with at most N intervening bases.
- **Mobilome Evidence Engine**: Inventories every parsed record and retains field-level mobilome evidence, replicon declarations, cautious component hypotheses, catalog provenance, and explicit external-analysis handoffs in text, JSON, or normalized TSV.
- **Cohort CLI (v0.9.3)**: Supports gzip/stdin input, lossless record selection, canonical-identity SQLite cohort indexes, snapshot-bound resumable batch runs, type-parity feature queries, deterministic JSONL/BED12/NCBI candidate exports, complete annotation diffs, evidence-preserving marker comparisons, and circular-aware operon proximity candidates.

---

## Installation & Setup

Requires Python 3.10+, `biopython>=1.80`, and `pyyaml>=6.0`.

```bash
# Clone the repository
git clone https://github.com/WhyAdr/genbank-parser.git
cd genbank-parser

# Install in editable mode
pip install -e .

# Or install with test dependencies
pip install -e .[test]

# Optional static neighborhood visualization
pip install -e .[viz]
```

---

## Unified CLI (`gbparse`) Usage

```bash
# 1. Structural & biological validation
gbparse validate input.gbff [--json]

# 2. Record metadata and contig summary
gbparse summary input.gbff

# 3. Export tab-delimited annotation TSV
gbparse extract input.gbff annotations.tsv

# 4. Search features by gene, product, KO, EC, Pfam, or regex
gbparse search input.gbff --gene ladA --format tsv

# 5. Single-locus qualifier deep-dive
gbparse locus input.gbff LOCUS_TAG

# 6. View genomic neighborhood window (+/- N genes, circular-aware)
gbparse neighborhood input.gbff LOCUS_TAG 5

# Emit schema-versioned data and an independent dna_features_viewer figure
gbparse neighborhood input.gbff LOCUS_TAG 5 \
  --format json --output neighborhood.json \
  --visualize --viz-output neighborhood.svg

# 7. Extract genomic sub-region with valid local coordinates
gbparse region input.gbff --locus LOCUS_TAG --flank-genes 5 --output region.gbk

# 8. Export protein FASTA from translations
gbparse fasta input.gbff proteins.faa

# 9. Extract genome FASTA (.fna) and CDS nucleotide slices (.ffn)
gbparse sequence input.gbff --fna genome.fna --ffn cds.ffn

# 10. Calculate codon usage, RSCU, and positional GC
gbparse codon input.gbff --output codon_usage.tsv

# 11. Functional profiling (COG distribution & pathway completeness)
gbparse functional input.gbff --format tsv

# 12. Scan annotation-supported mobilome islands
gbparse discover input.gbff --ruleset mobilome --cluster-gap 5000 --format text

# 13. Multi-genome comparative presence/absence matrix
gbparse compare ./genomes/ --targets "ladA,ssuD,K20938" --output matrix.tsv

# 14. Compare two annotation versions of the same genome
gbparse diff old.gbff new.gbff --format text

# 15. Extract annotation-based candidate phylogenetic markers
gbparse phylo input.gbff --markers all --min-length 50 --output-dir ./markers/

# 16. Scan CRISPR/Cas annotations and spatial proximity
gbparse crispr input.gbff --window 15000

# 17. Convert GenBank to standard GFF3
gbparse gff input.gbff output.gff3 --include-fasta

# 18. Parse Bakta multi-isolate summary tables
gbparse batch-summary ./isolates/ --csv summary.csv

# 19. Scan MEOR and petroleum-microbiology genomic potential
gbparse meor input.gbff --min-weight 2 --format json --output meor.json

# 20. Build a replicon-centric mobilome evidence report
gbparse mobilome input.gbff --format json --min-evidence 2 --output mobilome.json

# 21. Safe declarative feature query (no Python expression evaluation)
gbparse query input.gbff \
  --where 'type == "CDS" and (gene == "ladA" or ko == "K20938") and not pseudo' \
  --select record,locus_tag,gene,product,start,end,strand,ko \
  --format tsv

# 22. Unified interoperable exports
gbparse export input.gbff --format jsonl --output features.jsonl
gbparse export input.gbff --format bed12 --output features.bed --force
gbparse export input.gbff --format ncbi-table --output-dir table2asn-candidate

# 23. Circular-aware proximity-based operon candidates
gbparse operons input.gbff --max-gap 150 --min-gap -50 --format json --output operons.json

# 24. Inventory, filter, extract, or split whole records
gbparse records list input.gbff --format tsv
gbparse records filter input.gbff --topology circular --output circular.gbff
gbparse records split input.gbff --output-dir records/

# 25. Build and query a normalized annotation index
gbparse index build ./isolates cohort.gbidx --jobs 4 --report cohort-index-report.json
gbparse index query cohort.gbidx --where 'type == "CDS" and gene == "ladA"' --format jsonl
gbparse index status cohort.gbidx --format json

# 26. Apply a registered single-input command across a cohort
gbparse batch ./isolates --command validate --output-dir validate-run --jobs 4 -- --format json
gbparse batch ./isolates --command validate --output-dir validate-run --resume -- --format json
```

### Structured output and safety contracts

Machine-readable output is written as data only; diagnostics and warnings go
to stderr. Use `--force` to replace an existing output, and outputs that alias
an input are always rejected. Plain and gzip-compressed GenBank input are
equivalent, and `-` reads stdin when a command does not need to reread or
randomly access the source. JSON is deterministic UTF-8 with a final newline;
JSONL, TSV, CSV, FASTA, BED12, and NCBI table output is never silently
truncated.

`gbparse query` accepts `==`, `!=`, `<`, `<=`, `>`, `>=`, `contains`, `~=`,
`in`, `and`, `or`, `not`, parentheses, and explicit
`qualifier("name")` access. In PowerShell, use a double-quoted argument with
single-quoted string literals, for example `--where "type == 'CDS'"`; in Bash,
use single quotes around the expression. Windows CMD users can escape inner
double quotes, for example `--where "type == \"CDS\" and gene == \"ladA\""`.
Query results use the canonical `gbparse.feature.v1` projection and biological
`SeqFeature.extract()` semantics for FFN output.

`gbparse index` stores normalized annotation projections and raw source
fingerprints; it does not archive sequence or GenBank blobs. Source identity is
the normalized canonical path, while display spelling remains mutable metadata.
Indexing parses the same byte snapshot that is fingerprinted, and database plus
optional report publication is guarded as one outcome. Indexed results are
only as current as the verified source fingerprint at build/update time, and
indexed query is the same bounded declarative grammar, not arbitrary SQL.
`gbparse batch` runs registered single-input `gbparse` commands only. Its
manifest records canonical source identity, logical replay argv, input/output
hashes, stderr hashes, and child diagnostics; a durable sibling in-progress
tree makes interrupted runs resumable without publishing transient snapshot
paths. Required adapter outputs are checked before success, and resume keeps a
terminal outcome separate from its `resume_action`. Partial job success never
makes a partially failed run successful. Record selection selects source
records and does not rewrite or reconcile their annotations.

### v0.9.3 lifecycle contracts

Batch manifests use `pending`, `running`, `succeeded`, `threshold_failed`,
`failed`, and `not_requested` statuses. Resume actions are
`reused_unchanged`, `rerun_changed`, `replayed_running`,
`deferred_after_stop`, or `not_requested`. `--resume` reuses only a verified
source fingerprint, output contract, and stderr log; a durable
`.<output>.gbparse-inprogress` tree is the recovery point after interruption.
With `--on-error stop`, an unchanged failed or threshold-failed job establishes
an ordered barrier; changed jobs before that failure still run, while later
jobs remain deferred. Change the source and rerun successfully, or intentionally
restart with `--force`, to clear that state. Batch rejects custom `discover
--rules`, `meor --markers/--pathways`, and `mobilome --database-dir` resources
until auxiliary-resource snapshot provenance is implemented; direct commands
still accept them. Logical source labels and replay argv are retained in the
manifest while physical snapshots and disposable work-tree paths are removed
from published artifacts and stderr logs. New replay argv uses absolute source
and final-run output paths and includes overwrite permission; replay remains
valid only while the original source exists with the recorded SHA-256 and the
same command/environment contract.

If final batch publication fails, the complete `.<output>.gbparse-inprogress`
tree is retained as a recovery artifact. Resume bootstrapping copies a published
run into a unique sibling and atomically installs it, so an interrupted copy
cannot poison the canonical recovery path.

Indexes use schema revision 3. `gbparse index status` and `query` fail closed
on incomplete or incompatible old revisions; migrate them explicitly with
`gbparse index migrate cohort.gbidx`, or use `index update`, which performs
the supported in-place migration. Revision-1 and revision-2 migrations rebuild
the affected tables transactionally and match a fresh revision-3 schema;
current WAL mode is preserved by no-op inspection/migration. `--prune` treats
directory discovery as the authoritative cohort and removes no-longer-discovered
sources. Use `--report-force` when replacing an existing report. If a multi-file
publication cannot roll back automatically, the reported backup and staging
paths are retained for manual recovery.

Direct feature queries reject cohort-only `sample` and `sample_key` fields;
indexed queries expose those fields from the stored cohort identity.

The CLI returns 0 for success, 1 when a requested validation threshold is
reached after writing its report, 2 for argument errors, 3 for input/data
errors, and 4 for output or serialization errors. `gbparse --version` reports
the installed release.

GenBank coordinates are one-based inclusive. GFF3 uses the same feature
coordinate convention, BED12 uses zero-based half-open intervals, and PAF
uses zero-based half-open mappings. Origin-spanning circular features are
split into linked, bounded BED records with an explicit warning because BED
intervals cannot wrap.

`gbparse export --format ncbi-table` creates candidate `.fsa`, five-column
`.tbl`, and `.export.json` files. It does not invoke `table2asn`, create
submission metadata, or claim submission readiness. Operon results are
annotation-derived proximity candidates and are not evidence of transcription
or co-expression. The v0.8.5 exclusions are recorded in
[`gbparse-v0.8.5-deferred-scope.md`](gbparse-v0.8.5-deferred-scope.md).

### MEOR confidence and interpretation

`gbparse meor` uses the following evidence tiers (MEOR catalog `1.1`):

- **Weight 3**: active structured KEGG KO, fully specified structured EC number, or matching gene symbol.
- **Weight 2**: specific product-annotation regex.
- **Weight 1**: free-text note match, including KO or EC identifiers found only in `/note`.

Wildcard ECs are contextual catalog metadata and never produce Weight-3 evidence.

Marker hits indicate **annotation-supported genomic encoding potential**. They do not prove transcription, enzyme activity, hydrocarbon turnover, biosurfactant production, or field-scale enhanced oil recovery. See the [MEOR marker and scientific provenance reference](docs/meor_markers_reference.md) for the curated catalog and its limits.

### Mobilome evidence and interpretation

`gbparse discover --ruleset mobilome` remains the compact annotation-island scanner. `gbparse mobilome` is different: it inventories chromosomes, declared plasmids, and unresolved records; keeps each matching field/value/pattern; and emits a schema-versioned, replicon-centric report.

```bash
gbparse mobilome input.gbff --format text
gbparse mobilome input.gbff --include plasmid --min-evidence 2 --format tsv --output mobilome.tsv
gbparse mobilome input.gbff --database-dir ./reviewed-mobilome-catalog --format json
```

`--include` changes detailed presentation only: every parsed record remains in inventory and cross-record scanning. `--min-evidence` changes eligibility for configured hypotheses but retains lower-tier observations in JSON and TSV. The report treats generic Rep, pXO-numbered, phage-related, AMR-like, and virulence-associated annotations as bounded observations or candidates. It does not establish replication mechanism, plasmid ancestry, phenotype, transfer, co-transfer, expression, antimicrobial resistance, virulence, or phagemid identity. Spatial marker-cluster outputs are annotation-supported, spatially constrained candidate hypotheses; they are not autonomous counts of functional toxin-antitoxin systems. In particular, a retron RT-msr/msd core candidate does not require or establish a cognate effector. See the [mobilome evidence reference](docs/mobilome_evidence_reference.md) for catalog sources, thresholds, and external-tool handoffs.

### Neighborhood visualization

Visualization is optional: the core package does not install Matplotlib or a plotting backend. Install `genbank-parser[viz]` to enable static SVG, PNG, and PDF output. The default figure is a linear-unwrapped feature track, so neighborhoods crossing a circular origin remain continuous. Labels use `/gene`, then `/locus_tag`, then a shortened `/product`; the target CDS is highlighted without using color to encode strand.

`dna_features_viewer` is the supported single-neighborhood backend in this release. Comparative multi-genome synteny and a possible pyGenomeViz backend are intentionally deferred until explicit alignment or orthology evidence can support links between tracks.

Optional context overlays reuse canonical analysis logic:

```bash
gbparse neighborhood input.gbff LOCUS_TAG 8 \
  --visualize --viz-output locus-context.svg \
  --include-feature-types CDS,tRNA,rRNA,ncRNA,tmRNA,mobile_element \
  --color-by ruleset --ruleset mobilome \
  --show-operons --operon-gap 150
```

Ruleset colors mean only that an existing annotation matched a declared mobilome or xenobiotics term; the legend deliberately labels these as `annotation-rule match`. Operon links mean same-contig, adjacent, known-common-strand CDSs within the requested intervening-base threshold. Neither overlay demonstrates expression, horizontal transfer, operon co-transcription, or phenotype.

---

## Python API Usage

```python
from genbank_parser import read_genbank, extract_xrefs

# Load full document
doc = read_genbank("input.gbff")

# Iterate contigs/records
for rec in doc.records:
    print(f"Contig: {rec.id} ({rec.length:,} bp, topology: {rec.topology}, GC: {rec.gc_content:.1f}%)")
    
# Find a specific locus and extract biological sequence
match = doc.find_locus("ABC_00123")
if match:
    rec, feat = match
    print(feat.gene, feat.product, feat.length)
    nt_seq = feat.extract(rec.seq)
    xrefs = extract_xrefs(feat)
    print("KEGG KOs:", xrefs["kegg_kos"])
```

MEOR analysis is also available as a data-returning Python API:

```python
from genbank_parser.meor import analyze_meor

report = analyze_meor("input.gbff", min_weight=2, max_gap=200)
print(report.total_hits, report.pathways)
```

Mobilome analysis is available from its intentionally separate public subpackage:

```python
from genbank_parser.mobilome import analyze_mobilome

report = analyze_mobilome("input.gbff", include="all", min_evidence=2)
print(report.inventory, report.cross_record_hypotheses)
```

---

## Testing & Verification

Run the full pytest test suite:

```bash
pytest -v
```

---

## License

MIT License
