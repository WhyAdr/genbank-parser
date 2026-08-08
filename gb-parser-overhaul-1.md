# GenBank Parser Overhaul 1

**Repository:** `WhyAdr/genbank-parser`  
**Purpose:** Implementation handoff for a substantial correctness, architecture, testing, and feature overhaul of the current GenBank parser skill/tool suite.

---

## 1. Executive Summary

The codebase has already taken an important architectural step by migrating its core GenBank parsing logic from ad hoc regular-expression parsing to Biopython `Bio.SeqIO`. That was the right decision and solved several classes of problems around qualifier continuation, partial positions, strand handling, and compound locations.

However, the repository is currently in a halfway state:

- the **canonical parser** uses Biopython,
- but several **downstream scripts flatten or manually re-parse** sequence/metadata information,
- and some modules still operate as if the old regex parser were the source of truth.

The next overhaul should therefore prioritize **preserving biological semantics end-to-end**, rather than simply adding more analysis scripts.

The strongest long-term direction is to turn this project from a collection of useful GenBank scripts into a small, reusable **genome-annotation query and evidence engine** with:

1. a typed record/feature model,
2. faithful location/sequence semantics,
3. one canonical parser,
4. rigorous regression tests,
5. a unified CLI,
6. a documented Python API,
7. configurable rulesets,
8. optional external-tool evidence.

---

# 2. Priority Matrix

## P0 — Correctness Bugs and Biological Semantics

These should be fixed before adding major new features.

1. Fix `genbank_gff.py` compound-location detection.
2. Stop extracting CDS sequence by genomic `start:end` slicing.
3. Preserve and use full Biopython `SeqFeature.location`.
4. Use biological feature length rather than genomic span.
5. Handle strand `None` / `0` correctly.
6. Repair translation validation logic.
7. Implement correct GFF3 CDS phase.
8. Handle `codon_start`, `transl_table`, partial CDSs, pseudogenes, and compound CDSs.
9. Add regression tests for all of the above.

## P1 — Architecture / Maintainability

1. Introduce typed record and feature models.
2. Remove duplicate GenBank/ORIGIN/LOCUS parsing implementations.
3. Add a real pytest suite.
4. Add CI.
5. Add `pyproject.toml`.
6. Convert scripts into an installable package.
7. Add a unified CLI.
8. Reconcile `README.md` and `SKILL.md` with actual implementation.
9. Define a stable programmatic API.
10. Add `py.typed`.

## P2 — Scientific Depth

1. Sequence-aware comparative marker detection.
2. Reciprocal BLAST/MMseqs-style validation.
3. HMM-based phylogenomic marker extraction.
4. Configurable discovery rulesets.
5. Circular-aware genomic neighborhoods.
6. Annotation diffing.
7. Proper CRISPR sequence-analysis hooks.
8. Region/locus GenBank export.

## P3 — Advanced / Strategic Features

1. Evidence scoring framework.
2. Plugin-style external bioinformatics hooks.
3. Multi-isolate genome query engine.
4. Reusable pathway/ruleset definitions.
5. Provenance-rich JSON output.
6. Schema-versioned structured output.

---

# 3. Immediate Correctness Bugs

## 3.1 `genbank_gff.py`: `join_segments` condition is wrong

The canonical parser adds `join_segments` to **every parsed feature**:

```python
'join_segments': segments,
```

where ordinary features receive an empty list.

The GFF3 exporter currently checks:

```python
if 'join_segments' in f:
```

This condition is always true.

### Consequence

Ordinary simple features can be incorrectly routed through the compound-location branch.

### Required fix

Use:

```python
if f['join_segments']:
```

or preferably expose a property:

```python
feature.is_compound
```

### Required regression test

Create a fixture containing:

- one ordinary CDS,
- one joined CDS.

Assert:

- ordinary CDS produces the normal expected hierarchy,
- joined CDS produces the compound representation,
- ordinary features are not emitted as artificial parent/segment structures.

---

# 4. Preserve Biopython Location Semantics

## Problem

`parse_features()` currently normalizes a `SeqFeature` into a dictionary containing:

- `start`
- `end`
- `strand`
- `contig`
- `qualifiers`
- `partial_start`
- `partial_end`
- `join_segments`

This is convenient but loses important information contained in the original `SeqFeature.location`.

### Examples of lost or weakened semantics

- `CompoundLocation`
- ordering of compound parts
- exact fuzzy-position object types
- circular-origin-spanning representation
- extraction behavior
- strand uncertainty
- operator semantics such as `join` vs `order`

## Required architectural change

Retain the original location object.

For example:

```python
@dataclass(slots=True)
class GenBankFeature:
    record_id: str
    record_index: int
    feature_index: int
    type: str
    location: FeatureLocation | CompoundLocation
    qualifiers: dict[str, list[str]]
    record_length: int
    topology: str | None
```

Provide convenience properties:

```python
@property
def start(self) -> int:
    return int(self.location.start) + 1

@property
def end(self) -> int:
    return int(self.location.end)

@property
def strand(self) -> int | None:
    return self.location.strand

@property
def length(self) -> int:
    return len(self.location)

@property
def is_compound(self) -> bool:
    return isinstance(self.location, CompoundLocation)
```

And biological extraction:

```python
def extract(self, record_seq):
    return self.location.extract(record_seq)
```

or retain the original `SeqFeature` and call:

```python
seqfeature.extract(record.seq)
```

---

# 5. Stop Manually Slicing CDSs

## Current problem

Modules such as `genbank_sequence.py` and `genbank_codon.py` effectively do:

```python
seq[start:end]
```

and reverse-complement if needed.

This is only correct for a simple contiguous feature.

It is wrong for compound CDSs such as:

```text
join(100..200,300..400)
```

because the intervening sequence is incorrectly included.

## Required fix

Use Biopython:

```python
seqfeature.extract(record.seq)
```

or:

```python
feature.location.extract(record.seq)
```

for every biologically meaningful feature-sequence extraction.

## Modules that should be reviewed

At minimum:

- `genbank_sequence.py`
- `genbank_codon.py`
- GFF3 feature handling
- future region extraction
- validation
- any CDS nucleotide export logic

---

# 6. Feature Length Must Mean Biological Length

## Problem

Several modules calculate:

```python
end - start + 1
```

This is genomic span, not necessarily feature length.

For:

```text
join(100..200,300..400)
```

genomic span includes the gap.

## Required behavior

Expose both explicitly:

```python
feature.span_length
feature.sequence_length
```

For example:

```python
@property
def genomic_span(self):
    return self.end - self.start + 1

@property
def biological_length(self):
    return len(self.location)
```

Use `biological_length` for:

- CDS length,
- translation validation,
- codon statistics,
- exported feature sequence length.

Use genomic span only where spatial extent is intended.

---

# 7. Strand Handling

## Problem

Current logic effectively maps:

```python
strand = '-' if feat.location.strand == -1 else '+'
```

Therefore:

- `strand == 0`
- `strand is None`

are falsely converted to `+`.

## Required fix

Preserve the raw strand state.

Recommended internal representation:

```python
-1
+1
0
None
```

Only translate to symbols at rendering time:

```text
+1   -> +
-1   -> -
0    -> ?
None -> .
```

depending on output format.

---

# 8. Rename the Misleading `line` Field

## Problem

The parser increments:

```python
feat_counter += 1
```

and stores:

```python
'line': feat_counter
```

This is not a source-line number.

## Required change

Rename to:

```python
feature_index
```

Prefer a stable composite identity:

```python
(record_index, feature_index)
```

or:

```python
FeatureID(record_id, feature_index)
```

Do not describe it as a GenBank file line number in documentation.

---

# 9. Canonical Record Model

The current repository claims to have one canonical parser, but sequence and metadata are still independently re-parsed in several scripts.

Examples include custom parsing for:

- `ORIGIN`
- `LOCUS`
- metadata
- sequence extraction

This recreates multiple parsers.

## Proposed model

```text
GenBankDocument
├── Record
│   ├── id
│   ├── name
│   ├── description
│   ├── sequence
│   ├── length
│   ├── topology
│   ├── molecule_type
│   ├── annotations
│   ├── dbxrefs
│   └── features
├── Record
└── ...
```

Suggested dataclasses:

```python
@dataclass(slots=True)
class GenBankRecord:
    id: str
    name: str
    sequence: Seq
    annotations: dict
    features: list[GenBankFeature]

@dataclass(slots=True)
class GenBankDocument:
    path: Path
    records: list[GenBankRecord]
```

Then all modules should use:

```python
document = read_genbank(path)
```

No module should directly reparse `ORIGIN` or `LOCUS`.

---

# 10. Validator Overhaul

The current validator is useful but should become biologically aware.

## 10.1 Translation length validation

Current approximation resembles:

```python
nlen = end - start + 1
expected_aa = (nlen // 3) - 1
```

This breaks for:

- joined CDSs,
- `codon_start`,
- unusual genetic codes,
- partial CDSs,
- pseudogenes,
- `transl_except`,
- CDS without terminal stop codon,
- special recoding events.

## Required approach

For each CDS:

```text
location.extract(record.seq)
        ↓
codon_start
        ↓
transl_table
        ↓
partial_start / partial_end
        ↓
pseudo / pseudogene
        ↓
transl_except / exception
        ↓
translate
        ↓
compare against /translation
```

## 10.2 Validation severity

Introduce structured finding levels:

```text
ERROR
WARNING
INFO
```

Example:

```json
{
  "severity": "WARNING",
  "code": "CDS_TRANSLATION_MISMATCH",
  "record": "contig_1",
  "locus_tag": "ABC_00123",
  "message": "Computed translation differs from /translation"
}
```

## 10.3 Pseudogene-aware validation

Missing `/translation` should not automatically be an error for:

- `/pseudo`
- `/pseudogene`
- disrupted CDSs

## 10.4 Locus tag integrity

Current duplicate detection should be made feature-aware.

Normal:

```text
gene + CDS
```

with the same locus tag is expected.

Instead group by:

```text
(record, locus_tag, biological locus)
```

and flag only genuinely conflicting reuse.

## Additional useful QC checks

Add detection for:

- CDS length not divisible by three when not partial/pseudo
- internal stop codons
- missing `/product`
- missing `/locus_tag`
- duplicate `protein_id`
- conflicting `gene` names for one locus
- impossible coordinates
- zero-length features
- feature outside sequence bounds
- translation contains illegal residues
- `codon_start` outside 1–3
- unsupported `transl_table`
- malformed EC/KO/GO identifiers
- duplicate DB xrefs
- invalid feature hierarchy

---

# 11. GFF3 Exporter Overhaul

This module deserves a dedicated standards pass.

## 11.1 Fix compound-feature branch bug

As described above:

```python
if f['join_segments']:
```

not:

```python
if 'join_segments' in f:
```

## 11.2 Correct CDS phase

Do not assign:

```text
phase = 0
```

to every segment.

Calculate phase according to CDS segment order in biological 5' → 3' orientation.

Be especially careful for reverse-strand compound CDSs.

## 11.3 Correct `##sequence-region`

Use actual record length:

```text
##sequence-region CONTIG 1 RECORD_LENGTH
```

not min/max feature coordinates.

## 11.4 Preserve topology-aware semantics

For circular chromosomes/plasmids, support features crossing the origin.

## 11.5 GFF3 IDs

Guarantee globally unique IDs.

Avoid relying only on locus tags because:

- locus tags may be absent,
- some features share a locus tag,
- malformed files may contain duplicates.

Generate deterministic fallback IDs.

## 11.6 Parent-child hierarchy

Represent:

```text
gene
  └── mRNA/transcript if relevant
       └── CDS segments
```

appropriately.

For bacterial CDS-only records, avoid inventing biologically misleading transcript structure unless necessary for standards compliance.

## 11.7 Add validation tests

Fixture set:

- simple CDS
- reverse CDS
- joined forward CDS
- joined reverse CDS
- no locus tag
- partial CDS
- circular origin crossing

---

# 12. Codon Module Overhaul

The current RSCU utility is useful but should not yet be treated as expression-grade.

## Required improvements

### 12.1 Respect translation table

Read:

```text
/transl_table
```

per CDS.

Do not hard-code one table for all input.

### 12.2 Respect codon offset

Read:

```text
/codon_start
```

### 12.3 Compound features

Extract the biological coding sequence with `SeqFeature.extract()`.

### 12.4 Partial CDSs

Either:

- exclude by default,
- include with explicit flag,
- or report separately.

### 12.5 Pseudogenes

Exclude by default.

### 12.6 Stop codons

Do not blindly remove the final codon.

Determine whether a terminal stop exists.

## Optional advanced metrics

After correctness is solved, add:

- GC1
- GC2
- GC3
- GC3s
- ENC
- CAI
- per-gene CAI
- codon-pair bias
- amino-acid frequency
- synonymous codon frequency
- stop-codon usage

Optional reference-set selection:

```bash
gbparse codon genome.gbff --reference ribosomal
```

or:

```bash
gbparse codon genome.gbff --reference-loci ribosomal_tags.txt
```

---

# 13. Phylogenomic Marker Module

The current marker extractor is a useful annotation-based scout but should not be described as robust phylogenomics.

## Current strength

The logic that prioritizes:

1. exact `/gene`,
2. exact `/product`,
3. product fallback only when `/gene` is absent,

is sensible and helps suppress false positives.

## Main limitation

Annotation names do not prove orthology.

## Recommended split

### Mode 1 — quick annotation mode

```bash
gbparse phylo genome.gbff --mode annotation
```

Use current name-based behavior.

### Mode 2 — HMM mode

```bash
gbparse phylo genome.gbff --mode hmm
```

Use curated marker HMMs.

Possible marker sets:

- GTDB bacterial marker set
- GTDB archaeal marker set
- BUSCO lineage markers
- custom HMM directory
- custom protein FASTA references

## Important change to concatenation

Do not simply concatenate raw protein sequences and call the result alignment-ready.

Correct workflow:

```text
identify ortholog
   ↓
one sequence per genome per marker
   ↓
align marker independently
   ↓
trim
   ↓
concatenate aligned homologous columns
   ↓
partition file
   ↓
tree inference
```

The tool may export marker FASTAs, but actual concatenated phylogenomic alignment should be created only after per-marker alignment.

---

# 14. Comparative Module: Rename and Expand

`genbank_compare.py` is currently closer to a multi-genome annotation/marker scanner than a true synteny/orthology engine.

## Near-term documentation change

Describe it conservatively as:

> Multi-genome marker presence/absence scanner.

## Future evidence levels

Add progressive evidence modes:

```text
Level 1 — annotation / gene / xref
Level 2 — protein similarity
Level 3 — reciprocal best hit
Level 4 — profile HMM
Level 5 — orthogroup membership
Level 6 — conserved genomic context
```

## Suggested output

```text
marker: ladA
genome: isolate_01
annotation_match: true
ko_match: K20938
identity: 72.1
query_coverage: 96.4
subject_coverage: 94.2
reciprocal_hit: true
synteny_score: 0.82
confidence: HIGH
```

## External integrations

Optional hooks:

- MMseqs2
- DIAMOND
- BLAST+
- HMMER

Do not require all tools for the basic mode.

---

# 15. Discovery Engine → Declarative Rule Engine

`genbank_discover.py` currently embeds biological keyword dictionaries directly in Python.

This does not scale well.

## Move rules into external files

Suggested structure:

```text
rulesets/
├── mobilome.yaml
├── xenobiotics.yaml
├── metal_resistance.yaml
├── stress.yaml
├── secondary_metabolism.yaml
└── custom/
```

Example:

```yaml
id: long_chain_alkane_oxidation
name: Long-chain alkane oxidation

positive:
  gene:
    - ladA
  product:
    - long-chain alkane monooxygenase
    - alkane monooxygenase
  ko:
    - K20938

negative:
  product:
    - DNA repair protein AlkB

weights:
  exact_gene: 3
  exact_ko: 3
  product_phrase: 2
```

## Evidence-aware reporting

Output should explain the score:

```text
ABC_00123
score: 8

+3 exact gene: ladA
+3 KO match: K20938
+2 product phrase: long-chain alkane monooxygenase
```

This is substantially more useful than a black-box keyword hit.

## User-defined rules

Support:

```bash
gbparse discover genome.gbff --rules my_rules.yaml
```

---

# 16. Operon Module

Current logic is essentially:

> consecutive same-strand CDSs within a chosen gap.

This is useful but should be framed as a heuristic.

## Rename basic output

Use:

> co-directional proximity clusters

or:

> operon candidates

with an explicit confidence disclaimer.

## Allow overlapping genes

Current logic excludes negative intergenic gaps.

Add:

```bash
--min-gap -50
--max-gap 150
```

Small overlaps are common in bacterial operons.

## Future scoring

Optional score could consider:

- intergenic distance
- same strand
- overlap
- terminator
- promoter
- functional coherence
- shared pathway/ruleset
- transcriptional evidence if provided

---

# 17. CRISPR Module

Current behavior should be described as **annotation mining**, not sequence-level CRISPR detection.

## Suggested modes

```bash
gbparse crispr genome.gbff --mode annotation
```

and:

```bash
gbparse crispr genome.gbff --mode sequence
```

## Sequence mode may integrate

- MinCED
- CRISPRCasTyper
- other external tools

## Structured result fields

```text
array_id
record
start
end
repeat_consensus
repeat_length
spacer_count
nearest_cas_cluster
cas_subtype
distance
confidence
```

---

# 18. Circular Genome Awareness

Circular topology should become first-class.

## Why this matters

Circular chromosomes/plasmids can contain:

- CDS crossing coordinate 1
- neighborhood relationships across the origin
- operons spanning the coordinate break
- mobile elements near the origin
- region exports that wrap around

## Required infrastructure

Record model must preserve:

```text
topology = circular | linear | unknown
record_length
```

## Neighborhood API

For circular records:

```text
... gene_4998
gene_4999
gene_5000
gene_0001
gene_0002 ...
```

should be possible when the target lies near coordinate 1.

---

# 19. Region / Locus Export

This would be a very high-value addition.

## Examples

```bash
gbparse region genome.gbff LOCUS_TAG --flank-genes 10 \
    --output locus.gbk
```

and:

```bash
gbparse region genome.gbff \
    --record chromosome \
    --start 2150000 \
    --end 2180000 \
    --output region.gbk
```

Useful options:

```text
--flank-genes N
--flank-bp N
--rebase
--preserve-coordinates
--include-source
--output-format gbk|gff3|fna|faa|json
```

## Rebase mode

If region is:

```text
2,150,000..2,180,000
```

rebase to:

```text
1..30,001
```

and rewrite all feature coordinates consistently.

This will be useful for:

- clinker
- manual BGC inspection
- antiSMASH region extraction
- BLAST/HMM follow-up
- comparative loci

---

# 20. Generic Search / Query Command

Add a flexible query interface.

Examples:

```bash
gbparse search genome.gbff \
    --gene ladA
```

```bash
gbparse search genome.gbff \
    --product "monooxygenase"
```

```bash
gbparse search genome.gbff \
    --ko K20938
```

```bash
gbparse search genome.gbff \
    --gene-regex 'ladA|ssuD|alkB'
```

```bash
gbparse search genome.gbff \
    --feature CDS \
    --product-regex 'flavin.*monooxygenase'
```

Output:

```text
record
locus_tag
gene
product
start
end
strand
KO
COG
EC
Pfam
```

Formats:

```text
text
tsv
csv
json
```

---

# 21. Annotation Diff Mode

Add comparison between two annotation versions of the same genome.

Example:

```bash
gbparse diff old.gbff new.gbff
```

Potential comparisons:

- added CDS
- removed CDS
- coordinate changes
- product-name changes
- gene-name changes
- KO changes
- EC changes
- protein sequence changes
- pseudogene changes
- split/merged genes
- changed feature types

This would be highly useful for:

- Bakta version comparisons
- Bakta vs Prokka
- Bakta vs PGAP
- pre/post polishing annotations
- manuscript reproducibility

---

# 22. Cross-reference Extraction

The current central `extract_xrefs()` helper is a good direction.

Expand carefully.

Potential categories:

```text
GO
COG
KEGG KO
Pfam
Rfam
EC
UniProt
RefSeq
UniRef
InterPro
TIGRFAM
NCBI Protein
GeneID
BioCyc
MetaCyc
dbCAN/CAZy if encoded
antiSMASH-specific xrefs where relevant
```

## Important requirement

Never infer identifier semantics only from loose substring matching if a structured prefix exists.

Deduplicate results while preserving deterministic order.

---

# 23. Functional Pathway Module

The current hard-coded pathway completeness table is useful but should become externalized.

Move pathway definitions to YAML/JSON.

Example:

```yaml
id: glyoxylate_bypass
name: Glyoxylate bypass

steps:
  - name: Isocitrate lyase
    any_of:
      ko: [K01637]
      ec: ["4.1.3.1"]

  - name: Malate synthase
    any_of:
      ko: [K01638]
      ec: ["2.3.3.9"]
```

Then:

```bash
gbparse pathway genome.gbff --panel glyoxylate_bypass
```

or:

```bash
gbparse pathway genome.gbff --panel custom.yaml
```

Distinguish:

- complete
- nearly complete
- incomplete
- ambiguous
- annotation-insufficient

Do not imply metabolic capability solely from generic product-name matches.

---

# 24. Repository Packaging

Move from loose scripts toward:

```text
genbank-parser/
├── pyproject.toml
├── README.md
├── LICENSE
├── src/
│   └── genbank_parser/
│       ├── __init__.py
│       ├── cli.py
│       ├── io.py
│       ├── model.py
│       ├── xrefs.py
│       ├── validate.py
│       ├── sequence.py
│       ├── query.py
│       ├── compare.py
│       ├── discover.py
│       ├── codon.py
│       ├── gff.py
│       └── ...
├── tests/
│   ├── fixtures/
│   └── ...
└── rulesets/
```

## CLI

Expose:

```bash
gbparse
```

Subcommands:

```text
gbparse validate
gbparse summary
gbparse extract
gbparse search
gbparse locus
gbparse neighborhood
gbparse region
gbparse fasta
gbparse codon
gbparse functional
gbparse pathway
gbparse discover
gbparse compare
gbparse diff
gbparse phylo
gbparse crispr
gbparse gff
```

Retain thin compatibility wrappers for old script names initially if needed.

---

# 25. Programmatic API

Document a supported importable API.

Example:

```python
from genbank_parser import read_genbank

doc = read_genbank("genome.gbff")

for record in doc.records:
    print(record.id, record.length)

feature = doc.find_locus("ABC_00123")

print(feature.product)
print(feature.gene)
print(feature.length)

seq = feature.extract()
```

Useful query helpers:

```python
doc.find_locus(tag)
doc.search_gene(name)
doc.search_product(pattern)
doc.search_xref("K20938")
doc.features(type="CDS")
doc.neighborhood(tag, flank_genes=5)
```

Add:

```text
py.typed
```

and consistent type annotations.

---

# 26. Structured JSON Output and Schema Versioning

If JSON output is supported across commands, define a stable schema.

Example:

```json
{
  "schema_version": "1.0",
  "tool_version": "0.2.0",
  "input": {
    "path": "genome.gbff"
  },
  "records": []
}
```

`schema_version` must be independent from package version.

Document backward compatibility expectations.

---

# 27. Test Suite

This is a mandatory part of the overhaul.

## 27.1 Testing framework

Use:

```text
pytest
```

Recommended:

```text
pytest
pytest-cov
```

Optional:

```text
hypothesis
```

for property-based coordinate/location tests.

## 27.2 Core fixtures

Create minimal synthetic GenBank fixtures for:

1. simple bacterial CDS
2. reverse-strand CDS
3. partial `<start`
4. partial `>end`
5. `codon_start=2`
6. `codon_start=3`
7. `transl_table != 11`
8. pseudogene
9. `/pseudo`
10. joined CDS
11. reverse joined CDS
12. multi-record chromosome + plasmid
13. circular chromosome
14. origin-spanning feature
15. unknown strand
16. ambiguous nucleotide sequence
17. duplicate locus tag
18. CDS without locus tag
19. CDS without product
20. CDS without translation
21. `transl_except`
22. custom feature type
23. ncRNA
24. repeat region
25. antiSMASH-like region GBK

## 27.3 Real-world fixtures

Use small trimmed examples representing:

- Bakta
- Prokka
- NCBI RefSeq
- antiSMASH
- eukaryotic GenBank with joined CDS if practical

Do not commit huge genome files.

## 27.4 Key invariants

Examples:

```python
assert feature.extract(record.seq) == expected_sequence
```

```python
assert feature.biological_length == len(feature.location)
```

```python
assert simple_feature.biological_length == simple_feature.genomic_span
```

```python
assert compound_feature.biological_length < compound_feature.genomic_span
```

```python
assert ordinary_gff_feature_is_not_treated_as_compound
```

---

# 28. CI

Add GitHub Actions.

Suggested matrix:

```text
Python 3.10
Python 3.11
Python 3.12
Python 3.13
```

Run:

```bash
python -m pip install -e .[test]
pytest -q
```

Optional:

```bash
ruff check .
mypy src/
```

Use pinned minimum-supported Biopython in at least one job and current Biopython in another if desired.

---

# 29. Documentation Cleanup

## `README.md`

Update:

- actual module count,
- installation,
- unified CLI,
- supported input types,
- known limitations,
- API example,
- testing status.

## `SKILL.md`

Rewrite old archaeology inherited from the regex parser.

Do not claim:

- manual indentation parsing,
- custom multiline join skipping,
- warnings that no longer exist,
- unsupported constructs that Biopython already handles.

Describe the actual architecture:

```text
GenBank flatfile
     ↓
Bio.SeqIO
     ↓
SeqRecord / SeqFeature
     ↓
normalized typed model
     ↓
analysis commands
```

## Terminology

Prefer:

> GenBank flatfile parser

for `.gb/.gbk/.gbff`.

NCBI “feature table” is a different format.

Either:

1. fix terminology, or
2. add genuine NCBI `.tbl` support.

---

# 30. Add Genuine NCBI Feature Table Support Later

Potential normalized input stack:

```text
.gb/.gbk/.gbff → GenBank reader
.tbl           → NCBI feature-table reader
.gff/.gff3     → GFF reader
```

All feed the same internal feature model.

This is a good P3 target after the GenBank path is stable.

---

# 31. Bakta Multi-Isolate Summary

The newer `parse_bakta_summaries.py` functionality should be integrated into the package rather than living as an isolated utility.

Possible command:

```bash
gbparse batch-summary ./bakta-results/
```

Support:

```text
CSV
TSV
Markdown
JSON
```

Potential fields:

- genome size
- GC
- CDS
- tRNA
- tmRNA
- rRNA
- ncRNA
- pseudogenes
- hypothetical proteins
- UniRef-assigned
- RefSeq-assigned
- COG-assigned
- KO-assigned
- EC-assigned
- Pfam-assigned

Allow generic GenBank fallback when Bakta `.txt` is absent.

---

# 32. Additional High-Value Feature Ideas

## 32.1 Protein-domain result importer

Support importing external TSV output from:

- HMMER
- InterProScan
- pyhmmscan
- eggNOG-mapper
- KofamScan

Then attach these results to feature identities.

## 32.2 Neighborhood comparison

Given one locus in multiple genomes:

```bash
gbparse neighborhood-compare genomes/ --gene ladA
```

Produce compact synteny tables.

## 32.3 Annotation evidence provenance

For every functional claim retain:

```text
source
database
identifier
method
confidence
```

Example:

```json
{
  "type": "KO",
  "value": "K20938",
  "source": "Bakta",
  "qualifier": "db_xref"
}
```

## 32.4 Quality score per annotation

Potential future score from:

- named gene
- curated xref
- KO
- EC
- domain support
- sequence similarity
- reciprocal match
- neighborhood support
- contradictory evidence

---

# 33. Suggested Execution Order for the Makeover

## Phase 1 — Freeze current behavior

1. Create a branch.
2. Add tests capturing currently correct behavior.
3. Add a few real fixture files.
4. Run the entire existing CLI suite.

## Phase 2 — Fix the data model

1. Introduce `GenBankRecord`.
2. Introduce `GenBankFeature`.
3. Preserve `SeqFeature.location`.
4. Preserve record sequence/topology.
5. Expose convenience properties.

## Phase 3 — Migrate consumers

Migrate in this order:

1. `genbank_sequence.py`
2. `genbank_validate.py`
3. `genbank_gff.py`
4. `genbank_codon.py`
5. `genbank_neighborhood.py`
6. `genbank_operons.py`
7. `genbank_extract.py`
8. `genbank_locus.py`
9. `genbank_compare.py`
10. `genbank_discover.py`
11. `genbank_crispr.py`
12. `genbank_phylo.py`
13. Bakta summary utility

## Phase 4 — Delete duplicate parsers

Remove custom:

- ORIGIN parsing
- LOCUS parsing
- sequence parsing

where Biopython already provides the data.

## Phase 5 — Package

1. Add `pyproject.toml`.
2. Move source under `src/genbank_parser`.
3. Add CLI entry point.
4. Keep backward-compatible wrappers if practical.

## Phase 6 — Documentation

Update:

- README
- SKILL
- examples
- API docs
- known limitations

## Phase 7 — New features

Implement:

1. `search`
2. `region`
3. `diff`
4. configurable rules
5. improved compare
6. HMM phylo mode
7. optional CRISPR sequence mode

---

# 34. Acceptance Criteria

The overhaul should not be considered complete until all of these are true.

## Core parsing

- [ ] Multi-record GenBank files parse correctly.
- [ ] Compound locations are preserved.
- [ ] Partial positions are preserved.
- [ ] Strand `None` and `0` are not silently converted to `+`.
- [ ] Circular topology is retained.
- [ ] Feature identity is stable within one parse.

## Sequence extraction

- [ ] Simple CDS extraction equals Biopython `SeqFeature.extract()`.
- [ ] Reverse CDS extraction equals Biopython.
- [ ] Joined CDS extraction equals Biopython.
- [ ] Reverse joined CDS extraction equals Biopython.

## Validation

- [ ] Pseudogenes do not produce bogus missing-translation errors.
- [ ] `codon_start` is respected.
- [ ] `transl_table` is respected.
- [ ] Compound CDS length is handled correctly.
- [ ] Findings have severity and machine-readable codes.

## GFF3

- [ ] Ordinary features are not treated as compound.
- [ ] Compound CDSs are represented correctly.
- [ ] CDS phases are correct.
- [ ] sequence-region spans full records.
- [ ] IDs are unique.
- [ ] Circular-origin cases have defined behavior.

## Codon

- [ ] Uses extracted biological CDS sequences.
- [ ] Respects translation table.
- [ ] Respects codon start.
- [ ] Handles partial/pseudogene filtering explicitly.

## Testing

- [ ] pytest suite exists.
- [ ] GitHub Actions CI passes.
- [ ] Compound-location regression tests exist.
- [ ] At least one Bakta fixture exists.
- [ ] At least one Prokka/NCBI fixture exists.
- [ ] At least one antiSMASH-style fixture exists.

## Packaging

- [ ] `pip install -e .` works.
- [ ] `gbparse --help` works.
- [ ] Python API is importable.
- [ ] `py.typed` exists.
- [ ] package version is defined once.

## Documentation

- [ ] README matches current files.
- [ ] SKILL.md matches actual parser behavior.
- [ ] parser limitations are explicit.
- [ ] GenBank flatfile vs NCBI feature-table terminology is correct.

---

# 35. Recommended First Pull Request

Keep the first PR deliberately focused.

## PR 1 — Core correctness foundation

Implement only:

1. pytest infrastructure,
2. parser model refactor,
3. preserve `SeqFeature.location`,
4. preserve record sequence/topology,
5. rename `line` → `feature_index`,
6. proper strand state,
7. sequence extraction using Biopython,
8. compound-feature biological length,
9. fix GFF `join_segments` conditional,
10. regression fixtures.

Do **not** simultaneously implement every new feature.

This PR should establish the correct foundation.

---

# 36. Recommended Second Pull Request

## PR 2 — Validator + GFF3 correctness

Implement:

- translation-aware validator,
- pseudogene logic,
- `codon_start`,
- `transl_table`,
- GFF3 phase,
- full sequence-region,
- compound/reverse compound GFF tests,
- circular topology behavior.

---

# 37. Recommended Third Pull Request

## PR 3 — Package and CLI consolidation

Implement:

- `pyproject.toml`,
- `src/` package,
- `gbparse` console command,
- backwards-compatible script wrappers,
- typed API,
- `py.typed`,
- schema/version handling.

---

# 38. Recommended Fourth Pull Request

## PR 4 — Research workflow additions

Implement:

- `gbparse search`
- `gbparse region`
- `gbparse diff`
- externalized pathway definitions
- externalized discovery rules

---

# 39. Recommended Fifth Pull Request

## PR 5 — Sequence-aware comparative genomics

Implement:

- optional MMseqs2/BLAST/DIAMOND integration,
- RBH mode,
- HMM marker mode,
- confidence/evidence model,
- neighborhood conservation.

---

# 40. Strategic End State

The project should evolve from:

```text
GenBank file
   ↓
collection of scripts
   ↓
printed tables
```

toward:

```text
                     ┌──────────────┐
                     │ GenBank file │
                     └──────┬───────┘
                            ↓
                  ┌───────────────────┐
                  │ Biopython parser  │
                  └────────┬──────────┘
                           ↓
              ┌─────────────────────────┐
              │ Typed record/feature API│
              └────────────┬────────────┘
                           ↓
       ┌────────────────────────────────────┐
       │ Annotation + sequence evidence API │
       └──────────────┬─────────────────────┘
                      ↓
     ┌──────────────────────────────────────────┐
     │ query / validate / compare / discover    │
     │ region / diff / codon / GFF / phylo     │
     └──────────────────┬───────────────────────┘
                        ↓
       text / TSV / JSON / GenBank / GFF / FASTA
```

At that point the unique value of the repository is no longer simply:

> “parse GenBank files”

because Biopython already does that extremely well.

Its value becomes:

> **turn genome annotations into compact, traceable, evidence-aware biological analyses while preserving the full semantics of the source annotation.**

That is the direction the overhaul should optimize for.
