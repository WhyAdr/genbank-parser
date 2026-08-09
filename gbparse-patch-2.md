# gbparse Patch 2 — Correctness, Regression Hardening, and Roadmap Reconciliation

**Repository:** `WhyAdr/genbank-parser`<br>
**Target branch:** `main`<br>
**Recommended release:** `0.2.1`<br>
**Purpose:** Follow-up implementation plan after the `0.2.0` architectural overhaul. This patch should focus on correctness, regression hardening, test reliability, and aligning documentation/changelog claims with actual behavior before adding further scientific features.

---

## 1. Executive Summary

The `0.2.0` overhaul successfully implemented the most important architectural changes:

- a typed `GenBankDocument -> GenBankRecord -> GenBankFeature` model;
- preservation of Biopython `FeatureLocation` / `CompoundLocation`;
- biological sequence extraction via `location.extract(record.seq)`;
- distinction between biological feature length and genomic span;
- canonical GenBank parsing through `read_genbank()`;
- package restructuring under `src/genbank_parser/`;
- a unified `gbparse` CLI;
- `pyproject.toml`, `py.typed`, pytest scaffolding, and GitHub Actions configuration.

The core redesign is therefore sound and should **not** be replaced.

However, the current implementation still contains several correctness regressions and partially implemented roadmap items. The most important are:

1. pytest fixtures referenced by the test suite appear not to be committed because `.gitignore` excludes `*.gb`, `*.gbk`, and `*.gbff`;
2. `find_locus()` often returns the `gene` feature instead of the corresponding `CDS`, causing incorrect neighborhood and region behavior;
3. `region.py` can create internally inconsistent GenBank coordinates and destroys compound locations;
4. GFF3 phase handling still ignores `codon_start`, and synthetic hierarchy generation can duplicate IDs;
5. compatibility accessors reintroduce incorrect strand semantics for strand `0` / `None`;
6. codon analysis remains hard-coded to translation table 11;
7. `discover.py` advertises functionality that is not currently implemented (`operon_gap`, dark-matter analysis, TSV output);
8. comparative, phylogenomic, and CRISPR modules remain annotation-based rather than implementing the planned sequence/HMM layers;
9. several changelog and README claims are stronger than the implementation currently supports.

This patch should therefore be treated as a **correctness-and-truthfulness release**, not a feature-expansion release.

---

# 2. Patch Goals

## Required goals

1. Make a clean clone installable and testable.
2. Eliminate silent wrong-result behavior.
3. Preserve location semantics end-to-end.
4. Make GFF3 output biologically and structurally safer.
5. Reconcile public claims with actual functionality.
6. Add regression tests for every bug fixed here.
7. Avoid adding major new scientific dependencies in this patch.

## Explicit non-goals

The following should be deferred to a later `0.3.0` scientific-depth release:

- MMseqs2 / DIAMOND / BLAST sequence-aware comparative validation;
- reciprocal best-hit analysis;
- HMM-based GTDB/BUSCO phylogenomic marker recovery;
- MinCED / CRISPRCasTyper integration;
- full orthogroup inference;
- advanced evidence-weighted biological hypothesis scoring.

---

# 3. P0 — Test Suite Must Work From a Clean Clone

## Problem

`tests/conftest.py` references:

```text
tests/fixtures/simple_cds.gb
tests/fixtures/compound_joined.gb
tests/fixtures/special_cds.gb
tests/fixtures/multi_record_circular.gb
tests/fixtures/duplicate_locus.gb
```

but repository-wide `.gitignore` contains:

```gitignore
*.gbff
*.gbk
*.gb
```

This likely caused locally generated synthetic test fixtures to remain untracked.

A test suite that exists only in the developer working tree is not a regression suite.

## Required changes

### `.gitignore`

Add exceptions:

```gitignore
# Ignore general GenBank datasets
*.gbff
*.gbk
*.gb

# But always retain synthetic regression fixtures
!tests/fixtures/
!tests/fixtures/*.gb
!tests/fixtures/*.gbk
!tests/fixtures/*.gbff
```

### Commit fixtures explicitly

If they already exist locally:

```bash
git add -f tests/fixtures/
```

Confirm:

```bash
git ls-files tests/fixtures/
```

Expected output should include all five synthetic fixtures.

## Clean-clone verification

From a fresh temporary clone:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[test]
pytest -v
```

Do not use an existing editable install or developer environment for this verification.

## CI acceptance criterion

The GitHub Actions matrix must pass on all supported Python versions.

---

# 4. P0 — Fix `find_locus()` Ambiguity

## Problem

A normal bacterial GenBank annotation often contains both:

```text
gene
    /locus_tag="ABC_001"

CDS
    /locus_tag="ABC_001"
```

Current `GenBankRecord.find_locus()` returns the first matching feature.

In Biopython feature order, the `gene` feature frequently appears before the `CDS`.

Consequences:

- `gbparse neighborhood ABC_001` may receive a `gene` target but search a CDS-only list;
- the failure path currently falls back to `target_idx = 0`;
- this can silently show the neighborhood around the **first CDS on the contig**;
- `gbparse region --locus ABC_001 --flank-genes N` may silently fail to apply the requested gene window.

Silent wrong biological output is a P0 issue.

## Recommended model API

Add explicit locus resolution methods.

### `GenBankRecord`

```python
def find_locus_features(self, locus_tag: str) -> list[GenBankFeature]:
    return [
        f for f in self.features
        if f.locus_tag == locus_tag
    ]


def find_locus(
    self,
    locus_tag: str,
    prefer_type: str | None = "CDS",
) -> GenBankFeature | None:
    matches = self.find_locus_features(locus_tag)
    if not matches:
        return None

    if prefer_type is not None:
        for feat in matches:
            if feat.type.casefold() == prefer_type.casefold():
                return feat

    return matches[0]
```

### `GenBankDocument`

Mirror the same preference:

```python
def find_locus(
    self,
    locus_tag: str,
    prefer_type: str | None = "CDS",
) -> tuple[GenBankRecord, GenBankFeature] | None:
    ...
```

Optional convenience:

```python
def find_cds(self, locus_tag: str):
    return self.find_locus(locus_tag, prefer_type="CDS")
```

## Neighborhood behavior

Replace:

```python
except ValueError:
    target_idx = 0
```

with an explicit error.

Never silently substitute a different locus.

Suggested behavior:

```python
raise ValueError(
    f"Resolved target {locus_tag!r} is not present in CDS ordering"
)
```

or re-resolve the locus as a CDS before constructing the window.

## Required tests

Add fixture containing a normal paired `gene` + `CDS` with identical locus tag.

Tests:

```python
def test_find_locus_prefers_cds(...):
    ...

def test_neighborhood_target_is_requested_cds(...):
    ...

def test_neighborhood_never_falls_back_to_first_cds(...):
    ...

def test_region_flank_genes_uses_cds_for_paired_gene_feature(...):
    ...
```

The test should deliberately place the requested locus far from the first CDS so the old failure mode is obvious.

---

# 5. P0 — Repair Region Extraction Semantics

Affected file:

```text
src/genbank_parser/region.py
```

This module needs a focused rewrite.

## 5.1 Always keep extracted `SeqRecord` coordinates internally consistent

### Current issue

The returned `SeqRecord.seq` is always a sliced subsequence:

```python
sub_seq = full_seq[target_start - 1:target_end]
```

but when `rebase=False`, feature locations can remain approximately parent-genome coordinates.

That produces impossible records such as:

```text
sequence length = 20,000 bp
feature location = 2,005,000..2,006,000
```

A standalone GenBank `SeqRecord` must use coordinates relative to its own sequence.

### Required behavior

All returned region features should use local coordinates.

The extracted record should always be internally equivalent to:

```text
region nucleotide 1 == original target_start
```

Original coordinates should be preserved as metadata, not feature coordinates.

Possible qualifiers:

```text
/original_record="chromosome"
/original_location="2005000..2006000"
/original_start="2005000"
/original_end="2006000"
```

or record annotations:

```python
sub_record.annotations["parent_record"] = target_rec.id
sub_record.annotations["parent_start"] = target_start
sub_record.annotations["parent_end"] = target_end
```

### CLI recommendation

Keep `--rebase` temporarily for compatibility, but either:

1. deprecate it because local rebasing is required for valid standalone records, or
2. reinterpret it only as whether original-coordinate metadata is emitted.

Do **not** create a sliced `SeqRecord` with parent-coordinate feature locations.

## 5.2 Preserve `CompoundLocation`

### Current issue

Every retained feature is rebuilt as:

```python
FeatureLocation(new_st, new_en, strand=f.strand)
```

which turns:

```text
join(100..150,200..250)
```

into one continuous interval.

That violates the central architecture goal of preserving biological location semantics.

### Required implementation

Operate on the original Biopython location.

For every location part:

1. intersect it with the requested region;
2. clip only if necessary;
3. shift local coordinates by `-(target_start - 1)`;
4. preserve strand;
5. reconstruct a `CompoundLocation` if multiple parts survive.

Pseudo-interface:

```python
def _clip_and_shift_location(
    location,
    region_start_1based: int,
    region_end_1based: int,
):
    ...
```

The function should handle:

```text
FeatureLocation
CompoundLocation
BeforePosition
AfterPosition
```

as safely as practical.

If exact fuzzy-boundary preservation cannot be guaranteed after clipping, preserve the original location in metadata and document the limitation.

## 5.3 Circular region extraction

For circular records, support requested windows that cross origin.

Example:

```text
record length = 5,000,000
target near 4,999,000
flank = 5,000 bp
```

Expected sequence:

```text
4,994,000..5,000,000
+
1..4,000
```

This should produce a valid local region record.

Implement only if feasible in this patch; otherwise explicitly reject origin-crossing extraction with a clear error instead of silently truncating.

## Required tests

```python
def test_region_flank_genes_on_gene_cds_pair(...):
    ...

def test_region_coordinates_are_local_even_without_rebase_flag(...):
    ...

def test_region_preserves_compound_location(...):
    ...

def test_region_preserves_reverse_strand_compound_location(...):
    ...

def test_region_output_feature_end_never_exceeds_subrecord_length(...):
    ...
```

Optional:

```python
def test_region_circular_origin_wrap(...):
    ...
```

---

# 6. P0 — GFF3 Correctness Pass

Affected file:

```text
src/genbank_parser/gff.py
```

The original ordinary-feature `join_segments` bug is fixed. Preserve that fix.

Additional corrections are still required.

## 6.1 Initialize CDS phase from `codon_start`

Current logic starts every CDS with:

```python
current_phase = 0
```

and ordinary CDSs use:

```python
phase=0
```

This is incomplete for CDSs carrying:

```text
/codon_start=2
/codon_start=3
```

Use the feature's codon offset as the initial phase where appropriate:

```python
initial_phase = f.codon_start - 1
```

Then propagate phase over compound segments in biological 5' -> 3' order.

Be cautious with partial CDSs; document any assumptions.

## 6.2 Strengthen phase regression tests

Current joined-CDS fixture only verifies an easy case where phases happen to be:

```text
0, 0
```

Add cases that necessarily yield:

```text
0 -> 1
0 -> 2
1 -> ...
2 -> ...
```

for both strands.

## 6.3 Avoid duplicate synthetic `gene:` IDs

A normal GenBank file may already contain a `gene` feature with:

```text
/locus_tag="ABC_001"
```

If a compound CDS with the same locus tag is later encountered, synthetic generation of:

```text
gene:ABC_001
```

can duplicate the ID already emitted for the real gene.

### Recommended strategy

Pre-index features per locus tag:

```python
locus_features[tag] = [...]
```

When exporting a CDS:

1. if a corresponding real `gene` feature exists, reference its ID;
2. do not synthesize another gene;
3. synthesize transcript-level hierarchy only when needed and only with globally unique IDs.

Maintain a global emitted-ID set:

```python
emitted_ids: set[str]
```

and assert uniqueness before writing.

## 6.4 Circular-origin compound features

Do not claim full circular-origin GFF correctness unless it is actually tested.

At minimum:

- preserve segment-level output correctly;
- avoid representing an origin-spanning parent as a misleading whole-replicon span;
- add a test fixture for an origin-crossing compound feature on a circular record.

If a fully standards-compliant parent representation is deferred, document this limitation.

## Required tests

```python
def test_gff_phase_respects_codon_start_2(...):
    ...

def test_gff_phase_respects_codon_start_3(...):
    ...

def test_gff_compound_phase_nonzero_transition(...):
    ...

def test_gff_ids_are_globally_unique(...):
    ...

def test_gff_existing_gene_not_duplicated(...):
    ...

def test_gff_origin_spanning_compound_feature(...):
    ...
```

---

# 7. P0 — Fix Compatibility Strand Semantics

The typed model correctly maps:

```text
+1   -> +
-1   -> -
0    -> ?
None -> .
```

but compatibility methods currently use:

```python
'-' if self.strand == -1 else '+'
```

for:

```python
feature['strand']
feature.to_dict()['strand']
```

This reintroduces the old bug.

## Fix

Use:

```python
self.strand_symbol
```

everywhere.

Required test:

```python
def test_legacy_dict_preserves_unknown_and_unstranded_state(...):
    ...
```

Create synthetic features with strand `0` and `None` directly if GenBank fixture generation is awkward.

---

# 8. P0/P1 — Validator Refinements

Affected file:

```text
src/genbank_parser/validate.py
```

The validator is now substantially improved and should retain its present architecture.

## 8.1 Fix multi-record genome size reporting

Current report computes global feature coordinate span:

```python
max(feature.end) - min(feature.start) + 1
```

across all records.

For multiple contigs, this is not total genome length.

Use:

```python
doc.total_length
```

for total genome size.

Optionally report per-record span separately.

## 8.2 Exceptional translation qualifiers

Do not treat all CDSs as ordinary genetic-code translations.

Recognize presence of:

```text
/transl_except
/exception
```

and related exceptional biology.

Examples:

- selenocysteine;
- pyrrolysine;
- programmed ribosomal frameshifting;
- RNA editing;
- translational exceptions.

### Minimal patch behavior

If exceptional translation qualifiers are present:

- emit an `INFO` finding such as:

```text
EXCEPTIONAL_TRANSLATION_ANNOTATION
```

- do not escalate a computed mismatch or internal stop to normal severity unless explicitly justified.

Future releases can implement exact exception reconstruction.

## 8.3 Unknown translation table

Do not silently fall back to table 11 without reporting it.

Current behavior:

```python
except Exception:
    translate(table=11)
```

Instead:

```text
WARNING: UNKNOWN_TRANSLATION_TABLE
```

and skip strict comparison, or explicitly mark fallback behavior.

## Required tests

```python
def test_validator_total_length_multirecord(...):
    ...

def test_validator_transl_except_does_not_false_error(...):
    ...

def test_validator_unknown_translation_table_is_reported(...):
    ...
```

---

# 9. P1 — Make Codon Analysis Genetic-Code Aware

Affected file:

```text
src/genbank_parser/codon.py
```

Current model exposes `f.transl_table`, but codon analysis remains hard-coded to code 11.

This should be corrected for architectural consistency.

## Recommended implementation

Use Biopython codon tables:

```python
from Bio.Data import CodonTable
```

For each `transl_table`:

1. construct codon -> amino-acid mapping;
2. construct synonymous codon families;
3. count genes by translation table.

If a file contains multiple translation tables, either:

### Option A — pool per table and report separately

Best scientifically.

or:

### Option B — dominant table

Use the most common table and warn about excluded CDSs from other tables.

Do not silently reinterpret all CDSs as table 11.

## Stop-codon denominator

Current `total_codons` includes stops, but stop codons are not emitted in RSCU rows.

Choose one consistent definition:

- exclude terminal stops from `total_codons`, or
- explicitly output stop codon counts separately.

Preferred:

```text
sense_codons
terminal_stop_codons
ambiguous_codons
```

## Recommended additional output

```text
translation_tables_encountered
cds_excluded_partial
cds_excluded_pseudo
cds_excluded_nontriplet
ambiguous_codons
terminal_stops
```

---

# 10. P1 — Correct Coding Density Semantics

Current coding density sums all CDS biological lengths:

```python
sum(f.length for f in CDSs)
```

This double-counts overlapping CDSs.

## Recommended definition

Compute the union of CDS nucleotide intervals per record.

For compound CDSs:

- add each location part independently;
- merge overlapping/adjacent intervals;
- compute nonredundant covered nucleotides.

Expose:

```text
cds_feature_bp_sum
nonredundant_coding_bp
coding_density
```

Use nonredundant coverage for `coding_density`.

---

# 11. P1 — Repair `discover.py` Contract

Affected file:

```text
src/genbank_parser/discover.py
```

Current implementation and public description are inconsistent.

## 11.1 `operon_gap` is currently unused

The CLI accepts:

```bash
--operon-gap
```

but `discover_clusters()` does not use it.

Choose one:

### Preferred

Actually invoke the shared operon/proximity logic.

### Alternative

Remove `operon_gap` from this module until implemented.

Do not keep inert CLI parameters.

## 11.2 TSV output is advertised but not implemented

CLI accepts:

```bash
--format tsv
```

but the code falls through to text mode.

Implement a real TSV branch.

Suggested columns:

```text
record
cluster_id
cluster_start
cluster_end
locus_tag
gene
product
strand
rule_id
matched_term
weight
total_feature_score
```

## 11.3 “Dark-matter cluster” claim

The previous implementation contained explicit hypothetical/unknown-function cluster analysis.

Current module no longer performs that analysis.

Choose one:

1. restore it; or
2. stop calling this module a dark-matter detector.

For Patch 2, option 2 is acceptable.

## 11.4 Make packaged YAML the actual default

Move built-in rules into package resources:

```text
src/genbank_parser/rulesets/
    mobilome.yaml
    xenobiotics.yaml
```

Use:

```python
importlib.resources
```

to load them.

Avoid keeping a second copy of default mobilome rules hard-coded in Python.

CLI proposal:

```bash
gbparse discover input.gbff --ruleset mobilome
gbparse discover input.gbff --ruleset xenobiotics
gbparse discover input.gbff --rules custom.yaml
```

Allow multiple built-ins later.

## 11.5 Evidence-aware rule schema

Current flat schema:

```yaml
terms:
weight:
category:
```

is adequate for Patch 2.

Do not claim KO/domain-aware evidence unless implemented.

For `0.3.0`, evolve toward:

```yaml
id: alkane_monooxygenase
category: Alkane degradation
weight: 3

evidence:
  gene:
    - ladA
  product:
    - long-chain alkane monooxygenase
  ko:
    - K20938
  pfam:
    - ...

negative_patterns:
  - DNA repair AlkB
```

---

# 12. P1 — Restore Safe Comparative Matching

Affected file:

```text
src/genbank_parser/compare.py
```

Current product matching uses:

```python
t_low in prod
```

which can generate substring false positives.

Restore token/word-boundary logic.

Example helper:

```python
def _word_in(needle: str, haystack: str) -> bool:
    return bool(
        re.search(
            r'(?<![\w-])' + re.escape(needle) + r'(?![\w-])',
            haystack,
            re.IGNORECASE,
        )
    )
```

Typed xref matches should remain exact.

## Naming/documentation

Until sequence-based validation is implemented, describe this command as:

```text
annotation/xref marker scanner
```

not an ortholog inference engine.

---

# 13. P1 — Restore Phylogenetic Marker Safety

Affected file:

```text
src/genbank_parser/phylo.py
```

The current implementation silently keeps the first hit per marker.

This loses useful safety behavior from the previous version.

## Required fixes

For each marker collect **all** candidates:

```python
marker -> list[GenBankFeature]
```

Classify:

```text
0 hits  -> ABSENT
1 hit   -> SINGLE_COPY
>1 hits -> MULTI_COPY / REVIEW
```

Restore:

- minimum translation length;
- explicit multi-copy warnings;
- optional gene/product mismatch warnings.

Do not create a concatenated multi-marker “alignment” unless actual per-marker alignment is implemented.

## Terminology

Until HMM orthology is implemented:

```text
annotation-based phylogenetic marker candidates
```

is more accurate than “validated phylogenomic markers”.

---

# 14. P1 — CRISPR Annotation Scanner Cleanup

Affected file:

```text
src/genbank_parser/crispr.py
```

The module is useful as an annotation scanner, but not a sequence-level CRISPR detector.

## 14.1 Fix spatial filtering

Use one interval-distance helper:

```python
def interval_distance(a_start, a_end, b_start, b_end):
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0
```

Filter using the same value that is later printed.

## 14.2 Stabilize return schema

Current behavior can return:

```python
'colocalized': []
```

in one case and:

```python
'colocalized': 3
```

in another.

Use a stable schema:

```python
{
    "arrays": [...],
    "cas_genes": [...],
    "colocalized": [...],
    "colocalized_count": 3,
}
```

## Documentation

Rename wording to:

```text
CRISPR/Cas annotation scanner
```

until sequence-level detection is added.

---

# 15. P1 — Improve Annotation Diff Identity Matching

Affected file:

```text
src/genbank_parser/diff.py
```

Current matching identity:

```text
(record_id, start, end, strand)
```

means any boundary shift becomes:

```text
one removed CDS
+
one added CDS
```

rather than a coordinate shift.

## Recommended matching hierarchy

Use a staged matching approach:

### Tier 1 — stable identifiers

```text
same locus_tag
same protein_id
```

### Tier 2 — coordinate overlap

For unmatched features:

```text
same record
same strand
reciprocal overlap >= threshold
```

Suggested threshold:

```text
>= 80%
```

### Tier 3 — translated sequence hash

If translations are present:

```text
SHA-256(sequence)
```

can identify identical proteins despite coordinate/locus-tag changes.

## Classifications

Report:

```text
unchanged
product_changed
gene_changed
xref_changed
boundary_shifted
added
removed
possible_split
possible_merge
```

Full split/merge classification may be deferred.

## EC-number differences

The changelog currently claims KO/EC differences, but code only compares KOs.

Add EC comparison or weaken documentation.

---

# 16. P1 — Packaging and Version Hygiene

## 16.1 Python version

Current project metadata says:

```toml
requires-python = ">=3.8"
```

but source code uses PEP 604 union syntax:

```python
str | Path
int | None
```

This requires Python 3.10+.

Change:

```toml
requires-python = ">=3.10"
```

Update README badge and installation docs accordingly.

CI already begins at Python 3.10, which is appropriate.

## 16.2 Add actual MIT license file

Add:

```text
LICENSE
```

containing the MIT license.

Do not rely only on README badge / package metadata.

## 16.3 One package-version source of truth

Current version exists in:

```text
pyproject.toml
src/genbank_parser/__init__.py
```

Prefer metadata lookup:

```python
from importlib.metadata import version

__version__ = version("genbank-parser")
```

Keep package version only in `pyproject.toml`.

---

# 17. P1 — Documentation and Changelog Reconciliation

Several current claims should be narrowed until later roadmap stages are implemented.

## Replace or qualify

### Current-style claim

```text
100% backward compatibility
```

### Recommended

```text
Legacy script paths are preserved through lightweight wrappers.
Some command-line argument contracts changed in 0.2.0; users should prefer
the unified `gbparse` CLI.
```

### Current-style claim

```text
accurate GFF3 phase calculation
```

### Recommended after Patch 2

Only retain if `codon_start`, compound transitions, and reverse-strand cases are tested.

### Current-style claim

```text
phylogenomic marker extraction
```

### Recommended

```text
annotation-based candidate phylogenetic marker extraction
```

until HMM validation exists.

### Current-style claim

```text
CRISPR detector
```

### Recommended

```text
CRISPR/Cas annotation scanner
```

until sequence-level array detection exists.

### Current-style claim

```text
dark-matter discovery
```

Remove unless dark-matter clustering is restored.

### Current-style claim

```text
coordinate shifts
```

Only claim this after overlap/identifier-aware diff matching is implemented.

---

# 18. Test Matrix for Patch 2

The following test set should exist before tagging `0.2.1`.

## Parser/model

- ordinary forward CDS;
- ordinary reverse CDS;
- strand `None`;
- strand `0`;
- joined forward CDS;
- joined reverse CDS;
- partial start;
- partial end;
- `codon_start=2`;
- `codon_start=3`;
- alternate `transl_table`;
- pseudogene;
- exceptional translation qualifier;
- multi-record chromosome + plasmid;
- circular record.

## Locus/neighborhood

- paired `gene`+`CDS` locus prefers CDS;
- requested target is never replaced by index 0;
- circular neighborhood wraps correctly;
- linear neighborhood clamps correctly.

## Region

- locus extraction;
- ±N gene extraction;
- bp-flanked extraction;
- valid local coordinates;
- preserved compound location;
- reverse-strand compound location;
- all feature coordinates <= output sequence length;
- optional circular-origin extraction or explicit rejection.

## GFF3

- ordinary feature not split;
- compound feature split correctly;
- phase 0;
- phase 1;
- phase 2;
- `codon_start=2`;
- `codon_start=3`;
- reverse-strand joined CDS;
- globally unique IDs;
- no duplicate synthetic gene;
- complete `##sequence-region`;
- circular-origin feature fixture.

## Validator

- missing product;
- missing translation;
- pseudogene tolerance;
- translation mismatch;
- internal stop;
- alternative start;
- alternate translation table;
- exceptional translation;
- duplicate locus tag;
- out-of-bounds coordinate;
- total genome length across multiple records.

## Codon

- compound CDS biological extraction;
- reverse strand;
- translation table awareness;
- terminal stop accounting;
- ambiguous codon accounting;
- pseudogene exclusion.

## Discovery

- built-in YAML loads after package install;
- custom YAML loads;
- `--format text`;
- `--format json`;
- `--format tsv`;
- no unused/inert CLI arguments.

## Compare

- exact gene match;
- exact KO match;
- product word match;
- substring false-positive regression.

## Phylo

- zero-hit marker;
- one-hit marker;
- multi-copy marker;
- short-fragment exclusion.

## CRISPR

- interval-distance boundary;
- stable schema for zero hits;
- stable schema for positive hits.

## Diff

- identical feature;
- product rename;
- gene rename;
- KO change;
- EC change;
- 1–3 bp boundary shift recognized as shift rather than add/remove.

## CLI

Smoke-test every subcommand:

```bash
gbparse validate
gbparse summary
gbparse extract
gbparse search
gbparse locus
gbparse neighborhood
gbparse region
gbparse fasta
gbparse sequence
gbparse codon
gbparse functional
gbparse discover
gbparse compare
gbparse diff
gbparse phylo
gbparse crispr
gbparse gff
gbparse batch-summary
```

Also smoke-test retained scripts where backward wrappers are intentionally supported.

---

# 19. Suggested Implementation Order

## Commit 1 — test infrastructure

```text
test: commit GenBank fixtures and harden clean-clone CI
```

- `.gitignore` exceptions;
- add fixture files;
- verify pytest from fresh environment.

## Commit 2 — locus resolution

```text
fix: prefer CDS locus resolution and remove silent neighborhood fallback
```

- `model.py`;
- `neighborhood.py`;
- `region.py`;
- regression tests.

## Commit 3 — region correctness

```text
fix: preserve compound locations and valid local coordinates in region exports
```

## Commit 4 — GFF3

```text
fix: handle codon_start phases and unique feature hierarchy in GFF3
```

## Commit 5 — validator + codon

```text
fix: refine exceptional translation validation and genetic-code-aware codon analysis
```

## Commit 6 — discovery + comparison

```text
fix: reconcile discovery output contracts and restore safe marker matching
```

## Commit 7 — phylo + CRISPR + diff

```text
fix: restore marker-copy safeguards and stabilize CRISPR/diff semantics
```

## Commit 8 — packaging/docs

```text
chore: align package metadata, licensing, changelog, and documented capabilities
```

Then tag:

```text
v0.2.1
```

---

# 20. Release Acceptance Criteria

Do not tag the patch until all of the following are true:

- [ ] `pip install -e .[test]` succeeds from a fresh clone.
- [ ] All committed tests run without private/untracked fixtures.
- [ ] GitHub Actions passes on every supported Python version.
- [ ] `find_locus()` cannot silently resolve a gene when a CDS is required.
- [ ] neighborhood never silently falls back to the first CDS.
- [ ] region output always has internally valid coordinates.
- [ ] compound feature locations survive region extraction.
- [ ] GFF phase tests exercise 0, 1, and 2.
- [ ] GFF IDs are globally unique.
- [ ] compatibility strand output preserves `?` and `.`.
- [ ] validator reports total multi-record genome length correctly.
- [ ] codon analysis is translation-table aware or explicitly rejects unsupported mixtures.
- [ ] `discover --format tsv` actually emits TSV.
- [ ] all advertised discovery CLI options affect behavior.
- [ ] compare product matching is token/word safe.
- [ ] phylo reports multi-copy candidates rather than silently choosing the first.
- [ ] CRISPR return schema is stable.
- [ ] diff documentation matches actual shift/xref detection.
- [ ] `requires-python` is `>=3.10`.
- [ ] `LICENSE` exists.
- [ ] package version has one source of truth.
- [ ] README/SKILL/changelog no longer overclaim unimplemented P2 features.

---

# 21. Deferred `0.3.0` Scientific-Depth Roadmap

After `0.2.1` is stable, scientific expansion can proceed without destabilizing the core parser.

## Comparative genomics

Add optional:

```text
annotation evidence
+ protein similarity
+ reciprocal best hit
+ HMM profile
+ neighborhood conservation
```

Potential backends:

- MMseqs2;
- DIAMOND;
- BLAST+.

## HMM-based marker recovery

Add:

```bash
gbparse phylo --backend annotation
gbparse phylo --backend hmm
```

Potential marker sets:

- GTDB bacterial markers;
- GTDB archaeal markers;
- BUSCO lineage sets;
- user-provided HMM panels.

## Evidence-aware declarative discovery

Expand YAML schema to support:

```text
gene symbol
product
KO
COG
Pfam
HMM
negative evidence
neighborhood evidence
copy number
strand/orientation
```

Return an evidence trace explaining every score.

## CRISPR sequence mode

Potential hooks:

```text
MinCED
CRISPRCasTyper
```

Output:

```text
repeat consensus
repeat count
spacer count
Cas subtype
array–Cas distance
confidence
```

## Annotation-diff v2

Add sequence- and overlap-aware matching to classify:

```text
boundary shifts
splits
merges
renames
protein-preserving gene-model changes
```

---

# 22. Final Assessment

The `0.2.0` makeover should be preserved as the new architectural base.

The next patch should **not** be another rewrite.

The strongest path is:

```text
typed architecture
        ↓
adversarial regression tests
        ↓
correctness hardening
        ↓
truthful documentation
        ↓
0.2.1 stable core
        ↓
0.3.0 scientific-depth features
```

The central principle for every future module should remain:

> **Preserve the full biological semantics of the source annotation, and never silently replace uncertainty or failure with plausible-looking output.**
