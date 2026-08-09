# `gbparse-3-integration-plan.md`

## Integration plan: absorb `genbank-meor` into `gbparse`

**Target repository:** `WhyAdr/genbank-parser`<br>
**Source skill:** `WhyAdr/genbank-meor`<br>
**Target CLI:** `gbparse`<br>
**Recommended target release:** `genbank-parser` **v0.4.0**<br>
**Baseline inspected:**
- `genbank-parser` `main` at `f3f912ad494aa4f3c679080fbd76c46efe8d80aa` (`v0.3.0`, “deprecate legacy script wrappers”)
- `genbank-meor` `master` at `42d09105ca63a47dac78e0990767124ccde2dab9` (`v1.1.0` behavior plus note-only confidence fix)

---

## 1. Executive summary

The cleanest integration is **not** to copy `genbank_meor.py` into `genbank-parser`, and it is also **not** to reduce MEOR detection to another ordinary `gbparse discover --ruleset ...` YAML list.

`genbank-parser` has already crossed an architectural boundary in v0.3.0:

- one canonical Biopython-backed parser;
- typed `GenBankDocument`, `GenBankRecord`, and `GenBankFeature` models;
- a single supported `gbparse` CLI;
- packaged declarative rules for relatively simple annotation discovery;
- pytest coverage and Python 3.10–3.13 CI.

By contrast, `genbank-meor` remains a purpose-built standalone script containing:

- a vendored parser snapshot;
- 48 MEOR/hydrocarbon/biosurfactant marker definitions in Python literals;
- 9 biological categories;
- 7 pathway-completeness models;
- a three-tier evidence-confidence system;
- co-directional spatial clustering;
- feature-versus-marker multiplicity and physical-feature deduplication;
- domain-specific text, JSON, and TSV reporting;
- 50-kb density karyograms;
- a curated scientific provenance document.

The integration should therefore make MEOR a **first-class domain analysis package under `genbank_parser`**, exposed as:

```bash
gbparse meor INPUT.gbff
```

while reusing the canonical parser and shared typed models.

The recommended shape is:

```text
src/genbank_parser/
├── meor/
│   ├── __init__.py
│   ├── models.py
│   ├── database.py
│   ├── matcher.py
│   ├── scanner.py
│   ├── spatial.py
│   ├── pathways.py
│   ├── karyogram.py
│   └── report.py
└── data/
    └── meor/
        ├── markers.yaml
        ├── pathways.yaml
        └── provenance.yaml
```

The port should be **parity-first**. Before changing marker semantics, capture golden outputs from the current `genbank-meor` v1.1.0 and require the new `gbparse meor` implementation to match them. Only after parity should scientifically desirable changes—ambiguity resolution, EC wildcard semantics, broader source-aware evidence handling, improved BGC logic, or sequence-level confirmation—be introduced deliberately and versioned.

---

# 2. What exists today

## 2.1 `genbank-parser` v0.3.0

The current core has an installable `src/` package and a single entry point:

```toml
[project.scripts]
gbparse = "genbank_parser.cli:main"
```

with Python `>=3.10`, Biopython, PyYAML, pytest test extras, and packaged rulesets.

The canonical parsing path is:

```text
read_genbank()
    ↓
GenBankDocument
    ├── GenBankRecord
    │    └── GenBankFeature
    └── ...
```

Important existing invariants include:

1. **1-based inclusive presentation coordinates** through `GenBankFeature.start` and `.end`.
2. Preserved Biopython `FeatureLocation` / `CompoundLocation`.
3. Strand represented internally as `+1/-1/0/None`, with `strand_symbol`.
4. A globally unique `feature_index` across the parsed document.
5. Dictionary-style backward compatibility on `GenBankFeature`.
6. `parse_features()` already returning typed `GenBankFeature` objects.
7. `extract_xrefs()` extracting GO/COG/KEGG/Pfam/Rfam/EC evidence from `/db_xref`, `/note`, and `/EC_number`.

The current CLI has 18 analysis subcommands but no MEOR command.

The existing `discover` engine is declarative, but intentionally simple: rules consist primarily of terms and weights matched against gene/product/note text, followed by generic proximity clustering.

That mechanism is useful infrastructure, but its semantics are materially different from `genbank-meor`.

---

## 2.2 `genbank-meor` v1.1.0

`genbank-meor` currently provides:

### Marker system

- **48 marker definitions**
- **9 functional categories**
- evidence fields including:
  - gene-name regexes;
  - EC identifiers;
  - KEGG KO identifiers;
  - product regexes;
  - optional note regexes.

### Confidence system

Current matching semantics are:

| Weight | Evidence |
|---:|---|
| **3** | trusted structured KO, trusted EC, or matching gene symbol |
| **2** | specific product regex |
| **1** | note/product pattern or note-only KO/EC text |

The August 4 fix is especially important: a KO or EC merely written in `/note` is intentionally **not** treated as high-confidence structured evidence.

### Pathway system

Seven pathway-completeness models are evaluated as logical step sets. A step is considered present if **any** marker ID assigned to that step is found.

### Spatial clustering

Candidate MEOR operons/BGC-like clusters are created by:

- sorting hits by contig and coordinate;
- requiring same contig;
- requiring same strand;
- requiring hit-to-hit gap `<= --max-gap`;
- defaulting to `200 bp`;
- retaining clusters with at least two hit records;
- reporting both:
  - `hit_count`, representing marker evidence records;
  - `gene_count`, intended to represent distinct physical features.

### Reporting

Current outputs are:

```text
text
json
tsv
```

The text report contains:

1. category summary;
2. pathway completeness;
3. co-localized operon/BGC candidates;
4. 50-kb hit-density karyograms.

The JSON report returns:

```json
{
  "file": "...",
  "total_features": 0,
  "total_hits": 0,
  "total_clusters": 0,
  "hits": [],
  "clusters": [],
  "pathways": [],
  "karyograms": {}
}
```

The TSV report exposes one row per marker hit.

---

# 3. Integration goals

## 3.1 Primary goals

1. Make MEOR analysis a native `gbparse` capability:

   ```bash
   gbparse meor INPUT.gbff
   ```

2. Eliminate the duplicated/vendored `genbank_parser.py` from the MEOR implementation.

3. Use `read_genbank()` and typed `GenBankFeature` objects as the sole parsing source of truth.

4. Preserve `genbank-meor` v1.1.0 scientific behavior during the initial migration.

5. Preserve the three confidence tiers, especially the downgrade of note-only KO/EC evidence.

6. Preserve all nine categories, all 48 markers, all seven pathway models, and current reporting modes.

7. Move marker and pathway knowledge out of executable Python source into validated package data.

8. Add a proper automated regression suite for MEOR behavior.

9. Fold the specialized `genbank-meor` agent instructions into the main `genbank-parser` `SKILL.md`.

10. Leave room for future sequence-aware hydrocarbon annotation without making the first integration depend on HMMER, DIAMOND, BLAST, or external databases.

---

## 3.2 Non-goals for the first integration

The initial integration should **not** simultaneously attempt to:

- redesign the 48-marker biological database;
- replace annotation matching with sequence homology;
- add CANT-HYD HMM searches;
- infer actual petroleum-degradation phenotype;
- infer transcription/expression;
- make reservoir-scale MEOR performance claims;
- rewrite every `gbparse discover` abstraction;
- merge generic mobilome/xenobiotic discovery and MEOR into one universal rules language;
- preserve Python 3.8/3.9 support from the standalone skill.

The integrated command should inherit the parent package requirement of Python `>=3.10`.

---

# 4. Why MEOR should be a first-class `gbparse meor` command

A tempting design would be:

```bash
gbparse discover INPUT --ruleset meor
```

This should **not** be the primary interface.

The existing `discover` rules engine and MEOR engine differ in several fundamental ways.

| Capability | `gbparse discover` | `genbank-meor` |
|---|---|---|
| substring rules | yes | yes, plus regex |
| exact gene-pattern evidence | limited/general | yes |
| structured KEGG evidence | not first-class | yes |
| structured EC evidence | not first-class | yes |
| source-sensitive note downgrade | no | yes |
| evidence confidence tiers | generic summed weights | semantic tier 1/2/3 |
| biological categories | generic rules | 9 MEOR categories |
| pathway completeness | no | 7 pathway models |
| same-strand clustering | generic island clustering differs | required |
| marker-hit multiplicity | different | intentional |
| physical-feature dedup | not MEOR-specific | required |
| MEOR interpretation limits | no | yes |
| density karyogram | no | yes |
| domain report schema | no | yes |

Trying to encode all of this into the current generic rules format would either:

1. weaken MEOR behavior; or
2. force `discover` to become a much more complex expert-system framework.

For v0.4.0, a dedicated MEOR package is lower-risk and clearer.

Reusable primitives can still be extracted later where genuine commonality emerges.

---

# 5. Proposed package architecture

## 5.1 Target tree

```text
genbank-parser/
├── SKILL.md
├── README.md
├── pyproject.toml
├── src/
│   └── genbank_parser/
│       ├── __init__.py
│       ├── cli.py
│       ├── io.py
│       ├── model.py
│       ├── ...
│       ├── meor/
│       │   ├── __init__.py
│       │   ├── models.py
│       │   ├── database.py
│       │   ├── matcher.py
│       │   ├── scanner.py
│       │   ├── spatial.py
│       │   ├── pathways.py
│       │   ├── karyogram.py
│       │   └── report.py
│       └── data/
│           └── meor/
│               ├── markers.yaml
│               ├── pathways.yaml
│               └── provenance.yaml
├── docs/
│   └── meor_markers_reference.md
└── tests/
    ├── fixtures/
    │   └── meor_*.gb
    ├── test_meor_database.py
    ├── test_meor_matcher.py
    ├── test_meor_spatial.py
    ├── test_meor_pathways.py
    ├── test_meor_report.py
    └── test_meor_cli.py
```

A subpackage is preferable to one giant `meor.py` because the current standalone script already mixes at least six separate responsibilities.

---

# 6. Data model

Introduce explicit immutable dataclasses rather than passing anonymous dictionaries throughout the engine.

Suggested conceptual API:

```python
@dataclass(frozen=True)
class MeorMarker:
    id: str
    name: str
    category_id: str
    gene_patterns: tuple[str, ...]
    ecs: tuple[str, ...]
    kos: tuple[str, ...]
    product_patterns: tuple[str, ...]
    note_patterns: tuple[str, ...] = ()
    sources: tuple[str, ...] = ()

@dataclass(frozen=True)
class MeorCategory:
    id: str
    title: str
    role: str

@dataclass(frozen=True)
class MeorPathwayStep:
    name: str
    any_of: tuple[str, ...]

@dataclass(frozen=True)
class MeorPathway:
    id: str
    name: str
    category_id: str
    steps: tuple[MeorPathwayStep, ...]

@dataclass(frozen=True)
class MeorHit:
    feature_index: int
    contig: str
    locus_tag: str
    gene: str
    product: str
    start: int
    end: int
    strand: str
    category_id: str
    category_title: str
    marker_id: str
    marker_name: str
    weight: int
    evidence_type: str
    evidence_value: str

@dataclass(frozen=True)
class MeorCluster:
    cluster_id: str
    contig: str
    start: int
    end: int
    span_bp: int
    gene_count: int
    hit_count: int
    categories: tuple[str, ...]
    locus_tags: tuple[str, ...]
    genes: tuple[str, ...]
    members: tuple[MeorHit, ...]

@dataclass(frozen=True)
class MeorPathwayResult:
    pathway_id: str
    pathway: str
    completeness_pct: float
    steps_present: int
    total_steps: int
    steps: tuple[dict, ...]

@dataclass
class MeorReport:
    source_file: str
    total_features: int
    hits: list[MeorHit]
    clusters: list[MeorCluster]
    pathways: list[MeorPathwayResult]
    karyograms: dict[str, dict]
```

The public JSON representation can remain dictionary-based while internal logic gains stronger contracts.

---

# 7. Externalize the MEOR knowledge base

## 7.1 `markers.yaml`

Move all 48 marker definitions from `genbank_meor.py` into:

```text
src/genbank_parser/data/meor/markers.yaml
```

Example:

```yaml
schema_version: 1

categories:
  - id: cat3_long_alkane
    title: Long-Chain n-Alkanes & Heavy Paraffins (C18-C36+)
    role: Heavy paraffin wax degradation & pour point reduction

markers:
  - id: ladA
    category: cat3_long_alkane
    name: Long-chain Alkane Monooxygenase LadA/LadB
    genes:
      - '\blad[AB]\d?\b'
    ecs:
      - '1.14.14.28'
    kos:
      - K20938
    products:
      - 'long-chain alkane monooxygenase'
      - 'flavoprotein monooxygenase ladA'
    notes: []
    sources:
      - CANT-HYD
      - HMDB
```

### Loader validation

`database.py` should fail early if:

- category IDs are duplicated;
- marker IDs are duplicated;
- a marker references a missing category;
- a regex cannot compile;
- pathway marker references are unresolved;
- a weight/evidence field is malformed;
- expected schema version is unsupported.

Add integrity assertions:

```text
category count == 9
marker count == 48
pathway count == 7
```

These assertions are useful migration guards, even if counts change in a later curated release.

---

## 7.2 `pathways.yaml`

Represent pathway logic independently:

```yaml
schema_version: 1

pathways:
  - id: medium_alkane
    name: Medium Alkane Degradation Chain (C5-C16)
    category: cat2_medium_alkane
    steps:
      - name: Primary Alkane Hydroxylation
        any_of: [alkB, CYP153]
      - name: Electron Transfer (Rubredoxin System)
        any_of: [rubAB]
      - name: Alcohol Oxidation
        any_of: [alkJ]
      - name: Aldehyde Oxidation
        any_of: [alkH]
      - name: Acyl-CoA Activation
        any_of: [alkK]
      - name: Outer Membrane Alkane Import
        any_of: [alkL]
```

The key point is to encode the current semantics explicitly as `any_of`, not as an ambiguous generic list.

---

## 7.3 `provenance.yaml`

The existing `meor_markers_reference.md` contains useful scientific traceability that should not remain detached from the machine-readable marker database.

Suggested structure:

```yaml
sources:
  CANT-HYD:
    release: "1.0"
    accessed: "2026-08-04"
    note: "37 HMM profiles in source release"
  HMDB:
    release: "Hydrocarbon Monooxygenase Database"
    accessed: "2026-08-04"
  HADEG:
    release: "..."
    accessed: "2026-08-04"

markers:
  ladA:
    sources: [CANT-HYD, HMDB]
    interpretation: "Long-chain alkane oxidation ..."
```

Keep `docs/meor_markers_reference.md` as a human-facing document, but derive or validate it against these machine-readable records.

At minimum, add a test ensuring every marker has provenance.

---

# 8. Source-aware cross-reference handling

This is one of the most important integration points.

The current canonical `extract_xrefs()` intentionally searches both:

```text
/db_xref
/note
```

for convenience.

`genbank-meor`, however, requires a distinction between:

- structured annotation evidence; and
- free-text note evidence.

The standalone implementation currently works around this by cloning the qualifier dictionary, clearing `note`, and calling `extract_xrefs()` again.

That workaround should **not** be copied into the integrated code.

## Recommended change

Make source selection explicit while preserving backward compatibility:

```python
def extract_xrefs(
    feature: Any,
    *,
    include_notes: bool = True,
) -> dict[str, list[str]]:
    ...
```

Existing callers retain current behavior.

MEOR can then call:

```python
structured_xrefs = extract_xrefs(feature, include_notes=False)
all_xrefs = extract_xrefs(feature, include_notes=True)
```

An even stronger future API would attach provenance:

```python
XrefEvidence(
    kind="kegg_ko",
    value="K20938",
    source="db_xref",
)
```

but that is optional for v0.4.0.

### Acceptance rule

A KO/EC found **only in `/note`** must remain Weight 1, never Weight 3.

This deserves a dedicated regression test because it was explicitly fixed in `genbank-meor` v1.1.0.

---

# 9. Matching engine

Implement the current matching order exactly before scientific refinements.

Pseudo-code:

```python
def match_feature_to_marker(feature, marker) -> Match | None:
    structured = extract_xrefs(feature, include_notes=False)
    notes = feature.get_quals("note")
    gene = feature.gene
    product = feature.product

    # Weight 3
    if marker KO in structured KEGG:
        ...
    if marker EC in structured EC:
        ...
    if marker gene regex matches gene:
        ...

    # Weight 2
    if product regex matches product:
        ...

    # Weight 1
    if product/note regex matches free-text note:
        ...
    if marker KO/EC text appears only in note:
        ...

    return None
```

## Important: preserve precedence

A feature matching both an exact KO and a product regex should yield the highest-confidence reason according to the existing first-match ordering.

The matching function should return a structured object rather than `(bool, weight, reason)`:

```python
@dataclass(frozen=True)
class MarkerMatch:
    weight: int
    evidence_type: Literal[
        "kegg",
        "ec",
        "gene",
        "product",
        "note",
        "note_kegg",
        "note_ec",
    ]
    evidence_value: str
```

The report layer can recreate the legacy `reason` field for compatibility.

---

# 10. Scanner behavior

Use:

```python
doc = read_genbank(input_path)
```

and iterate records/features without reparsing the file.

Recommended scope:

```python
for record in doc.records:
    for feature in record.features:
        if feature.type not in {"CDS", "misc_feature"}:
            continue
```

Preserve current inclusion of `misc_feature` for parity.

Every hit should include `feature_index`.

That solves a weakness in the standalone code: physical-feature identity no longer needs to be reconstructed from `(contig, start, end, strand, locus_tag)`.

## Marker multiplicity

The current engine can assign one physical feature to more than one marker ID.

For the first integrated release:

> **Preserve this behavior.**

Do not silently collapse multi-marker assignments during the migration.

However:

```text
hit_count = number of marker-evidence records
gene_count = number of unique feature_index values
```

This gives the two counts unambiguous meanings.

---

# 11. Spatial clustering

MEOR clustering should remain separate from generic `discover` island clustering.

## Required parity semantics

For consecutive MEOR hits after coordinate sorting:

```text
same contig
AND same strand
AND next.start - previous.end <= max_gap
```

Default:

```text
max_gap = 200 bp
```

Retain only clusters containing at least two hit records.

## Recommended internal improvement

Deduplicate physical genes by:

```python
feature_index
```

rather than the standalone tuple:

```text
(contig, start, end, strand, locus_tag)
```

This is more reliable because:

- locus tags can be absent;
- duplicated annotation features can share coordinates;
- feature identity is already provided by the canonical parser.

## Do not accidentally substitute `find_operon_pairs()`

`genbank_parser.operons.find_operon_pairs()` is useful, but it is not automatically equivalent to legacy MEOR hit chaining.

The initial port should keep a dedicated `cluster_meor_hits()` function.

A generic shared spatial-clustering primitive may be extracted later after behavior is proven equivalent.

---

# 12. Pathway completeness

Port the seven pathway models without semantic change.

Current logic is:

```text
step is PRESENT if any target marker ID for that step exists anywhere among filtered hits
```

and:

```text
completeness = present_steps / total_steps × 100
```

This means pathway completeness is currently **genome-level marker completeness**, not:

- operon-local completeness;
- contig-local completeness;
- expression evidence;
- biochemical pathway flux.

The integrated documentation should say this explicitly.

## Optional future modes

Do not implement in the parity release, but leave the API open for:

```bash
--pathway-scope genome
--pathway-scope contig
--pathway-scope cluster
```

Genome-level should remain the default until a scientifically justified alternative is added.

---

# 13. Karyogram logic

Preserve the current lightweight karyogram output.

Default window:

```text
50,000 bp
```

For each contig:

```text
# = at least one MEOR hit in window
- = no MEOR hit
```

Expose the window size rather than hard-coding it internally:

```bash
gbparse meor INPUT --window-size 50000
```

The legacy default remains 50 kb.

Use `GenBankRecord.length` as the primary contig length rather than inferring contig span from the largest feature coordinate. This is an architectural improvement made possible by the typed parser and is biologically more faithful.

---

# 14. CLI design

Add a nineteenth first-class command to `src/genbank_parser/cli.py`.

Recommended interface:

```bash
gbparse meor INPUT.gbff
```

Full interface:

```bash
gbparse meor INPUT.gbff \
    --format text \
    --min-weight 1 \
    --max-gap 200 \
    --window-size 50000
```

Recommended arguments:

```text
input
--format {text,json,tsv}
--min-weight {1,2,3}
--max-gap INT
--window-size INT
--output PATH
```

### Optional but useful

```text
--markers PATH
--pathways PATH
```

These would allow experimental custom databases without editing package source.

They should be optional and default to packaged data.

### Do not add a new standalone wrapper

`genbank-parser` v0.3.0 explicitly removed legacy script wrappers.

Therefore do **not** add:

```text
scripts/genbank_meor.py
genbank_meor.py
```

inside the integrated repository.

The supported interface should be only:

```bash
gbparse meor ...
```

plus the Python API.

---

# 15. Python API

Expose a compact stable API:

```python
from genbank_parser.meor import analyze_meor

report = analyze_meor(
    "genome.gbff",
    min_weight=2,
    max_gap=200,
)
```

Possible signature:

```python
def analyze_meor(
    filepath: str | Path,
    *,
    min_weight: int = 1,
    max_gap: int = 200,
    window_size: int = 50_000,
    marker_database: MeorDatabase | None = None,
) -> MeorReport:
    ...
```

Critically:

> Core analysis functions should return data and not print.

Printing/serialization belongs in `report.py` or CLI dispatch.

This will make MEOR easier to test and reuse than the standalone implementation.

---

# 16. Output contracts

## 16.1 JSON

Preserve legacy fields while introducing explicit schema metadata.

Recommended:

```json
{
  "schema_version": "gbparse.meor.v1",
  "tool": "gbparse",
  "analysis": "meor",
  "file": "genome.gbff",
  "total_features": 0,
  "total_hits": 0,
  "total_clusters": 0,
  "parameters": {
    "min_weight": 1,
    "max_gap": 200,
    "window_size": 50000
  },
  "hits": [],
  "clusters": [],
  "pathways": [],
  "karyograms": {}
}
```

Keep the old major keys so downstream users can migrate with minimal breakage.

---

## 16.2 TSV

Preserve the existing header:

```text
contig
locus_tag
gene
start
end
strand
category_key
marker_id
marker_name
weight
reason
product
```

Optional additions such as:

```text
feature_index
evidence_type
```

should either:

1. be appended after the legacy columns; or
2. wait for a schema-versioned TSV v2.

For first-pass compatibility, append rather than reorder.

---

## 16.3 Text

Preserve the four-section mental model:

```text
[1] CATEGORY SUMMARY & HIT COUNTS
[2] MEOR PATHWAY COMPLETENESS PROFILE
[3] CO-LOCALIZED MEOR OPERONS & BGC CLUSTERS
[4] GENOMIC DENSITY KARYOGRAMS
```

Add the number of Weight-1 hits explicitly to the category table. The standalone text header currently emphasizes high and medium counts while total includes low-confidence evidence; this can be clearer without changing underlying results.

Suggested columns:

```text
High (W=3) | Medium (W=2) | Low (W=1) | Total
```

---

# 17. Agent-skill integration

The main `SKILL.md` should absorb the biologically specific routing language currently living in `genbank-meor`.

## 17.1 Frontmatter description

Extend the trigger description to include:

- MEOR;
- hydrocarbon degradation;
- petroleum microbiology;
- biosurfactants;
- bio-emulsifiers;
- `alkB`;
- `ladA`;
- `almA`;
- rhamnolipids;
- BTEX/PAH;
- anaerobic fumarate-addition markers.

The resulting skill should route these requests to:

```bash
gbparse meor
```

rather than instructing the agent to run a separate repository script.

---

## 17.2 CLI reference

Add:

| Subcommand | Purpose | Typical invocation |
|---|---|---|
| `meor` | MEOR, hydrocarbon-degradation, biosurfactant and bio-emulsifier genomic-potential analysis | `gbparse meor INPUT.gbff --min-weight 2 --format json` |

---

## 17.3 Interpretation discipline

Carry over this principle prominently:

> MEOR marker detection is evidence of **genomic encoding potential**, not proof of expression, enzyme activity, hydrocarbon turnover, biosurfactant production, or enhanced oil recovery performance.

Also distinguish:

- annotation-supported detection;
- sequence-confirmed family assignment;
- phenotype.

This becomes increasingly important if the engine later adds HMM/BLAST confirmation.

---

# 18. Documentation migration

## `README.md`

Add:

- MEOR to key capabilities;
- CLI examples;
- confidence-tier explanation;
- interpretation caveat;
- link to marker provenance reference.

## `docs/meor_markers_reference.md`

Move the current reference document into the parent repository.

Update wording from:

```text
implemented in genbank_meor.py
```

to:

```text
implemented by gbparse meor
```

## Standalone `genbank-meor`

After the integrated implementation is released and parity-tested:

1. retain the repository for historical provenance;
2. mark it deprecated/read-only in README;
3. point users to `genbank-parser >=0.4.0`;
4. optionally archive the repository after a transition period.

Do not delete its history—the standalone repository is useful provenance for the migration.

---

# 19. Testing strategy

MEOR currently has scientific logic dense enough that integration without dedicated tests would be risky.

The migration should start by writing tests **before** deleting the standalone parser.

---

## 19.1 Golden parity fixture

Create one synthetic GenBank fixture deliberately containing:

- a structured KO high-confidence hit;
- a structured EC high-confidence hit;
- an exact gene-name hit;
- a product-only Weight-2 hit;
- a note-only KO hit that must remain Weight 1;
- a note-only EC hit that must remain Weight 1;
- two same-strand nearby genes that should cluster;
- an opposite-strand gene that should split the cluster;
- two contigs;
- at least one feature matching more than one marker;
- one pathway that is complete;
- one pathway that is partial;
- one contig with no hits.

Run the legacy standalone engine and save:

```text
legacy.json
legacy.tsv
```

The integrated engine should match biologically meaningful fields exactly.

---

## 19.2 Matcher tests

`tests/test_meor_matcher.py`

Required cases:

### Weight 3

```text
/db_xref KEGG KO
/EC_number
exact gene regex
```

### Weight 2

```text
specific product regex
```

### Weight 1

```text
product keyword in /note
custom note regex
KO text in /note only
EC text in /note only
```

### Precedence

A feature matching KO + product should report Weight 3.

### Filtering

```text
min_weight=1
min_weight=2
min_weight=3
```

must produce deterministic subsets.

---

## 19.3 Database integrity tests

`tests/test_meor_database.py`

Assert:

```text
9 categories
48 markers
7 pathways
unique marker IDs
unique category IDs
all regexes compile
all markers map to valid categories
all pathway marker references resolve
all markers have provenance
```

---

## 19.4 Spatial tests

`tests/test_meor_spatial.py`

Test:

1. same contig + same strand + gap exactly 200 → clusters;
2. gap 201 → splits;
3. opposite strand → splits;
4. different contig → splits;
5. duplicate marker hits from one physical feature do not inflate `gene_count`;
6. duplicate marker hits do remain visible in `hit_count`;
7. missing locus tag does not break physical deduplication;
8. negative overlap behavior matches legacy semantics.

---

## 19.5 Pathway tests

`tests/test_meor_pathways.py`

Test:

```text
0/N steps
1/N steps
N/N steps
any_of alternatives
duplicate marker hits do not alter completeness
min_weight filtering changes pathway completeness predictably
```

---

## 19.6 Karyogram tests

Test:

- real record length is used;
- 50-kb default;
- custom window;
- boundary hit at exact window transition;
- multi-contig handling;
- empty-hit contigs.

---

## 19.7 Reporting tests

`tests/test_meor_report.py`

Test stable JSON keys:

```text
schema_version
file
total_features
total_hits
total_clusters
parameters
hits
clusters
pathways
karyograms
```

Test TSV header order.

Test that text output includes all four sections.

---

## 19.8 CLI tests

Extend `tests/test_cli.py`:

```python
assert main(["meor", str(meor_fixture)]) == 0
assert main(["meor", str(meor_fixture), "--min-weight", "2"]) == 0
assert main(["meor", str(meor_fixture), "--format", "json"]) == 0
```

Also test invalid values:

```text
--min-weight 0
--min-weight 4
--max-gap -1
--window-size 0
```

---

# 20. Scientific hardening issues discovered during integration

These should be documented now but handled **after parity**, unless a bug is severe enough to justify an intentional breaking change.

## 20.1 Ambiguous KO-to-marker mappings

Some current markers share or overlap functional identifiers.

A single annotation may therefore produce multiple marker hypotheses.

Example class of problem:

```text
one KO → more than one hydrocarbon-marker family
```

The current scanner preserves multiple hits.

For v0.4.0:

> preserve legacy marker multiplicity.

For a later release, add explicit ambiguity handling such as:

```text
candidate marker
supporting evidence
conflicting marker hypotheses
resolution status
```

Do not silently choose one family from an intrinsically non-specific KO.

---

## 20.2 Generic or wildcard EC identifiers

Definitions such as:

```text
1.14.13.-
1.14.14.-
```

are currently compared essentially as strings.

That does not implement true EC wildcard semantics against a specific annotation such as:

```text
1.14.13.7
```

A future matcher should support hierarchical EC matching explicitly.

However, enabling wildcard matching during the port would change hit counts and must therefore be a separately tested scientific change.

---

## 20.3 Broad annotation-derived false positives

Several markers—especially generic metabolic support genes, EPS export genes, phosphopantetheinyl transferases, or broad oxygenase descriptions—can be biologically real proteins without being specific evidence for an MEOR phenotype.

The confidence model mitigates this but does not eliminate it.

Future work should distinguish:

```text
core pathway-defining marker
supporting/accessory marker
contextual MEOR-relevant function
```

rather than treating all marker IDs as equally diagnostic.

---

## 20.4 Annotation evidence is not sequence confirmation

A Bakta/Prokka product named “alkane monooxygenase” is useful evidence but not equivalent to:

- HMM family confirmation;
- phylogenetic placement;
- catalytic residue validation;
- syntenic pathway evidence.

The integrated engine should retain conservative language such as:

```text
annotation-supported candidate
genomic potential
putative pathway
```

---

## 20.5 `gene_count` semantics

The standalone code describes feature deduplication as physical-gene deduplication but uses:

```text
contig + start + end + strand + locus_tag
```

The integrated parser already provides a better identity:

```text
feature_index
```

Use that internally.

This can be considered a non-breaking correctness improvement because it preserves the intended meaning of `gene_count`.

---

# 21. Implementation phases

## Phase 0 — Freeze behavior

**Goal:** make the current standalone engine testable as a reference.

Tasks:

- record exact source commit `42d0910`;
- create representative GenBank fixtures;
- generate golden JSON/TSV outputs;
- document expected counts;
- identify any known intentional inconsistencies.

Deliverable:

```text
tests/fixtures/meor_parity.gb
tests/golden/meor_v1_1_0.json
tests/golden/meor_v1_1_0.tsv
```

---

## Phase 1 — Import curated data

**Goal:** separate knowledge from code.

Tasks:

- create `data/meor/markers.yaml`;
- create `data/meor/pathways.yaml`;
- create `data/meor/provenance.yaml`;
- write validated loader;
- package data via `pyproject.toml`;
- add integrity tests.

No CLI yet.

---

## Phase 2 — Make xref evidence source-aware

**Goal:** remove the standalone qualifier-copy workaround.

Tasks:

- extend `extract_xrefs()` with `include_notes`;
- preserve default behavior for all existing commands;
- add parser regression tests;
- implement MEOR structured-vs-note evidence tests.

This phase must not change existing `gbparse search`, `functional`, etc. behavior.

---

## Phase 3 — Port matching and scanning

**Goal:** reproduce hit generation.

Tasks:

- implement `MeorHit`;
- implement marker matching;
- implement scan across typed features;
- preserve `CDS` + `misc_feature`;
- preserve marker multiplicity;
- implement `min_weight`;
- compare against golden hits.

---

## Phase 4 — Port spatial and pathway logic

Tasks:

- implement same-strand clustering;
- deduplicate by `feature_index`;
- retain `hit_count`;
- implement pathway completeness;
- implement record-length-aware karyograms;
- compare against legacy results.

---

## Phase 5 — Reporting and CLI

Tasks:

- implement text reporter;
- implement JSON serializer;
- implement TSV serializer;
- add `gbparse meor`;
- add `--output`;
- add `--window-size`;
- wire CLI dispatch;
- extend `tests/test_cli.py`.

---

## Phase 6 — Skill/documentation integration

Tasks:

- update `SKILL.md`;
- update `README.md`;
- migrate marker reference document;
- update package version to `0.4.0`;
- document source provenance from `genbank-meor v1.1.0`.

---

## Phase 7 — Deprecate standalone skill

Only after parity and release:

- update `genbank-meor` README with migration notice;
- stop maintaining its parser snapshot;
- point all new usage to `gbparse meor`;
- optionally archive the repository later.

---

# 22. Suggested commit sequence

Keep commits reviewable and bisectable.

```text
test(meor): add v1.1.0 parity fixtures and golden outputs

feat(meor): add validated MEOR marker and pathway data

refactor(io): support source-aware xref extraction

feat(meor): port tiered marker matching to typed GenBank features

feat(meor): add spatial clustering and pathway completeness

feat(meor): add karyogram and structured report models

feat(cli): add gbparse meor subcommand

test(meor): add matcher spatial pathway report and CLI regressions

docs(meor): integrate MEOR guidance and scientific provenance

chore: bump version to 0.4.0
```

Standalone repository, separately:

```text
docs: deprecate standalone genbank-meor in favor of gbparse >=0.4.0
```

---

# 23. Proposed CLI examples after integration

## Basic scan

```bash
gbparse meor genome.gbff
```

## Medium-or-higher confidence

```bash
gbparse meor genome.gbff --min-weight 2
```

## High-confidence only

```bash
gbparse meor genome.gbff --min-weight 3
```

## JSON

```bash
gbparse meor genome.gbff --format json --output meor.json
```

## TSV

```bash
gbparse meor genome.gbff --format tsv --output meor_hits.tsv
```

## Custom cluster gap

```bash
gbparse meor genome.gbff --max-gap 300
```

## Custom karyogram window

```bash
gbparse meor genome.gbff --window-size 100000
```

---

# 24. Expected end-state user workflow

Before:

```text
install/use genbank-parser
+
install/use genbank-meor
+
maintain duplicated parser snapshots
+
remember which script owns which analysis
```

After:

```text
pip install -e .
        ↓
     gbparse
        ↓
 ┌───────────────┐
 │ core parsing  │
 └───────┬───────┘
         │
 ┌───────┴────────────────────────────────────┐
 │ validate / search / region / functional... │
 │ meor                                        │
 └────────────────────────────────────────────┘
```

For the agent:

```text
user supplies .gbff/.gbk
        ↓
SKILL.md selects gbparse subcommand
        ↓
gbparse meor parses computationally
        ↓
agent reads compact text/JSON/TSV
        ↓
agent interprets genomic potential conservatively
```

No raw GenBank dump needs to be loaded into model context.

---

# 25. Acceptance criteria

The integration is ready to merge only when all of the following are true.

## Architecture

- [ ] No MEOR-local GenBank parser exists.
- [ ] MEOR uses `read_genbank()` / typed `GenBankFeature`.
- [ ] Marker knowledge is not hard-coded in one executable script.
- [ ] `gbparse` remains the sole supported CLI entry point.
- [ ] Package data installs correctly in editable and normal installs.

## Biological parity

- [ ] 9 categories load.
- [ ] 48 markers load.
- [ ] 7 pathway models load.
- [ ] note-only KO/EC evidence remains Weight 1.
- [ ] product-only evidence remains Weight 2.
- [ ] structured KO/EC and gene evidence remain Weight 3.
- [ ] `--min-weight` behavior matches standalone v1.1.0.
- [ ] same-strand clustering behavior matches legacy.
- [ ] opposite-strand hits are split.
- [ ] pathway completeness matches golden output.

## Reporting

- [ ] text output contains all four legacy analytical sections.
- [ ] JSON carries a schema version.
- [ ] legacy JSON keys are retained where practical.
- [ ] TSV preserves existing column order.
- [ ] `--output` works for JSON/TSV/text if supported.
- [ ] parameter values are recorded in JSON.

## Tests

- [ ] matcher unit tests pass.
- [ ] spatial tests pass.
- [ ] pathway tests pass.
- [ ] karyogram tests pass.
- [ ] report tests pass.
- [ ] CLI tests pass.
- [ ] all pre-existing `genbank-parser` tests pass.
- [ ] CI passes on Python 3.10–3.13.

## Documentation

- [ ] `README.md` documents `gbparse meor`.
- [ ] `SKILL.md` routes MEOR/hydrocarbon/biosurfactant requests correctly.
- [ ] genomic-potential interpretation limits are explicit.
- [ ] provenance reference is retained.
- [ ] standalone repository has a migration notice after release.

---

# 26. Post-integration roadmap

The important distinction is:

> **v0.4.0 should integrate the existing annotation-based skill faithfully.<br>
> Later versions can make it substantially smarter.**

Recommended later roadmap:

## v0.4.x — scientific correctness hardening

- EC wildcard/hierarchical matching;
- explicit marker ambiguity tracking;
- stricter role classes: diagnostic vs accessory;
- cluster-local pathway scoring;
- better handling of fragmented contigs;
- stronger tests using real Bakta and Prokka outputs.

## v0.5 — sequence-level confirmation layer

Optional confirmation modes:

```text
CANT-HYD HMM profiles
HMMER
DIAMOND/BLAST against curated proteins
domain architecture
phylogenetic placement for ambiguous oxygenases
```

Possible CLI:

```bash
gbparse meor genome.gbff --confirm hmm
```

The annotation scanner should remain usable without external databases.

## v0.6 — comparative MEOR mode

Potential interface:

```bash
gbparse meor-compare genomes/
```

or:

```bash
gbparse meor genomes/*.gbff --matrix meor_matrix.tsv
```

Outputs could compare:

- category potential;
- marker presence/absence;
- pathway completeness;
- high-confidence marker burden;
- biosurfactant classes;
- cluster architectures.

This would fit the comparative isolate-genomics use case better than repeatedly parsing single-genome reports manually.

## Future — evidence fusion

A mature MEOR framework could distinguish:

```text
annotation evidence
sequence-family evidence
phylogenetic evidence
syntenic/context evidence
experimental evidence
```

and report an evidence vector instead of a single integer weight.

That is a stronger long-term scientific model than continuously stretching the current 1/2/3 scale.

---

# 27. Recommended implementation decision

Proceed with the integration as a **new first-class `genbank_parser.meor` subpackage and `gbparse meor` command**, not as a copied script and not as a plain `discover` ruleset.

The critical architectural moves are:

1. delete the MEOR parser duplication;
2. reuse the typed `gbparse` core;
3. externalize marker/pathway/provenance data;
4. make cross-reference provenance explicit enough to preserve note-only confidence rules;
5. preserve v1.1.0 behavior through golden tests;
6. use `feature_index` for physical-feature identity;
7. keep domain-specific pathway, cluster, karyogram, and reporting semantics;
8. update the parent `SKILL.md` so the standalone skill is no longer needed;
9. defer biological-rule redesign until after parity.

This gives `gbparse` a much cleaner conceptual hierarchy:

```text
GenBank parsing
    ↓
typed genomic evidence
    ↓
general-purpose analyses
    +
domain-specific evidence engines
        └── MEOR / petroleum microbiology
```

That is a more scalable direction than maintaining independent GenBank-reading skills, each carrying its own parser snapshot.

---

# 28. Minimal implementation checklist for the coding agent

If this plan is handed directly to a coding agent, the shortest safe execution order is:

```text
1. Add MEOR parity fixture + golden outputs.
2. Add packaged marker/pathway/provenance YAML.
3. Add loader + validation tests.
4. Extend extract_xrefs(include_notes=True) without breaking callers.
5. Implement meor models.
6. Implement matcher.
7. Implement scanner.
8. Implement clustering.
9. Implement pathway evaluation.
10. Implement karyograms using GenBankRecord.length.
11. Implement serializers.
12. Add gbparse meor CLI.
13. Run golden parity tests.
14. Run complete pytest suite.
15. Update README + SKILL.md.
16. Bump to 0.4.0.
17. Only then deprecate standalone genbank-meor.
```

**Do not modify marker biology while the parity suite is still red.**

That separation between *migration* and *scientific improvement* is the most important guardrail in the entire integration.
