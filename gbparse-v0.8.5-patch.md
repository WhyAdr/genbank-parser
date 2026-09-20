# `genbank-parser` v0.8.5 CLI Interoperability and Query Patch

**Target repository:** `WhyAdr/genbank-parser`  
**Audited baseline:** `67456fd58b05d0240d2d99394033b5376887d0f7` (`main`, package v0.8.2)  
**Target package version:** `0.8.5`  
**Implementation audience:** Luna  
**Document status:** implementation plan; pseudocode blocks are intentionally non-literal  

## 1. Objective

Turn the existing collection of 20 `gbparse` subcommands into a more coherent,
pipe-safe, machine-readable CLI while adding three tightly scoped capabilities:

1. a safe declarative `gbparse query` command;
2. a unified `gbparse export` command for interoperable representations;
3. a standalone, structured `gbparse operons` command built on the existing
   circular-aware operon engine.

This release must also remove several downstream-workflow hazards in existing
commands: misleading format names, inconsistent output handling, validation
that cannot gate a pipeline, discarded comparison evidence, and silent
truncation of machine-readable annotation diffs.

The release remains a GenBank parsing and annotation-evidence tool. It must not
grow into a general workflow engine, sequence-search suite, pangenome package,
or annotation consensus caller.

The binding exclusions for this release are maintained separately in
`gbparse-v0.8.5-deferred-scope.md`.

## 2. Release invariants

These requirements apply to every workstream.

### 2.1 Biological semantics

- Continue using `read_genbank()` and the existing typed `GenBankDocument`,
  `GenBankRecord`, and `GenBankFeature` models as the canonical parser.
- Use `SeqFeature.extract()` semantics for biological nucleotide extraction.
- Preserve compound locations, biological segment order, partial positions,
  strand states `+1`, `-1`, `0`, and `None`, record topology, and
  origin-spanning features.
- Never infer function from a text or identifier match without reporting the
  exact field and matched value.
- Standalone operon results remain proximity-based **operon candidates**, not
  evidence of transcription or co-expression.
- NCBI export files are **candidate table2asn inputs**, not a claim that an
  assembly is submission-ready.

### 2.2 Input immutability and output safety

- No output path may resolve to the input GenBank file, even with `--force`.
- Existing output files require `--force` unless a legacy command already has
  a documented overwrite contract that must be retained temporarily.
- File output must use same-directory temporary files followed by atomic
  replacement where the platform permits it.
- A failed command must not leave a partial final output or a half-complete
  NCBI export directory.
- Data go to stdout; diagnostics, warnings, and progress go to stderr.

### 2.3 Determinism and structured data

- Stable sort keys must be explicit; never depend on set or filesystem order.
- JSON is UTF-8 with deterministic key order chosen by the serializer contract,
  a final newline, and no `NaN` or infinity.
- TSV/CSV use `newline=""` and fixed declared columns.
- Structured outputs use JSON `null` or empty strings according to the schema;
  never use display placeholders such as `"-"` as missing-data values.
- Human-readable truncation is permitted only when explicitly reported.
- JSON, JSONL, TSV, CSV, FASTA, BED, and submission-input outputs must never be
  silently truncated.

### 2.4 Backward compatibility

- Retain all current subcommand names and documented positional arguments.
- `metadata` remains an alias of `summary`.
- `extract`, `fasta`, `sequence`, and `gff` remain supported even though
  `export` provides a unified route.
- Preserve existing MEOR, mobilome, and neighborhood schema versions and
  golden files unless their payload genuinely changes.
- Do not silently change an existing JSON top-level type in a patch release.
  New schema envelopes apply to newly structured commands; legacy JSON arrays
  remain arrays in v0.8.5.
- Any deprecated flag must still work, emit one concise stderr warning, and be
  covered by a test.

## 3. Proposed v0.8.5 command surface

### 3.1 New commands

```text
gbparse query INPUT --where EXPR [--select FIELDS]
              [--format text|tsv|csv|json|jsonl]
              [--emit rows|faa|ffn] [--output PATH] [--force]

gbparse export INPUT --format annotations-tsv|jsonl|faa|ffn|fna|gff3|bed12|ncbi-table
               [--output PATH | --output-dir DIR]
               [--include-fasta] [--force]

gbparse operons INPUT [--max-gap 150] [--min-gap -50]
                [--min-genes 3] [--format text|tsv|json|gff3]
                [--output PATH] [--force]
```

### 3.2 Commonly normalized commands

```text
gbparse validate INPUT [--format text|json] [--output PATH]
                 [--fail-on never|error|warning|info] [--force]

gbparse summary INPUT [--format text|tsv|json] [--output PATH] [--force]
gbparse locus INPUT TARGET [--format text|tsv|json] [--output PATH] [--force]
gbparse crispr INPUT [--window BP] [--format text|tsv|json]
               [--output PATH] [--force]
gbparse phylo INPUT ... [--format text|tsv|json]
              [--output PATH] [--output-dir DIR] [--force]
gbparse functional INPUT [--format text|tsv|json] [--output PATH] [--force]
gbparse discover INPUT ... [--output PATH] [--force]
gbparse diff OLD NEW [--format text|json|jsonl|tsv]
             [--summary-only] [--max-display N]
             [--coordinate-map FILE.paf] [--output PATH] [--force]
gbparse compare INPUTS... --targets TARGETS
                [--mode count|presence|status]
                [--output MATRIX.tsv] [--evidence-output EVIDENCE.tsv]
                [--format text|tsv|json] [--force]
```

### 3.3 Root options

Add only universally meaningful root options:

```text
gbparse --version
gbparse --help
```

Do not add global `--threads`, `--jobs`, `--manifest`, or `--resume` in this
release.

## 4. Workstream F01 — shared I/O and renderer infrastructure

### Problem

Several command modules mix analysis, printing, file writing, and process
termination. This creates divergent behavior and makes Python API reuse hard.

### Required design

Introduce small internal utilities rather than a framework:

```text
src/genbank_parser/cli_io.py
src/genbank_parser/serializers.py
```

Suggested responsibilities:

```python
@dataclass(frozen=True)
class InputSource:
    label: str
    path: Path | None
    stream: TextIO
    compressed: bool

def open_input(source: str | Path) -> ContextManager[InputSource]: ...
def reject_input_output_collision(inputs, output): ...
def atomic_text_writer(path, *, force): ...
def write_text(rendered, output, *, force): ...
def eprint(message): ...
```

Renderer functions must accept data and return text; they must not parse files,
call `sys.exit()`, or decide the destination:

```python
def render_validation_text(report: ValidationReport) -> str: ...
def render_validation_json(report: ValidationReport) -> str: ...
```

### Rough pseudocode diff

```diff
- def analyze_functional(...):
-     ...
-     print(json.dumps(report, indent=2))
-     return report
+ def build_functional_report(...) -> FunctionalReport:
+     ...
+     return report
+
+ def render_functional(report, format_type: str) -> str:
+     if format_type == "json": ...
+     if format_type == "tsv": ...
+     return ...
```

### Acceptance criteria

- Core analyzers can be called under `capsys` without producing stdout/stderr.
- Renderers are deterministic and independently unit-tested.
- CLI orchestration owns printing, paths, and exit status.
- No library-level `sys.exit()` remains in functions imported by the unified
  CLI; standalone `main()` wrappers may translate exceptions to exit codes.

## 5. Workstream F02 — gzip and stdin parsing

### Required behavior

`read_genbank()` must accept:

- uncompressed `.gb`, `.gbk`, `.gbff`, and compatible text paths;
- gzip-compressed inputs, including `.gbff.gz`;
- `-` from the CLI for uncompressed or gzip-compressed stdin;
- a text stream through a new internal/API entry point.

Prefer magic-byte detection over suffix-only detection. Decode as UTF-8 with
the repository's existing replacement policy unless strict parsing is selected
in a later release.

### Proposed API

```python
def read_genbank(source: str | Path | TextIO) -> GenBankDocument: ...
def iter_genbank(source: str | Path | TextIO) -> Iterator[GenBankRecord]: ...
```

`read_genbank()` may materialize the returned records as it does today;
`iter_genbank()` establishes a reusable streaming path for exporters. Do not
redesign every analyzer around streaming in this release.

For stdin, `GenBankDocument.path` cannot truthfully be a real path. Either:

1. change the model to a source label that safely supports `Path | None`; or
2. retain the existing model and use a clearly synthetic `Path("-")` while
   structured reports separately declare `source_kind: "stdin"`.

Prefer option 1 if it does not cause broad API churn.

### Tests

- plain and gzipped fixture parse to equivalent typed documents;
- gzip is detected despite an unusual suffix;
- stdin plain/gzip parity;
- truncated gzip fails clearly;
- malformed input produces no output file;
- commands that require rereading or random access reject incompatible stdin
  usage with a concise message rather than a traceback.

## 6. Workstream F03 — CLI status, diagnostics, and output safety

### Required additions

```python
parser.add_argument(
    "--version",
    action="version",
    version=f"%(prog)s {__version__}",
)
```

Centralize domain-to-exit-code translation:

| Code | Meaning |
|---:|---|
| 0 | completed and policy threshold not reached |
| 1 | analysis completed but requested validation threshold reached |
| 2 | argparse usage error |
| 3 | input parse/data error |
| 4 | output serialization or publication error |

Do not broadly catch `Exception` and mislabel programming bugs as input errors.
Catch typed project exceptions, `OSError`, Biopython parsing errors where their
types are stable, and renderer validation failures.

### Output publication

```diff
- Path(output).write_text(rendered, encoding="utf-8")
+ with atomic_text_writer(output, force=args.force) as handle:
+     handle.write(rendered)
```

For multi-file export, build all files in a staging directory, validate them,
then publish the directory contents. If the destination exists, fail before
doing work unless `--force` is supplied.

### Tests

- output cannot overwrite any input path through symlink or path alias;
- existing output requires `--force`;
- serialization failure leaves the old output untouched;
- stdout contains only data when output is stdout;
- progress and warnings go to stderr;
- `gbparse --version` matches installed package metadata and wheel assertions.

## 7. Workstream F04 — normalize existing structured outputs

### Commands requiring structured modes

Implement renderer-independent result models for:

- `summary`;
- `locus`;
- `crispr`;
- `phylo`;
- `functional`.

Add `--output` and `--force` to `discover` and any normalized commands that
lack them.

### Functional TSV correction

`functional --format tsv` must emit actual TSV. Recommended row model:

```text
schema_version row_type pathway step status matched_refs steps_present
steps_total completeness_pct cog_id cog_count
```

Use separate `pathway`, `pathway_step`, and `cog` rows. Retain a `text` format
for the existing human report.

### Missing values

Human text may display `-`; structured serializers use empty TSV fields or JSON
`null`. Do not convert unknown strand (`0`) and unspecified strand (`None`) into
the same value.

### Compatibility

- Existing `validate --json` remains as a supported alias for
  `validate --format json`.
- Existing JSON top-level types remain unchanged in v0.8.5.
- Existing stdout defaults remain text unless the current command already
  defaults to another format.

## 8. Workstream F05 — validation as a pipeline gate

### Required CLI

```bash
gbparse validate INPUT --format json --output validation.json --fail-on error
```

Threshold semantics:

| `--fail-on` | Exit 1 when findings include |
|---|---|
| `never` | never; parse/output failures still use codes 3/4 |
| `error` | `ERROR` |
| `warning` | `ERROR` or `WARNING` |
| `info` | any finding |

Default to `never` for backward compatibility.

The report must always be fully written before returning exit 1. A validation
failure is a completed analysis result, not a serialization rollback.

Add summary counts to JSON without removing the existing findings array only
if this can be done compatibly. Otherwise keep `--format json` as the existing
array in v0.8.5 and offer summary counts in text/TSV; do not break consumers
merely to add an envelope.

### Additional safe validation checks

Within existing GenBank semantics, add only checks that do not require external
databases:

- feature start less than 1;
- every compound segment lies within record bounds;
- invalid or empty compound locations;
- CDS `/codon_start` outside 1–3 when explicitly present;
- duplicate CDS locus tags already handled, with deterministic ordering;
- translation length mismatch and internal stop behavior remain
  genetic-code-aware and exception-aware.

Do not add submission-policy checks that depend on organism-specific NCBI
rules; that belongs to `table2asn`.

## 9. Workstream F06 — complete annotation diffs

### Correctness fix

Remove `details[:200]` and per-category `[:50]` truncation from JSON and other
machine outputs.

```diff
- "changes": details[:200],
+ "changes": details,
+ "result_complete": True,
+ "change_count": len(details),
```

Human text may use `--max-display`, but must print both displayed and total
counts. `--summary-only` omits detail records by explicit request and records
that fact in structured metadata.

### TSV and JSONL

Add one row/object per classified feature pair with stable fields:

```text
record old_locus new_locus old_start old_end new_start new_end
old_strand new_strand match_tier classifications old_gene new_gene
old_product new_product old_kos new_kos old_ecs new_ecs
```

### Optional PAF coordinate map

Support importing a caller-supplied PAF file to compare annotations when Bakta
was rerun on a renamed, rotated, or modestly changed assembly:

```bash
gbparse diff old.gbff new.gbff --coordinate-map old-to-new.paf
```

Constraints:

- parse PAF only; do not invoke minimap2;
- require unique, deterministic usable mappings;
- report ambiguous/unmapped records instead of guessing;
- map strand and half-open PAF coordinates carefully into one-based inclusive
  GenBank feature coordinates;
- preserve both original and mapped coordinates in output;
- exact sequence-alignment identity remains external evidence;
- fail closed when a feature crosses a split mapping or an unsupported
  structural rearrangement.

If robust PAF mapping threatens the release schedule, defer only this optional
subfeature in a documented issue; full untruncated machine output is mandatory.

## 10. Workstream F07 — evidence-preserving `compare`

### Current gap

The current implementation gathers matching locus tags but emits only aggregate
counts. It also tests one target string simultaneously as gene, product token,
KO, EC, and COG, which is convenient but semantically ambiguous.

### Required outputs

Support:

```text
--mode count      integer copy counts
--mode presence   0/1 matrix
--mode status     ABSENT/SINGLE_COPY/MULTI_COPY
```

Add `--evidence-output` with one row per match:

```text
genome path target target_kind record feature_index locus_tag gene product
start end strand matched_field matched_value match_mode
```

The matrix remains concise; the evidence table makes every cell auditable.

### Typed target file

Retain comma-separated `--targets` for compatibility. Add an optional YAML/JSON
target definition:

```yaml
schema_version: 1
targets:
  - id: ladA_candidate
    any:
      - field: gene
        mode: casefold_exact
        value: ladA
      - field: ko
        mode: exact
        value: K20938
```

Suggested flag:

```text
--targets-file targets.yaml
```

Require exactly one of `--targets` and `--targets-file`. Reuse matching
primitives where possible, but do not couple this compact target schema to the
MEOR or mobilome catalogs.

### Tests

- matrix counts equal evidence-row counts;
- one feature matching two alternatives within one target counts once;
- one feature may legitimately match two different target IDs;
- gene exact match, product token boundaries, KO/EC/COG exactness;
- duplicate genome basenames cannot silently collapse rows—derive or require a
  unique sample ID;
- deterministic recursive file discovery;
- inaccessible or malformed files produce explicit per-input failure behavior.

## 11. Workstream F08 — safe declarative `query`

### Grammar

Implement a small parser; never use Python `eval`, `ast.literal_eval` as an
expression evaluator, SQL interpolation, or shell evaluation.

Minimum supported fields:

```text
record, record_index, feature_index, type, locus_tag, gene, product,
protein_id, start, end, strand, length, pseudo, partial,
ko, ec, cog, pfam, rfam, go, db_xref
```

Arbitrary qualifiers use an explicit accessor:

```text
qualifier("inference")
qualifier("note")
```

Minimum operators:

```text
==  !=  <  <=  >  >=  contains  ~=  in
and  or  not  ( ... )
```

Semantics:

- string equality is exact unless a field documents case-folded behavior;
- recommend explicit `casefold()` or a `casefold_exact` function rather than
  hidden case behavior;
- `~=` is regular-expression search with clear invalid-regex errors;
- scalar comparison against a multivalued field succeeds if any value matches;
- `not` operates on the complete predicate result;
- missing fields compare unequal except for an explicit null test;
- numeric comparisons reject strings and booleans;
- resource guards cap expression length, nesting depth, regex length, and
  compiled-regex count.

Example:

```bash
gbparse query isolate.gbff \
  --where 'type == "CDS" and (gene == "ladA" or ko == "K20938") and not pseudo' \
  --select record,locus_tag,gene,product,start,end,strand,ko \
  --format tsv
```

### Result model

Define one canonical `FeatureRow` or equivalent projection model shared by
`query` and `export`. It must retain:

- input/source label;
- record identity and index;
- feature index and type;
- all physical segments;
- start/end and strand state;
- partial/pseudo status;
- selected qualifiers;
- typed xrefs with their source field where available.

### Formats and emission

```text
--emit rows  + --format text|tsv|csv|json|jsonl
--emit faa
--emit ffn
```

- FAA uses annotated `/translation`; missing translations are either skipped
  with a reported count or rejected under `--strict`.
- FFN uses biological `SeqFeature.extract()` and requires sequence content.
- FASTA headers are deterministic and include unambiguous record/locus or
  feature-index identity.
- Do not implement queried GenBank window emission in this release; users can
  pass selected loci to the existing `region` command. This avoids undefined
  behavior for overlapping query windows.

### Query parser tests

- precedence and parentheses;
- repeated `not`;
- multivalued qualifiers;
- regex Unicode/case behavior;
- malformed tokens and unterminated strings;
- nesting/resource limits;
- no access to attributes, imports, callables, filesystem, or environment;
- compound/origin-spanning features;
- identical query results for plain and gzipped input;
- stable JSONL byte output.

## 12. Workstream F09 — unified `export`

### General behavior

`export` is a new coherent interface. Existing commands become thin aliases or
retain their current internals until safely migrated; do not remove them.

### 12.1 Annotation TSV

Equivalent to `extract`, but built from the canonical feature-row serializer.
Document the fixed columns and missing-value policy.

### 12.2 JSONL

One feature per line using schema `gbparse.feature.v1`:

```json
{"schema_version":"gbparse.feature.v1","record_id":"contig_1","feature_index":7}
```

The actual schema must include locations, segments, qualifiers, typed xrefs,
and biological flags. Validate representative JSONL objects in tests.

### 12.3 FASTA

Support:

- `faa`: annotated CDS translations;
- `ffn`: biological CDS nucleotide sequences;
- `fna`: record nucleotide sequences.

Reuse existing FASTA wrapping and extraction behavior. Define behavior for
missing sequence and missing translation explicitly.

### 12.4 GFF3

Delegate to the current GFF3 converter. `--include-fasta` applies only here.
Do not fork a second GFF implementation.

### 12.5 BED12

Coordinate conversion:

- GenBank: one-based inclusive;
- BED: zero-based half-open;
- `chromStart = start - 1`;
- `chromEnd = end`;
- compound parts become BED blocks;
- block starts are relative to `chromStart`;
- unknown strand becomes `.`; strand state `0` must not be mislabeled `+`;
- origin-spanning circular features require either a documented split-record
  representation or a fail-closed error; do not emit invalid negative or
  wrapped BED coordinates.

Choose and test one policy. Recommended v0.8.5 policy: split an
origin-spanning feature into linked BED records with a stable shared ID and a
warning field, because standard BED intervals cannot wrap.

### 12.6 NCBI table2asn input package

Generate:

```text
PREFIX.fsa
PREFIX.tbl
PREFIX.export.json
```

The `.fsa` and five-column `.tbl` files must have matching sequence IDs.
Preserve supported INSDC features and qualifiers conservatively. Unknown or
unsupported qualifiers must be reported in `PREFIX.export.json`; they must not
be silently dropped.

The JSON export report records:

- schema version;
- source path/label and SHA-256 when path-backed;
- gbparse version;
- record IDs;
- feature and qualifier counts;
- emitted, skipped, and unsupported items;
- warnings that `table2asn` validation is still required.

Do not generate `.sbt` author metadata, `.sqn`, or make submission-readiness
claims. Do not call external executables.

### Tests

- annotation TSV parity with existing `extract` for shared columns;
- FAA/FFN/FNA parity with existing exporters;
- GFF3 byte parity through delegation;
- JSONL schema and deterministic ordering;
- BED coordinate boundary and compound-feature fixtures;
- circular-origin BED policy;
- NCBI `.fsa`/`.tbl` ID agreement;
- multiline/join/complement feature-table cases;
- unsupported qualifier reporting;
- atomic multi-file publication rollback.

## 13. Workstream F10 — standalone `operons`

### Reuse rather than rewrite

The repository already contains `find_operon_pairs()`, `build_operon_result()`,
`OperonPair`, and `OperonCluster`. Promote this into the unified CLI and refactor
only enough to obtain renderer-independent, record-level structured output.

### Semantics

- consecutive CDS features only;
- same known strand (`+1` or `-1`);
- configured inclusive `min_gap <= gap <= max_gap`;
- record-local;
- circular last-to-first adjacency only for circular records;
- a cluster contains at least `--min-genes` CDSs, default 3;
- output wording always says `candidate` or `proximity candidate`.

### Required formats

- text: current report, cleaned and deterministic;
- TSV: one `pair` or `cluster_member` row per entity;
- JSON: schema-versioned full result;
- GFF3: candidate operon parents plus CDS references/children without rewriting
  the source annotation.

Suggested JSON schema identifier:

```text
gbparse.operons.v1
```

### Edge tests

- exact minimum and maximum gap;
- permitted overlaps;
- opposite and unknown strand rejection;
- two genes form a pair but not a three-gene cluster;
- plus and minus strand biological ordering;
- circular origin closure;
- linear record does not close last-to-first;
- compound CDS coordinates;
- duplicate locus tags do not collapse physical features.

## 14. Workstream F11 — schemas, documentation, and API boundaries

### Schemas

Add JSON Schemas for new complex contracts where practical:

```text
src/genbank_parser/data/schemas/feature-v1.schema.json
src/genbank_parser/data/schemas/operons-v1.schema.json
src/genbank_parser/data/schemas/query-v1.schema.json
```

Package them explicitly in `pyproject.toml` and verify them from the built
wheel. Do not relocate the existing mobilome schema in this patch.

### Documentation

Update:

- README CLI inventory and examples;
- changelog v0.8.5 section;
- structured-output contract;
- coordinate-system reference for GenBank, GFF3, BED, and PAF;
- query grammar and escaping examples for Bash, PowerShell, and Windows CMD;
- NCBI export limitations;
- operon candidate interpretation;
- migration notes for corrected `functional --format tsv` behavior;
- a reference to the companion `gbparse-v0.8.5-deferred-scope.md` document.

### Python API

Export stable builders/result classes deliberately from package modules. Do
not expose parser-internal query AST node classes as public API.

Retain compatibility wrappers that print only in CLI-specific functions. New
analysis functions should return typed immutable result objects where feasible.

## 15. Workstream F12 — testing and release gates

### Unit tests

Add focused suites:

```text
tests/test_cli_io.py
tests/test_query.py
tests/test_export.py
tests/test_operons_cli.py
tests/test_structured_outputs.py
tests/test_diff_complete.py
tests/test_compare_evidence.py
```

### CLI integration matrix

Test representative commands through `main([...])`, not only internal
functions:

| Scenario | Required assertion |
|---|---|
| stdout text | no diagnostic contamination |
| stdout TSV/JSON/JSONL | parseable and complete |
| file output | byte-equivalent to stdout data |
| existing path | rejected without `--force` |
| input/output collision | always rejected |
| parse error | expected nonzero code and no final output |
| validation threshold | report written, exit 1 |
| gzip/stdin | parity with ordinary input |

### Cross-version CI

Retain Python 3.10–3.13. Add wheel assertions for:

- `__version__ == "0.8.5"`;
- `gbparse --version`;
- packaged JSON schemas;
- one gzip parse;
- one `query --format jsonl` execution;
- one `export --format bed12` execution;
- one `operons --format json` execution;
- unchanged MEOR and mobilome resource versions unless independently bumped;
- optional visualization extra still imports and renders.

### Golden files

Use goldens only for stable public formats. Normalize line endings to LF.
Review golden changes manually; never regenerate and accept them blindly.

Suggested goldens:

```text
tests/golden/query_v1.jsonl
tests/golden/operons_v1.json
tests/golden/functional_v1.tsv
tests/golden/feature_v1.jsonl
tests/golden/export_v1.bed
```

### Property/invariant tests

- every query result points to a real source feature;
- compare matrix counts equal grouped evidence rows;
- all emitted BED intervals satisfy `0 <= start < end <= record_length`;
- all NCBI feature-table sequence IDs exist in the paired FASTA;
- JSONL object count equals emitted feature count;
- no machine serializer truncates results;
- plain/gzip/stdin parsing yields equivalent feature identities;
- atomic-output failures preserve pre-existing outputs.

## 16. Implementation order

Implement in this dependency order:

1. **F01–F03:** shared I/O, renderers, errors, output safety, version flag;
2. **F04–F05:** normalize existing outputs and validation gating;
3. **F06–F07:** complete `diff` and evidence-preserving `compare`;
4. **F08:** declarative query engine;
5. **F09:** unified exporters;
6. **F10:** standalone operons CLI;
7. **F11:** schemas and documentation;
8. **F12:** wheel/release gates and final regression pass.

Do not implement every workstream in one giant commit. Recommended commit
boundaries:

```text
refactor(cli): centralize input and atomic output handling
fix(cli): normalize structured serializers and validation exits
fix(diff): retain complete machine-readable change sets
feat(compare): emit typed marker evidence
feat(query): add safe declarative feature queries
feat(export): add unified interoperable exporters
feat(operons): expose structured candidate-operon CLI
docs(release): document v0.8.5 contracts and link deferred scope
chore(release): prepare genbank-parser 0.8.5
```

## 17. Release acceptance checklist

v0.8.5 is ready only when all items below are true:

- [ ] All old subcommands still resolve and their documented core use cases run.
- [ ] `gbparse --version` reports `0.8.5` from an installed wheel.
- [ ] Plain, gzip, and stdin fixtures produce equivalent parsed features.
- [ ] `functional --format tsv` is valid TSV rather than decorated text.
- [ ] `summary`, `locus`, `crispr`, and `phylo` have parseable structured output.
- [ ] `discover` and normalized commands can write safely to explicit outputs.
- [ ] `validate --fail-on` returns the documented code after writing its report.
- [ ] Machine-readable `diff` output contains every change unless
      `--summary-only` is explicitly requested.
- [ ] `compare` aggregate values reconcile exactly with its evidence table.
- [ ] Query expressions are parsed without `eval` and resource limits are tested.
- [ ] Query FAA/FFN emission preserves biological extraction semantics.
- [ ] JSONL, BED12, and NCBI table exporters pass coordinate/ID invariants.
- [ ] `operons` is circular-aware and labels every result as a candidate.
- [ ] Output collisions and partial publication are tested fail-closed.
- [ ] MEOR, mobilome, neighborhood, GFF3, and region regression tests remain green.
- [ ] Python 3.10, 3.11, 3.12, and 3.13 CI jobs pass.
- [ ] Source distribution and isolated wheel tests pass.
- [ ] README, changelog, and help snapshots match the final CLI.
- [ ] Features listed in `gbparse-v0.8.5-deferred-scope.md` are absent from the
      v0.8.5 command surface.

## 18. Definition of done

The release is complete when `gbparse` can be used safely in both interactive
and scripted settings without scraping decorated terminal output, losing
records to hidden limits, or writing partial artifacts; when users can express
nontrivial feature queries without writing Python; when common downstream
representations can be exported with explicit coordinate semantics; and when
the existing operon logic is available as a cautious, structured CLI analysis.

The guiding boundary is simple: v0.8.5 should make existing GenBank evidence
more searchable, portable, and reproducible—not attempt to perform every
downstream biological analysis itself.
