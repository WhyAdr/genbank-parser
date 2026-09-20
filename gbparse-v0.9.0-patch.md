# `genbank-parser` v0.9.0 Implementation Plan

**Target repository:** `WhyAdr/genbank-parser`  
**Audited baseline:** `e1138e24161d2936649feda4b362f556d7b564bc` (`main`, package v0.8.5)  
**Target release:** v0.9.0  
**Primary implementer:** Luna  
**Scope:** persistent cohort indexing, a universal single-input batch runner, and record-level listing/extraction/splitting/filtering  

## 1. Objective

v0.9.0 should turn `gbparse` from a strong single-file toolkit into a safe
cohort utility without weakening the biological and output guarantees established
in v0.8.5.

The release adds three connected capabilities:

1. `gbparse records` for inspecting and selecting whole GenBank records;
2. `gbparse index` for building and querying a persistent, normalized cohort
   index; and
3. `gbparse batch` for applying an existing single-input command across a
   deterministically discovered collection with resumability and provenance.

These features must reuse the canonical parser, feature projection, query
grammar, error taxonomy, and atomic publication helpers. They must not grow a
second parser, a second biological model, or an unbounded arbitrary-code/SQL
interface.

## 2. Baseline assessment

### 2.1 Verified state

At the audited v0.8.5 release commit:

- package metadata and `gbparse --version` report `0.8.5`;
- the full test suite completes with **228 passed, 8 skipped**;
- `read_genbank()` is the canonical parser for plain, gzip, stdin, and stream
  inputs;
- `FeatureRow` is the canonical schema-versioned feature projection;
- the bounded `query` grammar already handles scalar and multivalue fields;
- `write_text()` and `publish_directory()` provide atomic file and directory
  publication; and
- exit codes distinguish validation findings, usage errors, input failures,
  and output failures.

This is a suitable base. v0.9.0 should extend these contracts instead of
replacing them.

### 2.2 Gaps that affect this release

The current `batch-summary` implementation predates the v0.8.5 safety layer:

- discovery recognizes only uncompressed `.gb`, `.gbk`, and `.gbff` suffixes;
- output files are opened and overwritten directly;
- progress is printed to stdout;
- per-file exceptions are broadly caught and converted to warnings;
- the unified CLI returns success even when no file was found or every file
  failed; and
- its three outputs are not published transactionally as one artifact set.

Do not copy those behaviors into `gbparse batch` or `gbparse index`.

The typed `GenBankRecord` also lacks a retained source `SeqRecord`. Rebuilding a
whole GenBank record from the reduced typed fields risks losing record-level
references, database cross-references, keywords, taxonomy metadata, fuzzy
positions, and other Biopython-supported content. Whole-record operations must
therefore retain and serialize the original parsed `SeqRecord`.

Finally, `_dispatch()` is one large CLI function. The batch runner must not
duplicate every command's biological implementation. Introduce a small adapter
layer around existing command invocation and publication behavior.

## 3. Release invariants

### 3.1 Canonical biological semantics

- Parse GenBank only through `read_genbank()`/`iter_genbank()`.
- Preserve `FeatureLocation`, `CompoundLocation`, fuzzy positions, qualifiers,
  feature order, record order, topology, and biological strand states.
- Whole-record extraction is selection, not annotation correction.
- Index rows are projections of source annotations, not new biological calls.
- Batch execution must not alter the result of the wrapped command.

### 3.2 Safety and determinism

- Never overwrite an input, index, run directory, or output unless the relevant
  `--force` operation explicitly permits replacement.
- Build replacement indexes in a sibling temporary file and install them with
  `os.replace()` only after validation.
- Publish each batch job atomically. A failed job must not leave a plausible
  final artifact.
- Sort discovered files and emitted rows deterministically.
- Do not silently cap records, features, files, or query results.
- Send data to stdout only when stdout is the declared output. Progress,
  warnings, and summaries belong on stderr.

### 3.3 Bounded interfaces

- No `eval`, `exec`, shell command string, or arbitrary SQL execution.
- Batch command arguments remain an argv array and are validated against the
  selected command parser before any job starts.
- Index queries compile the existing bounded query AST into parameterized SQL.
- User-derived filenames are sanitized and collision-resolved deterministically.
- Regex limits must match the existing query engine's bounded behavior.

### 3.4 Compatibility

- Keep `batch-summary` available with its current name and default filenames.
- Keep every v0.8.5 command and output schema compatible unless this plan
  explicitly introduces a new versioned schema.
- SQLite is from the Python standard library; do not add DuckDB, Pandas,
  PyArrow, or an ORM as a core dependency.
- v0.9.0 indexes are schema-versioned artifacts and must fail clearly when
  opened by incompatible code.

## 4. Proposed command surface

### 4.1 `gbparse records`

```bash
# Inventory all records
gbparse records list assembly.gbff --format tsv

# Extract exact records into one GenBank file
gbparse records extract assembly.gbff \
  --record chromosome --record plasmid_p1 \
  --output selected.gbff

# Filter records using explicit record-level predicates
gbparse records filter assembly.gbff \
  --min-length 10000 --topology circular --has-feature CDS \
  --output circular_records.gbff

# Split selected records and emit an ID-to-file manifest
gbparse records split assembly.gbff \
  --output-dir split-records/ --format genbank
```

`records` has nested actions: `list`, `extract`, `filter`, and `split`. Do not
make an omitted action guess between read-only listing and file creation.

### 4.2 `gbparse index`

```bash
# Build a new cohort index
gbparse index build ./bakta-results cohort.gbidx \
  --jobs 4 --report cohort-index-report.json

# Add new files and replace changed sources; absent sources remain by default
gbparse index update cohort.gbidx ./new-results --jobs 4

# Remove indexed sources not present in the supplied discovery set
gbparse index update cohort.gbidx ./current-results --prune

# Query indexed feature annotations with the existing safe grammar
gbparse index query cohort.gbidx \
  --where 'type == "CDS" and (gene == "ladA" or ko == "K20938")' \
  --select source,record,locus_tag,gene,product,ko \
  --format tsv --output lada-candidates.tsv

# Inspect index metadata and cohort counts
gbparse index status cohort.gbidx --format json
```

The recommended extension is `.gbidx`, but the artifact is an ordinary SQLite
database. Do not depend on the suffix when validating an existing database.

### 4.3 `gbparse batch`

```bash
# Apply a single-input command to every discovered GenBank file
gbparse batch ./isolates \
  --command validate --output-dir validate-run --jobs 4 \
  -- --format json --fail-on warning

# Resume only jobs whose input or outputs no longer match the manifest
gbparse batch ./isolates \
  --command validate --output-dir validate-run --jobs 4 --resume \
  -- --format json --fail-on warning

# Batch a query; the runner inserts each input and owns each output path
gbparse batch ./isolates \
  --command query --output-dir ladA-query \
  -- --where 'gene == "ladA"' --format jsonl
```

Everything after `--` is a literal argument vector for the wrapped command.
The runner inserts the per-job input and output destinations. It must reject
runner-owned options in the tail, including `--output`, `--output-dir`,
`--force`, and positional attempts to replace the input.

## 5. Workstream F01 — shared cohort input discovery and identity

### Required design

Move generic discovery out of `bakta.py` into `discovery.py`.

Support:

- explicit `.gb`, `.gbk`, `.gbff`, `.gb.gz`, `.gbk.gz`, and `.gbff.gz` files;
- deterministic recursive directory discovery;
- overlapping input roots without duplicate physical paths;
- symlink-safe canonical identity while preserving the user-facing path;
- exclusion of the selected output directory/index and their staging siblings;
- stable ordering by normalized display path, not only basename; and
- raw-byte SHA-256, byte size, and modification metadata calculated by a shared
  fingerprint helper.

Do not recursively try to parse every arbitrary file. Explicit files with an
unusual suffix may still be accepted because `read_genbank()` detects gzip by
magic bytes; directory discovery should remain bounded to recognized suffixes.

Define immutable types such as:

```python
@dataclass(frozen=True)
class DiscoveredInput:
    requested_path: str
    resolved_path: Path
    display_path: str
    discovery_root: Path | None

@dataclass(frozen=True)
class FileFingerprint:
    sha256: str
    size_bytes: int
    mtime_ns: int
```

Sample/output keys must be readable and collision-safe:

1. remove the full recognized GenBank/gzip suffix;
2. sanitize to a portable basename;
3. use the readable stem when unique; and
4. append `-<digest8>` when two discovered inputs would receive the same key.

Never merge two files merely because their stems match. Exact content
duplicates remain separate sources unless the user explicitly asks for future
deduplication functionality.

### Rough pseudocode diff

```diff
+ src/genbank_parser/discovery.py
+   GENBANK_SUFFIXES = (...)
+   DiscoveredInput
+   FileFingerprint
+   discover_inputs(inputs, *, recursive=True, exclude=())
+   fingerprint_file(path, *, chunk_size=1 << 20)
+   assign_sample_keys(discovered)

  src/genbank_parser/bakta.py
- def discover_genbank_files(...): ...
+ from .discovery import discover_inputs
```

### Acceptance criteria

- Plain and recognized gzip files are discovered in the same deterministic
  order on Linux and Windows.
- Repeated files and overlapping directories produce one job per resolved path.
- Two `sample.gbff` files in different directories never collide.
- An output tree nested under an input root cannot be rediscovered during the
  same run.

## 6. Workstream F02 — lossless whole-record model and serialization

### Required model change

Retain the original Biopython `SeqRecord` on `GenBankRecord`:

```diff
  @dataclass
  class GenBankRecord:
      ...
+     raw_record: SeqRecord | None = field(default=None, repr=False, compare=False)
```

`iter_genbank()` must attach the source `SeqRecord`. This is not a second parser;
it is the canonical parser retaining the source object it already created.

Add a serializer that accepts selected `GenBankRecord` objects and emits
GenBank or FASTA without mutating the source records. For GenBank output, write
deep copies of `raw_record` so serialization helpers cannot alter the parsed
document.

If `raw_record` is absent for a programmatically constructed record, fail with
a precise error unless a separately tested lossless reconstruction path exists.
Do not quietly emit an impoverished GenBank record.

Support gzip output when the destination ends in `.gz`. Implement an atomic
binary writer or an atomic compressed-text wrapper; do not write compressed
output directly to its final path.

### Record projection

Add `RecordRow` with schema `gbparse.record.v1`. At minimum include:

- source;
- record index, ID, name, description;
- length, topology, molecule type, division, and date;
- accession(s), organism, taxonomy, strain/isolate;
- GC percentage over A/C/G/T bases;
- feature count and CDS/rRNA/tRNA/tmRNA/ncRNA/pseudogene counts; and
- whether nucleotide sequence is present.

Use JSON arrays for multivalue fields. TSV/CSV may join those fields with `;`
using the same empty-value policy as other serializers.

### Rough pseudocode diff

```diff
  src/genbank_parser/io.py
      yield GenBankRecord(
          ...,
+         raw_record=rec,
      )

+ src/genbank_parser/records.py
+   RecordRow
+   record_rows(document)
+   select_records(document, selector)
+   render_record_rows(rows, format_type)
+   render_selected_records(records, format_type)
+   split_record_artifacts(records, format_type)
```

### Fidelity tests

Round-trip a fixture containing:

- multiple records;
- record references and `dbxrefs`;
- taxonomy and structured comments;
- a circular record;
- compound and fuzzy feature locations;
- repeated qualifiers; and
- non-CDS features.

Parse the emitted GenBank again and compare semantic record/feature content,
not raw text formatting. Biopython is allowed to normalize whitespace and LOCUS
formatting.

## 7. Workstream F03 — `records` actions and selection semantics

### `records list`

```text
gbparse records list INPUT
  [--format text|tsv|csv|json|jsonl]
  [--output PATH] [--force]
  [selection predicates]
```

Default format is `text`; machine formats use `gbparse.record.v1`. An empty
selection is a successful empty listing with headers/valid JSON.

### `records extract`

```text
gbparse records extract INPUT
  --record ID [--record ID ...]
  [--format genbank|fasta]
  --output PATH [--force]
```

Rules:

- match exact record ID first, then exact name;
- if a selector matches multiple records, fail and report their record indices;
- if any requested selector is absent, fail before writing;
- preserve source order, even if selectors were supplied in another order;
- do not emit the same record twice; and
- require at least one selected record.

### `records filter`

```text
gbparse records filter INPUT
  [--record ID ...]
  [--record-regex REGEX]
  [--min-length N] [--max-length N]
  [--topology linear|circular|unknown]
  [--molecule-type TEXT]
  [--has-feature TYPE ...]
  [--has-locus LOCUS_TAG ...]
  [--invert]
  [--format genbank|fasta]
  --output PATH [--force]
```

Different predicate families combine with AND. Repeated values within one
family combine with OR, except repeated `--has-feature` and `--has-locus`, which
mean all requested values must be present. Apply `--invert` once to the final
predicate.

Require at least one predicate and at least one matching output record. This
prevents a typo from silently creating an empty or complete copy.

Bound record regex length and reject invalid expressions before parsing the
input. Do not reuse the feature query evaluator by pretending records are
features.

### `records split`

```text
gbparse records split INPUT
  [selection predicates]
  --output-dir DIR
  [--format genbank|fasta]
  [--force]
```

Publish the directory transactionally with `publish_directory()`. Include
`records.manifest.tsv` containing:

```text
schema_version  source  record_index  record_id  record_name  length  filename  sha256
```

Filenames are sanitized record IDs plus `.gbff` or `.fna`; collisions receive a
stable `-r<record_index>` suffix. The manifest, not a filename heuristic, is
the authoritative ID-to-file mapping.

### Rough CLI diff

```diff
  def create_parser():
+     p_records = subparsers.add_parser("records", ...)
+     record_actions = p_records.add_subparsers(dest="records_action", required=True)
+     add_records_list_parser(record_actions)
+     add_records_extract_parser(record_actions)
+     add_records_filter_parser(record_actions)
+     add_records_split_parser(record_actions)

  def _dispatch(args, parser):
+     if cmd == "records":
+         return dispatch_records(args)
```

## 8. Workstream F04 — SQLite cohort index schema

### Storage decision

Use normalized SQLite through the standard-library `sqlite3` module. Do not use
an ORM and do not add FTS in schema v1. Ordinary B-tree indexes cover exact
identifier queries, coordinates, types, and common equality predicates. FTS
would introduce tokenization and substring semantics that do not match the
existing `contains`/regex query language cleanly.

The index is an annotation/search cache, not an archive. Do not store full
genome sequence or serialized GenBank blobs. Store source fingerprints so
sequence-bearing follow-up commands can verify and reopen the original file.

### Schema v1

Use `PRAGMA user_version = 1`, `PRAGMA foreign_keys = ON`, and explicit schema
metadata. Recommended tables:

```sql
metadata(
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
)

sources(
  source_pk INTEGER PRIMARY KEY,
  display_path TEXT NOT NULL UNIQUE,
  resolved_path_at_index TEXT NOT NULL,
  sha256 TEXT NOT NULL,
  size_bytes INTEGER NOT NULL,
  mtime_ns INTEGER NOT NULL,
  compressed INTEGER NOT NULL CHECK (compressed IN (0,1)),
  record_count INTEGER NOT NULL,
  feature_count INTEGER NOT NULL
)

records(
  record_pk INTEGER PRIMARY KEY,
  source_pk INTEGER NOT NULL REFERENCES sources ON DELETE CASCADE,
  record_index INTEGER NOT NULL,
  record_id TEXT NOT NULL,
  record_name TEXT NOT NULL,
  description TEXT NOT NULL,
  length INTEGER NOT NULL,
  topology TEXT,
  molecule_type TEXT,
  organism TEXT,
  strain TEXT,
  gc_percent REAL,
  UNIQUE(source_pk, record_index)
)

features(
  feature_pk INTEGER PRIMARY KEY,
  record_pk INTEGER NOT NULL REFERENCES records ON DELETE CASCADE,
  feature_index INTEGER NOT NULL,
  type TEXT NOT NULL,
  locus_tag TEXT,
  gene TEXT,
  product TEXT,
  protein_id TEXT,
  start INTEGER NOT NULL,
  end INTEGER NOT NULL,
  strand INTEGER,
  biological_length INTEGER NOT NULL,
  partial INTEGER NOT NULL,
  pseudo INTEGER NOT NULL,
  UNIQUE(record_pk, feature_index)
)

segments(
  feature_pk INTEGER NOT NULL REFERENCES features ON DELETE CASCADE,
  ordinal INTEGER NOT NULL,
  start INTEGER NOT NULL,
  end INTEGER NOT NULL,
  PRIMARY KEY(feature_pk, ordinal)
)

qualifiers(
  feature_pk INTEGER NOT NULL REFERENCES features ON DELETE CASCADE,
  key TEXT NOT NULL,
  ordinal INTEGER NOT NULL,
  value TEXT NOT NULL,
  PRIMARY KEY(feature_pk, key, ordinal)
)

xrefs(
  feature_pk INTEGER NOT NULL REFERENCES features ON DELETE CASCADE,
  namespace TEXT NOT NULL,
  value TEXT NOT NULL,
  source_field TEXT NOT NULL,
  PRIMARY KEY(feature_pk, namespace, value, source_field)
)
```

Add indexes for:

- `sources(sha256)`;
- `records(record_id)`, `records(record_name)`, and `records(organism)`;
- `features(type)`, `features(locus_tag)`, `features(gene)`,
  `features(protein_id)`, and `features(record_pk, start, end)`;
- `xrefs(namespace, value)`; and
- `qualifiers(key, value)`.

Preserve qualifier ordinals and segment ordinals. Do not collapse multivalue
qualifiers into delimiter-separated strings.

Metadata must include at least:

- `schema_version = gbparse.index.v1`;
- `schema_revision = 1`;
- creating/updating `gbparse_version`;
- Python and Biopython versions;
- feature and record projection schema versions;
- query semantics version; and
- creation/update timestamps in UTC.

Timestamps are provenance, so byte-identical databases are not promised.
Logical row content and query order must still be deterministic.

## 9. Workstream F05 — index build, update, status, and verification

### Build

`index build INPUT... DATABASE` must:

1. discover and fingerprint every source deterministically;
2. reject an existing destination unless `--force`;
3. create a sibling temporary database;
4. parse sources through `read_genbank()`;
5. insert normalized rows inside transactions;
6. run `PRAGMA foreign_key_check` and `PRAGMA integrity_check`;
7. switch to a single-file journal state and close every connection;
8. write an optional schema-versioned report; and
9. atomically install the validated database.

Default `--on-error fail` means one bad input aborts publication. Optional
`--on-error skip` may publish the successful cohort only if the report and
stderr enumerate every skipped source and the command exits nonzero. Never
label such a build fully successful.

Parallel parsing may use worker processes, but only the parent process writes
SQLite. Consume worker results in discovery order, not completion order.

### Update

Updates must be crash-safe:

1. validate schema compatibility;
2. copy the current database to a sibling temporary database using SQLite's
   backup API, not a raw copy of a live WAL database;
3. leave unchanged sources untouched when their SHA-256 matches;
4. delete and reinsert a source transactionally when content changed;
5. add new sources;
6. retain absent sources unless `--prune` was explicitly supplied;
7. validate the temporary database; and
8. atomically replace the original.

`--prune` compares the supplied discovery set by stored display-path identity.
Report removed sources before replacement. It must not infer deletion merely
because one previously used input root was omitted accidentally; require the
update invocation to include `--prune` and at least one directory root.

### Status

`index status` reports:

- schema and tool versions;
- source, record, and feature counts;
- counts by feature type and topology;
- database size;
- duplicate source-content hashes;
- stored source paths that are currently absent or whose raw-byte digest has
  changed, only when `--verify-sources` is requested; and
- `quick_check`/foreign-key status when `--check` is requested.

Plain `status` must not hash every source unexpectedly.

### Rough pseudocode diff

```diff
+ src/genbank_parser/index/__init__.py
+ src/genbank_parser/index/schema.py
+ src/genbank_parser/index/build.py
+ src/genbank_parser/index/query.py
+ src/genbank_parser/index/report.py
+ src/genbank_parser/data/schemas/index-report-v1.schema.json

+ def build_index(inputs, destination, *, jobs, on_error): ...
+ def update_index(database, inputs, *, prune, jobs, on_error): ...
+ def inspect_index(database, *, check, verify_sources): ...
```

## 10. Workstream F06 — safe index query compilation

### Grammar reuse

Refactor the feature query parser so its AST and field registry are public
package internals shared by direct and indexed queries:

```diff
  src/genbank_parser/query.py
- class _QueryParser: ...
+ class QueryParser: ...
+ def parse_query_ast(expression) -> QueryNode: ...
+ def evaluate_query_ast(node, row) -> bool: ...
```

Do not create a second grammar. Add a compiler that maps the same AST to SQL
plus bound parameters.

Scalar fields map to columns. Multivalue cross-reference predicates compile to
`EXISTS` clauses against `xrefs`; `qualifier("name")` compiles to `EXISTS`
against `qualifiers`. Regex uses a connection-local bounded `REGEXP` function
with the existing regex length/count constraints.

All literal values remain bound parameters. Identifiers come only from a fixed
field mapping. No source expression fragment may be interpolated as SQL.

### Query parity

For every supported expression, direct:

```bash
gbparse query sample.gbff --where EXPR --format jsonl
```

and indexed:

```bash
gbparse index query sample.gbidx --where EXPR --format jsonl
```

must select the same physical features and expose the same canonical fields,
with indexed output adding cohort source identity.

`--select` accepts canonical `FeatureRow` fields plus `source` and sample key.
Formats are `text`, `tsv`, `csv`, `json`, and `jsonl`. Ordering is always
source display path, record index, feature index. `--limit` is optional and
explicit; there is no hidden default limit.

For schema v1, index query emits rows only. Do not offer `--emit faa|ffn`
unless the implementation reopens each source, verifies its digest, and uses
canonical biological extraction. Sequence emission is not required for v0.9.0.

### Compiler tests

Test:

- every scalar operator;
- `contains`, regex, `in`, boolean precedence, and nested negation;
- multivalue KO/EC/COG/Pfam/Rfam/GO fields;
- arbitrary qualifier access;
- `None`, empty strings, strand `0`, and unknown strand;
- Unicode qualifier values;
- quote/semicolon/comment-shaped injection strings;
- regex/query length and depth limits; and
- direct-versus-indexed result parity.

## 11. Workstream F07 — batch command adapter protocol

### Honest universality

The runner is universal across registered **single-input** `gbparse` commands,
not across arbitrary executables or intrinsically multi-input operations.

Create an explicit adapter registry. Each adapter declares:

```python
@dataclass(frozen=True)
class BatchCommandSpec:
    name: str
    input_arity: Literal[1]
    output_mode: Literal["file", "directory", "multi_file"]
    default_suffix: str | None
    batchable: bool
    forbidden_options: tuple[str, ...]
    build_job_argv: Callable[..., list[str]]
    collect_outputs: Callable[..., tuple[ExpectedOutput, ...]]
```

The registry must cover every v0.8.5 single-input command deliberately. At
minimum, v0.9.0 must support:

- single-file outputs: `validate`, `summary`, `extract`, `search`, `locus`,
  `neighborhood` data output, `region`, `fasta`, `codon`, `functional`,
  `discover`, `phylo` report output, `crispr`, `gff`, `meor`, `mobilome`,
  `query`, `operons`, and single-file `export` formats;
- multi-artifact outputs through dedicated adapters: `sequence`, `phylo
  --output-dir`, neighborhood data plus visualization, and `export
  --format ncbi-table`.

Explicitly reject:

- `compare` and `diff`, because their input cardinality is not one;
- `batch-summary`, because it already aggregates a cohort;
- `batch` and `index`, to prevent recursion;
- `records split`, because it expands one source into its own managed directory;
  and
- commands/options whose output topology has no registered adapter.

The error should name the unsupported command/option and explain the nearest
supported invocation. Never fall back to capturing unknown stdout and guessing
what it contains.

### Invocation strategy

Run each job in an isolated child Python process using an argv array, for
example:

```python
[sys.executable, "-m", "genbank_parser.cli", command, input_path, *tail]
```

This isolates stdout/stderr and optional plotting backends and allows safe
parallelism. Never set `shell=True`.

Pre-validate the target command once using its parser with a sentinel input and
staging output. A usage error aborts before any jobs start.

Provide a parser-construction function for a single command so validation does
not copy command definitions:

```diff
  src/genbank_parser/cli.py
+ def create_command_parser(command: str) -> argparse.ArgumentParser: ...
+ def create_parser() -> argparse.ArgumentParser:
+     # compose the same command parser definitions
```

## 12. Workstream F08 — batch execution, manifest, and resume

### Run layout

Recommended directory layout:

```text
RUN_DIR/
  batch-manifest.json
  outputs/
    <sample-key>/...
  logs/
    <sample-key>.stderr.log
```

No job writes directly to its final location. Each job receives a private
staging directory under the run's sibling work area. Publish the job output
only after the child exits and every expected artifact exists.

For a new run, assemble the complete run in a sibling staging directory and
publish it with `publish_directory()` semantics. A run with failed jobs may
still be published so its manifest and diagnostics survive, but:

- failed jobs have no final output artifacts;
- the manifest says `failed`;
- stderr clearly reports the failure count; and
- the batch process exits nonzero.

There must be no state in which a partial artifact is labeled successful.

### Manifest schema

Add `gbparse.batch.v1` with:

- schema version and gbparse/Python/Biopython/platform versions;
- UTC creation/update timestamps;
- normalized command and exact command-argument array;
- runner settings (`jobs`, error policy, discovery roots);
- one job per source with stable job ID and sample key;
- input requested/display/resolved paths, SHA-256, size, and mtime;
- status: `pending`, `running`, `succeeded`, `threshold_failed`, `failed`, or
  `skipped_unchanged`;
- child exit code and gbparse exit-class meaning;
- output relative paths, byte sizes, and SHA-256 hashes;
- stderr log relative path; and
- start/end timestamps and duration.

Do not store only a reconstructed shell string. The argv array is authoritative.

Write the manifest atomically after every completed job during resume-capable
execution. Normalize any stale `running` job to pending when resuming after an
interruption.

### Resume rules

`--resume` requires an existing run directory and compatible manifest. A job is
skipped only when all are true:

1. command and normalized arguments are unchanged;
2. gbparse version matches exactly;
3. input display identity and SHA-256 match;
4. previous status was `succeeded` or an explicitly accepted validation
   threshold state; and
5. every declared output exists and its SHA-256 matches.

Otherwise rerun the job. A changed input, truncated output, edited output,
missing log, or version/argument mismatch must never be silently accepted.

`--resume` and `--force` are mutually exclusive. Resume may add newly discovered
inputs. Inputs absent from the new discovery set remain in manifest history and
are marked `not_requested`; do not delete their outputs automatically.

### Exit aggregation

Use existing gbparse meanings with deterministic precedence:

1. output/publication failure -> 4;
2. input/data/job execution failure -> 3;
3. unexpected child usage failure -> 2;
4. validation threshold reached -> 1;
5. all requested jobs successful/skipped -> 0.

An `--on-error continue` policy continues remaining jobs but never changes the
final nonzero status. `--on-error stop` halts scheduling new jobs; already
running jobs may finish and be recorded.

### Rough pseudocode diff

```diff
+ src/genbank_parser/batch.py
+   BatchCommandSpec
+   BatchJob
+   BatchManifest
+   build_batch_plan(...)
+   run_job_subprocess(...)
+   validate_resume_job(...)
+   execute_batch(...)
+   aggregate_exit_code(...)

+ src/genbank_parser/data/schemas/batch-v1.schema.json
```

## 13. Workstream F09 — harden `batch-summary` on shared infrastructure

Keep `batch-summary` as the Bakta-focused cohort table. It is not replaced by
the universal runner because it is an aggregation command rather than one
independent artifact per input.

Refactor it to:

- use shared discovery, including gzip inputs;
- render CSV, TSV, and Markdown into memory;
- publish the three files as one atomic directory when `--output-dir` is used;
- retain the old `--csv`, `--tsv`, and `--md` options for compatibility, but
  preflight all three and write them atomically one by one without partial
  overwrite;
- add `--force`;
- keep diagnostics on stderr;
- return 3 when no usable inputs or no successful summaries exist;
- return nonzero for skipped parse failures unless an explicit policy accepts
  them; and
- stop catching `Exception` broadly around biological/programming errors.

Do not route `batch-summary` through `gbparse batch`; their output models are
different.

## 14. Workstream F10 — CLI integration and code organization

The v0.9.0 additions are large enough that `cli.py` should become routing code,
not acquire three more monolithic branches.

Recommended structure:

```text
src/genbank_parser/
  batch.py
  discovery.py
  records.py
  index/
    __init__.py
    schema.py
    build.py
    query.py
    report.py
```

Expose `dispatch_records()`, `dispatch_index()`, and `dispatch_batch()` from
their modules. Parser construction may remain centrally coordinated, but each
feature owns its action-specific argument additions.

Map expected failures into the existing exit taxonomy:

- invalid CLI combination -> argparse/2;
- malformed GenBank, incompatible/corrupt index, missing source -> 3;
- collision, permission failure, unsafe replacement, publication failure -> 4;
- wrapped validation threshold -> 1; and
- success -> 0.

Do not call `parser.error()` for runtime data failures because it converts them
to a usage exit and raises `SystemExit` inside the Python API.

## 15. Workstream F11 — schemas and Python API

Package the following schemas under `data/schemas/`:

- `record-v1.schema.json`;
- `index-report-v1.schema.json`;
- `index-status-v1.schema.json`; and
- `batch-v1.schema.json`.

Validate representative and golden outputs in tests with the optional
`jsonschema` dependency.

Recommended public Python API:

```python
from genbank_parser.records import (
    RecordRow,
    RecordSelector,
    record_rows,
    select_records,
)
from genbank_parser.index import build_index, update_index, query_index
from genbank_parser.batch import build_batch_plan, execute_batch
```

Keep low-level SQLite connections, mutable manifest internals, and CLI adapter
implementation private.

## 16. Workstream F12 — test plan

### Discovery tests

- nested directories and overlapping roots;
- plain and gzip suffixes;
- unusual explicit input suffix;
- symlink duplicate and symlink output collision;
- case-colliding and identical stems;
- deterministic ordering under shuffled directory enumeration; and
- nested output directory exclusion.

### Records tests

- complete `RecordRow` projection and schema validation;
- exact ID versus name precedence;
- duplicate/ambiguous IDs;
- missing requested record fails before output;
- all filter predicate families and `--invert`;
- empty list versus forbidden empty extraction;
- source-order preservation;
- split filename collision and manifest mapping;
- gzip output;
- output aliases input rejection; and
- semantic GenBank round-trip fidelity.

### Index tests

- schema creation, indexes, user version, and foreign keys;
- exact record/feature/qualifier/xref/segment counts;
- build failure leaves no destination;
- forced failed rebuild preserves the old destination;
- unchanged update creates no duplicate rows;
- changed input replaces exactly one source;
- explicit prune and non-pruning default;
- duplicate content at distinct paths remains distinct;
- SQLite backup-based update from a WAL-state fixture;
- corrupt and future-schema indexes fail clearly;
- `integrity_check` and `foreign_key_check`;
- deterministic query ordering;
- direct versus indexed query parity; and
- injection-shaped strings remain data.

### Batch tests

- argv prevalidation before job execution;
- forbidden runner-owned passthrough options;
- one registered example for every adapter class;
- explicit rejection for every non-batchable command;
- sample-key collision handling;
- `jobs=1` and `jobs>1` produce equivalent artifacts;
- child stdout/stderr isolation;
- one failed job does not leave a final artifact;
- continue and stop error policies;
- exit-code precedence;
- manifest schema and output hashes;
- interrupted `running` job recovery;
- unchanged resume skip;
- changed input, missing output, and corrupted output rerun;
- command/version mismatch rejects resume;
- new inputs append safely; and
- no shell interpretation of paths or arguments.

### Regression gates

- all v0.8.5 tests remain green;
- direct output of every wrapped command matches its batch-produced output;
- Windows and Linux CI cover filename, replacement, SQLite close/replace, and
  multiprocessing behavior;
- Python 3.10 through the project maximum CI version are covered; and
- optional visualization tests remain conditional on `[viz]`.

## 17. Workstream F13 — documentation and migration notes

Update:

- README command list and examples;
- structured-output and exit-code contracts;
- a dedicated index schema/query section;
- a batch manifest/resume section;
- record-selection semantics;
- Python API examples;
- `SKILL.md` command table;
- changelog with v0.9.0 compatibility notes; and
- package version and golden provenance values.

Document these limitations prominently:

- the index stores annotation projections, not genome sequence;
- indexed evidence is only as current as the source fingerprint at build/update;
- index query is not arbitrary SQL;
- batch runs only registered single-input gbparse commands;
- successful outputs from a partially failed batch do not make the run
  successful; and
- record filtering does not modify or reconcile annotations.

Replace the v0.8.5 deferred document's roadmap status only if desired, but do
not rewrite its historical statement of what was excluded from v0.8.5. Link
this v0.9.0 plan as the promotion of items 1, 2, and 5.

## 18. Explicitly out of scope for v0.9.0

Keep the following out of this release:

- multi-genome synteny/alignment;
- Bakta/PGAP/Prokka annotation reconciliation;
- automatic annotation-dialect auditing or correction;
- remote accession fetching;
- arbitrary external executable orchestration;
- Nextflow/Snakemake workflow generation;
- DuckDB, Parquet, or FTS indexes;
- a network service or long-running index daemon;
- storing complete sequences or GenBank blobs inside the index;
- arbitrary SQL supplied by users;
- sequence similarity, HMM, ORF, taxonomy, ANI, AMR, VF, or BGC callers; and
- biological rewriting of selected records.

These boundaries keep v0.9.0 focused on reliable cohort mechanics rather than
turning `gbparse` into a workflow manager or annotation platform.

## 19. Implementation order

1. **F01 discovery and fingerprinting** — foundation shared by all three goals.
2. **F02 raw-record retention and serializers** — establish lossless record
   handling before exposing `records`.
3. **F03 `records` actions** — deliver and test the smallest independent feature.
4. **F04 index schema** — freeze normalized storage and metadata contracts.
5. **F05 build/update/status** — make the database lifecycle crash-safe.
6. **F06 query compiler** — prove direct/indexed semantic parity.
7. **F07 adapter registry** — enumerate batchable output topologies explicitly.
8. **F08 executor/manifest/resume** — add isolation and concurrency after the
   adapter contract is stable.
9. **F09 `batch-summary` hardening** — reuse discovery/publication improvements.
10. **F10–F13 integration, schemas, tests, docs, and release gates**.

Suggested commits:

```text
feat(discovery): add deterministic GenBank cohort discovery
feat(records): retain and select lossless source records
feat(records): add list extract filter and split actions
feat(index): add normalized SQLite cohort schema and lifecycle
feat(index): compile safe feature queries to parameterized SQL
feat(batch): add registered batch adapters and run planning
feat(batch): add manifest-driven execution and verified resume
fix(batch-summary): adopt shared discovery and safe publication
docs: document v0.9.0 cohort workflows and contracts
test(release): add v0.9.0 cross-platform release gates
chore(release): prepare genbank-parser 0.9.0
```

## 20. Release acceptance checklist

- [ ] `gbparse records list|extract|filter|split` are implemented and documented.
- [ ] Whole-record GenBank round trips preserve semantic record and feature data.
- [ ] Split filenames cannot collide and every file is mapped in a manifest.
- [ ] `gbparse index build|update|query|status` are implemented.
- [ ] Index build/update are atomic and verified before replacement.
- [ ] Index schema is normalized, versioned, and foreign-key checked.
- [ ] Index query uses the existing grammar and parameterized SQL only.
- [ ] Direct and indexed query results agree on the same fixtures.
- [ ] `gbparse batch` supports every registered single-input output topology.
- [ ] Unsupported commands/options fail before any job is started.
- [ ] Batch manifests contain input and output hashes plus exact argv.
- [ ] Resume skips only hash-verified unchanged successful jobs.
- [ ] Partial batch failure exits nonzero and leaves no partial job artifact.
- [ ] `batch-summary` remains compatible but no longer silently succeeds on an
      empty or wholly failed cohort.
- [ ] All machine output remains data-only on stdout.
- [ ] All new JSON documents validate against packaged schemas.
- [ ] Linux and Windows CI cover SQLite replacement and batch concurrency.
- [ ] Existing v0.8.5 tests remain green.
- [ ] README, `SKILL.md`, changelog, version, and golden provenance are updated.

## 21. Definition of done

v0.9.0 is complete when a user can safely:

1. inspect, select, and split whole records without losing GenBank semantics;
2. build or update a persistent cohort index and obtain the same annotation
   matches as direct `gbparse query`;
3. run supported single-input commands across a cohort with collision-free
   outputs, explicit failure accounting, and a schema-versioned provenance
   manifest; and
4. resume an interrupted or partially failed run without trusting stale,
   missing, or modified artifacts.

The release is not done if success can be reported after all inputs failed,
if a changed input can be skipped during resume, if arbitrary SQL or shell text
is executed, if record extraction drops source semantics, or if an interrupted
index update can replace the last known-good database.
