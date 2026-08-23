# Mobilome Patch 0.1

Status: follow-up patch plan for the mobilome integration shipped in 0.6.0.
Derived from a commit-by-commit audit of the 11 implementation commits
(`d327ba1`..`a5cdcd5`) against `mobilome-parser-integration-plan.md`.

Audit date: 2026-08-23

Planning baseline: main at 589c27e (local release candidate, not pushed).

Audit verification state: full suite 103 passed / 6 optional-visualization
skips; scoped Ruff check and format clean on all new paths; wheel contains all
four data resources plus `py.typed`; `discover --ruleset mobilome` unchanged;
CLI error contract verified without tracebacks; atomic writer, protected-path
checks, schema validation, circular/linear distance math, and both
origin-wrapping strand directions verified by direct probing.

Scope: correctness defects plus plan-required gaps left unfinished by the
integrated implementation. Every item names its file and its acceptance check.

## P0 — Classification correctness (fix before trusting 0.6.0 output)

- [ ] **Fix the space-separated negative control in record tokens.**
  `src/genbank_parser/mobilome/replicons.py` `_RECORD_TOKEN_PATTERNS` blocks
  hyphenated forms (`plasmid-like`) but a DEFINITION/name containing
  `chromosome partition protein` still yields a tentative chromosome
  classification. Suppress token evidence when the token is part of a
  product-like phrase (e.g. `protein`/`proteins` within a bounded window after
  the token, or a curated negative-phrase lookahead). Acceptance: both declared
  negative controls (`plasmid-like protein`, `chromosome partition protein`)
  produce no classification evidence in name and description positions, while
  `Bacillus plasmid pXO1, complete sequence`-style statements remain tentative.
- [ ] **Add the missing negative-control fixture and regression.** No fixture
  or test exercises either declared negative control today (they exist only as
  `negative_controls` metadata in `provenance.yaml`). Extend
  `tests/fixtures/mobilome_inventory.gb` with product-phrase description
  records and assert `label == "unknown"`, `status == "insufficient"`.

## P1 — Required scientific regressions that were never written

- [ ] **Execute the forbidden-claim contract.** `FORBIDDEN_CLAIMS` (15
  phrases) in `tests/test_mobilome_fixtures.py` is dead code; only 6 phrases
  are asserted, and only against one record's hypothesis summaries. Assert the
  full table over the complete text, TSV, and JSON renderings of all three
  fixtures plus every configured `wording`/`interpretation`/`limitations`
  string in `markers.yaml` and `inference.yaml` (plan §16.6).
- [ ] **Input-error tests.** Empty file and malformed GenBank must raise
  `MobilomeInputError` from `analyze_mobilome()` and exit nonzero from the CLI
  without a traceback (plan §16.1, §16.8). Behavior verified by audit probe;
  currently untested.
- [ ] **pXO adversarial regressions.** Records with one, two, three,
  duplicated, and reordered `pXO[12]-NN` products must produce observations
  only and never an aggregate hypothesis (plan §16.5). The inference fixture
  carries a single pXO marker; the three-marker case was verified by probe but
  has no committed fixture or test.
- [ ] **Phage-module adversarial regressions.** Single and multiple distinct
  phage modules (integration/packaging/structure/lysis/replication), with and
  without source-declared plasmid plus eligible replication evidence, must
  never create a phagemid or phage-module-bearing-plasmid hypothesis (plan
  §16.5). No phage fixture beyond one Zot-like CDS exists today.
- [ ] **CRISPR comment-isolation regression.** A summary comment mentioning
  CRISPR must not create a hit; a `repeat_region` with `/rpt_type="CRISPR"`
  and a Cas CDS must remain separate observations with no combined hypothesis
  (plan §10.5, §16.5). The current fixture only proves the positive array+Cas
  path.
- [ ] **Pseudogenized-Rep-only regression.** A record whose only replication
  candidate is pseudogenized must report insufficient with no
  satellite/degenerate/evolutionary-transition wording (plan §16.5). The
  evidence fixture mixes `repS` pseudogene with eligible `repA`/`repB`, so the
  all-pseudo branch of `_infer_component_rule` is never exercised in
  isolation.
- [ ] **Mechanism non-selection regression.** RepA/RepB-only,
  replication-relaxation fusion, and pXO-bearing records must emit the
  unresolved-mechanism limitation and never theta/rolling-circle/strand-
  displacement wording (plan §16.5). Also note the catalog has no TPR marker,
  so the plan's TPR regression cannot run; either add the marker or record the
  deliberate omission in `docs/mobilome_evidence_reference.md`.
- [ ] **Database fail-closed matrix.** `tests/test_mobilome_database.py`
  covers only strength 4 and an extra custom file. Add rejection cases per
  plan §16.3: duplicate facet/marker/rule/component/source IDs; unknown
  facet/marker/component/source references; invalid regex; strength 0;
  field-ceiling violation; empty pattern/matcher/limitation lists; missing
  type-specific provenance fields; catalog/provenance/inference major-version
  disagreement; custom directory missing one resource file.
- [ ] **Scanner contract tests.** Add per plan §16.4: negative patterns
  suppress only their matcher; multi-value qualifiers yield traceable per-value
  reasons; exact duplicate reasons deduplicate while distinct reasons do not;
  exact/casefold_exact/regex behave distinctly; same text in gene/product/note
  does not inflate inference (only reasons per field).
- [ ] **Writer hardening tests.** Add per plan §16.7: symlink/hardlink alias
  rejection where the platform supports it; simulated interrupted write leaves
  no partial destination or temp file; source mutated between analysis and
  write is rejected by the SHA-256 comparison in `write_mobilome_report()`.
- [ ] **CLI matrix tests.** Add per plan §16.8: each `--include` value;
  invalid `--min-evidence`; `--database-dir` success and incomplete-directory
  failure; malformed and empty input without traceback.
- [ ] **Fixture completion.** Extend `tests/fixtures/mobilome_inference.gb`
  per plan §15.3 (pseudogenized-Rep-only, conflicting replication candidates,
  each incomplete mobility subset, complete record-local mobility pattern,
  helper+target on one record, incomplete helper, phage single/multi records,
  linear and circular TA pairs at the 2000 bp boundary) and
  `tests/fixtures/mobilome_evidence.gb` per plan §15.2 (lone toxin, lone
  antitoxin, AimR-like, efflux product, multi-value qualifiers, duplicate
  identical reasons, negative-pattern product, db_xref plus generic product).
  These fixtures gate the P1 regressions above; keep them synthetic and
  LF-normalized via the existing `.gitattributes` rules.

## P2 — Contract and provenance polish

- [ ] **Surface the unlocatable-feature limitation in TA pairing.** Plan §6.1
  requires unlocatable features to be excluded from distance rules "with an
  explicit limitation"; `_distance_between_hits()` currently drops the pair
  silently. Emit a record-level limitation when a same-record toxin/antitoxin
  pair is skipped for missing segments.
- [ ] **Fix annotation-provenance detection and normalization.**
  `_annotation_provenance()` in `mobilome/__init__.py` scans only string
  annotation values; Biopython delivers structured COMMENT blocks (Bakta's
  actual layout) as `annotations['structured_comment']` (a dict) and is
  skipped. Also normalize the captured database version (`v5.0.1` vs `1.9.4`
  currently keep and strip the `v` prefix inconsistently). Add a fixture with
  a realistic Bakta COMMENT; regenerate goldens if output changes.
- [ ] **Decide the null-pipeline entry.** Inputs with no detected pipeline
  emit an `annotation_pipelines` entry with `name: null` (visible in the
  golden). Either omit records with no detected pipeline or document the
  null-name semantics in `docs/mobilome_evidence_reference.md`; regenerate
  goldens if changed.
- [ ] **Add the MOB-suite literature anchor.** Plan §4.1 lists Robertson and
  Nash (2018) for MOB-suite database/threshold limitations; it is absent from
  `provenance.yaml` and the docs. Add the source record and cite it from the
  external-handoffs section of `docs/mobilome_evidence_reference.md`
  (handoffs carry no `source_ids` in schema v1, so this is a documentation and
  catalog fix, not a schema change).
- [ ] **Fix the broken README link.** Line 135: `the obilome evidence
  reference](docs/mobilome_evidence_reference.md)` is missing its opening
  bracket; renders as literal text with a dangling link.
- [ ] **Document or normalize ragged TSV rows.** `_tsv_values()` strips
  trailing empty cells, so data rows carry 15-53 columns against the fixed
  58-column header. Lossless and deterministic, but the deviation from the
  full-width contract is undocumented; state it in
  `docs/mobilome_evidence_reference.md` or emit full-width rows.
- [ ] **Exercise or annotate the record-annotation classification path.**
  Biopython never maps source `/plasmid` or `/chromosome` into
  `record.annotations`, so the `record_annotation_*` rules in
  `replicons.py` cannot fire from GenBank input. Either add a unit test with a
  synthetic `GenBankRecord` or mark the branch as defensive with a comment.
- [ ] **Record the sequencing deviations.** Phase 1 committed the final
  marker/provenance/inference catalog in full (later phase data scopes were
  no-ops), Phase 5 goldens pinned only a hand-constructed empty report, and
  the Phase 8 release commit modified `report.py` and two test files outside
  its staging scope. No code change required; keep the note in the changelog
  or audit trail so the phase-gate history is not misread as incremental.

## P3 — Deferred, separately authorized (plan §19)

- [ ] Standalone mobilome-parser comparison audit after patch verification.
- [ ] Optional real-input check with `PF_NNT_reoriented.gbff` (plan §16.9).
- [ ] Consider namespacing marker IDs that shadow facet IDs
  (`crispr_array`, `crispr_cas`) in a future catalog bump; also replace the
  hard-coded `pxo_numbered_product_annotation` whitelist in
  `database.py` with a declared retained-evidence token set.

## Gate for this patch set

Full pytest, scoped Ruff check and format for touched paths, `compileall`,
`pip check`, regenerated goldens where serializers change, and unchanged
`discover --ruleset mobilome` output. Catalog wording or threshold changes
bump `catalog_version`/`inference_version`; report-shape changes require
`gbparse.mobilome.v2`, so keep P2 provenance/TSV items additive or
documentation-only.

## Execution record — 2026-08-23

Implemented in the scoped follow-up patch:

- P0 record-token negative controls now suppress product-like phrases, with
  fixture and classification regressions for both declared controls.
- Structured Bakta comments are traversed through Biopython's typed
  `structured_comment` mapping, and `v`/`version` prefixes are normalized.
- Unlocatable toxin/antitoxin candidates now emit an explicit record-level
  spatial-pairing limitation.
- Curated provenance entries now require `supports`, and the MOB-suite
  literature anchor plus external-handoff documentation were added.
- Adversarial tests cover forbidden wording across all renderings and catalog
  text, malformed/empty input, pXO and phage aggregate negatives, CRISPR
  comment isolation, pseudogenized Rep-only output, mechanism non-selection,
  scanner matching contracts, CLI/database matrices, and safe-writer aliases,
  interruption cleanup, and source mutation.
- The v1 JSON golden was regenerated for the changed catalog-resource hash;
  the report schema and `discover --ruleset mobilome` implementation remain
  unchanged.

Validation completed: 45 mobilome tests pass; scoped Ruff, format, mypy,
`compileall`, `pip check`, and `git diff --check` pass. The full repository
suite passes with the writable `MPLCONFIGDIR`/`Agg` validation environment
(128 passed); the default local Tk/Matplotlib environment still exposes the
same six pre-existing visualization failures.
