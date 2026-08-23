# Mobilome Patch 0.2

Status: residual hardening plan after verifying the patch-0.1 implementation
(commit `c2a7b50`, "fix: harden mobilome patch contracts") against
`mobilome-patch-0.1.md`.

Re-audit date: 2026-08-23

Planning baseline: main at c2a7b50.

Re-audit verification state: full suite 122 passed / 6 optional-visualization
skips (45 mobilome tests, matching the execution record); scoped Ruff check
and format clean on all touched paths; P0 negative-control fix, Bakta
structured-comment provenance, and golden regeneration independently
re-probed and confirmed; live JSON render byte-matches the committed golden.

Scope: the remaining correctness and contract gaps. Patch 0.1 closed both P0
items and roughly two thirds of the P1 list; what remains is itemized below
with the same file-and-acceptance-check format.

## P1 — Correctness of the forbidden-claim contract itself

- [ ] **Extend the catalog forbidden-phrase scan to provenance.yaml.**
  `tests/test_mobilome_patch_01.py::
  test_forbidden_claims_are_absent_from_all_outputs_and_catalog_wording`
  scans only `markers.yaml` and `inference.yaml`. `provenance.yaml` currently
  contains a forbidden phrase at line 71 (erez-2017-arbitrium limitation:
  "One annotation does not prove active communication in a new genome.") —
  the exact phrase that was scrubbed from markers.yaml in this patch. Add
  `provenance.yaml` to the scanned names and reword that limitation (e.g.
  "communication activity", mirroring the markers.yaml fix). Acceptance: the
  test fails before the reword, passes after, over all three catalog files.

## P1 — Version discipline for changed catalogs

- [ ] **Bump catalog and provenance versions for content changes.**
  Commit `c2a7b50` changed markers.yaml (wording) and provenance.yaml (new
  Robertson and Nash source) while both versions remained "1.0.0"; the
  golden's resource hashes changed under unchanged version strings. Per plan
  §6.1 ("Marker/source content changes bump catalog_version") bump
  `catalog_version` to 1.0.1 and `provenance_version` to 1.0.1 (major
  versions must continue to agree), update the version assertions in
  `tests/test_mobilome_database.py` and the CI wheel gate, and regenerate
  `tests/golden/mobilome_v1.json` (and .tsv/.txt if versions surface there —
  they do in the JSON header block). Fold the patch-0.1 wording change and
  the P1 reword above into the same bump so the catalog converges on one
  consistent version. Acceptance: no resource hash changes under an
  unchanged version string in any future diff.

## P1 — Remaining plan §16 regressions still unwritten

- [ ] **Database fail-closed matrix.** Still only 3 rejection cases tested
  (strength 4, extra custom file, curated-without-supports). Add per plan
  §16.3: duplicate facet/marker/rule/component/source IDs; unknown
  facet/marker/component/rule references; invalid regex at load; strength 0;
  field-ceiling violation; empty pattern/matcher/limitations lists;
  local-rule sources missing mandatory fields; catalog/provenance/inference
  major-version disagreement; custom directory missing exactly one resource
  file. Acceptance: each case fails `load_mobilome_database()` with a
  distinct message.
- [ ] **CLI residual cases.** `--database-dir` pointing at a directory with
  one file missing (distinct from the extra-file case); invalid
  `--min-evidence 4` exits 2 without traceback. Acceptance: both error
  surfaces asserted in `tests/test_mobilome_patch_01.py` or
  `tests/test_mobilome_cli.py`.
- [ ] **TA boundary and lone-marker regressions.** Linear and circular
  toxin/antitoxin pairs at 1999/2000/2001 bp gaps (inclusive boundary
  behavior of `max_circular_gap_bp`); lone toxin and lone antitoxin produce
  observations but no pair hypothesis. The `_genbank()` builder in
  `tests/test_mobilome_patch_01.py` can be extended with coordinates and
  marker-gene parameters. Acceptance: boundary value 2000 pairs, 2001 does
  not; lone markers never emit a TA hypothesis.
- [ ] **Mobility component subsets.** Each incomplete subset (relaxase-only,
  T4CP-only, single-MPF, T4CP+one-MPF), the complete record-local pattern
  (oriT+relaxase on distinct features), and helper+target components on the
  same record (must NOT satisfy `distinct_helper_record`). Only the target
  side is currently exercised. Acceptance: exact `missing_components` tuples
  asserted per subset; same-record helper does not create a cross-record
  hypothesis.
- [ ] **Scanner residual contracts.** db_xref plus generic product reasons
  coexist on one hit; distinct reasons on the same field do not collapse;
  same annotation text in gene and product yields separate field-scoped
  reasons without eligibility inflation (the plan's "do not multiply
  evidence" clause). Acceptance in `tests/test_mobilome_scanner.py`.

## P2 — Documentation and trail hygiene

- [ ] **Fix the broken README link.** Line 135 still reads `the obilome
  evidence reference](docs/mobilome_evidence_reference.md)` — missing the
  opening bracket; patch 0.1 did not touch README.md. Acceptance: the link
  renders.
- [ ] **Add the patch changelog entry.** `changelogs.md` has no entry for
  c2a7b50, and the phase-gate history note from patch 0.1 (Phase 5
  placeholder goldens, Phase 8 out-of-scope serializer edit, Phase 1
  front-loaded catalog) was never recorded. One short 0.6.1 "Hardening"
  section covering both is enough; keep the sequencing note in the patch-0.1
  execution record if preferred. Acceptance: the release trail is
  reconstructible from changelogs.md alone.
- [ ] **Consider a v2 home for record-level inference limitations.** The
  unlocatable-TA limitation currently lives in
  `RepliconInventory.classification.limitations` — the only per-record
  limitations field available in schema v1 without a breaking change. Fine
  as additive v1 behavior, but it mixes classification caveats with spatial
  caveats. Record as a v2 design note (e.g. inventory-level
  `spatial_limitations`) rather than changing v1 now.

## P3 — Deferred, separately authorized (unchanged from patch 0.1)

- [ ] Standalone mobilome-parser comparison audit (plan §19).
- [ ] Optional real-input check with `PF_NNT_reoriented.gbff` (plan §16.9).
- [ ] v2 catalog hygiene: namespace marker IDs that shadow facet IDs
  (`crispr_array`, `crispr_cas`); replace the hard-coded
  `pxo_numbered_product_annotation` whitelist in `database.py` with a
  declared retained-evidence token set.

## Gate for this patch set

Identical to patch 0.1: full pytest, scoped Ruff check and format,
`compileall`, `pip check`, regenerated goldens where serializers or catalog
versions change, and unchanged `discover --ruleset mobilome` output. The
version-bump item makes golden regeneration mandatory this time. Hand the
catalog reword plus version bumps to one commit and the new tests to another
so the wording change and its enforcing scan land reviewably.

## Execution record — 2026-08-23

Implemented in the working tree:

- Bumped the mobilome catalog and provenance resources to 1.0.1, updated the
  packaged-wheel CI assertions, reworded the provenance forbidden phrase, and
  extended the forbidden-claim scan through all three catalog resources.
- Made custom YAML loading reject duplicate mapping keys, including duplicate
  component IDs, and added the requested fail-closed database matrix for
  references, regexes, strengths, empty lists, provenance fields, versions,
  and incomplete custom directories.
- Added CLI residual tests, inclusive 1999/2000/2001 bp linear/circular
  toxin-antitoxin boundaries, lone-marker negatives, mobility component
  subsets, complete local components, same-record helper exclusion, and
  field-scoped scanner reasons.
- Added the 0.6.1 changelog trail and documented the v2 home for spatial
  inference limitations. The README link was verified already corrected in
  the published 0.1 tree, so no README diff was necessary.

Validation completed: 74 mobilome tests and the full repository suite pass
under writable MPLCONFIGDIR/Agg (157 passed); scoped Ruff/format, mypy,
compileall, pip check, wheel install, discover, and diff checks pass.
