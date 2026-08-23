# Mobilome Parser Integration Plan

Status: implementation-ready plan; no implementation, commit, tag, archive, or push is authorized by this document.

Planning date: 2026-08-22

Target repository: genbank-feature-parser

Primary review: [mobilome-parser-review.md](../mobilome-parser/mobilome-parser-review.md)

Planning baseline: main and origin/main at fc480ebfe1d2df987857a5bd82ccb52e9388cd92 when this plan was drafted. Re-check the branch, worktree, tools, and tests before implementation because this is a point-in-time observation.

Target release: genbank-parser 0.6.0

Report schema: gbparse.mobilome.v1

Initial evidence-catalog version: 1.0.0

Initial inference-policy version: 1.0.0

## 1. Outcome

Integrate the useful requirements and scientifically supportable parts of the standalone mobilome-parser prototype into genbank-feature-parser as a native, typed, provenance-aware analysis:

~~~text
gbparse mobilome INPUT
    -> canonical read_genbank()
    -> complete typed record inventory
    -> field-aware, multi-label evidence scan
    -> conservative per-replicon assessments
    -> cautious cross-record hypotheses
    -> deterministic text, JSON, or TSV report
~~~

The implementation must use the existing Biopython-backed GenBankDocument and GenBankFeature model. The standalone regex parser is a requirements source and a defect oracle, not a code source and not a behavioral golden.

The completed feature will:

- inventory every record, including chromosomes, declared plasmids, and unresolved contigs;
- preserve compound and origin-spanning locations, unknown strands, partial features, and both GenBank pseudogene representations;
- retain every supported evidence facet for a feature instead of stopping at the first keyword category;
- separate observations, annotation-supported candidates, tentative hypotheses, conflicts, missing components, and external-tool recommendations;
- keep replication origin, relaxase, type IV coupling protein, and mating-pair-formation evidence distinct;
- emit schema-, tool-, catalog-, inference-, input-, and annotation-provenance;
- avoid claims of phenotype, transfer, co-transfer, replicon ancestry, replication mechanism, autonomy, antimicrobial resistance, virulence, or phagemid identity that the input does not establish;
- leave the existing gbparse discover --ruleset mobilome behavior available and unchanged.

## 2. Decisions and non-goals

### 2.1 Decisions

1. Add a native genbank_parser.mobilome subpackage and a gbparse mobilome subcommand.
2. Keep read_genbank() as the only GenBank parser.
3. Model mobilome evidence independently from the generic discover RuleMatch model because the new report needs exact field/value/pattern provenance and versioned semantics.
4. Keep the packaged rulesets/mobilome.yaml file and discover command intact. The discover command remains an annotation-island scanner; mobilome becomes a replicon-centric evidence assessment.
5. Use immutable typed report objects and pure serializers, following the newer MEOR subsystem.
6. Package three separately versioned scientific resources: markers.yaml, provenance.yaml, and inference.yaml.
7. Treat evidence strength as an ordinal measure of annotation specificity, never as probability, effect size, or biological certainty.
8. Retain all raw reasons in canonical JSON and TSV. The --min-evidence option changes eligibility for inference and the primary text view; it does not erase lower-tier observations.
9. Apply --include only to detailed presentation. Inventory and cross-record scanning always cover all parsed records.
10. Use one supported CLI executable, gbparse. Do not add another wrapper script.
11. Keep the public Python entry point at genbank_parser.mobilome rather than importing mobilome into the package root.
12. Use deterministic output with no generated-at timestamp. Input basename plus SHA-256 replaces an absolute path.

### 2.2 Non-goals for version 0.6.0

- Reimplementing or copying analyze_extrachromosomal.py.
- Parsing raw INSDC text with regular expressions.
- Calling MOB-suite, PlasmidFinder, AMRFinderPlus, VFDB searches, prophage callers, ICE callers, HMMER, or BLAST automatically.
- Inferring plasmid incompatibility group from product text.
- Selecting theta, rolling-circle, strand-displacement, or another replication mechanism from generic RepA/RepB, TPR, pXO, or replication-relaxation labels.
- Calling a record a satellite from absent markers or one pseudogenized replication candidate.
- Calling ICE/prophage boundaries or performing multi-isolate/cohort comparison
  in v1; deterministic JSON enables a separately designed future comparator.
- Calling mobilization, conjugation, obligate co-transfer, phenotype, expression, or experimentally active systems.
- Calling an AMR determinant or virulence factor from a generic keyword or domain alone.
- Calling an element a phagemid from one integrase, Zot-like, or other phage-related product.
- Treating a Bakta summary comment as the underlying CRISPR feature evidence.
- Changing read_genbank(), GenBankFeature, GenBankDocument, discover output, or the existing rulesets/mobilome.yaml in the first integration.
- Deleting or archiving the standalone mobilome-parser before the native implementation is validated.
- Committing real or private GBFF inputs.
- Pushing commits, tags, or releases without separate authorization.

## 3. Repository and baseline constraints

### 3.1 Existing components to reuse

| Concern | Existing authority | Required use |
|---|---|---|
| GenBank parsing | src/genbank_parser/io.py::read_genbank | Parse every input once; never add a second parser |
| Typed feature semantics | src/genbank_parser/model.py | Preserve segments, biological length, strand, partial state, and is_pseudo |
| Record identity | GenBankDocument record order plus feature_index | Use record_index to disambiguate duplicate record IDs |
| Pure subsystem pattern | src/genbank_parser/meor | Follow its database, models, orchestrator, and serializer split |
| Circular utilities | src/genbank_parser/spatial.py | Extend with named, tested segment-aware distance helpers only when needed |
| CLI registration | src/genbank_parser/cli.py | Add parser construction and dispatch without changing existing commands |
| Package resources | pyproject.toml | Add data/mobilome resources and verify them from an installed wheel |
| Output safety precedent | neighborhood/report code | Strengthen it with atomic writes and protected-path checks |

### 3.2 Files that remain unchanged through the initial phases

~~~text
src/genbank_parser/io.py
src/genbank_parser/model.py
src/genbank_parser/discover.py
src/genbank_parser/rulesets/mobilome.yaml
src/genbank_parser/__init__.py
~~~

If implementation reveals that one of these files must change, stop that phase, explain the contract change, add a focused regression first, and place the change in its own reviewable commit.

### 3.3 Recorded local baseline

The reliable Windows baseline at planning time was:

~~~powershell
$env:MPLBACKEND = "Agg"
$env:MPLCONFIGDIR = "C:\tmp\gbparse-mobilome-mpl"
python -m pytest -q -p no:cacheprovider
if ($LASTEXITCODE -ne 0) { throw "Recorded baseline pytest failed" }
python -m pip check
if ($LASTEXITCODE -ne 0) { throw "Recorded baseline dependency check failed" }
~~~

Observed result: 83 tests passed with one existing Biopython warning in the alternate-translation-table test; pip check reported no broken requirements.

Whole-repository Ruff and mypy are not currently valid release claims:

- Ruff 0.16.1 reports pre-existing findings outside this work.
- mypy 2.3.0 crashes internally in this checkout.

Therefore every phase blocks on Ruff check for new Python paths and any touched
legacy path after its named pre-existing diagnostic is resolved. Ruff format
check is blocking for new files and for legacy files already format-clean; it is
not run across the known format-dirty cli.py or tests/test_cli.py. The cli.py
diff is reviewed directly, tests/test_cli.py remains unchanged, and broad
formatting remains separate work. Mypy
becomes blocking only after Phase 0 records a working pinned invocation;
otherwise it remains an explicitly reported advisory check.

### 3.4 Worktree discipline

Before every phase:

~~~powershell
$ErrorActionPreference = "Stop"
$touchedPaths = @(
    # Populate from this phase's exact staging scope before running the gate.
)
if ($touchedPaths.Count -eq 0) { throw "Touched-path allowlist is empty" }

git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' status --short --branch
if ($LASTEXITCODE -ne 0) { throw "Could not inspect worktree status" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --check
if ($LASTEXITCODE -ne 0) { throw "Worktree diff check failed" }
$conflictHits = @(rg -n "^(<<<<<<<|=======|>>>>>>>)" -- $touchedPaths)
if ($LASTEXITCODE -gt 1) { throw "Conflict-marker scan failed" }
if ($conflictHits.Count -ne 0) {
    $conflictHits
    throw "Conflict markers found in touched paths"
}
~~~

Preserve the unrelated untracked paths observed when this plan was drafted:

~~~text
audit-reference/
gbparse-expandviz-plan.md
gbparse-meor-audit-assessment.txt
gbparse-meor-identifier-audit-plan.md
gbparse-patch-6-plan.md
scripts/inspect_audit_references.py
scripts/simulate_meor_audit.py
mobilome-parser-integration-plan.md
~~~

Stage only exact phase paths. Never use git add ., git add -A, or an unreviewed wildcard. When a future phase moves a file, stage its addition and deletion together.

## 4. Review findings converted into implementation requirements

The following table is the correspondence ledger between the standalone review, primary literature or authoritative documentation, the corrected implementation rule, and the regression that prevents reintroduction.

| Review or reference issue | Corrected implementation rule | Required regression |
|---|---|---|
| The prototype regex parser omitted feature classes and lost typed locations | Only read_genbank() may parse input; every FeatureRef is derived from GenBankFeature | Inventory counts every source, regulatory, rep_origin, repeat_region, CDS, and compound feature in synthetic fixtures |
| The prototype recognized /pseudo but not /pseudogene | Use GenBankFeature.is_pseudo, which covers pseudogene feature type, /pseudo, and /pseudogene | Separate fixtures for all three representations |
| all() over an empty replication-marker set produced a false satellite call | Require a nonempty eligible candidate set before any aggregate predicate; empty evidence is insufficient | Zero-marker record emits no satellite/autonomy/evolution claim |
| Generic RepA/RepB, TPR, pXO, and replication-relaxation labels were scored as mechanisms | Version 0.6.0 reports replication candidates and unresolved mechanism; it does not select a mechanism from those labels | Forbidden-mechanism assertions for each ambiguous marker and for the chromosome |
| A first-match category hid biological multiplicity | A feature may yield multiple MobilomeHit objects and every exact EvidenceReason is retained | Replication-relaxation fusion and Zot-like examples retain all applicable facets |
| Product keywords were promoted toward AMR/VF conclusions | Product and note patterns are candidates only; stronger AMR/VF status requires structured, provenance-bearing expert annotations or imported tool results | Generic beta-lactamase-related domain and efflux product remain candidates |
| Advertised oriT, CRISPR, phagemid, Inc, ICE, and cross-replicon behavior was not implemented | Report only implemented typed evidence; represent unsupported external analyses as not-run handoffs | Handoff rows say not_run with required input, database/version, and limitations |
| Mobility components were collapsed | Model oriT, relaxase, T4CP, and MPF/T4SS separately | Complete and incomplete component matrices list exact missing components |
| Cross-record evidence was called obligate co-mobilization | At most emit possible helper-dependent mobilization, with compatibility, cognate oriT, expression, and transfer listed as unresolved | Exact forbidden-phrase tests reject obligate co-transfer/co-transmission wording |
| ToxI was placed under Type I | Model ToxI as the RNA antitoxin of the Type III ToxIN system; allow non-CDS antitoxin evidence | ToxN CDS plus ToxI ncRNA/repeat evidence yields only a Type III pair candidate |
| HepT/MntA chemistry and organism correspondence were wrong | The reference states that the studied MntA used ATP-dependent poly-AMPylation of HepT Tyr104; the software label is a HEPN/MNT pair candidate and does not hard-code an unsupported TA type number | Documentation assertion plus marker/provenance test for the corrected source |
| Bakta was said not to annotate oriV/oriT | Treat support as version-sensitive; report observed typed features and captured Bakta version/database when present | Fixtures cover generic rep_origin separately from explicit oriV/oriT annotations |
| Multiple pXO1-xx products were called diagnostic | Version 0.6.0 emits traceable pXO-numbered product observations but no aggregate pXO-like classifier; a future rule requires explicit ordinal/synteny semantics plus positive and negative controls | One, two, three, duplicate, and reordered markers all remain observations without an aggregate identity claim |
| External typing was treated too categorically | Version 0.6.0 emits not-run handoffs that state the provenance a future imported result would require; it does not import or synthesize external results | Every handoff remains not_run and lists tool/database/threshold/taxonomic requirements |
| A lone phage-related feature approached a phagemid call | Version 0.6.0 preserves separate phage-module facets but disables aggregate phage-module-bearing-plasmid and phagemid hypotheses pending a validated operating threshold | Single and multiple module patterns remain observations and report the disabled-rule reason |
| AimR/AimP-like annotation was overinterpreted | Report a phage communication-system annotation candidate only; do not infer active communication, host response, cross-replicon signaling, or co-option | Forbidden-wording test for function and communication |
| Output could overwrite input and lacked a stable contract | Add schema-versioned pure serializers and an atomic protected-path writer | Input/custom catalog collision, hardlink/symlink where supported, interrupted write, and unchanged-input tests |

### 4.1 Scientific correspondence sources

The implementation documentation and provenance catalog must cite claims at assertion granularity rather than attach a general bibliography to an entire category:

- Yao et al. (2020), ATP-dependent polyadenylylation/AMPylation of HepT by MntA in the studied HEPN/MNT system: [Nucleic Acids Research article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/), DOI 10.1093/nar/gkaa855.
- Blower et al. (2012), ToxI as the RNA antitoxin of Type III ToxIN: [Nucleic Acids Research article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/), DOI 10.1093/nar/gks231.
- Smillie et al. (2010), mobilizable plasmids may use helper mating-pair-formation machinery but annotation co-occurrence does not establish obligate co-transfer: [Mobility of Plasmids](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/), DOI 10.1128/MMBR.00020-10.
- Akhtar and Khan (2012), two independent pXO1 replicons can each support full-length plasmid replication, so generic pXO-like annotations do not force one replication mechanism: [Plasmid article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4186245/), DOI 10.1016/j.plasmid.2011.12.012.
- Robertson and Nash (2018), MOB-suite performance and database/taxonomic limitations require versioned provenance and calibrated interpretation: [Microbial Genomics article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/), DOI 10.1099/mgen.0.000206.
- Current Bakta documentation lists origin, CRISPR, and pseudogene annotation capabilities; observed input features and annotation version take precedence over blanket tool claims: [Bakta documentation](https://bakta.readthedocs.io/en/latest/).
- Erez et al. (2017) supports the AimR/AimP/AimX phage communication correspondence, not an inference that any single AimR-like annotation proves active communication: [Nature article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5378303/), DOI 10.1038/nature21049.

Each entry in provenance.yaml must identify exactly which marker definition, wording choice, or limitation it supports. A publication supporting biochemical correspondence does not automatically validate an annotation pattern as functional in the user's input.

## 5. Proposed architecture

### 5.1 Data flow and separation of responsibilities

~~~text
Path
  |
  v
read_genbank() ------------------------------ canonical parsing boundary
  |
  v
GenBankDocument
  |\
  | +--> inventory_replicons() ------------- all records, metadata, counts
  |
  +----> scan_mobilome_features() ----------- all fields, all reasons, all facets
            |
            +--> eligible evidence ---------- threshold and pseudo/partial policy
                    |
                    +--> per-replicon inference
                    +--> cross-record inference
                    +--> external handoff descriptions
                              |
                              v
                         MobilomeReport
                           /    |    \
                        text   JSON   TSV ------------- pure serializers
                           \    |    /
                         safe atomic writer
~~~

Parsing, evidence matching, inference, serialization, and file writing must remain separate. Unit tests should construct typed objects without invoking the CLI, and CLI tests should only verify orchestration and error presentation.

### 5.2 File map

Illustrative target diff:

~~~text
 src/genbank_parser/
+  mobilome/
+    __init__.py
+    models.py
+    database.py
+    scanner.py
+    replicons.py
+    inference.py
+    report.py
+  data/mobilome/
+    markers.yaml
+    provenance.yaml
+    inference.yaml
+    report.schema.json
   cli.py
   spatial.py                    # only shared circular-distance helpers
   py.typed

+docs/mobilome_evidence_reference.md

 tests/fixtures/
+  mobilome_inventory.gb
+  mobilome_evidence.gb
+  mobilome_inference.gb
 tests/golden/
+  mobilome_v1.json
+  mobilome_v1.tsv
+  mobilome_v1.txt
+tests/test_mobilome_fixtures.py
+tests/test_mobilome_database.py
+tests/test_mobilome_replicons.py
+tests/test_mobilome_scanner.py
+tests/test_mobilome_inference.py
+tests/test_mobilome_spatial.py
+tests/test_mobilome_report.py
+tests/test_mobilome_api.py
+tests/test_mobilome_cli.py

+.gitattributes
 tests/conftest.py
 pyproject.toml
 .github/workflows/ci.yml
 README.md
 SKILL.md
 changelogs.md
~~~

Do not add a new parser file, standalone executable, duplicated feature model, or a copy of mobilome_biology.md.

Fence semantics are deliberate throughout this plan: only blocks fenced as
`diff` are applyable candidate patches and must pass `git apply --check` against
the recorded baseline. Diff-shaped blocks fenced as `text` are interface or
logic sketches; implementers must translate them into repository-native code
and must not apply them verbatim.

## 6. Typed contracts

The following is an interface sketch, not a patch to apply verbatim. Exact
imports and serialization helpers must follow repository style.

~~~text
diff --git a/src/genbank_parser/mobilome/models.py b/src/genbank_parser/mobilome/models.py
new file mode 100644
--- /dev/null
+++ b/src/genbank_parser/mobilome/models.py
@@
+from dataclasses import dataclass
+from typing import Literal
+
+SCHEMA_VERSION = "gbparse.mobilome.v1"
+
+RepliconClass = Literal["chromosome", "plasmid", "unknown"]
+RepliconTopology = Literal["circular", "linear", "unknown"]
+AssessmentStatus = Literal[
+    "supported", "tentative", "insufficient", "conflicting"
+]
+MatchMode = Literal["exact", "casefold_exact", "regex"]
+HandoffStatus = Literal["not_run"]
+ParticipantRole = Literal["subject", "target", "helper"]
+
+class MobilomeError(Exception):
+    """Base class for expected user-facing mobilome failures."""
+
+class MobilomeInputError(MobilomeError):
+    """The input cannot produce a typed GenBank document."""
+
+class MobilomeDatabaseError(MobilomeError):
+    """The packaged or custom evidence database is invalid."""
+
+class MobilomeParameterError(MobilomeError):
+    """A public API option is outside its declared closed set."""
+
+class MobilomeOutputError(MobilomeError):
+    """The requested report cannot be written safely."""
+
+@dataclass(frozen=True)
+class FeatureRef:
+    record_index: int
+    record_id: str
+    feature_index: int
+    feature_type: str
+    locus_tag: str | None
+    gene: str | None
+    product: str | None
+    start: int | None
+    end: int | None
+    strand: int | None
+    segments: tuple[tuple[int, int], ...]
+    location_operator: str | None
+    biological_length: int | None
+    wraps_origin: bool
+    is_partial: bool
+    is_pseudo: bool
+
+@dataclass(frozen=True)
+class RecordEvidence:
+    rule_id: str
+    source: Literal[
+        "source_qualifier",
+        "record_annotation",
+        "record_description",
+        "record_name",
+    ]
+    field: str
+    value: str
+    matched_text: str
+    pattern: str
+    source_feature_index: int | None
+    interpreted_as: RepliconClass
+    strength: int
+    source_ids: tuple[str, ...]
+
+@dataclass(frozen=True)
+class SourceQualifier:
+    key: str
+    values: tuple[str, ...]
+
+@dataclass(frozen=True)
+class CatalogResource:
+    name: str
+    sha256: str
+
+@dataclass(frozen=True)
+class AnalysisParameters:
+    include: Literal["all", "chromosome", "plasmid", "unknown"]
+    min_evidence: int
+    database_source: Literal["packaged", "custom"]
+
+@dataclass(frozen=True)
+class AnnotationProvenance:
+    record_indices: tuple[int, ...]
+    name: str | None
+    version: str | None
+    database_version: str | None
+    evidence_fields: tuple[str, ...]
+
+@dataclass(frozen=True)
+class EvidenceReason:
+    reason_id: str
+    field: str
+    match_mode: MatchMode
+    matched_value: str
+    matched_text: str
+    pattern: str
+    strength: int
+    eligible: bool
+    source_ids: tuple[str, ...]
+
+@dataclass(frozen=True)
+class MobilomeHit:
+    hit_id: str
+    marker_id: str
+    marker_label: str
+    facets: tuple[str, ...]
+    feature: FeatureRef
+    reasons: tuple[EvidenceReason, ...]
+    max_strength: int
+    eligible_for_inference: bool
+
+@dataclass(frozen=True)
+class RepliconClassification:
+    label: RepliconClass
+    status: AssessmentStatus
+    supporting_evidence: tuple[RecordEvidence, ...]
+    conflicting_evidence: tuple[RecordEvidence, ...]
+    limitations: tuple[str, ...]
+
+@dataclass(frozen=True)
+class HypothesisParticipant:
+    role: ParticipantRole
+    record_index: int
+    record_id: str
+
+@dataclass(frozen=True)
+class MobilomeHypothesis:
+    hypothesis_id: str
+    rule_id: str
+    kind: str
+    status: AssessmentStatus
+    summary: str
+    participants: tuple[HypothesisParticipant, ...]
+    supporting_hit_ids: tuple[str, ...]
+    missing_components: tuple[str, ...]
+    conflicting_hit_ids: tuple[str, ...]
+    limitations: tuple[str, ...]
+    source_ids: tuple[str, ...]
+
+@dataclass(frozen=True)
+class RepliconInventory:
+    record_index: int
+    record_id: str
+    name: str
+    description: str
+    length: int
+    topology: RepliconTopology
+    source_qualifiers: tuple[SourceQualifier, ...]
+    classification: RepliconClassification
+    feature_count: int
+    cds_count: int
+    pseudo_feature_count: int
+    pseudo_cds_count: int
+    pseudogene_feature_count: int
+
+@dataclass(frozen=True)
+class RepliconAssessment:
+    inventory: RepliconInventory
+    hits: tuple[MobilomeHit, ...]
+    hypotheses: tuple[MobilomeHypothesis, ...]
+
+@dataclass(frozen=True)
+class DisabledAggregateRule:
+    rule_id: str
+    reason: str
+    evidence_retained_as: tuple[str, ...]
+    enablement_requirements: tuple[str, ...]
+    source_ids: tuple[str, ...]
+
+@dataclass(frozen=True)
+class ExternalHandoff:
+    handoff_id: str
+    tool: str
+    purpose: str
+    status: HandoffStatus
+    reason: str
+    required_input: tuple[str, ...]
+    required_provenance_fields: tuple[str, ...]
+    limitations: tuple[str, ...]
+
+@dataclass(frozen=True)
+class MobilomeReport:
+    source_file: str
+    source_sha256: str
+    tool_version: str
+    catalog_version: str
+    inference_version: str
+    resources: tuple[CatalogResource, ...]
+    parameters: AnalysisParameters
+    annotation_provenance: tuple[AnnotationProvenance, ...]
+    inventory: tuple[RepliconInventory, ...]
+    replicons: tuple[RepliconAssessment, ...]
+    cross_record_hypotheses: tuple[MobilomeHypothesis, ...]
+    disabled_aggregate_rules: tuple[DisabledAggregateRule, ...]
+    handoffs: tuple[ExternalHandoff, ...]
+    limitations: tuple[str, ...]
~~~

### 6.1 Contract invariants

- JSON absence is null. TSV absence is an empty cell. Do not use sentinel strings such as -, NA, or NO_LOCUS.
- The strand is -1, 0, 1, or null, matching the typed parser. Do not coerce unknown to 0.
- start and end are one-based and inclusive, matching GenBankFeature.start/end.
  segments preserves ordered one-based inclusive compound parts and is the
  authority for origin-spanning features.
- segments contains one (start, end) pair for a simple located feature and all
  ordered join_segments for a compound feature. An unlocatable feature has an
  empty segment tuple and is retained as evidence but excluded from distance
  rules with an explicit limitation.
- When a raw location cannot be interpreted, report start, end, and
  biological_length are null rather than leaking the canonical model's defensive
  1/0 fallbacks as biological coordinates.
- wraps_origin is derived with record length and topology; min(start)/max(end) is not enough.
- A hit groups all reasons for one record index, feature index, and marker ID.
- Every record classification reason records its policy rule, matched
  field/value/text/pattern, source IDs, and source-feature index when the reason
  came from a source feature.
- record_index and feature_index inherit read_genbank() identity and are
  one-based; never renumber after --include filtering.
- A reason ID additionally contains field, matcher ordinal, and pattern ordinal. IDs are stable within one parsed input and catalog version, not across reannotation or catalog reorderings.
- No sequence or translation is copied into the report. A user can recover it from the source using record and feature identity.
- pseudo_feature_count counts every feature where is_pseudo is true;
  pseudo_cds_count is the subset whose feature type is CDS; and
  pseudogene_feature_count counts feature type pseudogene. These fields are
  explicit so a pseudogene-type observation cannot disappear behind a CDS-only
  count.
- Every tuple is deterministically sorted:
  - inventory by record_index;
  - features by feature_index;
  - hits by record_index, feature_index, marker_id;
  - reasons by field, matcher ordinal, pattern ordinal;
  - hypotheses by record index, rule_id, hypothesis_id;
  - source IDs and facets lexically.
- A record-local hypothesis has a subject participant. A cross-record hypothesis
  has explicit role-bearing participants such as target and helper. Consumers
  never parse opaque hit or hypothesis IDs to recover record participation.
- tool_version, catalog_version, inference_version, and all resource hashes are independent fields.
- A required/renamed field, enum semantic change, coordinate change, or
  scientific-status reinterpretation requires gbparse.mobilome.v2. Marker/source
  content changes bump catalog_version; thresholds, component logic, or
  configured wording changes bump inference_version, even when report schema
  remains v1.
- The report does not include wall-clock time, temporary paths, absolute paths, object repr strings, or platform-specific newlines.

## 7. Packaged database schemas

### 7.1 Evidence-strength semantics

Strength is a small ordinal describing annotation evidence specificity:

| Strength | Meaning | Examples | Prohibited interpretation |
|---|---|---|---|
| 3 | Explicit typed structure or recognized structured expert cross-reference | explicit oriT feature; curated AMRFinder identifier with tool/database provenance | Not 95 percent confidence; not proof of activity |
| 2 | Exact curated gene/product correspondence with bounded ambiguity | exact toxN gene plus matching product pattern | Not experimental validation |
| 1 | Ambiguous product/note keyword or broad family correspondence | Zot-like product; generic replication protein | Not determinant, mechanism, or phenotype |

Per-field ceilings are enforced:

- note cannot receive strength 3;
- generic product cannot receive strength 3;
- strength values do not add across repeated reasons;
- max_strength is a display convenience, not a statistical score;
- a pseudogene or partial feature remains observable but cannot satisfy a required functional component by default;
- a future structured-evidence adapter must not accept imported evidence without
  tool and database provenance; version 0.6.0 implements no such adapter.

### 7.2 markers.yaml

~~~yaml
schema_version: 1
catalog_version: "1.0.0"

facets:
  - id: replication_origin
    title: Annotated replication-origin evidence
  - id: replication_initiation
    title: Replication-initiation protein evidence
  - id: mobility_orit
    title: Origin-of-transfer evidence
  - id: mobility_relaxase
    title: Relaxase evidence
  - id: mobility_t4cp
    title: Type IV coupling-protein evidence
  - id: mobility_mpf
    title: Mating-pair-formation evidence
  - id: ta_toxin
    title: Toxin-antitoxin toxin evidence
  - id: ta_antitoxin
    title: Toxin-antitoxin antitoxin evidence
  - id: partition_maintenance
    title: Partition or plasmid-maintenance annotation evidence
  - id: recombination
    title: Recombination/integration annotation evidence
  - id: transposition
    title: Transposition annotation evidence
  - id: phage_integration
    title: Phage integration/excision annotation evidence
  - id: phage_packaging
    title: Phage packaging annotation evidence
  - id: phage_structure
    title: Phage structural annotation evidence
  - id: phage_lysis
    title: Phage lysis annotation evidence
  - id: phage_replication
    title: Phage replication annotation evidence
  - id: phage_regulation
    title: Phage regulation/communication annotation evidence
  - id: crispr_array
    title: CRISPR repeat-array annotation evidence
  - id: crispr_cas
    title: Cas annotation evidence
  - id: amr_candidate
    title: Antimicrobial-resistance annotation candidate
  - id: vf_candidate
    title: Virulence-associated annotation candidate

markers:
  - id: explicit_origin_of_transfer
    label: Explicitly annotated origin of transfer
    facets: [mobility_orit]
    feature_types: [oriT, regulatory]
    matchers:
      - field: feature_type
        mode: casefold_exact
        patterns: [oriT]
        negative_patterns: []
        strength: 3
      - field: regulatory_class
        mode: casefold_exact
        patterns: [oriT]
        negative_patterns: []
        strength: 3
    requires_non_pseudo: false
    requires_complete: false
    sources: [bakta-origin-docs, insdc-feature-table]
    interpretation: Supports the presence of an explicitly annotated oriT feature.
    limitations:
      - Annotation does not prove that the site is cognate with an observed relaxase.
      - Annotation does not demonstrate transfer.

  - id: generic_replication_origin
    label: Generic annotated replication origin
    facets: [replication_origin]
    feature_types: [rep_origin]
    matchers:
      - field: feature_type
        mode: casefold_exact
        patterns: [rep_origin]
        negative_patterns: []
        strength: 3
    requires_non_pseudo: false
    requires_complete: false
    sources: [insdc-feature-table]
    interpretation: Supports a generic annotated replication-origin feature.
    limitations:
      - Does not distinguish oriC from oriV unless an explicit qualifier does so.

  - id: toxn
    label: ToxN-family toxin annotation
    facets: [ta_toxin]
    feature_types: [CDS, pseudogene]
    matchers:
      - field: gene
        mode: casefold_exact
        patterns: [toxN]
        negative_patterns: []
        strength: 2
      - field: product
        mode: regex
        patterns:
          - "(?i)\\bToxN(?:-family)?\\b"
        negative_patterns: []
        strength: 2
    requires_non_pseudo: true
    requires_complete: true
    sources: [blower-2012-toxin]
    interpretation: Supports a ToxN-family toxin annotation.
    limitations:
      - Does not establish toxin activity.
~~~

Allowed matcher fields are explicit and closed:

~~~text
feature_type
gene
product
note
db_xref
inference
function
mobile_element_type
regulatory_class
rpt_type
~~~

Exact gene matching is whole-value and case-insensitive after trimming. Regex is permitted only where declared and compiled at database-load time. The database must never search one concatenated gene/product/note string.

### 7.3 provenance.yaml

~~~yaml
schema_version: 1
provenance_version: "1.0.0"

sources:
  - id: blower-2012-toxin
    type: primary_article
    title: Identification and classification of bacterial Type III toxin-antitoxin systems
    year: 2012
    doi: 10.1093/nar/gks231
    pmid: "22434880"
    pmcid: PMC3401426
    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/
    supports:
      - ToxI is the RNA antitoxin in the Type III ToxIN system.
    limitations:
      - Annotation correspondence does not establish activity in a new genome.

  - id: yao-2020-hepn-mnt
    type: primary_article
    title: HEPN/MNT toxin-antitoxin neutralization by polyadenylylation
    year: 2020
    doi: 10.1093/nar/gkaa855
    pmid: "33045733"
    pmcid: PMC7641770
    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/
    supports:
      - In the studied system MntA transfers three AMP moieties from ATP to HepT Tyr104.
    limitations:
      - The biochemical result must not be generalized to every HEPN/MNT annotation.

  - id: bakta-origin-docs
    type: official_documentation
    title: Bakta documentation
    url: https://bakta.readthedocs.io/en/latest/
    accessed: "2026-08-22"
    supports:
      - Bakta versions documented at access time can annotate oriC, oriV, oriT, CRISPR, and pseudogenes.
    limitations:
      - Capabilities and output representation are version-sensitive.

  - id: insdc-feature-table
    type: official_documentation
    title: INSDC feature table definitions
    url: https://www.ncbi.nlm.nih.gov/genbank/feature_table/
    accessed: "2026-08-22"
    supports:
      - Controlled source qualifiers and feature keys are record evidence.
    limitations:
      - A qualifier records the submitted annotation, not independent validation.

  - id: local-record-classification-policy
    type: local_rule
    rule_ids:
      - source_plasmid_qualifier
      - source_chromosome_qualifier
      - bounded_description_plasmid_token
      - bounded_description_chromosome_token
    rationale: Preserve controlled declarations and label bounded free-text clues as tentative.
    adjudicator: genbank-parser maintainers
    date: "2026-08-22"
    negative_controls:
      - plasmid-like protein
      - chromosome partition protein
    limitations:
      - Name/description token rules are annotation clues, not replicon prediction.

  - id: local-disabled-aggregate-policy
    type: local_rule
    rule_ids:
      - pxo_like_annotation_pattern_candidate
      - phage_module_bearing_plasmid_candidate
    rationale: Retain component observations while disabling unvalidated aggregate classifiers.
    adjudicator: genbank-parser maintainers
    date: "2026-08-22"
    negative_controls:
      - repeated product labels on one feature
      - unrelated records with isolated phage annotations
    limitations:
      - Enablement requires separately reviewed positive and negative controls.
~~~

Required provenance fields depend on source type:

- primary_article: title, year, at least one stable DOI/PMID/PMCID/URL identifier, supports, and limitations;
- official_documentation: title, URL, accessed date, supports, and limitations;
- curated_database: name, release or database version, URL, taxonomic scope, thresholds if applicable, and limitations;
- local_rule: a nonempty rule_ids list, rationale, adjudicator, date, negative
  controls, and limitations; every referenced classification/inference rule ID
  must exist and every local rule used by output must be referenced.

### 7.4 inference.yaml

inference.yaml declares policy, thresholds, wording, sources, and component mappings. Python implements typed predicates. Do not invent a general expression language or evaluate arbitrary YAML expressions.

~~~yaml
schema_version: 1
inference_version: "1.0.0"

components:
  replication_candidate:
    operator: any
    facets: [replication_initiation, replication_origin]
  mobilizable_core:
    operator: all
    facets: [mobility_orit, mobility_relaxase]
    min_distinct_features:
      mobility_orit: 1
      mobility_relaxase: 1
  helper_machinery:
    operator: all
    facets: [mobility_t4cp, mobility_mpf]
    min_distinct_features:
      mobility_t4cp: 1
      mobility_mpf: 2

rules:
  - id: possible_helper_dependent_mobilization
    kind: mobility
    required_components:
      target_record: [mobilizable_core]
      distinct_helper_record: [helper_machinery]
    distinct_features:
      target_record: true
      distinct_helper_record: true
    status: tentative
    wording: Possible helper-dependent mobilization
    limitations:
      - Compatibility is untested.
      - A cognate relaxase-oriT relationship is untested.
      - Record taxonomic and same-cell co-affiliation are untested.
      - Expression and transfer are untested.
      - Co-transfer of both records is not implied.
    sources: [smillie-2010-mobility]

  - id: type_iii_toxin_antitoxin_pair_candidate
    kind: toxin_antitoxin
    required_marker_ids: [toxn, toxi]
    same_record: true
    distinct_features: true
    max_circular_gap_bp: 2000
    status: tentative
    wording: Type III ToxIN annotation-pair candidate
    limitations:
      - Proximity is a versioned engineering heuristic, not proof of an operon.
      - Expression and toxin-antitoxin activity are untested.
    sources: [blower-2012-toxin]

disabled_aggregate_rules:
  - id: pxo_like_annotation_pattern_candidate
    reason: No v1 ordinal/synteny/null-model contract has been validated.
    evidence_retained_as: pxo_numbered_product_annotation
    enablement_requirements:
      - Curated marker-group and ordinal metadata.
      - Exact forward/reverse/circular-rotation synteny semantics.
      - Gap and interruption policy.
      - Committed positive and unrelated negative controls.
      - Sequence-homology provenance and reviewed operating threshold.
    sources: [anjum-2014-pxo1, local-disabled-aggregate-policy]

  - id: phage_module_bearing_plasmid_candidate
    reason: No v1 distinct-module operating threshold has been validated.
    evidence_retained_as:
      - phage_integration
      - phage_packaging
      - phage_structure
      - phage_lysis
      - phage_replication
      - phage_regulation
    enablement_requirements:
      - Committed positive and unrelated negative controls.
      - Versioned minimum distinct-module policy.
      - Source-declared plasmid and eligible replication evidence.
    sources: [local-disabled-aggregate-policy]
~~~

The loader rejects:

- duplicate facet, marker, rule, component, or source IDs;
- unknown fields, modes, facets, components, marker IDs, or source IDs;
- missing/unknown component operators; any means at least one listed facet and
  all means every listed facet, while each role's component-ID list is all-of;
- nonpositive or unknown-facet min_distinct_features entries, and configurations
  whose required count can be satisfied only by reusing the same feature;
- invalid regexes;
- strengths outside 1 through 3;
- missing or non-Boolean requires_non_pseudo/requires_complete policies;
- a protein marker that scans CDS but omits pseudogene feature type, because
  pseudogene observations must be retained even when ineligible;
- empty pattern, matcher, interpretation, limitations, or required-component lists;
- a matcher strength exceeding its field ceiling;
- a functional component that permits pseudogenes without an explicit override;
- an inference rule with no required evidence;
- a rule that allows the same hit to satisfy components declared distinct;
- a disabled aggregate rule duplicated in executable rules, or missing its
  reason, retained-evidence mapping, or enablement requirements;
- source-type records missing mandatory provenance;
- catalog, provenance, and inference major-version disagreement;
- custom markers without their matching provenance and inference resources;
- a packaged or custom resource whose exact bytes cannot be SHA-256 hashed
  deterministically.

### 7.5 Custom database policy

The CLI exposes one experimental --database-dir option rather than three independently mixable file options. The directory must contain exactly the expected resource names and compatible versions:

~~~text
markers.yaml
provenance.yaml
inference.yaml
~~~

All three are loaded and hashed together. This prevents a custom marker set from being silently paired with packaged citations or inference wording. MobilomeDatabase.source_paths is a tuple of concrete resolved Path objects for custom files only; it is empty for packaged importlib.resources Traversables. The input and all custom database files are protected from output collisions without pretending a zipped or abstract package resource is an operating-system path.

## 8. Replicon inventory and classification logic

Classification describes how the input annotates a record; it is not plasmid prediction.

### 8.1 Evidence order

1. Controlled source qualifiers:
   - /plasmid supports plasmid;
   - /chromosome supports chromosome.
2. Explicit record annotation fields, if a recognized controlled value is present.
3. A bounded whole-token record name or description statement is tentative.
4. Topology, record length, coverage, GC, gene content, pXO labels, and replication markers never classify a record.

If explicit evidence supports both chromosome and plasmid, classification is unknown/conflicting even if one source is stronger. Preserve both reasons. An undeclared circular contig is unknown. A source-declared linear plasmid remains plasmid. Missing topology is unknown rather than linear.

Controlled-qualifier reasons cite insdc-feature-table. Bounded
name/description-token reasons cite local-record-classification-policy and use
the exact declared rule ID and pattern. Substrings such as plasmid-like protein
or chromosome partition protein are negative controls, not record labels.

Illustrative logic:

~~~text
diff --git a/src/genbank_parser/mobilome/replicons.py b/src/genbank_parser/mobilome/replicons.py
new file mode 100644
--- /dev/null
+++ b/src/genbank_parser/mobilome/replicons.py
@@
+def classify_replicon(record: GenBankRecord) -> RepliconClassification:
+    reasons = collect_record_classification_evidence(record)
+    plasmid = tuple(r for r in reasons if r.interpreted_as == "plasmid")
+    chromosome = tuple(r for r in reasons if r.interpreted_as == "chromosome")
+
+    if plasmid and chromosome:
+        return RepliconClassification(
+            label="unknown",
+            status="conflicting",
+            supporting_evidence=(),
+            conflicting_evidence=plasmid + chromosome,
+            limitations=("Conflicting record-level declarations were retained.",),
+        )
+    if has_controlled_source_reason(plasmid):
+        return supported("plasmid", plasmid)
+    if has_controlled_source_reason(chromosome):
+        return supported("chromosome", chromosome)
+    if plasmid:
+        return tentative("plasmid", plasmid)
+    if chromosome:
+        return tentative("chromosome", chromosome)
+    return insufficient_unknown()
~~~

Inventory every source feature, but normalize source qualifiers into sorted tuples. Never retain raw BioPython objects in the public report.

## 9. Field-aware evidence scanning

### 9.1 Scanner algorithm

~~~text
diff --git a/src/genbank_parser/mobilome/scanner.py b/src/genbank_parser/mobilome/scanner.py
new file mode 100644
--- /dev/null
+++ b/src/genbank_parser/mobilome/scanner.py
@@
+def scan_feature(feature, database, min_evidence):
+    hits_by_marker = {}
+    for marker in database.markers:
+        if feature.type.casefold() not in marker.feature_types_casefold:
+            continue
+        for matcher_index, matcher in enumerate(marker.matchers):
+            for value in values_for_declared_field(feature, matcher.field):
+                for pattern_index, pattern in enumerate(matcher.patterns):
+                    matched_text = match_declared_value(
+                        value=value,
+                        mode=matcher.mode,
+                        pattern=pattern,
+                        negative_patterns=matcher.negative_patterns,
+                    )
+                    if matched_text is None:
+                        continue
+                    reason = make_reason(
+                        feature=feature,
+                        marker=marker,
+                        field=matcher.field,
+                        value=value,
+                        matched_text=matched_text,
+                        pattern=pattern,
+                        matcher_index=matcher_index,
+                        pattern_index=pattern_index,
+                        strength=matcher.strength,
+                        eligible=(
+                            matcher.strength >= min_evidence
+                            and (not marker.requires_non_pseudo or not feature.is_pseudo)
+                            and (not marker.requires_complete or not feature.is_partial)
+                        ),
+                    )
+                    hits_by_marker.setdefault(
+                        marker.id, HitBuilder(marker=marker, feature=feature)
+                    ).add_exact_reason(reason)
+    return finalize_sorted_hits(hits_by_marker)
~~~

### 9.2 Matching rules

- Evaluate each qualifier value separately; do not flatten a list into one string.
- Apply Unicode-safe trimming and casefolding for case-insensitive exact comparison.
- Gene symbols use casefold_exact, not substring matching.
- Product and note regexes must use bounded patterns and explicit negative patterns where ambiguity is known.
- Structured db_xref/inference evidence is evaluated before free text, but all unique reasons remain in the hit.
- A replication-relaxation fusion may support both replication and relaxase markers.
- One feature may support many facets. Catalog order never changes category membership.
- Deduplicate only an identical tuple of marker, feature, field, mode, matched value, matched text, pattern, strength, and source IDs.
- Scan all marker-declared feature types, including source, regulatory, rep_origin, oriT, repeat_region, ncRNA, mobile_element, misc_feature, and CDS.
- Keep pseudo and partial observations. A marker's requires_non_pseudo and
  requires_complete policy determines reason eligibility; inference still
  enforces component-level requirements.
- Capture every raw hit even below --min-evidence. Mark eligibility explicitly.
- Do not multiply evidence because the same annotation text appears in gene, product, and note.

### 9.3 Threshold behavior

--min-evidence has values 1, 2, or 3:

- it does not change parsing, inventory, replicon classification, resource provenance, or raw evidence retention;
- a reason is eligible when its strength is at least the threshold and the
  marker's requires_non_pseudo/requires_complete policy allows it;
- a hit is eligible when at least one reason is eligible;
- inference uses only eligible hits;
- canonical JSON and TSV retain both eligible and below-threshold reasons;
- text shows eligible evidence first and, when present, a clearly labeled below-threshold appendix;
- summary records raw_hits, eligible_hits, raw_reasons, and eligible_reasons.

## 10. Inference policy

### 10.1 General inference invariants

Every rule must:

1. require a nonempty eligible evidence set;
2. list supporting hit IDs;
3. list missing components;
4. list conflicting hit IDs or state that no configured conflict was observed;
5. carry limitations and source IDs;
6. state whether distinct features and/or the same record are required;
7. use only configured wording;
8. return insufficient instead of inventing a negative biological conclusion;
9. never treat absence of an annotation as evidence that a biological function is absent.

Pseudogenized candidates are reported as disruption observations. They do not establish that the record is non-autonomous, a satellite, degenerating, or at an evolutionary transition.

### 10.2 Replication assessment

Version 0.6.0 reports:

- explicit generic or specific annotated origin features;
- replication-initiation protein candidates;
- partition and maintenance candidates;
- eligible and pseudogenized candidates;
- mutually incompatible or ambiguous annotations as conflicts;
- unresolved mechanism when no validated mechanism-specific rule exists.

Version 0.6.0 does not emit theta, rolling-circle, strand-displacement, pXO-type mechanism, or Inc group solely from product names. A later version may add a mechanism rule only with sequence/domain requirements, reference-positive and reference-negative fixtures, explicit thresholds, and a schema/inference-version review.

Illustrative empty-evidence pseudocode after inference.py exists:

~~~text
diff --git a/src/genbank_parser/mobilome/inference.py b/src/genbank_parser/mobilome/inference.py
--- a/src/genbank_parser/mobilome/inference.py
+++ b/src/genbank_parser/mobilome/inference.py
@@ -1,2 +1,22 @@
-if all(feature.is_pseudo for feature in replication_features):
-    status = "degenerate/satellite"
+if eligible_replication_features:
+    status = assess_supported_replication_components(
+        eligible_replication_features
+    )
+elif raw_replication_features and all(
+    hit.feature.is_pseudo for hit in raw_replication_features
+):
+    status = "insufficient"
+    limitations.append(
+        "Observed replication candidates are pseudogenized; autonomy is unresolved."
+    )
+elif raw_replication_features:
+    status = "insufficient"
+    limitations.append(
+        "Annotated replication candidates were observed but are below the "
+        "configured eligibility policy; autonomy is unresolved."
+    )
+else:
+    status = "insufficient"
+    limitations.append(
+        "No annotated replication candidate was observed; autonomy is unresolved."
+    )
~~~

The pseudogene branch deliberately inspects raw candidates because pseudogenes
are not eligible functional candidates by default. Raw evidence can explain why
the result is insufficient without satisfying a positive component.

### 10.2.1 Partition, maintenance, recombination, and transposition

Version 0.6.0 preserves these as observation-only facets:

- partition/maintenance markers from bounded gene/product correspondences;
- recombination/integration markers from typed mobile_element features and
  bounded CDS annotations;
- transposition markers from typed mobile_element qualifiers and bounded
  transposase/insertion-sequence annotations.

An integrase can retain both recombination and a phage-integration facet. A
mobile_element feature can retain recombination and transposition reasons.
Observations do not establish a complete partition system, stable inheritance,
element boundaries, mobility, integration activity, or transposition activity.
No v1 aggregate rule uses these facets unless the rule explicitly names and
tests the required distinct components.

### 10.3 Mobility assessment

Keep these components separate:

| Component | Evidence object | What it supports | What it does not support |
|---|---|---|---|
| oriT | explicit annotated oriT | candidate transfer origin | cognate relaxase, transfer |
| relaxase | curated relaxase annotation/domain | candidate DNA processing protein | active nicking, compatibility |
| T4CP | coupling-protein annotation | candidate substrate-coupling component | complete transfer apparatus |
| MPF/T4SS | multiple distinct structural components | candidate mating-pair formation machinery | complete/expressed/conjugative system |

Per-record output can say:

- annotated mobility components observed;
- mobilizable-core annotation candidate when distinct oriT and relaxase evidence meet the configured rule;
- conjugation-machinery annotation candidate when the configured minimum distinct T4CP/MPF components are present;
- incomplete component pattern with missing components.

It cannot say transfer-capable, conjugative plasmid, demonstrated mobilization, or successful transfer.

Cross-record output can say possible helper-dependent mobilization only when:

- the target record has eligible oriT and relaxase components on distinct features, unless a curated fusion rule explicitly permits one feature;
- a different record has the configured T4CP and MPF component set;
- record identities are unambiguous by record_index.

It must list compatibility, cognate oriT-relaxase relationship, record taxonomic
and assembly affiliation, physical coexistence in the same cell, expression,
transfer, and co-transfer as untested. Merely sharing one multi-record file is
not co-residence evidence.

### 10.4 Toxin-antitoxin neighborhoods

- Do not reuse find_operon_pairs(), because that helper deliberately assumes adjacent same-strand CDS features and ToxI can be noncoding RNA.
- Require curated toxin and antitoxin marker pairs on distinct features in the same record.
- Use segment-aware circular distance, not min/max bounds.
- Record orientation but do not require same strand unless a specific source-backed rule does.
- A default max_circular_gap_bp of 2000 is a versioned engineering heuristic and appears in output provenance.
- ToxN plus ToxI may yield a tentative Type III ToxIN annotation-pair candidate.
- HepT plus MntA may yield a tentative HEPN/MNT annotation-pair candidate. Do not hard-code an unsupported type number.
- A lone toxin or antitoxin is an observation, not a functional TA system.

### 10.5 CRISPR evidence

- CRISPR-array evidence comes from an actual repeat_region or other catalog-declared typed feature with explicit CRISPR qualifiers.
- Cas evidence comes from typed CDS/gene/product correspondences.
- A summary count in comments may be preserved only as annotation provenance/context and cannot create a CRISPR hit.
- Version 0.6.0 reports array and Cas observations separately. It has no combined
  CRISPR/Cas proximity hypothesis because no source-backed distance rule has been
  declared.
- The output does not claim active immunity, locus completeness, or plasmid
  targeting.

### 10.6 pXO-numbered annotation observations

Version 0.6.0 retains each bounded pXO-numbered product match with its field,
pattern, feature, and literature limitation, but emits no aggregate pXO-like
hypothesis. Counts alone do not define marker identity, ordinal extraction,
forward/reverse order, circular rotation, allowable gaps, or a null model.

The disabled-rule metadata records the controls needed for a later proposal.
One, two, three, duplicate, reordered, or clustered product labels all remain
observations. They are not diagnostic, a homology result, an ancestry assignment,
an incompatibility group, or a replication mechanism.

### 10.7 Phage-related, Zot-like, AimR-like, AMR, and VF evidence

- Assign phage-related observations to distinct integration, packaging,
  structure, lysis, replication, or regulation facets rather than one undifferentiated
  phage label.
- Version 0.6.0 emits no aggregate phage-module-bearing plasmid or phagemid
  hypothesis because no validated distinct-module operating threshold is
  available. The disabled-rule reason is reportable provenance.
- Zot-like product evidence is ambiguous and can be a phage-related and VF-candidate facet simultaneously. Product text alone never establishes a virulence factor.
- AimR/AimP/AimX correspondence may support an annotation candidate. It does not establish active arbitrium communication or cross-replicon regulation.
- A generic beta-lactamase-related domain, transporter, or efflux product is an AMR candidate only.
- Stronger AMR/VF evidence would require a future validated adapter for a
  structured expert result with tool version, database version, identifier,
  thresholds, and taxonomic limitations.
- Version 0.6.0 does not import such results. A later imported result would still
  be an annotation/tool result, not a phenotype.

### 10.8 External handoffs

The report may recommend but never silently run:

| Tool class | Purpose | Required provenance |
|---|---|---|
| MOB-suite or PlasmidFinder | plasmid typing/mobility comparison | tool version, database release, identity/coverage thresholds, taxonomic scope |
| AMRFinderPlus | structured AMR annotation | tool version, database version, organism option, method/threshold provenance |
| VFDB or curated VF workflow | virulence-associated comparison | database release, search method, thresholds, coverage, interpretation limits |
| Prophage/ICE caller | dedicated element-boundary analysis | tool/database version, caller thresholds, boundary uncertainty |
| HMM/domain/sequence workflow | mechanism-specific evidence | model/accession/version, cutoffs, alignment coverage |

Default status is not_run. Imported results require a separate validated adapter and must never be inferred from the existence of a local executable.

## 11. Circular spatial logic

Add small pure helpers to src/genbank_parser/spatial.py only when the first consumer is implemented:

~~~text
diff --git a/src/genbank_parser/spatial.py b/src/genbank_parser/spatial.py
--- a/src/genbank_parser/spatial.py
+++ b/src/genbank_parser/spatial.py
@@
+def intervening_gap_bp(
+    left_segments: Sequence[tuple[int, int]],
+    right_segments: Sequence[tuple[int, int]],
+) -> int:
+    """Return the non-negative number of bases strictly between features."""
+
+def circular_feature_distance_bp(
+    first_segments: Sequence[tuple[int, int]],
+    second_segments: Sequence[tuple[int, int]],
+    *,
+    record_length: int,
+) -> int:
+    """Return the minimum intervening gap in either direction on a circle."""
~~~

Required semantics:

- public start/end and segment coordinates are interpreted consistently with the
  canonical one-based, inclusive GenBankFeature properties;
- for two non-overlapping linear intervals, the gap is
  max(0, right_start - left_end - 1); adjacent features therefore have a gap of
  zero;
- overlapping features have a gap of zero;
- ordered compound segments are retained;
- origin-spanning join(271..300,1..30)-style features are not represented as a record-wide linear span;
- wraps_origin is detected by strand-aware modulo unwrapping of the ordered
  CompoundLocation parts, not by requiring a segment to touch coordinate 1 or
  record_length; join(900..950,20..50) on a 1,000-bp circular record is a
  required positive fixture;
- circular pair distance is the minimum nonnegative clockwise intervening gap
  across segment boundaries in both directions after overlap is checked;
- linear and circular records have separate tested behavior;
- missing or invalid record length is an explicit error;
- no TA inference uses crispr.interval_distance() unchanged.

## 12. Output contracts

### 12.1 Canonical JSON

~~~json
{
  "schema_version": "gbparse.mobilome.v1",
  "tool": "gbparse",
  "analysis": "mobilome",
  "tool_version": "0.6.0",
  "catalog_version": "1.0.0",
  "inference_version": "1.0.0",
  "input": {
    "file": "input.gbff",
    "sha256": "..."
  },
  "parameters": {
    "include": "all",
    "min_evidence": 1,
    "database_source": "packaged"
  },
  "provenance": {
    "resources": [
      {"name": "markers.yaml", "sha256": "..."},
      {"name": "provenance.yaml", "sha256": "..."},
      {"name": "inference.yaml", "sha256": "..."},
      {"name": "report.schema.json", "sha256": "..."}
    ],
    "annotation_pipelines": [
      {
        "record_indices": [1, 2, 3, 4, 5],
        "name": "Bakta",
        "version": null,
        "database_version": null,
        "evidence_fields": ["record_annotations"]
      }
    ]
  },
  "summary": {
    "records_scanned": 5,
    "records_reported": 4,
    "raw_hits": 0,
    "eligible_hits": 0,
    "raw_reasons": 0,
    "eligible_reasons": 0,
    "hypotheses": 0
  },
  "inventory": [],
  "replicons": [],
  "cross_record_hypotheses": [],
  "disabled_aggregate_rules": [],
  "handoffs": [],
  "limitations": []
}
~~~

Representative nested v1 shapes:

~~~json
{
  "inventory_member": {
    "record_index": 1,
    "record_id": "contig_1",
    "name": "contig_1",
    "description": "example",
    "length": 1000,
    "topology": "circular",
    "source_qualifiers": [
      {"key": "plasmid", "values": ["unnamed"]}
    ],
    "classification": {
      "label": "plasmid",
      "status": "supported",
      "supporting_evidence": [
        {
          "rule_id": "source_plasmid_qualifier",
          "source": "source_qualifier",
          "field": "plasmid",
          "value": "unnamed",
          "matched_text": "unnamed",
          "pattern": "qualifier-present",
          "source_feature_index": 1,
          "interpreted_as": "plasmid",
          "strength": 3,
          "source_ids": ["insdc-feature-table"]
        }
      ],
      "conflicting_evidence": [],
      "limitations": []
    },
    "feature_count": 12,
    "cds_count": 8,
    "pseudo_feature_count": 1,
    "pseudo_cds_count": 1,
    "pseudogene_feature_count": 0
  },
  "hit": {
    "hit_id": "r1:f4:mtoxn",
    "marker_id": "toxn",
    "marker_label": "ToxN-family toxin annotation",
    "facets": ["ta_toxin"],
    "feature": {
      "record_index": 1,
      "record_id": "contig_1",
      "feature_index": 4,
      "feature_type": "CDS",
      "locus_tag": "EX_0004",
      "gene": "toxN",
      "product": "ToxN-family toxin",
      "start": 101,
      "end": 450,
      "strand": 1,
      "segments": [[101, 450]],
      "location_operator": null,
      "biological_length": 350,
      "wraps_origin": false,
      "partial": false,
      "pseudo": false
    },
    "reasons": [
      {
        "reason_id": "r1:f4:mtoxn:gene:0:0",
        "field": "gene",
        "match_mode": "casefold_exact",
        "matched_value": "toxN",
        "matched_text": "toxN",
        "pattern": "toxN",
        "strength": 2,
        "eligible": true,
        "source_ids": ["blower-2012-toxin"]
      }
    ],
    "max_strength": 2,
    "eligible_for_inference": true
  },
  "cross_record_hypothesis": {
    "hypothesis_id": "h:possible_helper_dependent_mobilization:1:2",
    "rule_id": "possible_helper_dependent_mobilization",
    "kind": "mobility",
    "status": "tentative",
    "summary": "Possible helper-dependent mobilization",
    "participants": [
      {"role": "target", "record_index": 1, "record_id": "contig_1"},
      {"role": "helper", "record_index": 2, "record_id": "contig_2"}
    ],
    "supporting_hit_ids": [
      "r1:f6:mexplicit_origin_of_transfer",
      "r1:f7:mrelaxase",
      "r2:f10:mt4cp",
      "r2:f11:mmpf_component",
      "r2:f12:mmpf_component"
    ],
    "missing_components": [],
    "conflicting_hit_ids": [],
    "limitations": [
      "Record taxonomic and same-cell co-affiliation are untested."
    ],
    "source_ids": ["smillie-2010-mobility"]
  }
}
~~~

The labels inventory_member, hit, and cross_record_hypothesis above identify
member shapes; they are not extra top-level keys in a real report. The serializer
maps FeatureRef.is_partial/is_pseudo to JSON partial/pseudo and applies the same
names in TSV.

The packaged report.schema.json is the normative structural contract. It must:

- use JSON Schema draft 2020-12;
- set additionalProperties to false on normative objects;
- define required fields and nullable values explicitly;
- constrain strength to integers 1 through 3;
- constrain status, classification, topology, match mode, and handoff status enums;
- constrain participant roles to subject, target, or helper and require
  role-bearing participant objects; require at least two participant objects for
  a cross-record hypothesis;
- require each disabled aggregate rule to expose retained evidence, reason,
  source IDs, and nonempty enablement requirements;
- define unique stable IDs where practical;
- permit additive content only through a future schema version, not silent unversioned fields.

Runtime code does not add a mandatory jsonschema dependency. Phase 1 adds
jsonschema to the test extra, uses it to validate every golden and CLI JSON
payload against the packaged schema, and also keeps direct structural assertions
for readable failures. A pure validate_mobilome_report_semantics() function
enforces constraints JSON Schema cannot express generically, including unique
participant roles, distinct target/helper record indices, hit-ID references,
and consistency between summary counts and arrays; serializers call it before
emitting output.

### 12.2 Long-form TSV

A single wide summary row cannot represent inventory, evidence reasons, hypotheses, and handoffs without losing provenance. Use one normalized table with row_type:

~~~text
inventory
evidence
hypothesis
cross_record_hypothesis
disabled_aggregate_rule
handoff
~~~

Required columns, in fixed order:

~~~text
schema_version
row_type
record_index
record_id
record_name
record_description
record_length
topology
replicon_class
classification_status
feature_count
cds_count
pseudo_feature_count
pseudo_cds_count
pseudogene_feature_count
feature_index
feature_type
locus_tag
gene
product
start
end
strand
segments_json
wraps_origin
partial
pseudo
hit_id
marker_id
facets_json
reason_id
evidence_field
match_mode
matched_value
matched_text
pattern
strength
eligible
source_ids_json
hypothesis_id
hypothesis_rule_id
hypothesis_kind
hypothesis_status
participants_json
supporting_hit_ids_json
missing_components_json
conflicting_hit_ids_json
summary
limitations_json
disabled_rule_id
disabled_rule_reason
disabled_evidence_retained_as_json
disabled_enablement_requirements_json
handoff_id
handoff_tool
handoff_status
handoff_reason
handoff_provenance_json
~~~

Emit:

- one inventory row for every input record, even with zero hits;
- one evidence row per reason, not merely per hit;
- one row per replicon hypothesis;
- one row per cross-record hypothesis, with role-bearing record_index/record_id
  objects in participants_json;
- one row per disabled aggregate rule with its retained evidence, reason, sources,
  and enablement requirements;
- one row per handoff.

Collection cells use compact JSON, not comma or semicolon conventions. Use the csv module with tab delimiter, correct quoting for tabs/newlines/quotes, UTF-8, and newline="". Empty values remain empty cells.

### 12.3 Text

The text report is deterministic and human-oriented:

~~~text
Mobilome evidence report
Input and provenance
Record inventory
Replicon assessments
  Observed eligible evidence
  Tentative hypotheses
  Conflicts and missing components
  Below-threshold observations
Cross-record hypotheses
Disabled aggregate rules and retained evidence
External handoffs not run
Interpretation limits
~~~

Avoid decorative scores and categorical verdicts. Every hypothesis includes rule ID, support, missing evidence, and limitations. Zero-hit records remain visible.

### 12.4 Safe output writer

Default output is stdout. Supplying --output writes the file and emits no duplicate report to stdout.

~~~text
diff --git a/src/genbank_parser/mobilome/report.py b/src/genbank_parser/mobilome/report.py
new file mode 100644
--- /dev/null
+++ b/src/genbank_parser/mobilome/report.py
@@
+def write_mobilome_report(
+    report,
+    *,
+    source_path,
+    database_paths,
+    format_type,
+    output_path=None,
+):
+    verify_source_sha256(source_path, report.source_sha256)
+    rendered = serialize_mobilome_report(report, format_type)
+    if output_path is None:
+        return rendered
+    assert_not_same_path(output_path, (source_path, *database_paths))
+    atomic_write_text(output_path, rendered)
+    return ""
~~~

Protection requirements:

1. source_path and database_paths are mandatory keyword arguments; the writer
   has no unsafe empty protected-path default.
2. Verify that source_path still has the SHA-256 recorded in the report.
3. Resolve absolute normalized paths without requiring the output to exist.
4. For existing paths, use os.path.samefile where supported to catch symlink and hardlink aliases.
5. Protect the input and every custom database resource.
6. Create the parent directory if requested.
7. Create a temporary file in the destination directory.
8. Write UTF-8 with LF newlines, flush, and fsync.
9. Replace atomically with os.replace.
10. Remove the temporary file on failure.
11. Never create a default sidecar.
12. Verify input bytes are unchanged after every output-path test.

## 13. Public API and CLI

### 13.1 Python API

~~~text
diff --git a/src/genbank_parser/mobilome/__init__.py b/src/genbank_parser/mobilome/__init__.py
new file mode 100644
--- /dev/null
+++ b/src/genbank_parser/mobilome/__init__.py
@@
+def load_mobilome_database(
+    database_dir: str | Path | None = None,
+) -> MobilomeDatabase: ...
+
+def inventory_replicons(
+    document: GenBankDocument,
+) -> tuple[RepliconInventory, ...]: ...
+
+def scan_mobilome_features(
+    features: Iterable[GenBankFeature],
+    database: MobilomeDatabase,
+    *,
+    min_evidence: int = 1,
+) -> tuple[MobilomeHit, ...]: ...
+
+def infer_replicon_hypotheses(
+    replicon: RepliconInventory,
+    hits: Sequence[MobilomeHit],
+    database: MobilomeDatabase,
+) -> tuple[MobilomeHypothesis, ...]: ...
+
+def infer_cross_record_hypotheses(
+    assessments: Sequence[RepliconAssessment],
+    database: MobilomeDatabase,
+) -> tuple[MobilomeHypothesis, ...]: ...
+
+def analyze_mobilome(
+    filepath: str | Path,
+    *,
+    include: Literal["all", "chromosome", "plasmid", "unknown"] = "all",
+    min_evidence: int = 1,
+    database: MobilomeDatabase | None = None,
+) -> MobilomeReport:
+    """Return a pure annotation-supported mobilome evidence report."""
~~~

analyze_mobilome() must:

- validate include and min_evidence against their closed runtime sets before
  reading; Literal annotations alone are not runtime validation;
- parse once with read_genbank();
- hash the exact input bytes before and after parsing and fail if they differ, so
  the reported SHA-256 identifies the content actually analyzed;
- raise MobilomeInputError when no records are parsed;
- inventory and scan all records;
- infer all record-local and cross-record hypotheses;
- apply include only to the detailed replicons tuple;
- keep the complete inventory and records_scanned count;
- return data without printing or writing.

scan_mobilome_features() calls the same min_evidence runtime validator before
iterating its feature iterable, so direct API callers cannot bypass the 1/2/3
contract.

serialize_mobilome_report() and write_mobilome_report() likewise validate
format_type before serialization or filesystem work. Invalid API parameters
raise MobilomeParameterError; the CLI's argparse choices remain the earlier
user-facing guard.

Document:

~~~python
from genbank_parser.mobilome import analyze_mobilome

report = analyze_mobilome("sample.gbff", include="all", min_evidence=2)
~~~

Do not add a root-level import in genbank_parser/__init__.py. This avoids expanding the root API and avoids circular import risk around package metadata.

### 13.2 CLI

~~~diff
diff --git a/src/genbank_parser/cli.py b/src/genbank_parser/cli.py
--- a/src/genbank_parser/cli.py
+++ b/src/genbank_parser/cli.py
@@ -20,6 +20,8 @@
 from .locus import inspect_locus
-from .metadata import extract_metadata
 from .meor import analyze_meor
 from .meor.database import load_meor_database
 from .meor.report import serialize_report
+from .metadata import extract_metadata
+from .mobilome import MobilomeError, analyze_mobilome, load_mobilome_database
+from .mobilome.report import write_mobilome_report
 from .neighborhood import build_neighborhood, write_neighborhood
@@ -312,5 +314,28 @@
     p_meor.add_argument("--output", help="Output path (default: stdout)")
     p_meor.add_argument("--markers", help="Experimental custom marker YAML")
     p_meor.add_argument("--pathways", help="Experimental custom pathway YAML")

+    # 20. mobilome
+    p_mobilome = subparsers.add_parser(
+        "mobilome",
+        help="Build an annotation-supported replicon and mobilome evidence report",
+    )
+    p_mobilome.add_argument("input")
+    p_mobilome.add_argument(
+        "--format", choices=["text", "json", "tsv"], default="text"
+    )
+    p_mobilome.add_argument(
+        "--include",
+        choices=["all", "chromosome", "plasmid", "unknown"],
+        default="all",
+    )
+    p_mobilome.add_argument(
+        "--min-evidence", type=int, choices=[1, 2, 3], default=1
+    )
+    p_mobilome.add_argument(
+        "--database-dir",
+        help="Experimental directory containing a complete mobilome database",
+    )
+    p_mobilome.add_argument("--output")
+
     return parser
@@ -470,7 +495,28 @@
             if input_path == output_path:
                 parser.error("--output must not overwrite the input GenBank file")
             Path(args.output).write_text(rendered, encoding="utf-8", newline="")
         else:
             sys.stdout.write(rendered)

+    elif cmd == "mobilome":
+        try:
+            database = load_mobilome_database(args.database_dir)
+            report = analyze_mobilome(
+                args.input,
+                include=args.include,
+                min_evidence=args.min_evidence,
+                database=database,
+            )
+            rendered = write_mobilome_report(
+                report,
+                source_path=args.input,
+                database_paths=database.source_paths,
+                format_type=args.format,
+                output_path=args.output,
+            )
+        except MobilomeError as exc:
+            parser.error(str(exc))
+        if args.output is None:
+            sys.stdout.write(rendered)
+
     return 0
~~~

Error contract:

- user/input/database/output errors return a nonzero CLI status without a traceback;
- diagnostics go to stderr;
- stdout contains only a requested report;
- malformed GenBank, empty parsed document, invalid catalog, invalid threshold, and output collision are distinguishable messages;
- public database, analysis, and writer functions translate expected
  OSError/ValueError/parser failures into the appropriate MobilomeError subtype
  while preserving the original exception as the cause;
- existing commands and their help remain unchanged.

## 14. Package and CI diffs

### 14.1 Package data

~~~diff
diff --git a/pyproject.toml b/pyproject.toml
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -20,4 +20,5 @@
 [project.optional-dependencies]
 test = [
     "pytest>=8.0",
+    "jsonschema>=4.23,<5",
 ]
@@ -35,4 +36,9 @@
 [tool.setuptools.package-data]
-genbank_parser = ["rulesets/*.yaml", "data/meor/*.yaml"]
+genbank_parser = [
+    "rulesets/*.yaml",
+    "data/meor/*.yaml",
+    "data/mobilome/*.yaml",
+    "data/mobilome/*.json",
+]
-
+
 [tool.pytest.ini_options]
~~~

No new required runtime dependency is expected beyond Biopython and PyYAML.

### 14.2 Installed-wheel gate

Illustrative CI additions:

~~~diff
diff --git a/.github/workflows/ci.yml b/.github/workflows/ci.yml
--- a/.github/workflows/ci.yml
+++ b/.github/workflows/ci.yml
@@ -53,7 +53,11 @@
     - name: Verify installed API, data, CLI, and dependencies
       run: |
-        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; assert __version__ == '0.5.0'; assert load_meor_database().catalog_version == '1.1'; assert files('genbank_parser').joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file()"
+        .wheel-venv/bin/python -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.meor import load_meor_database; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.0'; assert load_meor_database().catalog_version == '1.1'; assert db.catalog_version == '1.0.0'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/meor/markers.yaml').is_file(); assert files('genbank_parser.rulesets').joinpath('mobilome.yaml').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
         .wheel-venv/bin/gbparse --help
         .wheel-venv/bin/gbparse meor tests/fixtures/meor_parity.gb --format json --output wheel-meor.json
-        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-meor.json')); assert p['tool_version'] == '0.5.0'; assert p['catalog_version'] == '1.1'"
+        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-meor.json')); assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.1'"
+        .wheel-venv/bin/gbparse mobilome --help
+        .wheel-venv/bin/gbparse mobilome tests/fixtures/mobilome_evidence.gb --format json --output wheel-mobilome.json
+        .wheel-venv/bin/python -c "import json; p=json.load(open('wheel-mobilome.json', encoding='utf-8')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.0'; assert p['catalog_version'] == '1.0.0'"
+        .wheel-venv/bin/gbparse discover tests/fixtures/mobilome_evidence.gb --ruleset mobilome
         .wheel-venv/bin/python -m pip check
~~~

Update every hard-coded 0.5.0 assertion in CI and tests in the release phase, including tests/test_meor_report.py. Retain Python 3.10 through 3.13 matrix coverage.

## 15. Synthetic fixture design

The fixtures are small, hand-adjudicated GenBank records. They must be readable and reviewable in Git, contain no private sequence, and validate through the existing parser before mobilome code exists.

### 15.1 mobilome_inventory.gb

Include separate records for:

- controlled /chromosome;
- controlled /plasmid on a circular record;
- definition-only tentative plasmid;
- undeclared circular unknown;
- source-declared linear plasmid;
- missing topology;
- contradictory /plasmid and /chromosome evidence;
- duplicate record IDs distinguished by record_index;
- multiple source features;
- zero features;
- generic rep_origin;
- explicit oriV representation;
- explicit oriT representation;
- regulatory and source features that the old parser omitted.

### 15.2 mobilome_evidence.gb

Include:

- /pseudo CDS;
- /pseudogene="unitary" CDS;
- pseudogene feature type;
- partial CDS;
- unknown strand;
- compound joined feature;
- origin-spanning compound feature;
- replication-relaxation fusion;
- exact gene and ambiguous product variants;
- multiple qualifier values;
- duplicate identical reasons;
- negative-pattern products;
- structured db_xref plus generic product;
- repeat_region with explicit CRISPR qualifier;
- Cas CDS;
- ToxN CDS and ToxI noncoding feature;
- HepT and MntA pair;
- lone toxin and lone antitoxin;
- Zot-like, integrase, AimR-like, generic beta-lactamase-domain, and efflux examples.

### 15.3 mobilome_inference.gb

Use multiple records to cover:

- empty marker set;
- one pseudogenized Rep;
- generic RepA/RepB only;
- TPR only;
- one, two, and three distinct numbered pXO-like markers;
- duplicated pXO marker on one feature;
- conflicting replication candidates;
- each incomplete mobility component subset;
- complete record-local mobility annotation pattern;
- helper machinery and candidate target on distinct records;
- helper and target components on the same record;
- incomplete helper;
- phage single marker;
- multiple distinct phage modules without plasmid evidence;
- multiple phage modules with source-declared plasmid and replication evidence;
- linear and circular TA pairs at distance boundaries.

Do not use the ignored real PF_NNT_reoriented.gbff as a committed fixture. It remains an optional local integration check.

## 16. Detailed testing matrix

### 16.1 Parser and inventory

- Every synthetic fixture parses with read_genbank().
- Feature counts include source, regulatory, rep_origin, oriT, repeat_region, ncRNA, mobile_element, misc_feature, and CDS.
- /pseudo, /pseudogene, and pseudogene feature type all set is_pseudo.
- Inventory separately asserts all-pseudo, pseudo-CDS, and pseudogene-feature
  counts so each representation is visible.
- Compound segment order and biological length are preserved.
- Origin-spanning features set wraps_origin correctly.
- Unknown strand remains null; strand 0 remains 0.
- An unlocatable feature retains evidence with null start/end/biological_length,
  empty segments, and no spatial inference.
- Duplicate record IDs do not collide.
- Empty parsed input raises MobilomeInputError.
- Malformed GenBank exits cleanly in CLI.

### 16.2 Replicon classification

- Controlled chromosome and plasmid declarations are supported.
- Definition/name-only evidence is tentative.
- Undeclared circular remains unknown.
- Source-declared linear plasmid remains plasmid.
- Missing topology remains unknown.
- Conflicting declarations are unknown/conflicting and retain all reasons.
- Classification reasons retain rule/pattern/source provenance and the source
  feature index where applicable.
- Topology, size, Rep labels, and pXO labels never classify a record.
- --include filters details but not inventory or scanning.

### 16.3 Database validation

- Packaged resources load through importlib.resources from source and wheel.
- Valid versions and exact-byte SHA-256 hashes are deterministic.
- Duplicate IDs fail.
- Unknown facets, fields, match modes, sources, markers, rules, and components fail.
- Invalid regex fails at load time.
- Strength 0 and 4 fail.
- Field-strength ceiling violation fails.
- Empty pattern/matcher/interpretation/limitation/rule requirement fails.
- A rule that reuses one hit where distinct features are required fails or is rejected at evaluation.
- Missing or mismatched custom resource files fail closed.
- Source records missing type-specific provenance fail.
- Unsupported schema major version fails with a clear message.

### 16.4 Scanner

- Exact, casefold exact, and regex modes behave distinctly.
- Gene matching is whole-value, not substring.
- Field evidence reports exact qualifier field and original value.
- matched_text and configured pattern are both preserved.
- Multi-value qualifiers produce traceable reasons.
- Structured xref and product reasons coexist.
- One fusion can have several markers/facets.
- No first-match break.
- Exact duplicate reasons deduplicate; distinct reasons do not.
- Negative patterns suppress only the intended matcher.
- Non-CDS markers are scanned.
- Typed mobile_element and bounded CDS annotations can retain partition,
  recombination, and transposition facets; a multi-role integrase keeps every
  applicable reason.
- Pseudo/partial state appears in FeatureRef.
- A pseudogene feature-type product can match a protein marker as raw evidence
  but is ineligible when requires_non_pseudo is true.
- Pseudogene hit is raw evidence but not functional support by default.
- min-evidence changes eligibility, not raw hit retention.
- Repeated scans are byte-for-byte deterministic after serialization.

### 16.5 Inference adversarial tests

- Empty evidence is insufficient.
- One pseudogenized Rep does not yield satellite, non-autonomous, degenerate, or evolutionary-transition wording.
- RepA, RepB, TPR, pXO, and replication-relaxation alone do not choose a mechanism.
- Conflicting replication observations remain visible.
- Generic rep_origin is not silently called oriV.
- Each mobility component subset lists exact missing components.
- One MPF hit is insufficient for helper_machinery; the initial policy requires
  one T4CP feature and at least two distinct MPF-feature markers.
- Required distinct components cannot be satisfied by duplicated reasons on one feature.
- A curated fusion exception is explicit and rule-specific.
- Cross-record helper rule uses different record indices.
- Cross-record target/helper roles are explicit structured participants; no
  consumer must decode IDs to identify records.
- Invalid participant roles, duplicate roles, and duplicate target/helper record
  indices are rejected.
- Helper evidence says possible helper-dependent mobilization only.
- Same-file records do not satisfy taxonomic, assembly-affiliation, or same-cell
  co-residence evidence.
- No obligate co-transfer/co-transmission wording.
- ToxI can be a non-CDS antitoxin feature.
- TA pairing uses circular distance and distinct features.
- Same-strand adjacency is not required by a generic TA rule.
- CRISPR comment alone does not create a hit.
- Explicit repeat-region/Cas features create annotation candidates only.
- One, two, three, duplicate, reordered, and clustered pXO-numbered products do
  not yield an aggregate pXO-like hypothesis in v1.
- Single or multiple phage-module observations do not yield an aggregate
  phage-module-bearing-plasmid or phagemid hypothesis in v1.
- Disabled aggregate rules and their enablement requirements remain traceable
  in report provenance.
- Generic beta-lactamase-domain and efflux products are not confirmed AMR.
- Zot-like product is not confirmed VF.
- AimR-like annotation does not prove communication.

### 16.6 Scientific wording tests

Maintain an explicit forbidden-claim list over all text output and configured inference wording:

~~~text
obligate co-mobilisation
obligate co-mobilization
obligate co-transfer
guaranteed transfer
confirmed conjugative
diagnostic pXO
proves pXO ancestry
confirmed phagemid
confirmed resistance
confirmed virulence
active communication
satellite plasmid
evolutionary transition
theta replication
rolling-circle replication
~~~

The last two are forbidden unless a future, independently reviewed mechanism rule supplies the required evidence. Tests should assert acceptable exact wording too, not merely absence of prohibited phrases.

### 16.7 Serialization and writer

- JSON key order, schema version, enums, nulls, source hashes, and no absolute path.
- Disabled aggregate rules round-trip through JSON and TSV with reason, retained
  evidence, sources, and enablement requirements.
- JSON UTF-8 and final newline.
- TSV exact columns and row types.
- One inventory row per zero-hit record.
- One evidence row per reason.
- Collection cells are valid compact JSON.
- Embedded tab, quote, newline, and non-ASCII data round-trip.
- Text sections and cautious wording match golden semantics.
- Repeated serialization is byte-identical.
- Different --include values do not alter records_scanned or complete inventory.
- Different --min-evidence values retain raw counts and reasons.
- Input basename and SHA-256 are stable.
- Output cannot equal input or a custom database file.
- Symlink/hardlink aliases are rejected where the platform supports them.
- Atomic replace leaves no partial destination or temporary file after simulated failure.
- Input bytes remain unchanged.
- Source mutation between analysis and write is rejected by SHA-256 comparison.
- Input mutation during parse is rejected rather than producing mismatched
  provenance.
- --output keeps stdout empty.

### 16.8 CLI/API regressions

- gbparse mobilome --help.
- Default text stdout.
- JSON and TSV stdout.
- Each --include value.
- API-invalid include and format_type fail before input read or output creation.
- Direct scanner API rejects min_evidence outside 1, 2, or 3 before consuming
  the iterable.
- Each --min-evidence value and invalid values.
- Packaged and complete custom database.
- Incomplete custom database failure.
- Output path success and collision failure.
- Malformed and empty input with no traceback.
- Public import from genbank_parser.mobilome.
- Root package import remains unchanged.
- Existing gbparse discover --ruleset mobilome behavior remains unchanged.
- Existing parser, model, query, region, neighborhood, operon, CRISPR, MEOR, GFF, and visualization tests remain green.

### 16.9 Optional real-input checks

When the ignored PF_NNT_reoriented.gbff is locally available, a non-CI script or opt-in test may assert:

- five records are inventoried;
- four records have explicit source-level plasmid evidence;
- all 21 canonical pseudogene CDS features are retained;
- the chromosome's generic rep_origin is not called an oriV or a replication mechanism;
- no record is classified from circularity alone;
- no generic product produces a categorical AMR, VF, phagemid, satellite, or transfer claim.

The optional check must skip clearly when the file is absent and must never stage the file.

### 16.10 Line-ending reproducibility

The checkout currently uses core.autocrlf=true and has no repository
.gitattributes. Phase 0 adds narrow LF rules so resource SHA-256 values and
byte-identical goldens do not vary by platform:

~~~gitattributes
src/genbank_parser/data/mobilome/*.yaml text eol=lf
src/genbank_parser/data/mobilome/*.json text eol=lf
tests/fixtures/mobilome_*.gb text eol=lf
tests/golden/mobilome_v1.* text eol=lf
~~~

Do not add a repository-wide line-ending rule and do not renormalize unrelated
files. Test both Git blob bytes and working-tree resource bytes where the
contract requires exact hashes.

## 17. Reusable phase gate

Every implementation phase has a focused test list below. Before its commit, it also runs:

~~~powershell
$ErrorActionPreference = "Stop"
$env:MPLBACKEND = "Agg"
$env:MPLCONFIGDIR = "C:\tmp\gbparse-mobilome-mpl"

python -m pytest -q -p no:cacheprovider
if ($LASTEXITCODE -ne 0) { throw "Full pytest gate failed" }
python -m compileall -q src scripts tests
if ($LASTEXITCODE -ne 0) { throw "compileall gate failed" }
python -m pip check
if ($LASTEXITCODE -ne 0) { throw "Dependency check failed" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --check
if ($LASTEXITCODE -ne 0) { throw "Worktree diff check failed" }
~~~

Focused quality checks use new Python paths plus explicitly cleaned legacy paths:

~~~powershell
python -m ruff check --no-cache <new-paths-and-named-cleaned-legacy-paths>
if ($LASTEXITCODE -ne 0) { throw "Scoped Ruff check failed" }
python -m ruff format --check --no-cache <new-paths-and-format-clean-legacy-paths>
if ($LASTEXITCODE -ne 0) { throw "Scoped Ruff format check failed" }
~~~

If Phase 0 establishes a working mypy version:

~~~powershell
python -m mypy --cache-dir C:\tmp\gbparse-mobilome-mypy src\genbank_parser\mobilome
if ($LASTEXITCODE -ne 0) { throw "Pinned mypy gate failed" }
~~~

Otherwise record the mypy crash as advisory and do not claim type-check success.

Before each commit:

~~~powershell
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' status --short --branch
if ($LASTEXITCODE -ne 0) { throw "Could not inspect pre-staging status" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --cached --quiet
if ($LASTEXITCODE -ne 0) {
    throw "Index is not empty before phase staging; inspect existing staged work"
}

$allowedPhasePaths = @(
    # Populate only with the exact paths listed in this phase's staging scope.
)
if ($allowedPhasePaths.Count -eq 0) { throw "Phase allowlist is empty" }

git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' add -- $allowedPhasePaths
if ($LASTEXITCODE -ne 0) { throw "Exact-path staging failed" }
$cachedPaths = @(
    git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --cached --name-only
)
if ($LASTEXITCODE -ne 0) { throw "Could not enumerate staged paths" }
if ($cachedPaths.Count -eq 0) { throw "Phase staged no changes" }
$unexpectedPaths = @($cachedPaths | Where-Object { $_ -notin $allowedPhasePaths })
if ($unexpectedPaths.Count -ne 0) {
    throw "Unexpected staged paths: $($unexpectedPaths -join ', ')"
}

$lfControlledPaths = @(
    $cachedPaths | Where-Object {
        $_ -match '^(src/genbank_parser/data/mobilome/.*\.(yaml|json)|tests/fixtures/mobilome_.*\.gb|tests/golden/mobilome_v1\..*)$'
    }
)
foreach ($controlledPath in $lfControlledPaths) {
    $attributes = @(git check-attr text eol -- $controlledPath)
    if ($LASTEXITCODE -ne 0) { throw "Could not read attributes for $controlledPath" }
    if (-not ($attributes -match ': text: set$')) {
        throw "text attribute is not set for $controlledPath"
    }
    if (-not ($attributes -match ': eol: lf$')) {
        throw "eol=lf is not set for $controlledPath"
    }
    $indexHash = (git rev-parse --verify ":$controlledPath").Trim()
    if ($LASTEXITCODE -ne 0) { throw "No staged blob for $controlledPath" }
    $worktreeHash = (git hash-object --no-filters -- $controlledPath).Trim()
    if ($LASTEXITCODE -ne 0 -or $indexHash -ne $worktreeHash) {
        throw "Index/worktree byte mismatch for $controlledPath"
    }
    $controlledBytes = [IO.File]::ReadAllBytes(
        (Join-Path (Get-Location) $controlledPath)
    )
    if (13 -in $controlledBytes) {
        throw "CR byte found in LF-controlled path $controlledPath"
    }
}

$cachedPaths
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --cached --check
if ($LASTEXITCODE -ne 0) { throw "Staged diff check failed" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --cached
if ($LASTEXITCODE -ne 0) { throw "Could not render staged diff for review" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' commit -m "<phase message>"
if ($LASTEXITCODE -ne 0) { throw "Phase commit failed" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' diff --cached --quiet
if ($LASTEXITCODE -ne 0) { throw "Index is not empty after phase commit" }
git -c safe.directory='D:/W/Skills Claude/genbank-feature-parser' status --short --branch
if ($LASTEXITCODE -ne 0) { throw "Could not inspect post-commit status" }
~~~

A phase is not complete merely because its focused tests pass. It is complete only when its focused tests, reusable full gate, scoped lint, diff review, and exact-path staged review all pass. Do not push as part of any phase.

## 18. Implementation phases

### Phase 0 — Freeze adjudicated fixtures and contracts

Goal: make corrected parser behavior and scientific wording requirements executable before adding production code.

Tasks:

- Re-check branch, baseline SHA, dirty paths, Python versions, and current full-suite result.
- Add the three small synthetic fixtures described above.
- Add narrow .gitattributes LF rules for only the new mobilome resources,
  fixtures, and goldens; do not renormalize existing paths.
- Add tests/test_mobilome_fixtures.py that tests only current read_genbank()/model behavior.
- Add an in-test contract table for forbidden wording and required corrected labels.
- Record the intended report, catalog, and inference versions in test constants.
- Confirm no fixture copies a private or ignored biological dataset.
- Reproduce the mypy 2.3.0 crash. Test a supported pinned version in isolation; record a working invocation or explicitly keep mypy advisory.
- Do not add mobilome production imports yet.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_fixtures.py tests\test_parser.py tests\test_model.py
if ($LASTEXITCODE -ne 0) { throw "Phase 0 focused pytest gate failed" }
~~~

Commit:

~~~text
test(mobilome): add adjudicated typed-parser fixtures
~~~

Exact staging scope:

~~~text
tests/fixtures/mobilome_inventory.gb
tests/fixtures/mobilome_evidence.gb
tests/fixtures/mobilome_inference.gb
tests/test_mobilome_fixtures.py
.gitattributes
~~~

Exit criteria:

- fixtures demonstrate both pseudogene qualifier forms, compound/origin locations, non-CDS evidence, ambiguous negatives, conflicts, and empty evidence;
- all expectations are canonical-parser facts or explicit future contracts, not parity with the regex prototype;
- the commit is green without production mobilome code.

### Phase 1 — Add typed models and a validated packaged database

Goal: establish versioned data contracts and fail-closed scientific resources before matching or inference.

Tasks:

- Add mobilome/models.py with immutable contracts and explicit enums.
- Add mobilome/database.py with importlib.resources loading and custom-directory support.
- Add a minimal mobilome/__init__.py so the subpackage is explicit; export only
  database/model names needed by Phase 1 and expand it when orchestration lands.
- Add markers.yaml, provenance.yaml, inference.yaml, and report.schema.json.
- Add jsonschema>=4.23,<5 to the test extra only and validate packaged schema
  examples without changing runtime dependencies.
- Refresh the implementation environment with python -m pip install -e ".[test]"
  after editing pyproject.toml and before running schema tests. This is an
  environment prerequisite, not a staged production artifact.
- Implement SHA-256 hashes over the exact packaged or custom resource bytes.
- Validate all IDs, references, regexes, strengths, field ceilings, source-type fields, policy constraints, and cross-file versions.
- Add package-data patterns.
- Add tests/test_mobilome_database.py.
- Add no CLI entry and no inference behavior.

Focused gate:

~~~powershell
python -m pip install -e ".[test]"
if ($LASTEXITCODE -ne 0) { throw "Could not refresh the test extra" }
python -m pytest -q -p no:cacheprovider tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 1 focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 1 scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 1 scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add validated evidence catalog and models
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/models.py
src/genbank_parser/mobilome/database.py
src/genbank_parser/mobilome/__init__.py
src/genbank_parser/data/mobilome/markers.yaml
src/genbank_parser/data/mobilome/provenance.yaml
src/genbank_parser/data/mobilome/inference.yaml
src/genbank_parser/data/mobilome/report.schema.json
tests/test_mobilome_database.py
pyproject.toml
~~~

Exit criteria:

- packaged and complete custom resources load deterministically;
- every scientific assertion is linked to a source or declared local policy;
- invalid databases fail before an input is analyzed;
- resource loading works without adding a runtime dependency.

### Phase 2 — Add inventory, classification, and field-aware scanning

Goal: produce complete typed observations with no biological synthesis.

Tasks:

- Add mobilome/replicons.py.
- Add mobilome/scanner.py.
- Build FeatureRef only from GenBankFeature.
- Implement conservative record classification and topology normalization.
- Scan declared fields independently and preserve all reasons and facets.
- Implement deterministic IDs and sorting.
- Keep raw and eligible evidence.
- Add tests/test_mobilome_replicons.py and tests/test_mobilome_scanner.py.
- Add fixture factories to tests/conftest.py only if they are shared.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_replicons.py tests\test_mobilome_scanner.py tests\test_model.py tests\test_parser.py
if ($LASTEXITCODE -ne 0) { throw "Phase 2 focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_replicons.py tests\test_mobilome_scanner.py
if ($LASTEXITCODE -ne 0) { throw "Phase 2 scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_replicons.py tests\test_mobilome_scanner.py
if ($LASTEXITCODE -ne 0) { throw "Phase 2 scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): classify replicons and typed feature evidence
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/replicons.py
src/genbank_parser/mobilome/scanner.py
tests/conftest.py
tests/test_mobilome_replicons.py
tests/test_mobilome_scanner.py
~~~

Omit `tests/conftest.py` from the Phase 2 allowlist when that optional shared
fixture file is unchanged.

Exit criteria:

- all records and declared feature types are represented;
- one feature can retain multiple labels and multiple exact reasons;
- no topology or marker-content plasmid prediction appears;
- pseudogene and compound-location behavior is inherited, not reimplemented.

### Phase 3 — Add conservative record-local inference

Goal: assess replication and mobility components without making mechanism, autonomy, transfer, AMR, VF, or phenotype claims.

Tasks:

- Add mobilome/inference.py.
- Implement explicit typed predicates for replication observations and mobility component matrices.
- Require nonempty evidence before all aggregate logic.
- Enforce non-pseudo and distinct-feature rules.
- Preserve ambiguous alternatives, conflicts, and missing components.
- Emit unresolved replication mechanism in v1.
- Add tests/test_mobilome_inference.py covering all adversarial cases.
- Add exact required and forbidden wording assertions.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_inference.py tests\test_mobilome_scanner.py
if ($LASTEXITCODE -ne 0) { throw "Phase 3 focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 3 scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 3 scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add cautious replication and mobility assessments
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/inference.py
tests/test_mobilome_inference.py
~~~

Exit criteria:

- empty evidence, one pseudo Rep, generic RepA/TPR/pXO, and conflicts have calibrated output;
- mobility components remain separate;
- inference rows name support, missing evidence, conflicts, sources, and limitations;
- forbidden scientific claims are absent.

### Phase 4A — Add circular TA logic and typed CRISPR observations

Goal: add only the record-local features that need spatial semantics, keeping
CRISPR observations separate in schema v1.

Tasks:

- Add intervening_gap_bp() and circular_feature_distance_bp() with explicit
  one-based inclusive coordinate semantics.
- Add tests/test_mobilome_spatial.py for linear, circular, adjacent,
  overlapping, origin-spanning, compound, and unlocatable cases.
- Add TA-pair rules that allow non-CDS ToxI and do not require same-strand CDS
  adjacency.
- Add actual repeat-region and Cas observation markers without treating summary
  comments as feature evidence or creating a combined CRISPR/Cas hypothesis.
- Replace the existing successive-pair zip() in spatial.py with
  itertools.pairwise() while that file is touched, resolving its pre-existing
  RUF007 diagnostic without formatting or rewriting unrelated code.
- Extend tests/test_mobilome_inference.py for TA and CRISPR negative/positive
  observations.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_spatial.py tests\test_mobilome_inference.py tests\test_neighborhood.py tests\test_operons.py tests\test_patch2.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4A focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome src\genbank_parser\spatial.py tests\test_mobilome_spatial.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4A scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome src\genbank_parser\spatial.py tests\test_mobilome_spatial.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4A scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add circular TA and CRISPR evidence
~~~

Exact staging scope:

~~~text
src/genbank_parser/spatial.py
src/genbank_parser/mobilome/inference.py
src/genbank_parser/data/mobilome/inference.yaml
src/genbank_parser/data/mobilome/markers.yaml
src/genbank_parser/data/mobilome/provenance.yaml
tests/test_mobilome_spatial.py
tests/test_mobilome_inference.py
~~~

Exit criteria:

- circular boundary cases pass with canonical coordinates;
- TA pairing uses distinct typed features and reports the heuristic gap;
- CRISPR arrays and Cas features are separate observations;
- comments alone do not create CRISPR evidence;
- the touched spatial.py is Ruff-check and Ruff-format clean without broad churn.

### Phase 4B — Add cross-record mobility and external handoffs

Goal: add the only v1 cross-record hypothesis and make external analyses explicit
not-run recommendations.

Tasks:

- Add possible helper-dependent mobilization across distinct target/helper
  records using the declared component operators and minimum distinct features.
- Add structured target/helper participants to every cross-record hypothesis.
- Add missing compatibility, cognate oriT/relaxase, expression, coexistence,
  transfer, and co-transfer limitations.
- Add not-run ExternalHandoff records with required future provenance fields.
- Extend inference and schema tests for participant roles, incomplete helpers,
  feature reuse, and forbidden wording.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_inference.py tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4B focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_inference.py tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4B scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_inference.py tests\test_mobilome_database.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4B scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add cross-record mobility handoffs
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/models.py
src/genbank_parser/mobilome/inference.py
src/genbank_parser/data/mobilome/inference.yaml
src/genbank_parser/data/mobilome/provenance.yaml
src/genbank_parser/data/mobilome/report.schema.json
tests/test_mobilome_database.py
tests/test_mobilome_inference.py
~~~

Exit criteria:

- target/helper roles are structured and unambiguous;
- one T4CP plus fewer than two distinct MPF features is incomplete;
- the only cross-record wording is possible helper-dependent mobilization;
- handoffs are not_run and list the provenance a future adapter would require.

### Phase 4C — Add provenance-bounded AMR/VF and phage observations

Goal: preserve useful ambiguous observations without enabling unsupported
aggregate classifiers.

Tasks:

- Add separate phage integration, packaging, structure, lysis, replication, and
  regulation facets.
- Add observation-only partition/maintenance, recombination, and transposition
  facets over typed mobile_element and bounded CDS annotations.
- Add product-level Zot, AimR/AimP/AimX, generic AMR, and VF candidate markers
  with field-specific strength ceilings and source limitations.
- Add pXO-numbered product observations.
- Add disabled-rule metadata for pXO-like and phage-module-bearing-plasmid
  aggregates, including exact future enablement requirements.
- Assert that one, two, three, duplicate, reordered, and clustered pXO markers
  never create an aggregate hypothesis in v1.
- Assert that one or many phage-module observations never create a phagemid or
  phage-module-bearing-plasmid hypothesis in v1.
- Add no imported-results adapter.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_database.py tests\test_mobilome_scanner.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4C focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_database.py tests\test_mobilome_scanner.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4C scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_database.py tests\test_mobilome_scanner.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 4C scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add provenance-bounded annotation candidates
~~~

Exact staging scope:

~~~text
src/genbank_parser/data/mobilome/markers.yaml
src/genbank_parser/data/mobilome/provenance.yaml
src/genbank_parser/data/mobilome/inference.yaml
tests/test_mobilome_database.py
tests/test_mobilome_scanner.py
tests/test_mobilome_inference.py
~~~

Exit criteria:

- every ambiguous observation retains exact reasons and sources;
- pXO/phage aggregate rules are visibly disabled and never evaluated;
- generic AMR/VF, Zot, and AimR observations remain candidates only;
- no phenotype, ancestry, communication, mechanism, or element-identity claim is
  emitted.

### Phase 5 — Add deterministic serializers and atomic output

Goal: make the report contract lossless, stable, and safe before exposing a command.

Tasks:

- Add mobilome/report.py.
- Implement canonical JSON, normalized TSV, and human text serializers.
- Implement validate_mobilome_report_semantics() and call it from every
  serializer before rendering.
- Add atomic writer and protected-path collision logic.
- Add the three goldens.
- Generate Phase 5 goldens with the then-current package metadata version
  (0.5.0 on the recorded baseline); do not hard-code the future release version
  inside production serializers.
- Add tests/test_mobilome_report.py.
- Test nulls, Unicode, tabs/newlines, zero-hit inventory, deterministic hashes, and stdout behavior.
- Validate output against report.schema.json without adding a required runtime dependency.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_report.py tests\test_mobilome_inference.py
if ($LASTEXITCODE -ne 0) { throw "Phase 5 focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\mobilome tests\test_mobilome_report.py
if ($LASTEXITCODE -ne 0) { throw "Phase 5 scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_report.py
if ($LASTEXITCODE -ne 0) { throw "Phase 5 scoped Ruff format check failed" }
~~~

Commit:

~~~text
feat(mobilome): add versioned deterministic reports
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/report.py
tests/golden/mobilome_v1.json
tests/golden/mobilome_v1.tsv
tests/golden/mobilome_v1.txt
tests/test_mobilome_report.py
~~~

Exit criteria:

- repeated renders are byte-identical;
- JSON is the lossless canonical representation;
- TSV retains one row per reason and every zero-hit record;
- no output can overwrite input or a custom database resource;
- an output path causes no duplicate stdout.

### Phase 6 — Add orchestrator, public subpackage API, and gbparse mobilome

Goal: expose the tested pipeline through one native API and CLI without disturbing existing commands.

Tasks:

- Extend mobilome/__init__.py with analyze_mobilome() and curated exports.
- Register the mobilome parser and dispatch in cli.py.
- Add --format, --include, --min-evidence, --database-dir, and --output.
- Normalize only the touched CLI import block, resolving its pre-existing I001
  while adding the mobilome imports; do not run broad Ruff formatting on cli.py.
- Add dedicated tests/test_mobilome_api.py and tests/test_mobilome_cli.py so
  legacy tests/test_cli.py remains an unchanged regression.
- Preserve discover --ruleset mobilome exactly.
- Do not add a root package export or wrapper executable.

Focused gate:

~~~powershell
python -m pytest -q -p no:cacheprovider tests\test_mobilome_api.py tests\test_mobilome_cli.py tests\test_cli.py tests\test_mobilome_fixtures.py tests\test_mobilome_database.py tests\test_mobilome_replicons.py tests\test_mobilome_scanner.py tests\test_mobilome_inference.py tests\test_mobilome_spatial.py tests\test_mobilome_report.py tests\test_query_region_diff.py
if ($LASTEXITCODE -ne 0) { throw "Phase 6 focused pytest gate failed" }
python -m ruff check --no-cache src\genbank_parser\cli.py src\genbank_parser\mobilome tests\test_mobilome_api.py tests\test_mobilome_cli.py
if ($LASTEXITCODE -ne 0) { throw "Phase 6 scoped Ruff check failed" }
python -m ruff format --check --no-cache src\genbank_parser\mobilome tests\test_mobilome_api.py tests\test_mobilome_cli.py
if ($LASTEXITCODE -ne 0) { throw "Phase 6 scoped Ruff format check failed" }
python -m genbank_parser.cli mobilome --help
if ($LASTEXITCODE -ne 0) { throw "Phase 6 mobilome help smoke failed" }
~~~

Commit:

~~~text
feat(cli): add gbparse mobilome
~~~

Exact staging scope:

~~~text
src/genbank_parser/mobilome/__init__.py
src/genbank_parser/cli.py
tests/test_mobilome_api.py
tests/test_mobilome_cli.py
~~~

Exit criteria:

- API is pure and importable;
- CLI stdout/stderr/error contracts pass;
- include does not truncate analysis context;
- no existing command output changes.

### Phase 7 — Add scientific reference, user documentation, and skill guidance

Goal: make interpretation limits and provenance discoverable at the same time as the feature.

Tasks:

- Add docs/mobilome_evidence_reference.md from the corrected correspondence ledger.
- Document every evidence-strength tier, hypothesis status, source, threshold, and limitation.
- Add README CLI/API examples and explain discover versus mobilome.
- Update SKILL.md routing, CLI table, output guidance, and epistemic boundaries.
- Rewrite the already long SKILL frontmatter description compactly rather than appending beyond the validator's 1,024-character limit.
- Add a migration note: the standalone prototype remains available for audit comparison but is not the native parser.
- Do not copy the old erroneous biological prose.

Focused gate:

~~~powershell
python "C:\Users\LIHTR-UA\.codex\skills\.system\skill-creator\scripts\quick_validate.py" .
if ($LASTEXITCODE -ne 0) { throw "Phase 7 skill validation failed" }
$wordingHits = @(
    rg -n "obligate|diagnostic|confirmed|satellite|theta|rolling-circle|ToxI|HepT|MntA|oriV|oriT" README.md SKILL.md docs\mobilome_evidence_reference.md
)
if ($LASTEXITCODE -gt 1) { throw "Phase 7 wording scan failed" }
$wordingHits
python -m pytest -q -p no:cacheprovider tests\test_mobilome_fixtures.py tests\test_mobilome_database.py tests\test_mobilome_replicons.py tests\test_mobilome_scanner.py tests\test_mobilome_inference.py tests\test_mobilome_spatial.py tests\test_mobilome_report.py tests\test_mobilome_api.py tests\test_mobilome_cli.py tests\test_cli.py
if ($LASTEXITCODE -ne 0) { throw "Phase 7 focused pytest gate failed" }
~~~

Every search hit must be manually adjudicated; a blanket word deletion is not the goal.

Commit:

~~~text
docs(mobilome): document provenance and interpretation bounds
~~~

Exact staging scope:

~~~text
docs/mobilome_evidence_reference.md
README.md
SKILL.md
~~~

Exit criteria:

- docs match implemented names and exact output;
- scientific corrections are source-linked;
- limitations are centralized but evidence-specific boundaries remain beside each claim;
- skill validation passes.

### Phase 8 — Prepare and validate the local 0.6.0 release candidate

Goal: prove the feature from fresh local artifacts and prepare, but do not
publish, a 0.6.0 release-candidate commit. The established CI split remains:
editable-install tests on Python 3.10 through 3.13 and installed-wheel validation
on Python 3.12.

Tasks:

- Bump project version from 0.5.0 to 0.6.0.
- Refresh the editable test installation after the version bump so
  importlib.metadata reports 0.6.0 before golden/full-suite checks.
- Update tests/test_meor_report.py and CI hard-coded version expectations in the same commit.
- Regenerate the JSON, TSV, and text mobilome goldens so their tool_version
  fields move from the Phase 5 development baseline to 0.6.0.
- Add the final dated 0.6.0 entry to changelogs.md; Phase 8 is the sole changelog
  owner.
- Extend the wheel job for mobilome resources, public API, py.typed, CLI help, JSON generation, schema fields, discover compatibility, and pip check.
- Preserve editable-install CI on Python 3.10 through 3.13 and the installed
  wheel job on Python 3.12; do not claim wheel coverage on every interpreter.
- Build a clean wheel and sdist.
- Install the wheel in an isolated environment, outside the checkout import path.
- Check wheel contents for all four mobilome resources.
- Run installed gbparse mobilome and gbparse discover smokes.
- Preserve the existing separate optional-visualization wheel job.
- Inspect the final diff, staged paths, and worktree.

Focused gate:

This phase's focused gate is the complete local release gate:

~~~powershell
$ErrorActionPreference = "Stop"
$releaseToken = [Guid]::NewGuid().ToString("N")
$releaseRoot = Join-Path "C:\tmp" ("gbparse-mobilome-release-" + $releaseToken)
$artifactDir = Join-Path $releaseRoot "artifacts"
$wheelVenv = Join-Path $releaseRoot "wheel-venv"
$smokeDir = Join-Path $releaseRoot "smoke"
$smokeOutput = Join-Path $smokeDir "mobilome.json"
$meorOutput = Join-Path $smokeDir "meor.json"
$env:MPLBACKEND = "Agg"
$env:MPLCONFIGDIR = Join-Path $releaseRoot "mpl"
$fixturePath = (Resolve-Path "tests\fixtures\mobilome_evidence.gb").Path
$meorFixturePath = (Resolve-Path "tests\fixtures\meor_parity.gb").Path

$null = New-Item -ItemType Directory -Path $artifactDir, $smokeDir, $env:MPLCONFIGDIR

python -m pip install -e ".[test]"
if ($LASTEXITCODE -ne 0) { throw "Could not refresh 0.6.0 editable metadata/test dependencies" }
python -m pytest -q -p no:cacheprovider
if ($LASTEXITCODE -ne 0) { throw "Full pytest gate failed" }
python -m compileall -q src scripts tests
if ($LASTEXITCODE -ne 0) { throw "compileall gate failed" }
python -m pip check
if ($LASTEXITCODE -ne 0) { throw "Local dependency check failed" }
python -m build --outdir $artifactDir
if ($LASTEXITCODE -ne 0) { throw "Artifact build failed" }

$wheelFiles = @(Get-ChildItem -LiteralPath $artifactDir -Filter "genbank_parser-0.6.0-*.whl")
$sdistFiles = @(Get-ChildItem -LiteralPath $artifactDir -Filter "genbank_parser-0.6.0.tar.gz")
if ($wheelFiles.Count -ne 1) { throw "Expected exactly one 0.6.0 wheel" }
if ($sdistFiles.Count -ne 1) { throw "Expected exactly one 0.6.0 sdist" }

python -m venv $wheelVenv
if ($LASTEXITCODE -ne 0) { throw "Fresh wheel environment creation failed" }
$wheelPython = Join-Path $wheelVenv "Scripts\python.exe"
$wheelGbparse = Join-Path $wheelVenv "Scripts\gbparse.exe"
& $wheelPython -m pip install --force-reinstall $wheelFiles[0].FullName
if ($LASTEXITCODE -ne 0) { throw "Wheel installation failed" }

Push-Location $smokeDir
try {
    & $wheelPython -c "from importlib.resources import files; from genbank_parser import __version__; from genbank_parser.mobilome import load_mobilome_database; db=load_mobilome_database(); root=files('genbank_parser'); assert __version__ == '0.6.0'; assert db.catalog_version == '1.0.0'; assert root.joinpath('py.typed').is_file(); assert root.joinpath('data/mobilome/markers.yaml').is_file(); assert root.joinpath('data/mobilome/provenance.yaml').is_file(); assert root.joinpath('data/mobilome/inference.yaml').is_file(); assert root.joinpath('data/mobilome/report.schema.json').is_file()"
    if ($LASTEXITCODE -ne 0) { throw "Installed API/resource smoke failed" }
    & $wheelGbparse --help
    if ($LASTEXITCODE -ne 0) { throw "Installed CLI help failed" }
    & $wheelGbparse mobilome $fixturePath --format json --output $smokeOutput
    if ($LASTEXITCODE -ne 0) { throw "Installed mobilome smoke failed" }
    & $wheelPython -c "import json, pathlib; p=json.loads(pathlib.Path('mobilome.json').read_text(encoding='utf-8')); assert p['schema_version'] == 'gbparse.mobilome.v1'; assert p['tool_version'] == '0.6.0'"
    if ($LASTEXITCODE -ne 0) { throw "Installed mobilome JSON check failed" }
    $discoverJson = & $wheelGbparse discover $fixturePath --ruleset mobilome --format json
    if ($LASTEXITCODE -ne 0 -or -not $discoverJson) { throw "Installed discover smoke failed" }
    & $wheelGbparse meor $meorFixturePath --format json --output $meorOutput
    if ($LASTEXITCODE -ne 0) { throw "Installed MEOR smoke failed" }
    & $wheelPython -c "import json, pathlib; p=json.loads(pathlib.Path('meor.json').read_text(encoding='utf-8')); assert p['schema_version'] == 'gbparse.meor.v1'; assert p['tool_version'] == '0.6.0'"
    if ($LASTEXITCODE -ne 0) { throw "Installed MEOR JSON check failed" }
    & $wheelPython -m pip check
    if ($LASTEXITCODE -ne 0) { throw "Installed dependency check failed" }
}
finally {
    Pop-Location
}
~~~

Commit:

~~~text
chore(release): prepare genbank-parser 0.6.0
~~~

Exact staging scope:

~~~text
pyproject.toml
tests/test_meor_report.py
.github/workflows/ci.yml
tests/golden/mobilome_v1.json
tests/golden/mobilome_v1.tsv
tests/golden/mobilome_v1.txt
changelogs.md
~~~

Exit criteria:

- the local full suite and fresh artifact checks pass;
- installed-wheel import and resource checks pass outside the source tree;
- installed CLI emits valid gbparse.mobilome.v1 JSON;
- discover and MEOR wheel checks remain green;
- no unrelated or real biological input is staged;
- the repository is a local release candidate but not pushed or tagged.

Remote CI is a separate gate after explicit publication authorization. If the
release-candidate commit is then pushed, require the Python 3.10–3.13 editable
matrix, Python 3.12 wheel job, and optional visualization wheel job to pass
before any tag or release. A local Phase 8 commit does not claim remote CI
success.

## 19. Post-release migration, separately authorized

Only after 0.6.0 is validated against synthetic fixtures and optional real inputs:

1. Compare native and standalone outputs as an audit, not equality target.
2. Confirm every useful standalone category has a native evidence facet, a deliberate handoff, or a documented rejection.
3. Document that discrepancies caused by prototype defects are expected.
4. Decide whether to mark mobilome-parser deprecated, archive it, or retain it as historical reference.
5. Perform any standalone-repository change in that repository and in a separate commit.
6. Request explicit authorization before pushing, tagging, releasing, archiving, or deleting anything.

## 20. Acceptance criteria

### 20.1 Architecture

- [ ] read_genbank() is the sole parser.
- [ ] No regex parsing of INSDC structure is introduced.
- [ ] Existing GenBankFeature location and pseudogene semantics are reused.
- [ ] discover --ruleset mobilome remains compatible.
- [ ] The root package API remains intentionally small.
- [ ] New logic is split among models, database, scanner, replicons, inference, report, and orchestration.

### 20.2 Scientific and provenance behavior

- [ ] Every hit records exact field, value, matched text, pattern, strength, sources, and feature identity.
- [ ] A feature can retain multiple markers/facets.
- [ ] Empty evidence yields insufficient.
- [ ] /pseudo, /pseudogene, and pseudogene feature type are honored.
- [ ] Generic Rep/pXO/TPR text does not choose a mechanism.
- [ ] Generic rep_origin is not silently called oriV.
- [ ] Mobility components remain separate.
- [ ] Partition/maintenance, recombination, and transposition observations are
  retained without activity, completeness, or boundary claims.
- [ ] Cross-record evidence uses only possible helper-dependent mobilization wording.
- [ ] ToxI is represented as the RNA antitoxin of Type III ToxIN.
- [ ] HEPN/MNT correspondence reflects the cited ATP-dependent HepT Tyr104 modification and is not generalized beyond annotation correspondence.
- [ ] CRISPR evidence comes from typed features, not summary comments.
- [ ] pXO-numbered and phage-module evidence remains traceable, while both
  aggregate classifiers are explicitly disabled in v1.
- [ ] Zot, AimR, AMR, VF, and phage observations are not promoted to phenotype or identity claims.
- [ ] External handoffs are not silently run and carry version/database/threshold limitations.

### 20.3 Output and safety

- [ ] JSON validates as gbparse.mobilome.v1.
- [ ] Tool, catalog, inference, annotation, input, and resource provenance are present.
- [ ] JSON null and TSV empty values replace sentinels.
- [ ] TSV is normalized and lossless at reason level.
- [ ] Every record appears in inventory, including zero-hit records.
- [ ] --include is a presentation filter only.
- [ ] --min-evidence does not delete raw evidence.
- [ ] Serialization is deterministic and platform-independent.
- [ ] Input/custom database collisions are rejected.
- [ ] Writes are atomic and input bytes remain unchanged.
- [ ] --output produces no duplicate stdout.

### 20.4 Testing and distribution

- [ ] All focused phase tests pass.
- [ ] Full pytest passes under Agg with a writable Matplotlib config.
- [ ] New Python paths are Ruff-check and Ruff-format clean; touched legacy
  paths pass Ruff check after named diagnostics are resolved, with pre-existing
  full-file format debt neither hidden nor broadly rewritten.
- [ ] mypy status is truthfully reported and only blocking with a working pinned toolchain.
- [ ] compileall, pip check, and git diff --check pass.
- [ ] Skill validator passes.
- [ ] Wheel contains py.typed and all mobilome resources.
- [ ] Installed API and CLI work outside the checkout.
- [ ] After a separately authorized push, editable-install CI passes on Python
  3.10, 3.11, 3.12, and 3.13; the installed-wheel and optional-visualization
  jobs pass on Python 3.12.

### 20.5 Git and scope

- [ ] Each implementation phase is one independently green commit.
- [ ] Only exact phase paths are staged.
- [ ] Conflict markers are checked before staging.
- [ ] Existing unrelated untracked files are preserved.
- [ ] Real/ignored GBFF data is not committed.
- [ ] No standalone deletion is mixed into target integration.
- [ ] No push, tag, archive, or release occurs without separate authorization.

## 21. Final implementation sequence

The required order is:

~~~text
Baseline re-check
  -> Phase 0 fixtures/contracts
  -> Phase 1 models/database
  -> Phase 2 inventory/scanner
  -> Phase 3 record-local inference
  -> Phase 4A circular TA/CRISPR observations
  -> Phase 4B cross-record mobility/handoffs
  -> Phase 4C bounded AMR/VF/phage observations
  -> Phase 5 serializers/writer
  -> Phase 6 API/CLI
  -> Phase 7 documentation/skill
  -> Phase 8 local release-candidate validation
  -> separately authorized remote CI
  -> separately authorized migration/tag/release/publication
~~~

Do not collapse the phases to save commits. Later logic depends on earlier contracts, and each gate is designed to make scientific and software regressions independently reviewable.
