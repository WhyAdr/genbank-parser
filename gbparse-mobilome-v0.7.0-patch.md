# Mobilome Patch 0.7.0 — Deferred Review Items

Status: implementation plan derived from `mobilome-deferred-review.md`. This
patch resolves the actionable v0.7.0 scope from the priority matrix:
catalog-level fixes, marker expansions, inference wording fixes, and
validator hygiene.

Audit date: 2026-09-06

Planning baseline: main at e7cdd11 ("docs(mobilome): intake and document
paywalled reference papers for integrons, CRISPR-Cas, and relaxases").

Fence semantics (inherited from the integration plan): only blocks fenced as
`diff` below are applyable candidate patches. Blocks fenced as `text`,
`python`, or `yaml` are discussion sketches, not patches.

Scope exclusions (deferred to future milestones):
- Finding C1 (data-driven rule engine) — architectural refactor, high effort.
- Finding C3 (strand-aware TA pairing) — requires new YAML field, inference
  branch, and test fixtures.
- Finding C4 (O(N×M) TA pair filtering) — algorithm refactor.
- Finding C5 (tripartite TA systems) — requires N-marker cluster spatial
  engine.
- TenpIN / CptIN markers — need validated regex patterns and provenance.

## Patch scope summary

| Item | Deferred Finding | Effort | Impact |
|---|---|---|---|
| 1 | C2: Topology-aware TA gap wording | Low | Medium |
| 2 | A4: Retained-evidence validator hardcode removal | Low | Medium |
| 3 | C7.2: Subtype-suffixed `cas` gene regex | Low | High |
| 4 | C7.1: `virB` / `trb` MPF gene coverage | Medium | High |
| 5 | C7.3: `repX` / `repL` replication initiators | Medium | Medium |
| 6 | A3: Handoff `source_ids` in schema and model | Medium | Medium |
| 7 | A5: Standardize `missing_components` with nested `missing_facets` | Medium | High |
| 8 | Spatial limitations home in `RepliconInventory` | Low | Medium |
| 9 | Catalog ID hygiene: disambiguate shadowed CRISPR markers | Low | Medium |
| 10 | Version bumps, golden regeneration | Low | — |

---

## 1. C2: Topology-aware TA gap wording

### 1.1 Finding

In `inference.py` line 296, the TA pair hypothesis limitation text hardcodes
"maximum circular gap" regardless of the actual record topology. On a linear
chromosome, this reads as "Configured maximum circular gap: 2000 bp" which is
misleading.

### 1.2 Diff — `src/genbank_parser/mobilome/inference.py`

```diff
diff --git a/src/genbank_parser/mobilome/inference.py b/src/genbank_parser/mobilome/inference.py
--- a/src/genbank_parser/mobilome/inference.py
+++ b/src/genbank_parser/mobilome/inference.py
@@ -293,7 +293,8 @@ def _infer_toxin_antitoxin_pairs(
                         supporting=(first, second),
                         limitations=tuple(rule.limitations)
                         + (
-                            f"Configured maximum circular gap: {rule.max_circular_gap_bp} bp; observed gap: {gap} bp.",
+                            f"Configured maximum gap ({inventory.topology}): "
+                            f"{rule.max_circular_gap_bp} bp; observed gap: {gap} bp.",
                         ),
                     )
                 )
```

### 1.3 Golden impact

The golden file `tests/golden/mobilome_v1.json` line 1160 changes from:
```text
"Configured maximum circular gap: 2000 bp; observed gap: 9 bp."
```
to:
```text
"Configured maximum gap (circular): 2000 bp; observed gap: 9 bp."
```

The same string appears in `tests/golden/mobilome_v1.tsv` inside the
`limitations_json` column for the TA hypothesis row.

---

## 2. A4: Retained-evidence validator hardcode removal

### 2.1 Finding

In `database.py` line 691, the disabled-rule validator subtracts the hardcoded
string `"pxo_numbered_product_annotation"` from the facet-id set. This is
because `pxo_numbered_product_annotation` is a marker ID, not a facet ID,
but the validator was built assuming all retained-evidence tokens are facet IDs.

### 2.2 Diff — `src/genbank_parser/mobilome/database.py`

This diff changes `_parse_disabled_rules` to accept `marker_ids` as an
additional parameter and validate retained-evidence tokens against the union
of `facet_ids | marker_ids`.

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -655,7 +655,10 @@ def _parse_component_support(
 
 
 def _parse_disabled_rules(
-    inference_data: dict[str, Any], facet_ids: set[str], rule_ids: set[str]
+    inference_data: dict[str, Any],
+    facet_ids: set[str],
+    marker_ids: set[str],
+    rule_ids: set[str],
 ) -> tuple[DisabledAggregateRule, ...]:
     raw_rules = inference_data.get("disabled_aggregate_rules")
     if not isinstance(raw_rules, list) or not all(
@@ -687,8 +690,9 @@ def _parse_disabled_rules(
             raw.get("evidence_retained_as"),
             f"disabled rule {rule_id}.evidence_retained_as",
         )
+        allowed = facet_ids | marker_ids
         unknown = sorted(
-            set(retained) - facet_ids - {"pxo_numbered_product_annotation"}
+            set(retained) - allowed
         )
         if unknown:
             raise MobilomeDatabaseError(
```

### 2.3 Diff — call site in `load_mobilome_database`

The call site in `load_mobilome_database` must pass `marker_ids`:

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -855,7 +858,7 @@ def load_mobilome_database(database_dir: str | Path | None = None) -> MobilomeDa
         {marker.id for marker in markers},
     )
     disabled = _parse_disabled_rules(
-        inference_data, facet_ids, {rule.id for rule in rules}
+        inference_data, facet_ids, {marker.id for marker in markers}, {rule.id for rule in rules}
     )
     handoffs = _parse_handoffs(inference_data)
     _validate_source_references(markers, rules, disabled, sources)
```

---

## 3. C7.2: Subtype-suffixed `cas` gene regex

### 3.1 Finding

The current `crispr_cas` gene regex `(?i)^cas(?:1|2|...|14|A|...|F)$` rejects
subtype-suffixed names like `cas8f`, `cas12a`, `cas13d`, which are standard
Cas nomenclature (Makarova 2020). It also risks a Python `re.error` when
inline flags follow alternation in certain patterns.

### 3.2 Diff — `src/genbank_parser/data/mobilome/markers.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -423,7 +423,7 @@
       - field: gene
         mode: regex
         patterns:
-          - "(?i)^cas(?:1|2|3|4|5|6|7|8|9|10|12|13|14|A|B|C|D|E|F)$"
+          - "(?i)^cas(?:(?:[1-9]|1[0-4])[a-z]?|[A-F])$"
         negative_patterns: []
         strength: 2
```

Validated matches: `cas1` through `cas14`, `cas8f`, `cas12a`, `cas13d`,
`casA` through `casF`.

Validated rejections: `casein`, `cascade`, `cas15`, `cas123`, `cas0`.

### 3.3 Product regex: no change needed

The product regex `(?i)\b(?:CRISPR-associated|Cas(?:1|...|14))\b` matches
free-text product descriptions, not gene symbols. Subtype suffixes in product
text are already captured by the `CRISPR-associated` branch. No change needed.

---

## 4. C7.1: `virB` / `trb` MPF gene coverage

### 4.1 Finding

Currently, `mpf_component` gene regex only matches IncF-style `tra` genes.
Agrobacterium Ti-plasmid MPF (`virB1`-`virB11`) and IncP plasmid MPF
(`trbB`-`trbL`) are standard conjugation components that should be detected.

### 4.2 Diff — `src/genbank_parser/data/mobilome/markers.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -188,6 +188,8 @@
       - field: gene
         mode: regex
         patterns:
           - "(?i)^tra(?:A|B|C|E|F|K|L|N|U|W)$"
+          - "(?i)^virB(?:[1-9]|1[01])$"
+          - "(?i)^trb[B-L]$"
         negative_patterns: []
         strength: 2
```

### 4.3 Provenance: no new sources needed

The `virB` and `trb` patterns are covered by `christie-2025-t4ss` and
`garcillan-barcia-2025-extended-mobility`, both already in the marker's
`sources` array. No provenance.yaml change needed.

---

## 5. C7.3: `repX` / `repL` replication initiators

### 5.1 Finding

Replication-initiation markers currently match `repA`, `repB`, `repS`, `repF`,
`repFR` by casefold-exact gene. The tubulin-like `repX` on pXO1 and `repL` are
established replication initiators that should be detected.

### 5.2 Diff — `src/genbank_parser/data/mobilome/markers.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -118,7 +118,7 @@
       - field: gene
         mode: casefold_exact
-        patterns: [repA, repB, repS, repF, repFR]
+        patterns: [repA, repB, repL, repS, repF, repFR, repX]
         negative_patterns: []
         strength: 2
```

### 5.3 Provenance: new source entries

Two new provenance entries are required for `repX` and `repL`.

```diff
diff --git a/src/genbank_parser/data/mobilome/provenance.yaml b/src/genbank_parser/data/mobilome/provenance.yaml
--- a/src/genbank_parser/data/mobilome/provenance.yaml
+++ b/src/genbank_parser/data/mobilome/provenance.yaml
@@ -63,6 +63,30 @@
     supports:
       - pXO annotations alone do not specify a single replication mechanism.
     limitations:
       - Product labels are not a sequence-homology or ancestry result.
+  - id: anand-2008-repx
+    type: primary_article
+    title: "GTP-dependent polymerization of the tubulin-like RepX replication protein encoded by the pXO1 plasmid of Bacillus anthracis"
+    year: 2008
+    doi: 10.1111/j.1365-2958.2007.06100.x
+    pmid: "18179418"
+    url: https://doi.org/10.1111/j.1365-2958.2007.06100.x
+    supports:
+      - RepX is a tubulin-like replication initiator encoded on pXO1 of Bacillus anthracis.
+    limitations:
+      - Annotation alone does not identify the replication mechanism (rolling-circle, theta, or strand-displacement).
+  - id: tinsley-2006-repx
+    type: primary_article
+    title: "A novel FtsZ-like protein is involved in replication of the anthrax toxin-encoding pXO1 plasmid in Bacillus anthracis"
+    year: 2006
+    doi: 10.1128/JB.188.8.2829-2835.2006
+    pmid: "16585744"
+    pmcid: PMC1447009
+    url: https://pmc.ncbi.nlm.nih.gov/articles/PMC1447009/
+    supports:
+      - An FtsZ-like protein (RepX) is required for pXO1 replication in Bacillus anthracis.
+    limitations:
+      - An FtsZ-like annotation does not establish a replication mechanism without further analysis.
```

### 5.4 Marker source update

The `replication_initiation_candidate` marker must cite the new provenance
sources alongside the existing local policy:

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -130,7 +130,7 @@
     requires_non_pseudo: true
     requires_complete: true
-    sources: [local-annotation-candidate-policy]
+    sources: [anand-2008-repx, local-annotation-candidate-policy, tinsley-2006-repx]
     interpretation: Retains a replication-initiation annotation candidate.
     limitations:
       - Does not identify a replication mechanism or incompatibility group.
```

---

## 6. A3: Handoff `source_ids`

### 6.1 Finding

In Schema v1, handoff objects lack a `source_ids` array. Unlike hits and
hypotheses, which cite specific literature keys in `provenance.yaml`,
handoffs cannot carry citations directly in the JSON without triggering
schema validation errors.

### 6.2 Diff — `src/genbank_parser/mobilome/models.py` (ExternalHandoff)

```diff
diff --git a/src/genbank_parser/mobilome/models.py b/src/genbank_parser/mobilome/models.py
--- a/src/genbank_parser/mobilome/models.py
+++ b/src/genbank_parser/mobilome/models.py
@@ -222,6 +222,7 @@ class ExternalHandoff:
     required_input: tuple[str, ...]
     required_provenance_fields: tuple[str, ...]
     limitations: tuple[str, ...]
+    source_ids: tuple[str, ...] = ()
```

### 6.3 Diff — `src/genbank_parser/data/mobilome/report.schema.json` (handoff)

```diff
diff --git a/src/genbank_parser/data/mobilome/report.schema.json b/src/genbank_parser/data/mobilome/report.schema.json
--- a/src/genbank_parser/data/mobilome/report.schema.json
+++ b/src/genbank_parser/data/mobilome/report.schema.json
@@ -168,8 +168,8 @@
     "handoff": {
       "type": "object", "additionalProperties": false,
-      "required": ["handoff_id", "tool", "purpose", "status", "reason", "required_input", "required_provenance_fields", "limitations"],
-      "properties": {"handoff_id": {"type": "string"}, "tool": {"type": "string"}, "purpose": {"type": "string"}, "status": {"const": "not_run"}, "reason": {"type": "string"}, "required_input": {"$ref": "#/$defs/stringArray"}, "required_provenance_fields": {"$ref": "#/$defs/stringArray"}, "limitations": {"$ref": "#/$defs/stringArray"}}
+      "required": ["handoff_id", "tool", "purpose", "status", "reason", "required_input", "required_provenance_fields", "limitations", "source_ids"],
+      "properties": {"handoff_id": {"type": "string"}, "tool": {"type": "string"}, "purpose": {"type": "string"}, "status": {"const": "not_run"}, "reason": {"type": "string"}, "required_input": {"$ref": "#/$defs/stringArray"}, "required_provenance_fields": {"$ref": "#/$defs/stringArray"}, "limitations": {"$ref": "#/$defs/stringArray"}, "source_ids": {"$ref": "#/$defs/stringArray"}}
     }
```

### 6.4 Diff — `src/genbank_parser/data/mobilome/inference.yaml` (handoff sources)

Each handoff should now declare its `sources` key so the database loader can
populate `source_ids`:

```diff
diff --git a/src/genbank_parser/data/mobilome/inference.yaml b/src/genbank_parser/data/mobilome/inference.yaml
--- a/src/genbank_parser/data/mobilome/inference.yaml
+++ b/src/genbank_parser/data/mobilome/inference.yaml
@@ -125,30 +125,35 @@
     limitations:
       - A future tool result remains dependent on database composition and thresholds.
+    sources: [carattoli-2014-plasmidfinder, robertson-2018-mob-suite]
   - id: amr_annotation
     tool: AMRFinderPlus
     purpose: Structured antimicrobial-resistance annotation
     reason: Generic product text is retained only as an annotation candidate.
     required_input: [protein or nucleotide sequence]
     required_provenance_fields: [tool_version, database_version, organism_option, method, thresholds]
     limitations:
       - A tool annotation is not a phenotype measurement.
+    sources: [feldgarden-2021-amrfinderplus]
   - id: vf_annotation
     tool: VFDB or curated virulence workflow
     purpose: Virulence-associated comparison
     reason: Product text is retained only as an annotation candidate.
     required_input: [protein or nucleotide sequence]
     required_provenance_fields: [database_release, search_method, thresholds, coverage, interpretation_limits]
     limitations:
       - A similarity result is not virulence phenotype evidence.
+    sources: [chen-2005-vfdb, liu-2022-vfdb]
   - id: element_boundary
     tool: Prophage or ICE caller
     purpose: Dedicated element-boundary analysis
     reason: Component annotations do not establish element boundaries.
     required_input: [assembly sequence and annotation]
     required_provenance_fields: [tool_version, database_version, caller_thresholds, boundary_uncertainty]
     limitations:
       - Predicted boundaries require method-specific review.
+    sources: [arndt-2016-phaster, camargo-2024-genomad]
   - id: mechanism_evidence
     tool: HMM, domain, or sequence workflow
     purpose: Mechanism-specific evidence
     reason: Generic Rep and pXO labels cannot identify a replication mechanism.
     required_input: [protein sequence]
     required_provenance_fields: [model_accession, model_version, cutoffs, alignment_coverage]
     limitations:
       - Domain evidence still requires a reviewed mechanism rule.
+    sources: [local-annotation-candidate-policy]
```

### 6.5 Diff — `src/genbank_parser/mobilome/database.py` (handoff parser)

The `_parse_handoffs` function must parse the new optional `sources` key and
populate `source_ids`. The `_expect_keys` set must include `"sources"`.

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -733,6 +733,7 @@ def _parse_handoffs(inference_data: dict[str, Any]) -> tuple[ExternalHandoff, ..
                 "reason",
                 "required_input",
                 "required_provenance_fields",
                 "limitations",
+                "sources",
             },
         )
@@ -756,6 +757,12 @@ def _parse_handoffs(inference_data: dict[str, Any]) -> tuple[ExternalHandoff, ..
                 limitations=_string_list(
                     raw.get("limitations"), f"handoff {handoff_id}.limitations"
                 ),
+                source_ids=tuple(
+                    sorted(
+                        _string_list(
+                            raw.get("sources", []), f"handoff {handoff_id}.sources"
+                        )
+                    )
+                ),
             )
         )
```

### 6.6 Diff — `src/genbank_parser/mobilome/database.py` (`_validate_source_references`)

Add handoffs to the source-reference cross-check:

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -765,6 +772,7 @@ def _validate_source_references(
 def _validate_source_references(
     markers: tuple[MobilomeMarker, ...],
     rules: tuple[InferenceRule, ...],
     disabled: tuple[DisabledAggregateRule, ...],
+    handoffs: tuple[ExternalHandoff, ...],
     sources: tuple[ProvenanceSource, ...],
 ) -> None:
     known_sources = {source.id for source in sources}
@@ -775,6 +783,10 @@ def _validate_source_references(
             (f"disabled aggregate rule {rule.rule_id}", rule.source_ids)
             for rule in disabled
         ),
+        *(
+            (f"handoff {handoff.handoff_id}", handoff.source_ids)
+            for handoff in handoffs
+        ),
     ]:
```

### 6.7 Diff — call site for `_validate_source_references`

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -862,7 +872,7 @@ def load_mobilome_database(database_dir: str | Path | None = None) -> MobilomeDa
     handoffs = _parse_handoffs(inference_data)
-    _validate_source_references(markers, rules, disabled, sources)
+    _validate_source_references(markers, rules, disabled, handoffs, sources)
```

### 6.8 Diff — `src/genbank_parser/mobilome/report.py` (handoff JSON renderer)

The JSON renderer must include `source_ids`:

```diff
diff --git a/src/genbank_parser/mobilome/report.py b/src/genbank_parser/mobilome/report.py
--- a/src/genbank_parser/mobilome/report.py
+++ b/src/genbank_parser/mobilome/report.py
@@ -237,6 +237,7 @@ def _handoff_dict(handoff: ExternalHandoff) -> dict[str, object]:
         "required_provenance_fields": _list(handoff.required_provenance_fields),
         "limitations": _list(handoff.limitations),
+        "source_ids": _list(handoff.source_ids),
     }
```

### 6.9 Diff — `src/genbank_parser/mobilome/report.py` (TSV handoff rows)

The TSV handoff rows should include `source_ids_json`:

```diff
diff --git a/src/genbank_parser/mobilome/report.py b/src/genbank_parser/mobilome/report.py
--- a/src/genbank_parser/mobilome/report.py
+++ b/src/genbank_parser/mobilome/report.py
@@ -625,6 +626,7 @@ def render_tsv(report: MobilomeReport) -> str:
                     }
                 ),
+                "source_ids_json": _compact_json(handoff.source_ids),
             }
         )
```

---

## 7. A5: Standardize `missing_components` with nested `missing_facets`

### 7.1 Finding

In `inference.py`, `_infer_component_rule` currently emits a mixed vocabulary in `MobilomeHypothesis.missing_components`:
- When mobility helper evaluation fails, it emits Facet IDs (e.g. `("mobility_mpf",)`).
- When replication evaluation fails with no hits, it emits a Component ID (`("replication_candidate",)`).

Furthermore, `missing_components` was defined as `tuple[str, ...]`, losing the structural relationship between which component was missing and which facets within that component were unsatisfied.

### 7.2 Architecture & Model

We introduce a structured dataclass `MissingComponent` in `models.py`:
```python
@dataclass(frozen=True)
class MissingComponent:
    """One unsatisfied component requirement with its missing facets."""

    component_id: str
    missing_facets: tuple[str, ...] = ()
```

`MobilomeHypothesis.missing_components` becomes `tuple[MissingComponent, ...]`. In the JSON report schema, `missing_components` is an array of objects:
```json
[
  {
    "component_id": "helper_machinery",
    "missing_facets": ["mobility_mpf"]
  }
]
```
For fully satisfied hypotheses, `missing_components` remains empty `[]`.

### 7.3 Diff — `src/genbank_parser/mobilome/models.py`

```diff
diff --git a/src/genbank_parser/mobilome/models.py b/src/genbank_parser/mobilome/models.py
--- a/src/genbank_parser/mobilome/models.py
+++ b/src/genbank_parser/mobilome/models.py
@@ -165,6 +165,14 @@ class MobilomeHit:
 
 
+@dataclass(frozen=True)
+class MissingComponent:
+    """One unsatisfied component requirement with its missing facets."""
+
+    component_id: str
+    missing_facets: tuple[str, ...] = ()
+
+
 @dataclass(frozen=True)
 class MobilomeHypothesis:
     """Cautious record- or inter-record-level inference from retained hits."""
@@ -173,7 +181,7 @@ class MobilomeHypothesis:
     status: HypothesisStatus
     summary: str
     participants: tuple[HypothesisParticipant, ...]
     supporting_hit_ids: tuple[str, ...]
-    missing_components: tuple[str, ...]
+    missing_components: tuple[MissingComponent, ...]
     conflicting_hit_ids: tuple[str, ...]
     limitations: tuple[str, ...]
     source_ids: tuple[str, ...]
@@ -375,6 +383,7 @@ __all__ = [
     "EvidenceReason",
     "FeatureRef",
     "HypothesisParticipant",
+    "MissingComponent",
     "MobilomeHit",
     "MobilomeHypothesis",
     "MobilomeMarker",
```

### 7.4 Diff — `src/genbank_parser/mobilome/inference.py`

```diff
diff --git a/src/genbank_parser/mobilome/inference.py b/src/genbank_parser/mobilome/inference.py
--- a/src/genbank_parser/mobilome/inference.py
+++ b/src/genbank_parser/mobilome/inference.py
@@ -14,6 +14,7 @@ from .models import (
     HypothesisParticipant,
     InferenceComponent,
     InferenceRule,
+    MissingComponent,
     MobilomeDatabase,
     MobilomeHit,
     MobilomeHypothesis,
@@ -90,16 +91,18 @@ def _component_support(
         return None, missing or component.facets
     if missing:
         return None, missing
+    missing_facets: list[str] = []
     selected_hits: list[MobilomeHit] = []
     used = set(forbidden)
     for facet in component.facets:
         chosen = _choose_hits(
             candidates_by_facet[facet],
             minimums.get(facet, 1),
             forbidden_features=used if distinct_features else set(),
         )
         if chosen is None:
-            return None, (facet,)
+            missing_facets.append(facet)
+            continue
         selected_hits.extend(chosen)
         if distinct_features:
             used.update(_feature_key(hit) for hit in chosen)
+    if missing_facets:
+        return None, tuple(missing_facets)
     return tuple(selected_hits), ()
@@ -112,13 +115,13 @@ def _rule_component_support(
     database: MobilomeDatabase,
     *,
     role: str,
-) -> tuple[tuple[MobilomeHit, ...] | None, tuple[str, ...]]:
+) -> tuple[tuple[MobilomeHit, ...] | None, tuple[MissingComponent, ...]]:
     """Evaluate every component required for one role in a rule."""
 
     component_ids = rule.components_by_role.get(role, ())
     selected: list[MobilomeHit] = []
     used: set[tuple[int, int]] = set()
-    missing: list[str] = []
+    missing: list[MissingComponent] = []
     for component_id in component_ids:
         component = database.component_map[component_id]
         component_hits, component_missing = _component_support(
@@ -128,13 +131,18 @@ def _rule_component_support(
             distinct_features=rule.distinct_features,
         )
         if component_hits is None:
-            missing.extend(component_missing or (component_id,))
+            missing.append(
+                MissingComponent(
+                    component_id=component_id,
+                    missing_facets=component_missing,
+                )
+            )
             continue
         selected.extend(component_hits)
         if rule.distinct_features:
             used.update(_feature_key(hit) for hit in component_hits)
     if missing:
-        return None, tuple(sorted(set(missing)))
+        return None, tuple(missing)
     return tuple(selected), ()
@@ -154,7 +162,7 @@ def _hypothesis(
     rule: InferenceRule,
     inventory: RepliconInventory,
     supporting: Sequence[MobilomeHit],
-    missing: Sequence[str] = (),
+    missing: Sequence[MissingComponent] = (),
     summary: str | None = None,
     limitations: Sequence[str] | None = None,
 ) -> MobilomeHypothesis:
@@ -166,7 +174,7 @@ def _hypothesis(
         summary=summary or rule.wording,
         participants=(_participant(inventory),),
         supporting_hit_ids=tuple(sorted({hit.hit_id for hit in supporting})),
-        missing_components=tuple(sorted(set(missing))),
+        missing_components=tuple(missing),
         conflicting_hit_ids=(),
         limitations=tuple(limitations if limitations is not None else rule.limitations),
         source_ids=tuple(sorted(rule.sources)),
@@ -220,10 +228,17 @@ def _infer_component_rule(
         extra_limitations.append(
             "No annotated replication candidate was observed; autonomy is unresolved."
         )
+    default_missing = tuple(
+        MissingComponent(
+            component_id=cid,
+            missing_facets=database.component_map[cid].facets,
+        )
+        for cid in component_ids
+    )
     return _hypothesis(
         hypothesis_id=f"h:{rule.id}:{inventory.record_index}",
         rule=rule,
         inventory=inventory,
         supporting=raw,
-        missing=missing or ("replication_candidate",),
+        missing=missing or default_missing,
         summary="Replication evidence insufficient; mechanism unresolved"
```

### 7.5 Diff — `src/genbank_parser/data/mobilome/report.schema.json`

```diff
diff --git a/src/genbank_parser/data/mobilome/report.schema.json b/src/genbank_parser/data/mobilome/report.schema.json
--- a/src/genbank_parser/data/mobilome/report.schema.json
+++ b/src/genbank_parser/data/mobilome/report.schema.json
@@ -150,8 +150,13 @@
     "participant": {
       "type": "object", "additionalProperties": false,
       "required": ["role", "record_index", "record_id"],
       "properties": {"role": {"enum": ["subject", "target", "helper"]}, "record_index": {"type": "integer", "minimum": 1}, "record_id": {"type": "string"}}
     },
+    "missingComponent": {
+      "type": "object", "additionalProperties": false,
+      "required": ["component_id", "missing_facets"],
+      "properties": {"component_id": {"type": "string"}, "missing_facets": {"$ref": "#/$defs/stringArray"}}
+    },
     "hypothesis": {
       "type": "object", "additionalProperties": false,
       "required": ["hypothesis_id", "rule_id", "kind", "status", "summary", "participants", "supporting_hit_ids", "missing_components", "conflicting_hit_ids", "limitations", "source_ids"],
@@ -155,5 +160,5 @@
-        "supporting_hit_ids": {"$ref": "#/$defs/stringArray"}, "missing_components": {"$ref": "#/$defs/stringArray"}, "conflicting_hit_ids": {"$ref": "#/$defs/stringArray"}, "limitations": {"$ref": "#/$defs/stringArray"}, "source_ids": {"$ref": "#/$defs/stringArray"}
+        "supporting_hit_ids": {"$ref": "#/$defs/stringArray"}, "missing_components": {"type": "array", "items": {"$ref": "#/$defs/missingComponent"}}, "conflicting_hit_ids": {"$ref": "#/$defs/stringArray"}, "limitations": {"$ref": "#/$defs/stringArray"}, "source_ids": {"$ref": "#/$defs/stringArray"}
       }
     },
```

### 7.6 Diff — `src/genbank_parser/mobilome/report.py`

```diff
diff --git a/src/genbank_parser/mobilome/report.py b/src/genbank_parser/mobilome/report.py
--- a/src/genbank_parser/mobilome/report.py
+++ b/src/genbank_parser/mobilome/report.py
@@ -210,7 +210,13 @@ def _hypothesis_dict(hypothesis: MobilomeHypothesis) -> dict[str, object]:
         "summary": hypothesis.summary,
         "participants": [_participant_dict(item) for item in hypothesis.participants],
         "supporting_hit_ids": _list(hypothesis.supporting_hit_ids),
-        "missing_components": _list(hypothesis.missing_components),
+        "missing_components": [
+            {
+                "component_id": item.component_id,
+                "missing_facets": _list(item.missing_facets),
+            }
+            for item in hypothesis.missing_components
+        ],
         "conflicting_hit_ids": _list(hypothesis.conflicting_hit_ids),
         "limitations": _list(hypothesis.limitations),
@@ -546,7 +552,15 @@ def render_tsv(report: MobilomeReport) -> str:
             "participants_json": _compact_json(
                 [_participant_dict(item) for item in hypothesis.participants]
             ),
             "supporting_hit_ids_json": _compact_json(hypothesis.supporting_hit_ids),
-            "missing_components_json": _compact_json(hypothesis.missing_components),
+            "missing_components_json": _compact_json(
+                [
+                    {
+                        "component_id": item.component_id,
+                        "missing_facets": _list(item.missing_facets),
+                    }
+                    for item in hypothesis.missing_components
+                ]
+            ),
             "conflicting_hit_ids_json": _compact_json(hypothesis.conflicting_hit_ids),
@@ -683,5 +697,9 @@ def render_text(report: MobilomeReport) -> str:
         if insufficient:
             print("    Conflicts and missing components", file=output)
             for item in insufficient:
-                missing = ", ".join(item.missing_components) or "none configured"
+                missing = ", ".join(
+                    f"{comp.component_id} (missing {', '.join(comp.missing_facets)})"
+                    if comp.missing_facets
+                    else comp.component_id
+                    for comp in item.missing_components
+                ) or "none configured"
                 print(f"      {item.summary}; missing: {missing}", file=output)
```

### 7.7 Diff — `src/genbank_parser/mobilome/__init__.py`

```diff
diff --git a/src/genbank_parser/mobilome/__init__.py b/src/genbank_parser/mobilome/__init__.py
--- a/src/genbank_parser/mobilome/__init__.py
+++ b/src/genbank_parser/mobilome/__init__.py
@@ -27,6 +27,7 @@ from .models import (
     ExternalHandoff,
     FeatureRef,
     HypothesisParticipant,
+    MissingComponent,
     MobilomeAnalysisError,
     MobilomeDatabase,
     MobilomeDatabaseError,
@@ -252,6 +253,7 @@ __all__ = [
     "FeatureRef",
     "HypothesisParticipant",
     "MatchMode",
+    "MissingComponent",
     "MobilomeAnalysisError",
     "MobilomeDatabase",
     "MobilomeDatabaseError",
```

### 7.8 Diff — `tests/test_mobilome_patch_01.py`

```diff
diff --git a/tests/test_mobilome_patch_01.py b/tests/test_mobilome_patch_01.py
--- a/tests/test_mobilome_patch_01.py
+++ b/tests/test_mobilome_patch_01.py
@@ -17,6 +17,7 @@ from genbank_parser.mobilome import (
     infer_replicon_hypotheses,
     inventory_replicons,
     load_mobilome_database,
+    MissingComponent,
     record_inference_limitations,
     scan_mobilome_features,
 )
@@ -419,7 +420,9 @@ def test_partial_conjugation_helper_marks_missing_components_clearly() -> None:
                 ),
             ),
             "helper_machinery_annotation_candidate",
-            ("mobility_mpf",),
+            (
+                MissingComponent("helper_machinery", ("mobility_mpf",)),
+            ),
         ),
         "single_mpf": (
             (
@@ -431,7 +434,9 @@ def test_partial_conjugation_helper_marks_missing_components_clearly() -> None:
                 ),
             ),
             "helper_machinery_annotation_candidate",
-            ("mobility_mpf", "mobility_t4cp"),
+            (
+                MissingComponent("helper_machinery", ("mobility_t4cp", "mobility_mpf")),
+            ),
         ),
         "t4cp_plus_one_mpf": (
             (
@@ -449,7 +454,9 @@ def test_partial_conjugation_helper_marks_missing_components_clearly() -> None:
                 ),
             ),
             "helper_machinery_annotation_candidate",
-            ("mobility_mpf",),
+            (
+                MissingComponent("helper_machinery", ("mobility_mpf",)),
+            ),
         ),
     }
```

---

## 8. Spatial limitations home in `RepliconInventory`

### 8.1 Finding

In Schema v1, spatial pairing caveats from `record_inference_limitations` (such as unlocatable features skipping spatial pairing) are pushed into `RepliconInventory.classification.limitations`. This mixes replicon-level plasmid/chromosome classification caveats (e.g. `"Source plasmid declaration missing"`) with locus-level feature coordinate limitations.

We introduce `RepliconInventory.spatial_limitations: tuple[str, ...]` and route `record_inference_limitations` results directly into this field.

### 8.2 Diff — `src/genbank_parser/mobilome/models.py`

```diff
diff --git a/src/genbank_parser/mobilome/models.py b/src/genbank_parser/mobilome/models.py
--- a/src/genbank_parser/mobilome/models.py
+++ b/src/genbank_parser/mobilome/models.py
@@ -194,6 +194,7 @@ class RepliconInventory:
     pseudo_feature_count: int
     pseudo_cds_count: int
     pseudogene_feature_count: int
+    spatial_limitations: tuple[str, ...] = ()
```

### 8.3 Diff — `src/genbank_parser/mobilome/replicons.py`

```diff
diff --git a/src/genbank_parser/mobilome/replicons.py b/src/genbank_parser/mobilome/replicons.py
--- a/src/genbank_parser/mobilome/replicons.py
+++ b/src/genbank_parser/mobilome/replicons.py
@@ -218,6 +218,7 @@ def inventory_replicons(document: GenBankDocument) -> tuple[RepliconInventory,
             pseudo_feature_count=pseudo_feature_count,
             pseudo_cds_count=pseudo_cds_count,
             pseudogene_feature_count=pseudogene_feature_count,
+            spatial_limitations=(),
         )
```

### 8.4 Diff — `src/genbank_parser/mobilome/__init__.py`

```diff
diff --git a/src/genbank_parser/mobilome/__init__.py b/src/genbank_parser/mobilome/__init__.py
--- a/src/genbank_parser/mobilome/__init__.py
+++ b/src/genbank_parser/mobilome/__init__.py
@@ -192,10 +192,7 @@ def analyze_mobilome(
         record_hits = tuple(
             hit for hit in hits if hit.feature.record_index == member.record_index
         )
         limitations = record_inference_limitations(member, record_hits, active_database)
         if limitations:
             member = replace(
                 member,
-                classification=replace(
-                    member.classification,
-                    limitations=member.classification.limitations + limitations,
-                ),
+                spatial_limitations=member.spatial_limitations + limitations,
             )
         all_assessments.append(
```

### 8.5 Diff — `src/genbank_parser/mobilome/report.py`

```diff
diff --git a/src/genbank_parser/mobilome/report.py b/src/genbank_parser/mobilome/report.py
--- a/src/genbank_parser/mobilome/report.py
+++ b/src/genbank_parser/mobilome/report.py
@@ -138,6 +138,7 @@ def _inventory_dict(inventory: RepliconInventory) -> dict[str, object]:
                 "limitations": _list(classification.limitations),
             },
+            "spatial_limitations": _list(inventory.spatial_limitations),
             "feature_count": inventory.feature_count,
             "cds_count": inventory.cds_count,
```

### 8.6 Diff — `src/genbank_parser/data/mobilome/report.schema.json`

```diff
diff --git a/src/genbank_parser/data/mobilome/report.schema.json b/src/genbank_parser/data/mobilome/report.schema.json
--- a/src/genbank_parser/data/mobilome/report.schema.json
+++ b/src/genbank_parser/data/mobilome/report.schema.json
@@ -108,7 +108,7 @@
     "inventory": {
       "type": "object", "additionalProperties": false,
-      "required": ["record_index", "record_id", "name", "description", "length", "topology", "source_qualifiers", "classification", "feature_count", "cds_count", "pseudo_feature_count", "pseudo_cds_count", "pseudogene_feature_count"],
+      "required": ["record_index", "record_id", "name", "description", "length", "topology", "source_qualifiers", "classification", "spatial_limitations", "feature_count", "cds_count", "pseudo_feature_count", "pseudo_cds_count", "pseudogene_feature_count"],
       "properties": {
         "record_index": {"type": "integer", "minimum": 1}, "record_id": {"type": "string"}, "name": {"type": "string"}, "description": {"type": "string"}, "length": {"type": "integer", "minimum": 0},
         "topology": {"enum": ["circular", "linear", "unknown"]},
         "source_qualifiers": {"type": "array", "items": {"$ref": "#/$defs/sourceQualifier"}},
         "classification": {"$ref": "#/$defs/classification"},
+        "spatial_limitations": {"$ref": "#/$defs/stringArray"},
         "feature_count": {"type": "integer", "minimum": 0}, "cds_count": {"type": "integer", "minimum": 0}, "pseudo_feature_count": {"type": "integer", "minimum": 0}, "pseudo_cds_count": {"type": "integer", "minimum": 0}, "pseudogene_feature_count": {"type": "integer", "minimum": 0}
       }
     },
```

---

## 9. Catalog ID hygiene: disambiguate shadowed CRISPR markers

### 9.1 Finding

In `markers.yaml`, markers `crispr_array` and `crispr_cas` share identical IDs with their enclosing facets `crispr_array` and `crispr_cas`. This shadowing causes ambiguity in rule declarations and documentation.

We rename:
- Marker `crispr_array` -> `crispr_repeat_array`
- Marker `crispr_cas` -> `crispr_cas_candidate`

And in `database.py`, we add an explicit validation assertion in `_parse_markers` ensuring that no marker ID ever collides with or shadows an existing facet ID.

### 9.2 Diff — `src/genbank_parser/data/mobilome/markers.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -392,7 +392,7 @@ markers:
     limitations:
       - Annotation does not establish active immunity or target specificity.
 
-  - id: crispr_array
+  - id: crispr_repeat_array
     label: CRISPR repeat-array annotation
     facets: [crispr_array]
     feature_types: [repeat_region, misc_feature, regulatory]
@@ -415,7 +415,7 @@ markers:
     limitations:
       - Annotation does not establish a complete or active CRISPR system.
 
-  - id: crispr_cas
+  - id: crispr_cas_candidate
     label: Cas annotation candidate
     facets: [crispr_cas]
     feature_types: [CDS, pseudogene]
```

### 9.3 Diff — `src/genbank_parser/mobilome/database.py`

```diff
diff --git a/src/genbank_parser/mobilome/database.py b/src/genbank_parser/mobilome/database.py
--- a/src/genbank_parser/mobilome/database.py
+++ b/src/genbank_parser/mobilome/database.py
@@ -291,6 +291,10 @@ def _parse_markers(
     for raw_marker in raw_markers:
         marker_id = _string(raw_marker.get("id"), "marker.id")
+        if marker_id in facet_ids:
+            raise MobilomeDatabaseError(
+                f"Marker ID '{marker_id}' shadows an existing facet ID"
+            )
         _expect_keys(raw_marker, f"marker {marker_id}", allowed)
         facets = _string_list(raw_marker.get("facets"), f"marker {marker_id}.facets")
```

### 9.4 Diff — Test suites (`test_mobilome_scanner.py`, `test_mobilome_patch_01.py`, `test_mobilome_database.py`)

```diff
diff --git a/tests/test_mobilome_scanner.py b/tests/test_mobilome_scanner.py
--- a/tests/test_mobilome_scanner.py
+++ b/tests/test_mobilome_scanner.py
@@ -17,8 +17,8 @@ def test_mobilome_scanner_finds_all_synthetic_categories() -> None:
     assert by_marker["explicit_origin_of_transfer"].reasons[0].field == "regulatory_class"
-    assert by_marker["crispr_array"].reasons[0].field == "note"
-    assert {reason.field for reason in by_marker["crispr_array"].reasons} == {
+    assert by_marker["crispr_repeat_array"].reasons[0].field == "note"
+    assert {reason.field for reason in by_marker["crispr_repeat_array"].reasons} == {
         "note",
         "rpt_type",
     }
@@ -34,7 +34,7 @@ def test_mobilome_scanner_finds_all_synthetic_categories() -> None:
     assert by_marker["toxn"].reasons[0].field == "gene"
     assert by_marker["toxi"].reasons[0].field == "gene"
-    cas = next(hit for hit in hits if hit.marker_id == "crispr_cas")
+    cas = next(hit for hit in hits if hit.marker_id == "crispr_cas_candidate")
     assert {reason.field for reason in cas.reasons} == {"gene", "product"}
diff --git a/tests/test_mobilome_patch_01.py b/tests/test_mobilome_patch_01.py
--- a/tests/test_mobilome_patch_01.py
+++ b/tests/test_mobilome_patch_01.py
@@ -306,8 +306,8 @@ def test_database_loads_all_expected_facets_and_markers() -> None:
         "phage_structural_candidate",
-        "crispr_array",
-        "crispr_cas",
+        "crispr_repeat_array",
+        "crispr_cas_candidate",
         "pxo_numbered_product_annotation",
diff --git a/tests/test_mobilome_database.py b/tests/test_mobilome_database.py
--- a/tests/test_mobilome_database.py
+++ b/tests/test_mobilome_database.py
@@ -257,7 +257,7 @@ def test_database_loads_packaged_data() -> None:
     with open(Path(__file__).parents[1] / "src" / "genbank_parser" / "data" / "mobilome" / "markers.yaml") as f:
         p = yaml.safe_load(f)
         m = next(
-            marker for marker in p["markers"] if marker["id"] == "crispr_array"
+            marker for marker in p["markers"] if marker["id"] == "crispr_repeat_array"
         )
         assert m["requires_complete"] is False
```

### 9.5 Golden impact

In golden files (`tests/golden/mobilome_v1.json`, `.tsv`, `.txt`):
- Hit IDs change from `r1:f13:mcrispr_array` to `r1:f13:mcrispr_repeat_array`
- Marker ID changes from `crispr_array` to `crispr_repeat_array`
- Hit IDs change from `r1:f14:mcrispr_cas` to `r1:f14:mcrispr_cas_candidate`
- Marker ID changes from `crispr_cas` to `crispr_cas_candidate`
- Enclosing Facet IDs remain `["crispr_array"]` and `["crispr_cas"]`

---

## 10. Version bumps

### 10.1 `pyproject.toml`

```diff
diff --git a/pyproject.toml b/pyproject.toml
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -6,7 +6,7 @@
 name = "genbank-parser"
-version = "0.6.5"
+version = "0.7.0"
```

### 10.2 `markers.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/markers.yaml b/src/genbank_parser/data/mobilome/markers.yaml
--- a/src/genbank_parser/data/mobilome/markers.yaml
+++ b/src/genbank_parser/data/mobilome/markers.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-catalog_version: "1.2.0"
+catalog_version: "1.3.0"
```

### 10.3 `provenance.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/provenance.yaml b/src/genbank_parser/data/mobilome/provenance.yaml
--- a/src/genbank_parser/data/mobilome/provenance.yaml
+++ b/src/genbank_parser/data/mobilome/provenance.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-provenance_version: "1.2.0"
+provenance_version: "1.3.0"
```

### 10.4 `inference.yaml`

```diff
diff --git a/src/genbank_parser/data/mobilome/inference.yaml b/src/genbank_parser/data/mobilome/inference.yaml
--- a/src/genbank_parser/data/mobilome/inference.yaml
+++ b/src/genbank_parser/data/mobilome/inference.yaml
@@ -1,5 +1,5 @@
 schema_version: 1
-inference_version: "1.2.0"
+inference_version: "1.3.0"
```

---

## 11. Verification plan

### 11.1 Commands

```bash
# 1. Regenerate golden files
python -m pytest tests/ -x -k "mobilome" --tb=short 2>&1 | head -100

# 2. Full test suite
python -m pytest tests/ -x --tb=short

# 3. Lint
python -m ruff check src/ tests/
python -m ruff format --check src/ tests/

# 4. Type check
python -m mypy src/genbank_parser/mobilome/ --strict
```

### 11.2 Manual verification checklist

- [ ] `cas12a`, `cas8f`, `cas13d` gene names match the updated regex
- [ ] `cas15`, `casein`, `cascade` are rejected by the updated regex
- [ ] `virB1` through `virB11` match the new MPF gene pattern
- [ ] `trbB` through `trbL` match the new MPF gene pattern
- [ ] `trbA`, `trbW`, `virB0`, `virB12` are rejected
- [ ] `repX`, `repL` match replication-initiation casefold-exact list
- [ ] TA pair hypothesis on circular record reads "maximum gap (circular)"
- [ ] TA pair hypothesis on linear record reads "maximum gap (linear)"
- [ ] Handoff JSON objects include `source_ids` arrays
- [ ] Disabled-rule validator accepts `pxo_numbered_product_annotation` without
  the hardcoded string subtraction
- [ ] `missing_components` emits `MissingComponent` instances containing Component IDs and nested `missing_facets`
- [ ] `RepliconInventory.spatial_limitations` decouples locus spatial-pairing warnings from replicon classification limitations
- [ ] `crispr_repeat_array` and `crispr_cas_candidate` disambiguate marker IDs from facet IDs
- [ ] `database.py` rejects marker IDs that shadow facet IDs with `MobilomeDatabaseError`
- [ ] Schema validation passes with `jsonschema` on regenerated golden JSON
- [ ] Golden TSV handoff rows include `source_ids_json` column data

### 11.3 Golden regeneration

After applying all diffs, the golden files must be regenerated:

```bash
python -c "
from genbank_parser.mobilome import analyze_mobilome
from genbank_parser.mobilome.report import render_json, render_tsv
report = analyze_mobilome('tests/fixtures/mobilome_synthetic.gbff')
open('tests/golden/mobilome_v1.json', 'w').write(render_json(report))
open('tests/golden/mobilome_v1.tsv', 'w').write(render_tsv(report))
"
```

---

## 12. Commit plan

Commits should be atomic and ordered by dependency:

1. **`feat(mobilome): add repX and repL provenance entries`** — provenance.yaml only
2. **`feat(mobilome): expand replication initiators, MPF gene coverage, and cas regex`** — markers.yaml
3. **`fix(mobilome): resolve catalog ID hygiene and prevent marker/facet shadowing`** — markers.yaml, database.py, tests
4. **`fix(mobilome): remove hardcoded pxo whitelist from disabled-rule validator`** — database.py
5. **`feat(mobilome): add source_ids to handoff schema, model, and renderer`** — models.py, database.py, report.py, report.schema.json, inference.yaml
6. **`feat(mobilome): decouple spatial limitations into RepliconInventory`** — models.py, replicons.py, __init__.py, report.py, report.schema.json
7. **`feat(mobilome): standardize missing_components with nested missing_facets`** — models.py, inference.py, report.py, report.schema.json, __init__.py, tests
8. **`fix(mobilome): use topology-aware wording in TA gap limitation text`** — inference.py
9. **`chore(mobilome): bump versions to 1.3.0 / 0.7.0 and regenerate goldens`** — pyproject.toml, markers.yaml, provenance.yaml, inference.yaml, goldens

---

## 13. Cross-reference to deferred review

| Deferred Finding | Patch Section | Status |
|---|---|---|
| A3: Handoff `source_ids` | §6 | Addressed |
| A4: Validator hardcode | §2 | Addressed |
| A5: Mixed `missing_components` vocab | §7 | Addressed |
| Spatial limitations home | §8 | Addressed |
| Catalog ID hygiene | §9 | Addressed |
| C1: Data-driven rule engine | — | Deferred to v2 architecture |
| C2: TA gap wording | §1 | Addressed |
| C3: Strand-aware TA pairing | — | Deferred |
| C4: O(N×M) TA filtering | — | Deferred |
| C5: Tripartite TA systems | — | Deferred |
| C6: Annotation pipeline scope | — | No action needed |
| C7.1: virB / trb MPF | §4 | Addressed |
| C7.2: cas subtype regex | §3 | Addressed |
| C7.3: repX / repL | §5 | Addressed |
| TenpIN / CptIN markers | — | Deferred |

---

## 14. Execution Record

- **Executed**: 2026-09-06
- **Status**: Complete & fully verified
- **Changes Applied**:
  - §1: Topology-aware TA gap wording in `src/genbank_parser/mobilome/inference.py`.
  - §2: Removed hardcoded `pxo_numbered_product_annotation` whitelist from `_parse_disabled_rules` in `src/genbank_parser/mobilome/database.py`.
  - §3: Expanded Cas gene regex to `(?i)^cas(?:(?:[1-9]|1[0-4])[a-z]?|[A-F])$` in `src/genbank_parser/data/mobilome/markers.yaml`.
  - §4: Expanded MPF gene patterns with Ti-plasmid `virB1-11` and IncP `trbB-L` in `src/genbank_parser/data/mobilome/markers.yaml`.
  - §5: Added `repX` and `repL` replication initiators with verified provenance sources `anand-2008-repx` and `tinsley-2006-repx`.
  - §6: Added `source_ids` to `ExternalHandoff` model, schema, `inference.yaml`, and JSON/TSV report serializers.
  - §7: Introduced `MissingComponent` dataclass with nested `missing_facets`, updating model, inference, report, schema, and tests.
  - §8: Decoupled locus-level spatial pairing warnings into `RepliconInventory.spatial_limitations`.
  - §9: Renamed shadowed markers to `crispr_repeat_array` and `crispr_cas_candidate`, added validation against facet-shadowing in `database.py`, and updated test suites.
  - §10: Bumped package to 0.7.0, catalog/provenance/inference to 1.3.0, and updated CI wheel assertions and changelogs.
- **Verification Summary**:
  - Full pytest suite: 161 passed (1 known Biopython locus warning, 0 failures).
  - Scoped Ruff check and format: all checks passed cleanly on mobilome modules and test files.
  - Python compileall: passed cleanly across `src/` and `tests/`.
  - Pip check: no broken requirements found.
  - Git diff check: pure ASCII, zero conflict markers.
  - Goldens: regenerated with LF formatting for `mobilome_v1.json`, `.tsv`, and `.txt`.
