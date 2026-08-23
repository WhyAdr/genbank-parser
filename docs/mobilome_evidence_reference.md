# Mobilome evidence reference

`gbparse mobilome` is an annotation-supported, replicon-centric evidence report. It uses the canonical typed GenBank parser, retains exact field/value/pattern matches, and records the catalog, inference-policy, input, and resource hashes used for each result.

## Evidence tiers

- Strength 3 is an explicit typed structure or recognized structured annotation, such as an explicitly annotated origin-of-transfer feature.
- Strength 2 is an exact curated gene or bounded product correspondence.
- Strength 1 is an ambiguous product or note correspondence.

Strength is an ordinal annotation-specificity label, not a probability, activity measurement, or phenotype result. Notes and generic products cannot receive Strength 3. Pseudogene and partial observations remain visible, but do not satisfy functional components by default.

## Scope and cautious hypotheses

The command inventories every input record. Source `/plasmid` and `/chromosome` declarations support the inventory class; bounded record-name or description tokens are tentative, and a circular record without a declaration remains unresolved. Topology, length, GC content, Rep labels, and pXO labels do not classify a record.

The v1 inference policy can emit only configured annotation candidates, including a ToxN/ToxI proximity pair, a HepT/MntA proximity pair, and possible helper-dependent mobilization when separate records satisfy the required oriT/relaxase and T4CP/MPF component patterns. The report lists untested compatibility, cognate oriT-relaxase relationship, same-cell affiliation, expression, transfer, and co-transfer conditions.

Generic Rep, pXO-numbered products, and phage-module observations do not select a replication mechanism or element identity. The pXO-like and phage-module aggregate rules are intentionally disabled in v1, with their future enablement requirements carried in every report. Generic AMR-like, virulence-associated, Zot-like, and AimR/AimP/AimX annotations remain candidates only.

## Source correspondence

- ToxI is represented as the RNA antitoxin in a Type III ToxIN correspondence, following [Blower et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/).
- The HEPN/MNT reference retains the reported MntA/HepT biochemical correspondence without generalizing it to every annotation match; see [Yao et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/).
- Helper-dependent mobilization wording is limited by the distinction between annotated component co-occurrence and demonstrated transfer described by [Smillie et al. (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/).
- The external MOB-suite or PlasmidFinder handoff requires versioned database and threshold provenance; see [Robertson and Nash (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/) for MOB-suite performance and database/taxonomic limitations.
- AimR/AimP/AimX labels are retained as annotation correspondence only; [Erez et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5378303/) does not make a single new-genome annotation functional evidence.
- INSDC feature and source-qualifier meanings are taken from the [INSDC feature table](https://www.ncbi.nlm.nih.gov/genbank/feature_table/); observed input feature representation remains authoritative.
- Bakta capabilities and output forms are version-sensitive; use the observed typed record and feature annotations alongside the [Bakta documentation](https://bakta.readthedocs.io/en/latest/).

## External handoffs

The report never runs external programs. It records `not_run` handoffs for MOB-suite or PlasmidFinder, AMRFinderPlus, a curated virulence workflow, prophage/ICE callers, and HMM/domain or sequence workflows. A future imported result must retain its tool/database release, method, thresholds, coverage, taxonomic scope, and limitations.

Records without a detected annotation pipeline retain a `null` pipeline name in
the schema-versioned `annotation_pipelines` array; this explicitly records that
the input did not declare a recognized pipeline rather than inferring one.

TSV uses the fixed 58-column header and lossless empty-cell semantics. To keep
rows readable, trailing empty cells are omitted from individual data rows;
consumers should map fields by the header and tolerate these deterministic
ragged rows.

The catalog intentionally has no TPR-specific marker in v1. Generic RepA/RepB,
replication-relaxation, and pXO annotations therefore remain observations and
do not select theta, rolling-circle, or another replication mechanism.
