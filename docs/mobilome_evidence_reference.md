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

- ToxI is represented as the RNA antitoxin in a Type III ToxIN correspondence, following [Blower et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/); the pair-rule wording also carries the general toxin-antitoxin classification context of [Qiu et al. (2022)](https://doi.org/10.1016/j.micres.2022.127159).
- The HEPN/MNT reference retains the reported MntA/HepT biochemical correspondence without generalizing it to every annotation match; see [Yao et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7641770/). Current classification reviews and the fused MNT-HEPN study of [Yao et al. (2026)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13000451/) document HEPN/MNT-type systems among Type VII toxin-antitoxin systems; the rule text deliberately assigns no type number to an annotation match.
- Helper-dependent mobilization wording is limited by the distinction between annotated component co-occurrence and demonstrated transfer described by [Smillie et al. (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2937521/); the modern extended-mobility framing of [Garcillán-Barcia et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12282955/) anchors the same caution without narrowing mobility to canonical conjugation.
- The external MOB-suite or PlasmidFinder handoff requires versioned database and threshold provenance; see [Robertson and Nash (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/) for MOB-suite performance and database/taxonomic limitations.
- AimR/AimP/AimX labels are retained as annotation correspondence only; [Erez et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5378303/) does not make a single new-genome annotation functional evidence.
- Integrative and conjugative element context follows [Johnson and Grossman (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5180612/); integration annotations in one record do not delimit an ICE.
- Integron correspondence follows [Mazel (2006)](https://doi.org/10.1038/nrmicro1462); an integron-like annotation is not a cassette-array inventory without sequence analysis.
- Antimicrobial-resistance gene mobility context follows [Partridge et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6148190/); association in that review does not establish that a generic resistance-like annotation is mobile.
- Virulence-factor determination context follows [Chen et al. (2005)](https://academic.oup.com/nar/article/33/suppl_1/D325/2505203) and the VFDB 2022 classification scheme of [Liu et al. (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728188/); version provenance is required for any comparison.
- Bakta output and database-version context is cited from [Schwengers et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544/) alongside the [Bakta documentation](https://bakta.readthedocs.io/en/latest/); observed typed input features take precedence over general tool claims.
- INSDC feature and source-qualifier meanings are taken from the [INSDC feature table](https://www.ncbi.nlm.nih.gov/genbank/feature_table/); observed input feature representation remains authoritative.

## External handoffs

The report never runs external programs. It records `not_run` handoffs for MOB-suite or PlasmidFinder, AMRFinderPlus, a curated virulence workflow, prophage/ICE callers, and HMM/domain or sequence workflows. A future imported result must retain its tool/database release, method, thresholds, coverage, taxonomic scope, and limitations.

Primary descriptions for the handoff tools carried in provenance: MOB-suite ([Robertson and Nash, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/)), PlasmidFinder ([Carattoli et al., 2014](https://journals.asm.org/doi/10.1128/aac.02412-14)), AMRFinderPlus ([Feldgarden et al., 2021](https://www.nature.com/articles/s41598-021-91456-0)), VFDB ([Chen et al., 2005](https://academic.oup.com/nar/article/33/suppl_1/D325/2505203); [Liu et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728188/)), prophage and element boundary callers such as PHASTER ([Arndt et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4987931/)) and geNomad ([Camargo et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/)), and sequence-level mobile-element callers such as ISEScan ([Xie and Tang, 2017](https://doi.org/10.1093/bioinformatics/btx433)) and MobileElementFinder ([Johansson et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7729385/)). External handoffs declare their `source_ids` linking directly to versioned provenance.

Records without a detected annotation pipeline retain a `null` pipeline name in
the schema-versioned `annotation_pipelines` array; this explicitly records that
the input did not declare a recognized pipeline rather than inferring one.

TSV uses the fixed 59-column header and lossless empty-cell semantics. To keep
rows readable, trailing empty cells are omitted from individual data rows;
consumers should map fields by the header and tolerate these deterministic
ragged rows.

The catalog intentionally has no TPR-specific marker in v1. Generic RepA/RepB,
replication-relaxation, and pXO annotations therefore remain observations and
do not select theta, rolling-circle, or another replication mechanism.

Record-level spatial pairing warnings (such as unlocatable features skipping
spatial pairing) are decoupled into `RepliconInventory.spatial_limitations`
rather than mixing classification and pairing caveats. They render in the text
report and are carried in the TSV `spatial_limitations_json` column, so all
three formats expose the same record-level pairing caveats.

Subtype-suffixed Cas gene names follow the updated evolutionary classification
of CRISPR-Cas systems ([Makarova et al., 2025](https://doi.org/10.1038/s41564-025-02180-8));
the gene and product matchers of `crispr_cas_candidate` both recognize the
suffixed nomenclature, while an annotation match still does not establish a
CRISPR-Cas subtype.

The repX-specific citations (Anand et al. 2008; Tinsley and Khan 2006) anchor
the repX pattern only; the remaining generic rep gene names (repA, repB,
repS, repF, repFR) remain local-policy annotation candidates.

Type III toxin-antitoxin superfamilies beyond ToxIN follow Blower et al. 2012
([NAR, PMC3401426](https://pmc.ncbi.nlm.nih.gov/articles/PMC3401426/)): the
three sequence-similarity families are ToxIN, TenpIN (named for the
*Photorhabdus luminescens* locus), and CptIN (named for the *Coprococcus
catus* GD/7 locus). `tenpn`/`tenpi` and `cptn`/`cpti` retain the toxin and RNA
antitoxin annotations of the latter two families, and their pair rules mirror
the ToxIN rule. An annotation match still does not establish Type III
activity, abortive infection, or phage resistance.

Tripartite (N-marker) toxin-antitoxin modules are now supported by spatial
cluster rules: toxin-antitoxin-chaperone (TAC) systems with HigBA plus a
SecB-like chaperone ([Bordes et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3100995/);
[Mets et al. 2024](https://doi.org/10.1016/j.chom.2024.05.003);
[Nakamoto et al. 2026](https://doi.org/10.1016/j.celrep.2026.117024)),
the omega-epsilon-zeta addiction module of pSM19035
([Camacho et al. 2002](https://doi.org/10.1515/BC.2002.191);
[Brzozowska et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3953265/);
[Dmowski et al. 2016](https://doi.org/10.1099/mic.0.000308)), and the MqsRAC
tripartite system, MqsR (RNase toxin) / MqsA (antitoxin) / MqsC (SecB-type
chaperone) ([Fernández-García et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10783111/)).
Retron antiphage modules pair the retron-encoded reverse transcriptase with
the msDNA region; the effector toxin is system-specific (RcaT in Retron-Sen2)
and is deliberately not required by the annotation rule
([Bobonis et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC11938430/)).

Toxin-antitoxin spatial inference is cluster-based rather than combinatorial:
eligible hits of one rule's markers are joined by single-linkage proximity
edges (different markers, within the configured same-record gap), and one
aggregate hypothesis is emitted per cluster that contains every required
marker. A tandem array therefore collapses into one cluster hypothesis with
all supporting hits instead of one hypothesis per cross-product pair, an
isolated toxin without a cognate emits nothing, and two-feature clusters keep
the previous pair wording and observed-gap limitation.

The repL gene-name pattern now has dedicated curation
(Catchpole and Dyke 1992, FEMS Microbiol Lett, PMID 1577254: repL disruption
abolishes pT48 replication and functional repL complements in trans;
Catchpole and Dyke 1991, Mol Microbiol, PMID 1906970; Projan et al. 1987,
J Bacteriol, PMID 2822666; Khan and Novick 1982, J Bacteriol, PMID 7056699;
Kwong et al. 2017, [Front Microbiol](https://pmc.ncbi.nlm.nih.gov/articles/PMC5703833/),
family-level PF05732/pSN2 anchor; Šprincová et al. 2005, Plasmid, PMID
15907537, RepL family outside the staphylococci). The remaining generic rep
gene names (repA, repB, repS, repF, repFR) stay local-policy annotation
candidates, and no rep annotation selects a replication mechanism.

The unified type IV secretion nomenclature review
([Christie et al.](https://doi.org/10.1093/femsre/fuaf069)) was released
online-first on 2025-12-31 and is assigned to the 2026 volume 50 print issue;
the provenance record carries the 2026 print year while keeping its
`christie-2025-t4ss` stable identifier.

An optional real-genome benchmark (`scripts/benchmark_mobilome.py`) runs the
packaged analysis over a local, gitignored `PF_NNT_reoriented.gbff` and
reports runtime, rendered output volume, and the section 16.9 calibration
expectations; `tests/test_mobilome_real_genome.py` performs the same checks
as an opt-in test that skips clearly when the file is absent. Neither ever
stages the real input. The standalone mobilome-parser prototype remains an
audit-time comparison artifact in its own repository; the native
`genbank_parser.mobilome` engine is canonical, and any archival of the
prototype happens in that repository in a separate commit.

Paywalled full-text papers held for curation live in the gitignored
`reference_pdfs/` directory; they are intentionally outside version control,
and the catalog cites the published records rather than repository PDFs.
