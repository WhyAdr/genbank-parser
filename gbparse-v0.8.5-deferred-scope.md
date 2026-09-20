# `genbank-parser` v0.8.5 Deferred Scope

**Target repository:** `WhyAdr/genbank-parser`  
**Audited baseline:** `67456fd58b05d0240d2d99394033b5376887d0f7` (`main`, package v0.8.2)  
**Companion implementation plan:** `gbparse-v0.8.5-patch.md`  
**Status:** binding exclusions for the v0.8.5 implementation  

## Explicitly deferred and out of scope

Do **not** implement, scaffold, advertise as imminent, or add placeholder CLI
commands for any of the following in v0.8.5:

1. `gbparse index` or any persistent cohort database;
2. a universal batch runner, resumable executor, or universal run-provenance
   manifest; keep the existing Bakta-focused `batch-summary` command;
3. `gbparse synteny` or multi-genome neighborhood alignment;
4. `gbparse reconcile` or automatic Bakta/PGAP/Prokka reconciliation;
5. `gbparse records` list/extract/split/filter commands;
6. `gbparse annotation-audit` or annotation-dialect detection;
7. beyond-parsing optional extras, including HMM search, Smith-Waterman
   verification, ORF calling, translation-start prediction, CRISPR sequence
   discovery, taxonomy, ANI, AMR/VF/BGC calling, or invocation of external
   analytical tools.

Also excluded:

- no SQLite, DuckDB, Parquet, or FTS index;
- no Nextflow/Snakemake executor;
- no remote accession fetcher;
- no call to `table2asn`; v0.8.5 may only generate candidate input files;
- no automatic biological correction or rewriting of source annotations;
- no change to MEOR or mobilome biological inference catalogs unless required
  to preserve compatibility with the shared I/O refactor.

## Scope-control rule

If implementation of an accepted v0.8.5 workstream appears to require one of
the deferred capabilities above, Luna must stop at the integration boundary,
document the dependency or follow-up opportunity, and keep the deferred
capability out of the patch. Deferred work must not be introduced indirectly
through an undocumented helper, hidden experimental flag, placeholder parser,
or bundled external-tool invocation.

These exclusions do not reject the ideas permanently. They separate future
roadmap candidates from the bounded v0.8.5 delivery surface.
