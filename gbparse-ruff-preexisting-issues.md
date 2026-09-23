# Pre-existing Ruff findings

Status: remediation complete; retained as an audit note
Audit date: 2026-09-23
Observed commit: `dae5267d5e7f5e6200bfdce725ac105065bba96c`
Comparison base: `1ef9b9d50aae52813b28bec412093c8b9f0a3f1f`
Resolution commit: `9df0de1` (`chore(quality): clear pre-existing Ruff findings`)

## Scope and evidence

The repository-wide check was run with:

```text
ruff check src tests --output-format concise
```

It reports **16 findings**: eight `I001` import-order findings, seven `F401`
unused-import findings, and one `F841` unused-local finding. Ruff reports that
15 are automatically fixable and one requires a manual edit.

These are pre-existing relative to the v0.9.3 patch: none of the affected
files appears in the v0.9.3 diff from the comparison base to the observed
commit. This establishes that the v0.9.3 implementation did not introduce
these findings; it does not claim that they originated in one particular older
commit.

## Resolution status

The proposed import cleanup and unused-symbol removal were applied in
`9df0de1`. The manual `categories` assignment was removed from the MEOR
scanner, and the import-only changes were kept separate from the v0.9.3
feature-hardening commits.

Current repository-wide verification:

```text
ruff check src tests
```

Result: **0 findings**.

## Finding inventory

| # | Location | Rule | Finding | Proposed treatment |
|---:|---|---|---|---|
| 1 | `src/genbank_parser/extract.py:2` | `I001` | Standard-library imports are not Ruff/isort ordered. | Reorder the import block. |
| 2 | `src/genbank_parser/extract.py:9` | `F401` | `.io.get_qual` is imported but unused. | Remove the import after confirming no dynamic use. |
| 3 | `src/genbank_parser/extract.py:10` | `F401` | `.model.GenBankDocument` is imported but unused. | Remove the import. |
| 4 | `src/genbank_parser/fasta.py:2` | `I001` | Standard-library imports are not Ruff/isort ordered. | Reorder the import block. |
| 5 | `src/genbank_parser/fasta.py:9` | `F401` | `.model.GenBankDocument` is imported but unused. | Remove the import. |
| 6 | `src/genbank_parser/meor/database.py:2` | `I001` | `re` is ordered after `from ...` standard-library imports. | Reorder the import block. |
| 7 | `src/genbank_parser/meor/report.py:2` | `I001` | The import block has Ruff/isort formatting drift. | Apply the formatter-only import correction. |
| 8 | `src/genbank_parser/meor/scanner.py:19` | `F841` | Local `categories = database.category_map` is assigned but never read. | Remove the dead assignment; retain the existing iteration over `database.categories`. |
| 9 | `tests/conftest.py:2` | `I001` | A blank line is missing between standard-library and third-party imports. | Reformat the import block. |
| 10 | `tests/test_meor_database.py:2` | `I001` | `re` is ordered after `from pathlib import Path`. | Reorder the import block. |
| 11 | `tests/test_meor_parity.py:2` | `I001` | The import block has Ruff/isort formatting drift. | Apply the formatter-only import correction. |
| 12 | `tests/test_model.py:5` | `F401` | `GenBankDocument` is imported but unused. | Remove the import. |
| 13 | `tests/test_model.py:5` | `F401` | `GenBankFeature` is imported but unused. | Remove the import. |
| 14 | `tests/test_model.py:5` | `F401` | `GenBankRecord` is imported but unused. | Remove the import. |
| 15 | `tests/test_neighborhood.py:3` | `I001` | Third-party and first-party imports are not ordered. | Reorder the import blocks without changing test behavior. |
| 16 | `tests/test_sequence_codon.py:6` | `F401` | `genbank_parser.sequence.parse_sequences` is imported but unused. | Remove the import. |

## Recommended remediation

Keep this cleanup separate from feature or release-hardening work. The safe
mechanical pass is:

```text
ruff check src tests --fix
```

Then review the resulting diff, especially the unused-import removals and the
manual `categories` deletion. The expected source behavior is unchanged: the
MEOR scanner already iterates over `database.categories`, while the unused
lookup is never consulted.

Suggested validation after the cleanup:

```text
ruff check src tests
python -m compileall -q src tests
python -m pytest -q tests/test_meor_database.py tests/test_meor_parity.py tests/test_model.py tests/test_neighborhood.py tests/test_sequence_codon.py tests/test_v085_cli.py
python -m pytest -q
git diff --check
```

The focused tests cover the modules whose imports or dead assignment are being
cleaned; the full suite protects the broader parser, MEOR, CLI, and release
contracts. Acceptance is zero Ruff findings, no test regressions, and a diff
limited to the eleven files represented in the inventory above.

Suggested commit boundary:

```text
chore(quality): clear pre-existing Ruff findings
```
