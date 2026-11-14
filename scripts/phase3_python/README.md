# Phase 3 scripts (Python route)

Phase 3 adds the fifth evidence layer (single cell, from CELLxGENE),
runs the pre-committed time split check against DepMap 20Q2, and
produces the two deliverables the project is designed to end at: one
reason card per shortlist gene and the validation experiment plan.

Run in order, from the repository root, with the project venv:

```
.venv/bin/python scripts/phase3_python/00_download_phase3.py
.venv/bin/python scripts/phase3_python/01_singlecell.py
.venv/bin/python scripts/phase3_python/02_scoring_v3.py
.venv/bin/python scripts/phase3_python/03_shortlist_v3.py
.venv/bin/python scripts/phase3_python/04_timesplit.py
.venv/bin/python scripts/phase3_python/05_reason_cards.py
.venv/bin/python scripts/phase3_python/06_validation_plan.py
```

| Script | What it does | Main outputs |
| --- | --- | --- |
| 00_download_phase3.py | Downloads the CELLxGENE LUAD single cell H5AD and the DepMap Public 20Q2 files (figshare, size and md5 verified before the file lands) | data/raw/cellxgene, data/raw/depmap20q2 (gitignored) |
| 01_singlecell.py | Builds the malignant cell share layer from the CELLxGENE dataset: author labels mapped to malignant / immune / stromal / other, LUAD cells only, light QC, counts summed per compartment in chunks | singlecell_layer.csv, singlecell_qc.json, singlecell_compartments.csv |
| 02_scoring_v3.py | Five-layer ranking with the spec weights; Cox outcome layer re-selected on the v3 preliminary ranking (reuses Phase 2 values, runs only newly selected genes); sensitivity analysis (equal and expression-first weights); answer key v3; figure | target_scores_phase3.csv, sensitivity_phase3.csv, cox_phase3.csv, unscoreable_genes_phase3.csv, answer_key_check.md, top20_phase3.png |
| 03_shortlist_v3.py | Decision shortlist: top 20 by main score among genes with dependency support, robustness recomputed within the subset; KM figure for the top shortlist gene; stacked compartment shares figure | shortlist_dependency_supported_phase3.csv, km_shortlist_top_phase3.png, singlecell_shortlist_shares.png |
| 04_timesplit.py | Time split check: dependency layer rebuilt from DepMap 20Q2 CERES scores, frozen four-layer ranking (single cell layer excluded, it cannot be frozen), comparison with the Phase 2 ranking | dependency_layer_20q2.csv, depmap_line_summary_20q2.json, target_scores_frozen.csv, timesplit_comparison.csv, timesplit_check.md |
| 05_reason_cards.py | One plain-language reason card per v3 shortlist gene, generated from the result tables only | results/reason_cards/<GENE>.md (20 cards) |
| 06_validation_plan.py | Picks the test and control cell lines for the top 5 shortlist genes from DepMap 24Q4 (mutation status and per-line Chronos) and writes the validation plan | validation_lines.csv, docs/validation-plan.md |

Deviations from docs/project-spec.md are recorded in
results/run_log.md. The R templates in scripts/ are kept as design
references; the Python route is what actually ran, because the R /
Bioconductor stack is not installed in this environment.
