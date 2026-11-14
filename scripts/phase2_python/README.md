# Phase 2 scripts (Python route)

Phase 2 turns from the full LUAD cohort to the STK11/KEAP1-mutant
subgroup and adds the two layers Phase 1 was missing: DepMap CRISPR
dependency and GTEx safety. It also fixes the Phase 1 fold change
artifact with LFC shrinkage.

Run in order, from the repository root, with the project venv:

```
.venv/bin/python scripts/phase2_python/01_download_phase2.py
.venv/bin/python scripts/phase2_python/02_de_shrinkage_subgroup.py
.venv/bin/python scripts/phase2_python/02b_subgroup_de.py   # only if 02 is interrupted after Part A
.venv/bin/python scripts/phase2_python/03_depmap.py
.venv/bin/python scripts/phase2_python/04_gtex.py
.venv/bin/python scripts/phase2_python/05_scoring.py
.venv/bin/python scripts/phase2_python/06_dependency_shortlist.py
```

| Script | What it does | Main outputs |
| --- | --- | --- |
| 01_download_phase2.py | Downloads cBioPortal LUAD mutations (STK11, KEAP1, KRAS), GTEx v8 median TPM, DepMap Public 24Q4 files (figshare, size and md5 verified before the file lands) | data/raw/cbioportal, data/raw/gtex, data/raw/depmap24q4 (gitignored) |
| 02_de_shrinkage_subgroup.py | Tumor vs normal DE with apeglm LFC shrinkage; rebuilds the Phase 1 ranking with shrunk fold changes; mutant vs wild type DE within tumors (combined and per-gene, groups of 20+ only) | de_tumor_vs_normal_shrunk.csv, target_scores_phase1_shrunk.csv, de_mutant_vs_wt.csv, de_stk11_vs_wt.csv, de_keap1_vs_wt.csv, de_both_vs_wt.csv, mutation_groups.json |
| 02b_subgroup_de.py | Companion to 02: runs Part B (subgroup DE) alone with identical code, for use if 02 is interrupted after Part A | de_mutant_vs_wt.csv, de_stk11_vs_wt.csv, de_keap1_vs_wt.csv, de_both_vs_wt.csv |
| 03_depmap.py | Selective essentiality of every gene in STK11/KEAP1-mutant lung lines vs other lung lines; common essential flags | dependency_layer.csv, depmap_line_summary.json |
| 04_gtex.py | Safety window per gene from GTEx heart, brain, liver expression | safety_layer.csv |
| 05_scoring.py | Four-layer ranking with the spec weights, Cox outcome layer for the top 500 preliminary candidates in the mutant subgroup, sensitivity analysis (equal and expression-first weights), answer key v2, figures | target_scores_phase2.csv, sensitivity_phase2.csv, cox_phase2.csv, unscoreable_genes.csv, answer_key_check.md, 4 figures |
| 06_dependency_shortlist.py | Decision shortlist: top 20 by main score among genes with dependency support, robustness recomputed within the subset, KM figure for the top shortlist gene | shortlist_dependency_supported.csv, km_shortlist_top.png |

Deviations from docs/project-spec.md and the R templates are recorded
in results/run_log.md. The R templates in scripts/ are kept as design
references; the Python route is what actually ran, because the R /
Bioconductor stack is not installed in this environment.
