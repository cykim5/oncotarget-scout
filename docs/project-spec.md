# Project spec

## Goal

Produce a ranked, evidence backed target shortlist for STK11/KEAP1-mutant lung adenocarcinoma, with reason cards and a validation experiment plan, from public data only, in a repository a third party can rerun.

## Working rules

- Every number in a result file must be produced by a committed script. No hand edited result tables.
- No raw data in git. Download steps are scripted or documented in `data/README.md`.
- Claims are labeled by layer: expression claim, dependency claim, single cell claim, safety claim, outcome claim. A conclusion may combine layers only when each supporting layer is named.
- The grading checks in the README were written before the analysis. If a check fails, the result is reported and the scoring is revised in the open, with the change logged.

## Phases

### Phase 1, MVP (TCGA only)

Deliverables:
- `scripts/01_download_tcga.R` running end to end: LUAD RNA-seq counts, mutations, clinical data saved under `data/processed/`.
- Differential expression: tumor vs normal (DESeq2), with volcano and PCA figures in `results/figures/`.
- Pathway analysis: GSEA on Hallmark and Reactome, table and figure.
- Survival: Kaplan-Meier by median expression for the top expressed candidates, plus a Cox model table adjusted for stage and age.
- Full cohort target table v1: expression and survival layers only, clearly labeled as a partial score.
- Answer key check v1: rank positions of EGFR, KRAS, CD274 reported in `results/answer_key_check.md`.

Success criteria: the pipeline reruns from a clean clone with the documented steps, and the answer key genes are reported wherever they land, high or low.

### Phase 2, subgroup and function

Deliverables:
- Mutation stratification: STK11 mutant, KEAP1 mutant, both, neither; group sizes reported, small groups merged with the rule stated in advance if any group has under 20 tumors.
- Differential expression: mutant vs wild type within LUAD.
- DepMap layer: CRISPR Chronos scores joined by gene; dependency contrast between STK11/KEAP1-mutant lung lines and other lung lines; common essential genes flagged and penalized, because a gene every cell needs is not a subgroup target.
- Safety layer: GTEx median expression across vital tissues, summarized as a safety window value per gene.
- Subgroup target table v2 with the full scoring from `docs/target-score.md`.

Success criteria: the subgroup ranking is demonstrably different from the full cohort ranking, and the differences are explainable gene by gene.

### Phase 3, single cell, grading, and the lab plan

Deliverables:
- Single cell layer: one LUAD CELLxGENE dataset processed with Scanpy (QC, clustering, cell type annotation); per gene share of expression in malignant vs immune vs stromal compartments.
- Target table v3, final score, top 20 reason cards in `results/reason_cards/`, one short markdown file per gene generated from the result tables.
- Time split check: evidence restricted to a pre-2021 snapshot where the data source allows it, ranking rebuilt, later validated targets located in it. Where a source cannot be time frozen, the limitation is stated and the check is run on the layers that can.
- Validation experiment plan in `docs/validation-plan.md`: for the top 5 candidates, the DepMap cell lines carrying STK11/KEAP1 mutations, the proposed perturbation (CRISPR knockout or available inhibitor), the readout, and the confirm and reject outcomes.

### Phase 4, optional interface

- A small LLM question interface over the result tables: ask why a gene ranks where it does, receive an answer that cites the table values. Built last, and only over computed results, so the interface can never invent evidence.

## Time estimate

Phase 1 is a focused weekend to a two week part time build. Phase 2 is a similar block. Phase 3 is the longest, mostly the single cell processing. The project is presentable after Phase 1 and strong after Phase 2.

## Risks and mitigations

- Subgroup size: STK11 and KEAP1 groups in TCGA LUAD are each on the order of 80 to 100 tumors, workable for differential expression with effect size reporting, and the merged mutant group is larger. Mitigation: report effect sizes and uncertainty, rely on the dependency layer for confirmation.
- DepMap lung lines with both annotations and the right mutations are limited. Mitigation: use all lung lines with mutation calls, and report the line count behind every dependency claim.
- Single cell datasets vary in annotation quality. Mitigation: choose a CELLxGENE dataset with author cell type labels, and treat the layer as a location check, not a discovery step.
