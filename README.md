# OncoTarget Scout

**Finding drug targets for the lung cancer patients immunotherapy leaves behind.**

A reproducible, multi-modal target prioritization project for one specific group of lung adenocarcinoma (LUAD) patients: tumors with STK11 and/or KEAP1 mutations. These tumors are immunologically cold, respond poorly to PD-1/PD-L1 immunotherapy and to chemotherapy, have worse outcomes across treatments, and have no approved therapy that targets the alteration itself. They are a real, named unmet need in the clinic, not a textbook example.

Status: Phase 2 complete (2026-11-12). The STK11/KEAP1 subgroup ranking across four layers (expression incl. mutant vs wild type, DepMap dependency, GTEx safety, subgroup survival) is in results/, written up in results/phase2-report.md, with methods and deviations in results/run_log.md. The decision shortlist (dependency-supported top 20, led by INTS6 and NFE2L2) is in results/tables/shortlist_dependency_supported.csv. Phase 1's fold change artifact is fixed by LFC shrinkage; see the dated correction note in results/phase1-report.md. Phase 3 (single cell layer, time split check, reason cards, validation plan) is next.

## Research question

For STK11/KEAP1-mutant LUAD, which genes are the best drug targets, and can a fully public-data pipeline prove its ranking is trustworthy before anyone spends money testing it in a lab?

## What the project produces

1. A ranked shortlist of about 20 candidate target genes for this patient subgroup.
2. A reason card for each gene, showing the evidence behind its rank in plain language.
3. A ready to run validation experiment plan for the top candidates (which cell lines to test, which knockout or drug assay, what result would confirm or kill the candidate).

Most target projects stop at item 1. This one is designed to end at item 3.

## Why this subgroup

- STK11 is mutated in roughly 15 to 30 percent of non-small cell lung cancer, KEAP1 in roughly 17 to 23 percent, and they often occur together, frequently with KRAS.
- These mutations are linked to low PD-L1 expression, low CD8 T cell infiltration, and primary resistance to immune checkpoint inhibitors.
- Testing for STK11 and KEAP1 at diagnosis is recommended in current NCCN guidance for prognostic assessment and trial enrollment, which means the subgroup is already identified in real clinical workflows.
- There are currently no FDA approved therapies that directly target STK11 or KEAP1 alterations. New approaches under study focus on the metabolic and oxidative stress dependencies these tumors develop. That dependency angle is exactly what this project searches for.

## The five evidence layers

A candidate has to earn its rank across five independent layers. A gene that looks good in only one layer does not rank.

1. **Tumor activity (TCGA, bulk RNA-seq, R).** Is the gene much more active in LUAD tumors than in healthy lung, and is it more active in STK11/KEAP1-mutant tumors than in wild type tumors of the same cancer?
2. **Functional dependency (DepMap, CRISPR).** When this gene is knocked out in lung cancer cell lines, do cells with STK11/KEAP1 mutations actually die or stop growing, while other cells are less affected? This layer turns correlation into evidence the cancer depends on the gene. Chronos gene effect scores are used, with more negative meaning more essential.
3. **Right cell (CZ CELLxGENE, single cell, Python).** Is the gene active in the malignant cells themselves, and not only in immune or stromal cells around them?
4. **Safety window (GTEx and the Human Cell Atlas).** Is the gene quiet in vital normal tissue (heart, brain, liver, normal lung cell types)? A target that is active everywhere is a safety problem waiting to happen.
5. **Patient outcome (TCGA clinical).** Is high activity of the gene linked to shorter survival in this subgroup? That ties the ranking to outcomes patients care about.

Pathway analysis (GSEA, Hallmark and Reactome) runs across the ranked list to show what the winning genes do as a group. In this subgroup, a meaningful result should speak the language of metabolism and oxidative stress, because that is the biology STK11 and KEAP1 loss creates.

## How the ranking is graded

Three checks, all decided before the ranking is run:

1. **Answer key check.** The pipeline should recover established LUAD biology in the full cohort ranking (for example EGFR, KRAS, and CD274, which encodes PD-L1) near the top, without being told to look for them.
2. **Time split check, in the style of a time machine evaluation.** Freeze the evidence at an earlier date, rank, then ask whether later clinical developments (for example the KRAS G12C inhibitor approvals that began in 2021) appear high in the frozen ranking. A method that would have found tomorrow's targets yesterday is worth trusting on new names.
3. **Dependency check.** Top candidates must be supported by the DepMap CRISPR layer, not by expression alone. Expression only candidates are reported separately and labeled as such.

## What is different here

- **Subgroup first.** Large resources such as the Open Targets Platform score targets for a disease as a whole. This project scores targets for one mutation defined subgroup with a specific treatment failure, which is how treatment decisions are actually made.
- **Function, not just expression.** Many public TCGA projects end with a list of differentially expressed genes. Differential expression is layer 1 of 5 here. The CRISPR dependency layer can remove a gene that is highly expressed and completely dispensable, and it can surface a metabolic dependency gene that expression alone would rank in the middle.
- **An experiment is the deliverable.** The output includes a validation plan mapped to real DepMap cell lines that carry STK11/KEAP1 mutations. The computational work is designed to close the loop into a physical experiment, which is the working style of autonomous and automated labs.
- **Pre-committed grading.** The three checks above are written down before the analysis runs, so the result can fail publicly. A ranking that cannot fail is a story, not evidence.

## Repository structure

```
oncotarget-scout/
  README.md                  this file
  docs/
    novelty-and-rationale.md why this project, and what is different
    project-spec.md          phases, deliverables, success criteria
    target-score.md          the scoring system, in detail
  data/
    README.md                data sources and download instructions (no raw data is committed)
  scripts/
    01_download_tcga.R       TCGA LUAD download via TCGAbiolinks (template)
    02_differential_expression.R  DESeq2 tumor vs normal and mutant vs wild type (template)
    03_survival_and_scoring.py    survival association and score assembly (template)
    README.md                how the pipeline is meant to run
    phase1_python/           the Python route that actually ran Phase 1
    phase2_python/           the Python route that ran Phase 2
  results/                   generated tables, figures, and reports (Phase 1 and 2 complete)
  environment.yml            Python environment for the single cell and scoring steps
```

## Data sources (all public)

- TCGA LUAD, via the Genomic Data Commons and TCGAbiolinks: RNA-seq counts, somatic mutations, clinical outcomes.
- DepMap (Broad Institute): CRISPR gene effect (Chronos) and cell line mutation annotations.
- CZ CELLxGENE Discover: lung tumor single cell RNA-seq datasets.
- GTEx: normal tissue expression for the safety layer.
- Human Cell Atlas lung references: normal lung cell type expression.

Raw data files are not committed to this repository. `data/README.md` records the exact source and download step for each dataset so any reader can rebuild the inputs.

## Roadmap

- Phase 1, MVP: TCGA download, differential expression, pathway analysis, survival association, first full cohort ranking, answer key check. R plus Python.
- Phase 2: mutant vs wild type analysis for the STK11/KEAP1 subgroup, DepMap dependency layer, safety layer, subgroup ranking.
- Phase 3: single cell layer, time split check, reason cards, validation experiment plan.
- Phase 4 (optional): a small LLM interface that answers a biologist's question about any ranked gene by citing the project's own result tables.

## Limits, stated plainly

- Public cohort data is retrospective and mostly from one population mix. Findings are hypotheses for experiments, not clinical claims.
- TCGA treatment data is incomplete, so immunotherapy resistance is represented by the mutation subgroup, which is the clinically used proxy, and not by measured response in this cohort.
- Cell line dependency is evidence about cell lines. The validation plan exists because this evidence has to be re-earned in a lab.
