# OncoTarget Scout

**Finding drug targets for the lung cancer patients immunotherapy leaves behind.**

A reproducible, multi-modal target prioritization project for one specific group of lung adenocarcinoma (LUAD) patients: tumors with STK11 and/or KEAP1 mutations. These tumors are immunologically cold, respond poorly to PD-1/PD-L1 immunotherapy and to chemotherapy, have worse outcomes across treatments, and have no approved therapy that targets the alteration itself. They are a real, named unmet need in the clinic, not a textbook example.

Status: Phases 1 to 5 complete. All five evidence layers are in: the five-layer subgroup ranking is in results/, written up in results/phase3-report.md, with methods and deviations in results/run_log.md. The decision shortlist (dependency-supported top 20, led by INTS6, NFE2L2, and UXS1) is in results/tables/shortlist_dependency_supported_phase3.csv, with one plain-language reason card per gene in results/reason_cards/ and the validation experiment plan in docs/validation-plan.md. The time split check against DepMap 20Q2 is in results/timesplit_check.md. Phase 2's report and shortlist stay in place as the four-layer record. Phase 4 (2026-11-17) added a druggability check on the top 5, AutoDock Vina docking of UXS1, and a virtual screen of every approved drug against it. Phase 5 (2026-11-17) fixed the headline ranking, which had been scoring genes with missing layers above genes with complete evidence.

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

### Check 1 failed, and fixing it changed the headline

The answer key check is the one that failed. KRAS sat at 4,104 of 28,328 and EGFR at 5,182, and the top 20 was non-coding RNAs and pseudogenes.

The cause was in how the layers were combined. A gene's score was the weighted mean over the layers it actually had, with the weights renormalized each time, so that a missing layer was never imputed as zero. That is correct for describing one gene and wrong for ranking genes against each other, because the fewer layers a gene has, the higher its ceiling. A gene with expression evidence alone can reach 1.0; a gene with all five almost never can. Every one of the top 20 was missing the dependency layer.

The fix is a rule, not a tuning. A gene is ranked against other genes only if it has the layer that carries the most weight, which is dependency. Genes without DepMap coverage are not deleted. They are ranked separately in `results/tables/no_dependency_layer_phase3.csv` and reported with their coverage, so nothing is hidden.

| Gene | Before | After |
| --- | --- | --- |
| KRAS | 4,104 of 28,328 | 157 of 16,711, top 1 percent |
| EGFR | 5,182 of 28,328 | 588 of 16,711, top 4 percent |
| CD274, which encodes PD-L1 | 9,694 of 28,328 | 4,656 of 16,711 |

KRAS and EGFR now land where established biology says they should, without being told to look for them. CD274 stays far down the list, which is the right answer, because it predicts who responds to immunotherapy rather than naming something this subgroup depends on.

The decision output never moved. The shortlist already required dependency support, so the twenty genes anyone would have acted on are the same before and after. The headline was wrong while the decision was right, which is worth saying plainly, because the headline is what a reader sees first. The rule was chosen because the old one rewarded missing data. Where it puts KRAS is the check on that rule, not the reason for it.

Written up in `results/answer_key_check_corrected.md`, reproduced by `scripts/phase3_python/07_corrected_ranking.py`.

## What is different here

- **Subgroup first.** Large resources such as the Open Targets Platform score targets for a disease as a whole. This project scores targets for one mutation defined subgroup with a specific treatment failure, which is how treatment decisions are actually made.
- **Function, not just expression.** Many public TCGA projects end with a list of differentially expressed genes. Differential expression is layer 1 of 5 here. The CRISPR dependency layer can remove a gene that is highly expressed and completely dispensable, and it can surface a metabolic dependency gene that expression alone would rank in the middle.
- **An experiment is the deliverable.** The output includes a validation plan mapped to real DepMap cell lines that carry STK11/KEAP1 mutations. The computational work is designed to close the loop into a physical experiment, which is the working style of autonomous and automated labs.
- **Pre-committed grading.** The three checks above are written down before the analysis runs, so the result can fail publicly. A ranking that cannot fail is a story, not evidence.

## From a ranked gene to a molecule

Phase 4 took the third-ranked gene, UXS1, as far as computation goes. The result is negative, and it is the most useful thing in the repository.

UXS1 is an enzyme with two solved crystal structures and no known inhibitor, and an independent 2025 Cancer Research paper (PMID 40966319) found it to be a selective dependency in KEAP1-mutant lung cancer. The druggability check on the top 5 made it the most actionable candidate on structure.

| Where | What is in it |
| --- | --- |
| `docs/druggability-top5.md` | Open Targets tractability and ChEMBL inhibitor search for the top 5 |
| `docking-uxs1/` | AutoDock Vina docking of the natural ligands and NAD, holo and apo receptor |
| `screening-uxs1/` | A 100 compound pilot, with a measured correction for scoring function size bias |
| `screening-full-uxs1/` | Every approved drug in ChEMBL: 2,757 compounds, 2,065 distinct structures after salt forms collapse, scored in 6.7 hours |
| `docs/screen-full-findings.md` | What the screen found, and why the ranked output is not a hit list |
| `docs/next-steps-after-the-screen.md` | What to do instead |

The screen failed its own positive control. Docked in the same box as the library, the enzyme's own substrate scored below two thirds of random approved drugs, so a compound that beats the substrate has demonstrated nothing. Three repairs were tried: a corrected search box, a second scoring function, and the right protonation state for a ligand that carries about three negative charges. Each helped, and the best combination still left the substrate at the 79th percentile against a bar of the top few percent that was set before the screen ran. The one lead that looked real, three sulfonylureas ranking together as a chemical series, did not survive the corrected box.

So the target stands and the docking route closes. The published biology is that losing UXS1 is what sensitizes these tumors, which points at a degrader or a knockdown rather than an active site inhibitor, and at a CRISPR experiment that needs no compound at all.

## Checking the pipeline's own output

The Phase 5 defect survived three phases of my own review, so Phase 6 asks whether an automated layer would have caught it. Three specialists with separate remits run in parallel, and a fourth reads their separate outputs and looks for one cause behind them.

```bash
python3 scripts/run_verification_experiment.py
```

Against the ranking as it stood before the fix, it raises six critical findings, including the one that matters: 0 of the top 20 carry the dependency layer against 59 percent of the population, p = 1.8e-08. Against the corrected ranking, one. That one is a false positive, it is left in deliberately, and `docs/verification-layer.md` explains why silencing it would be tuning the checker to the answer.

The architecture result is the more interesting half, and it is not the flattering one. At rule speed, running the specialists in parallel is **slower** than running them in sequence, 0.94x, because thread setup costs more than the work. Once a check costs anything, where you cut the work matters more than whether you cut it: splitting by agent caps at 1.80x, which is exactly 9 judgments over a busiest specialist of 5, while splitting by judgment reaches 8.62x on the same checks, the same backend and the same threads.

Full numbers, the GIL caveat, and what is simulated rather than measured: `docs/verification-layer.md`.

## Tests

```bash
python3 -m pytest tests/
```

The scoring rules live in `scripts/scoring_core.py` so that the pipeline and the tests run the same code. `tests/test_scoring_core.py` covers the 0 to 1 scaling, the weighted score, the coverage count and the essentiality gate, and includes a regression test for the ranking bug above: a gene with two layers must not outrank a gene with five.

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
  tests/
    test_scoring_core.py     tests for the scoring rules, including the ranking regression
  scripts/
    scoring_core.py          the scoring rules, shared by the pipeline and the tests
    01_download_tcga.R       TCGA LUAD download via TCGAbiolinks (template)
    02_differential_expression.R  DESeq2 tumor vs normal and mutant vs wild type (template)
    03_survival_and_scoring.py    survival association and score assembly (template)
    README.md                how the pipeline is meant to run
    phase1_python/           the Python route that actually ran Phase 1
    phase2_python/           the Python route that ran Phase 2
    phase3_python/           the Python route that ran Phase 3
  docking-uxs1/              AutoDock Vina docking of UXS1 (Phase 4)
  screening-uxs1/            100 compound pilot screen (Phase 4)
  screening-full-uxs1/       full approved drug screen and its controls (Phase 4)
  results/                   generated tables, figures, and reports
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
- Phase 3: single cell layer, time split check, reason cards, validation experiment plan. Complete.
- Phase 4 (done, November 2026): druggability check on the top 5, AutoDock Vina docking of UXS1, and a virtual screen of every approved drug against it. The roadmap had reserved Phase 4 for an LLM interface that answers a biologist's question by citing the project's own tables. That is still a possible next step. The chemistry took the slot because the druggability check made UXS1 the first candidate worth testing against a real structure.
- Phase 5 (done, November 2026): corrected the headline ranking so that genes with missing layers no longer outrank genes with complete evidence, and added tests for the scoring rules.
- Phase 6 (done, November 2026): a verification layer over the pipeline's own output, with specialists that run in parallel and cross check each other, measured against a single pass on the defect Phase 5 fixed. `docs/verification-layer.md`.
- Next: the CRISPR experiment in `docs/next-steps-after-the-screen.md`. It tests the claim the whole project rests on, and it needs no compound.

## Limits, stated plainly

- Public cohort data is retrospective and mostly from one population mix. Findings are hypotheses for experiments, not clinical claims.
- TCGA treatment data is incomplete, so immunotherapy resistance is represented by the mutation subgroup, which is the clinically used proxy, and not by measured response in this cohort.
- Cell line dependency is evidence about cell lines. The validation plan exists because this evidence has to be re-earned in a lab.
