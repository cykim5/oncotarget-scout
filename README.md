# OncoTarget Scout

**Finding drug targets for the lung cancer patients immunotherapy leaves behind.**

A reproducible target prioritization pipeline for one group of lung adenocarcinoma (LUAD) patients: tumors carrying STK11 or KEAP1 mutations, or both. These tumors are immunologically cold. They respond poorly to PD-1/PD-L1 immunotherapy and to chemotherapy, they carry worse outcomes across treatments, and no approved therapy targets the alteration itself. This is a named unmet need in the clinic, not a textbook example.

## Status

Phases 1 to 6 complete. All five evidence layers are in.

| Output | Where |
| --- | --- |
| Five-layer subgroup ranking | `results/`, written up in `results/phase3-report.md` |
| Decision shortlist, 20 genes led by INTS6, NFE2L2 and UXS1 | `results/tables/shortlist_dependency_supported_phase3.csv` |
| One plain-language reason card per gene | `results/reason_cards/` |
| Validation experiment plan | `docs/validation-plan.md` |
| Time split check against DepMap 20Q2 | `results/timesplit_check.md` |
| Methods and deviations | `results/run_log.md` |

Phase 2's report and shortlist stay in place as the four-layer record.

## Research question

For STK11/KEAP1-mutant LUAD, which genes are the best drug targets? And can a pipeline built only on public data show its ranking is trustworthy before anyone spends money testing it?

## What the project produces

1. A ranked shortlist of about 20 candidate target genes for this subgroup.
2. A reason card for each gene, giving the evidence behind its rank in plain language.
3. A validation experiment plan for the top candidates: which cell lines to test, which assay to run, and what result would confirm or kill the candidate.

Most target projects stop at item 1. This one is built to end at item 3.

## Why this subgroup

- STK11 is mutated in roughly 15 to 30 percent of non-small cell lung cancer and KEAP1 in roughly 17 to 23 percent. They often occur together and often alongside KRAS.
- Both are linked to low PD-L1 expression, low CD8 T cell infiltration and primary resistance to checkpoint inhibitors.
- Current NCCN guidance recommends testing for STK11 and KEAP1 at diagnosis. The subgroup is already identified in real clinical workflows.
- No FDA approved therapy targets STK11 or KEAP1 alterations directly. Work under study focuses on the metabolic and oxidative stress dependencies these tumors develop. That dependency angle is what this project searches for.

## The five evidence layers

A candidate has to earn its rank across five independent layers. A gene that looks good in one layer does not rank.

1. **Tumor activity.** TCGA bulk RNA-seq. Is the gene more active in LUAD tumors than in healthy lung, and more active in mutant tumors than in wild type tumors of the same cancer?
2. **Functional dependency.** DepMap CRISPR screens. Knock the gene out in lung cancer cell lines. Do mutant cells die or stop growing while other cells carry on? This turns correlation into evidence of dependency. Scores are Chronos gene effect, more negative meaning more essential.
3. **Right cell.** CZ CELLxGENE single cell data. Is the gene active in the malignant cells themselves rather than in the immune or stromal cells around them?
4. **Safety window.** GTEx and the Human Cell Atlas. Is the gene quiet in heart, brain, liver and normal lung? A target active everywhere is a safety problem waiting to happen.
5. **Patient outcome.** TCGA clinical data. Does high activity track with shorter survival in this subgroup?

Pathway analysis (GSEA, Hallmark and Reactome) runs across the ranked list to show what the winning genes do as a group. For this subgroup a meaningful result should speak the language of metabolism and oxidative stress, because that is the biology STK11 and KEAP1 loss creates.

## How the ranking is graded

Three checks, all fixed before the ranking runs.

1. **Answer key check.** The pipeline should recover established LUAD biology near the top of the full cohort ranking. EGFR, KRAS and CD274, which encodes PD-L1. Without being told to look for them.
2. **Time split check.** Freeze the evidence at an earlier date and rank. Then ask whether later developments show up high in the frozen ranking, such as the KRAS G12C inhibitors approved from 2021. A method that would have found tomorrow's targets yesterday is worth trusting on new names.
3. **Dependency check.** Top candidates must be supported by the CRISPR layer, not by expression alone. Expression-only candidates are reported separately and labeled.

### Check 1 failed, and fixing it changed the headline

The answer key check is the one that failed. KRAS sat at 4,104 of 28,328 and EGFR at 5,182. The top 20 was non-coding RNAs and pseudogenes.

The cause was in how the layers were combined. A gene's score was the weighted mean over the layers it actually had, renormalized each time so a missing layer was never imputed as zero. That is right for describing one gene and wrong for ranking genes against each other. The fewer layers a gene has, the higher its ceiling. A gene with expression evidence alone can reach 1.0. A gene with all five almost never can. Every one of the top 20 was missing the dependency layer.

The fix is a rule, not a tuning. A gene is ranked against other genes only if it carries the layer that holds the most weight, which is dependency. Genes without DepMap coverage are not deleted. They are ranked separately in `results/tables/no_dependency_layer_phase3.csv` and reported with their coverage, so nothing is hidden.

| Gene | Before | After |
| --- | --- | --- |
| KRAS | 4,104 of 28,328 | 157 of 16,711, top 1 percent |
| EGFR | 5,182 of 28,328 | 588 of 16,711, top 4 percent |
| CD274, which encodes PD-L1 | 9,694 of 28,328 | 4,656 of 16,711 |

KRAS and EGFR now land where established biology says they should. CD274 stays far down the list, which is the right answer. It predicts who responds to immunotherapy rather than naming something this subgroup depends on.

The decision output never moved. The shortlist already required dependency support, so the twenty genes anyone would have acted on are the same before and after. The headline was wrong while the decision was right. Worth saying plainly, because the headline is what a reader sees first.

Written up in `results/answer_key_check_corrected.md` and reproduced by `scripts/phase3_python/07_corrected_ranking.py`.

## What is different here

- **Subgroup first.** Resources like the Open Targets Platform score targets for a disease as a whole. This project scores them for one mutation-defined subgroup with a specific treatment failure. That is how treatment decisions are actually made.
- **Function, not just expression.** Many public TCGA projects end with a list of differentially expressed genes. Here that is layer 1 of 5. The CRISPR layer can drop a gene that is highly expressed and completely dispensable, and it can surface a dependency that expression alone would rank in the middle.
- **An experiment is the deliverable.** The validation plan maps to real DepMap cell lines carrying STK11/KEAP1 mutations. The computational work is designed to close the loop into a physical experiment.
- **Pre-committed grading.** The three checks are written down before the analysis runs, so the result can fail publicly. A ranking that cannot fail is a story, not evidence.

## From a ranked gene to a molecule

Phase 4 took the third-ranked gene, UXS1, as far as computation goes.

UXS1 is an enzyme with two solved crystal structures and no known inhibitor. An independent 2025 Cancer Research paper (PMID 40966319) found it to be a selective dependency in KEAP1-mutant lung cancer. The druggability check on the top 5 made it the most actionable candidate on structure.

| Where | What is in it |
| --- | --- |
| `docs/druggability-top5.md` | Open Targets tractability and ChEMBL inhibitor search for the top 5 |
| `docking-uxs1/` | AutoDock Vina docking of the natural ligands and NAD, holo and apo receptor |
| `screening-uxs1/` | A 100 compound pilot with a measured correction for scoring function size bias |
| `screening-full-uxs1/` | Every approved drug in ChEMBL against 4GLL: 2,757 compounds in 6.7 hours |
| `screening-2b69/` | The same library against 2B69, the structure that works |
| `docs/control-2b69-result.md` | Why the receptor was the problem |
| `docs/screen-2b69-findings.md` | What the corrected screen found |

### The screen failed its own control

Docked in the same box as the library, the enzyme's own substrate scored below two thirds of random approved drugs. A compound that beats the substrate has shown nothing.

Three repairs were tried: a corrected search box, a second scoring function, and the right protonation state for a ligand carrying about three negative charges. Each helped. The best combination still left the substrate at the 79th percentile against a bar set at the top few percent before the screen ran. The one lead that looked real, three sulfonylureas ranking together as a series, did not survive the corrected box.

### The control was wrong too, and the cause was mine

PDB 4GLL has no UDP bound. With the substrate site empty, ARG272 had relaxed into it, backbone included. Every run above docked into a site the receptor was occupying.

Re-running against PDB 2B69, the same enzyme at 1.21 A with UDP actually bound:

| Measurement | 4GLL | 2B69 |
| --- | --- | --- |
| Natural ligand percentiles | 30th, 55th, 58th | 99.5th, 99.0th, 88.4th |
| Score against molecular size, Spearman rho | -0.67 | +0.11 |

A blocked pocket cannot form real contacts, so the only thing separating compounds was bulk. That is why the ranking was substantially a ranking by size. Full comparison in `docs/control-2b69-result.md`.

### The rescreen

2,738 compounds against the open structure in 6 hours. Inside the control's size window the score no longer tracks molecular size, rho +0.098 against the +0.11 the control measured on a different set of molecules. 4.7 percent of the library now beats the substrate, where two thirds of it did before.

Of the 59 best scoring structures, 5 carry a negative group and contact all four arginines that grip the substrate's phosphate: betamethasone phosphate, carbenicillin indanyl, alvimopan, fluorescein and pemirolast. Vina scores hydrophobic contact and hydrogen bonds and has no electrostatics term at all. Recovering the site's own chemistry from a library never designed for it is what a working setup looks like.

None of the five is a lead. Two are prodrugs whose anionic group exists to be cleaved off. One is built not to cross membranes and UXS1 sits in the Golgi. One is a dye that binds most things.

**The approved drug library holds no credible binder for this site.** That is unsurprising, because oral drugs are rarely polyanions. The finding is that the protocol can now tell the difference and the 4GLL version demonstrably could not. Details in `docs/screen-2b69-findings.md`.

The target stands either way, because it was never selected by docking. The published biology is that losing UXS1 sensitizes these tumors, which points at a degrader or a knockdown rather than an active site inhibitor. That points in turn at a CRISPR experiment needing no compound at all.

## Checking the pipeline's own output

The Phase 5 defect survived three phases of my own review. Phase 6 asks whether an automated layer would have caught it. Three specialists with separate remits run in parallel. A fourth reads their outputs and looks for one cause behind them.

```bash
python3 scripts/run_verification_experiment.py
```

Against the ranking as it stood before the fix it raises six critical findings, including the one that matters: 0 of the top 20 carry the dependency layer against 59 percent of the population, p = 1.8e-08. Against the corrected ranking it raises one. That one is a false positive and it is left in deliberately. `docs/verification-layer.md` explains why silencing it would be tuning the checker to the answer.

The architecture result is the more interesting half and it is not the flattering one. At rule speed, running the specialists in parallel is **slower** than running them in sequence. 0.94x, because thread setup costs more than the work.

Once a check costs anything, where you cut the work matters more than whether you cut it. Splitting by agent caps at 1.80x, which is 9 judgments over a busiest specialist of 5. Splitting by judgment reaches 8.62x on the same checks, the same backend and the same threads.

Full numbers, the GIL caveat and what is simulated rather than measured: `docs/verification-layer.md`.

## Tests

```bash
python3 -m pytest tests/
```

The scoring rules live in `scripts/scoring_core.py` so the pipeline and the tests run the same code. `tests/test_scoring_core.py` covers the 0 to 1 scaling, the weighted score, the coverage count and the essentiality gate. It includes a regression test for the ranking bug above: a gene with two layers must not outrank a gene with five.

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
  screening-full-uxs1/       full approved drug screen against 4GLL and its controls
  screening-2b69/            the same screen against the open structure
  control-2b69/              the positive control that settled the receptor question
  results/                   generated tables, figures, and reports
  environment.yml            Python environment for the single cell and scoring steps
```

## Data sources (all public)

- TCGA LUAD, via the Genomic Data Commons and TCGAbiolinks: RNA-seq counts, somatic mutations, clinical outcomes.
- DepMap (Broad Institute): CRISPR gene effect (Chronos) and cell line mutation annotations.
- CZ CELLxGENE Discover: lung tumor single cell RNA-seq datasets.
- GTEx: normal tissue expression for the safety layer.
- Human Cell Atlas lung references: normal lung cell type expression.

Raw data files are not committed. `data/README.md` records the exact source and download step for each dataset so any reader can rebuild the inputs.

## Roadmap

- **Phase 1.** TCGA download, differential expression, pathway analysis, survival association, first full cohort ranking, answer key check. R plus Python.
- **Phase 2.** Mutant vs wild type analysis for the subgroup, DepMap dependency layer, safety layer, subgroup ranking.
- **Phase 3.** Single cell layer, time split check, reason cards, validation experiment plan.
- **Phase 4.** Druggability check on the top 5, AutoDock Vina docking of UXS1, and a virtual screen of every approved drug against it. The roadmap had reserved this slot for an LLM interface answering a biologist's question from the project's own tables. That is still a possible next step. The chemistry took the slot because the druggability check made UXS1 the first candidate worth testing against a real structure.
- **Phase 5.** Corrected the headline ranking so genes with missing layers no longer outrank genes with complete evidence. Added tests for the scoring rules.
- **Phase 6.** A verification layer over the pipeline's own output, measured against a single pass on the defect Phase 5 fixed. Then the receptor correction and the rescreen against 2B69.
- **Next.** The CRISPR experiment. It tests the claim the whole project rests on and it needs no compound.

## Limits, stated plainly

- Public cohort data is retrospective and mostly from one population mix. Findings are hypotheses for experiments, not clinical claims.
- TCGA treatment data is incomplete. Immunotherapy resistance is represented by the mutation subgroup, which is the clinically used proxy, rather than by measured response in this cohort.
- Cell line dependency is evidence about cell lines. The validation plan exists because this evidence has to be re-earned in a lab.
