# Novelty and rationale

This note answers one question in full sentences, for interviews, applications, and the portfolio page: what is new and important about this project, compared with what already exists?

## The problem is real and specific

Lung cancer treatment has two big success stories, targeted therapy for patients with an actionable driver mutation, and immunotherapy for patients whose tumors engage the immune system. STK11/KEAP1-mutant lung adenocarcinoma sits outside both stories. These tumors tend to be immunologically cold, with low PD-L1 expression and low CD8 T cell infiltration, they respond poorly to checkpoint inhibitors and chemotherapy, and there is no approved therapy aimed at the alteration itself. Guidelines already recommend testing for these mutations at diagnosis, for prognosis and trial enrollment. So the patients are identified, the need is acknowledged, and the target list for them is short. That is a meaningful gap to work on.

## What others have already done, honestly

- Large platforms, above all Open Targets, integrate genetics, expression, pathway, drug, and literature evidence and score every target for a disease as a whole. That resource is excellent and this project uses the same evidence philosophy. It does not answer the subgroup question, and it does not produce an experiment.
- Thousands of TCGA analyses, papers, and student projects compute differentially expressed genes for a cancer and run enrichment analysis. The honest summary of that genre: a gene list, a pathway plot, and no way to tell which entries matter, no functional evidence, and no test of whether the method would have found anything real.
- DepMap based studies rank CRISPR dependencies, usually across a whole lineage, and show that selective dependencies exist. Fewer connect dependency back to a mutation defined patient subgroup, a safety window in normal tissue, patient survival, and a validation plan in one pipeline.
- Industry target discovery teams do all of this, behind closed doors, with proprietary data. Their work is not reproducible by a reader.

## The novel combination

No single step here is new. The combination, the subgroup, and the grading are the contribution:

1. **One decision, one subgroup.** Every layer answers the same narrow question: what should be tested next for STK11/KEAP1-mutant LUAD? Narrow scope is what makes the output usable.
2. **Expression must survive contact with function.** Requiring DepMap CRISPR support separates "the gene is on" from "the cancer needs the gene." Those are different claims, and drug programs fail when they are confused. This distinction is the scientific center of the project.
3. **Safety is scored at the same time, not after.** Normal tissue and normal lung cell type expression enter the same score, so the ranking prefers targets a drug could plausibly hit without hitting the patient everywhere.
4. **The method grades itself against the future.** The time split check borrows the evaluation style used in the Gene-ius project: build the ranking with yesterday's evidence, then see whether it contained the targets the field later validated. Most published rankings are never tested this way, because the test can fail.
5. **The last output is a lab plan.** Each top candidate is paired with DepMap cell lines that carry the relevant mutations, a knockout or inhibition assay, and the result that would confirm or reject the candidate. This is where an Automated Sciences background changes the project: the computational result is designed from the start to be executed and measured by an automated experiment, closing the loop between prediction and physical evidence.

## Why it is worth doing, for a career

The Boehringer Ingelheim Scientist III posting that prompted this project lists, as daily work: multi-modal omics analysis, processing public and internal datasets into biological insights, developing innovative analysis approaches when standard methods are not adequate, interpreting and presenting results to biologists, and FAIR, reproducible, documented analysis. This project is that job description performed in public, on public data:

- bulk and single cell NGS analysis, in R and Python, on Unix based tooling
- public oncology datasets used by working teams: TCGA, DepMap, CELLxGENE, GTEx
- pathway analysis as an interpretation step, not a decoration
- a scoring approach designed for the question, with its weights and its failure tests written down
- a GitHub repository a hiring manager can clone, read, and rerun

## What success looks like

- The full cohort ranking recovers known LUAD targets without being steered to them.
- The subgroup ranking is enriched for metabolism and oxidative stress biology, which would show the method is reading the subgroup's actual biology rather than recycling generic cancer genes.
- At least a few top subgroup candidates have strong DepMap support in STK11/KEAP1-mutant lines, a clean safety window, and no existing approved drug, which makes them genuinely new leads worth an experiment.
- Equally valuable and fully reportable: a clean negative result, for example the finding that expression based ranking adds nothing over dependency alone in this subgroup. The pre-committed grading makes a negative result a finding, not a failure of effort.

## One sentence version

For the lung cancer subgroup that immunotherapy fails, this project asks CRISPR screens, tumor genomes, single cell atlases, and normal tissue data to vote together on the next drug targets, proves its own ranking against targets the field validated later, and hands the lab a plan to test the winners.
