# Target score design

One score per gene, built from five layer values. Each layer value is scaled from 0 to 1 across all scored genes, so the layers can be combined and a reader can see exactly what moved a gene up or down.

## Layers

1. **Tumor specificity (weight 0.25).** Combines the tumor vs normal log fold change and the mutant vs wild type log fold change from DESeq2, with adjusted p values used as a gate rather than a ranking input. High expression everywhere does not score; enrichment in the subgroup does.
2. **Functional dependency (weight 0.30, the highest weight, on purpose).** From DepMap Chronos scores: how essential the gene is in STK11/KEAP1-mutant lung lines, minus how essential it is across all lines. A strongly negative selective value scores high. Pan-essential genes are penalized hard, because they offer no therapeutic window.
3. **Malignant cell share (weight 0.15).** From the single cell layer: the share of the gene's total tumor expression that sits in malignant cells. Immune compartment genes are interesting for other strategies and are labeled, but they do not score as direct tumor targets here.
4. **Safety window (weight 0.20).** From GTEx and normal lung cell references: 1 minus the gene's normalized expression in vital normal tissues. Quiet in heart, brain, and liver scores high.
5. **Outcome link (weight 0.10, the lowest weight, on purpose).** From the Cox model in the subgroup: association between high expression and shorter survival. This layer is supporting evidence only, because survival association is the layer most easily confounded.

Missing layers are never imputed as zero. The score reports a coverage column (how many of the 5 layers support the value), and genes are compared within the same coverage where a decision depends on it.

## Why these weights

The ranking philosophy, stated before any results: a target is a gene the cancer needs, that a drug can reach in the tumor, and that the patient can tolerate losing elsewhere. Dependency is the closest measurement of "needs," so it carries the most weight. Expression is necessary and noisy, so it is second. Safety decides real programs, so it is third. Cell location refines. Survival association decorates, so it is last.

## Sensitivity analysis

The final report includes the ranking rebuilt under two alternate weightings (equal weights, and expression first) plus a rank correlation between them. Candidates that stay in the top 20 under all three are labeled robust. Candidates that appear only under one weighting are labeled fragile. This is the same honesty standard the MCED cost-effectiveness project applies to ICERs: show which assumptions move the answer.

## Exclusions, decided in advance

- Genes with no DepMap coverage and no expression signal are excluded from the ranked list and listed separately as unscoreable.
- Known approved drug targets are never excluded. They are the answer key and the benchmark; they are flagged as known so a reader can separate rediscovery from novelty at a glance.
