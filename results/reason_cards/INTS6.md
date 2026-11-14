# INTS6

Shortlist rank 1 of 20 (overall rank 144 of 28,328). Score 0.652. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by dependency (0.924), safety (0.877). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.000. Log2 fold change, tumor vs healthy lung (shrunk): 0.088. Log2 fold change, mutant vs wild type tumors: 0.036.
  Mutant vs wild type adjusted p value: 6.44e-01.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.924. Selective essentiality: -0.319 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -1.047.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.713 (raw share of the gene's total tumor expression in malignant cells: 0.713; immune 0.154, stromal 0.058, other 0.075). Expressed in 27.6% of malignant cells, 7.7% of immune cells, 16.4% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.877 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.928. Cox hazard ratio per standard deviation of expression: 0.996 (log hazard ratio -0.004), p value 0.974, n=132 mutant subgroup patients.

## Watch out

- Expression does not distinguish mutant tumors from wild type tumors; the rank rests on the other layers.
- No significant survival link in the mutant subgroup, so this is a dependency hypothesis, not a prognosis story.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.