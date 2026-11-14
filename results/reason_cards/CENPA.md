# CENPA

Shortlist rank 12 of 20 (overall rank 485 of 28,328). Score 0.579. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by safety (0.979), single cell location (0.832). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.181. Log2 fold change, tumor vs healthy lung (shrunk): 3.606. Log2 fold change, mutant vs wild type tumors: 0.420.
  Mutant vs wild type adjusted p value: 2.99e-03.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.394. Selective essentiality: -0.136 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -1.013.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.832 (raw share of the gene's total tumor expression in malignant cells: 0.832; immune 0.134, stromal 0.025, other 0.008). Expressed in 7.0% of malignant cells, 1.2% of immune cells, 2.0% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.979 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.954. Cox hazard ratio per standard deviation of expression: 1.221 (log hazard ratio 0.199), p value 0.115, n=132 mutant subgroup patients.

## Watch out

- No significant survival link in the mutant subgroup, so this is a dependency hypothesis, not a prognosis story.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.