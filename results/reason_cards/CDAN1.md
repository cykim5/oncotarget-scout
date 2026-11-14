# CDAN1

Shortlist rank 17 of 20 (overall rank 684 of 28,328). Score 0.520. It is in the top 20 under the main weighting and some but not all alternate weightings.

## What the evidence says

The rank is driven most by dependency (0.663), safety (0.716). Coverage is 4 of 5 layers (expression+dependency+safety+singlecell).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.011. Log2 fold change, tumor vs healthy lung (shrunk): 0.292. Log2 fold change, mutant vs wild type tumors: 0.071.
  Mutant vs wild type adjusted p value: 2.38e-01.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.663. Selective essentiality: -0.229 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -0.907.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.821 (raw share of the gene's total tumor expression in malignant cells: 0.821; immune 0.089, stromal 0.030, other 0.060). Expressed in 18.0% of malignant cells, 1.9% of immune cells, 4.6% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.716 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Survival was not modeled for this gene (it was outside the Cox set), so this layer does not contribute to the score.

## Watch out

- No single weak point stands out in the tables; the usual caveats still apply (40 mutant lines behind the dependency claim, one single cell dataset, survival from 132 patients).

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.