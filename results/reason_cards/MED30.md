# MED30

Shortlist rank 18 of 20 (overall rank 799 of 28,328). Score 0.514. It is in the top 20 under the main weighting and some but not all alternate weightings.

## What the evidence says

The rank is driven most by dependency (0.639), safety (0.742). Coverage is 4 of 5 layers (expression+dependency+safety+singlecell).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.000. Log2 fold change, tumor vs healthy lung (shrunk): 0.071. Log2 fold change, mutant vs wild type tumors: -0.034.
  Mutant vs wild type adjusted p value: 6.82e-01.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.639. Selective essentiality: -0.220 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -2.096.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.815 (raw share of the gene's total tumor expression in malignant cells: 0.815; immune 0.105, stromal 0.042, other 0.038). Expressed in 14.7% of malignant cells, 2.4% of immune cells, 5.3% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.742 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Survival was not modeled for this gene (it was outside the Cox set), so this layer does not contribute to the score.

## Watch out

- The time split check did not keep it in the frozen pre-2021 shortlist, so its case is weaker in the older data.
- Expression does not distinguish mutant tumors from wild type tumors; the rank rests on the other layers.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.