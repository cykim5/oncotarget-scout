# ARF4

Shortlist rank 10 of 20 (overall rank 483 of 28,328). Score 0.582. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by dependency (0.808), safety (0.612). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.017. Log2 fold change, tumor vs healthy lung (shrunk): 0.427. Log2 fold change, mutant vs wild type tumors: -0.019.
  Mutant vs wild type adjusted p value: 7.68e-01.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.808. Selective essentiality: -0.279 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -0.829.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.770 (raw share of the gene's total tumor expression in malignant cells: 0.770; immune 0.130, stromal 0.057, other 0.043). Expressed in 32.8% of malignant cells, 10.6% of immune cells, 20.4% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.612 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.977. Cox hazard ratio per standard deviation of expression: 1.464 (log hazard ratio 0.381), p value 0.00423, n=132 mutant subgroup patients.

## Watch out

- The time split check did not keep it in the frozen pre-2021 shortlist, so its case is weaker in the older data.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.