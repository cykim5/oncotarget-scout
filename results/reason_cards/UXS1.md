# UXS1

UXS1 is an enzyme in sugar metabolism.

Shortlist rank 3 of 20 (overall rank 216 of 28,328). Score 0.634. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by dependency (0.921), safety (0.722). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.000. Log2 fold change, tumor vs healthy lung (shrunk): -0.356. Log2 fold change, mutant vs wild type tumors: -0.148.
  Mutant vs wild type adjusted p value: 1.35e-02.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.921. Selective essentiality: -0.318 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -0.671.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.834 (raw share of the gene's total tumor expression in malignant cells: 0.834; immune 0.095, stromal 0.023, other 0.048). Expressed in 31.6% of malignant cells, 7.1% of immune cells, 8.6% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.722 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.887. Cox hazard ratio per standard deviation of expression: 0.717 (log hazard ratio -0.333), p value 0.0138, n=132 mutant subgroup patients.

## Watch out

- Expression does not distinguish mutant tumors from wild type tumors; the rank rests on the other layers.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.