# GMPS

Shortlist rank 6 of 20 (overall rank 338 of 28,328). Score 0.615. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by dependency (0.789), safety (0.731). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.026. Log2 fold change, tumor vs healthy lung (shrunk): 0.663. Log2 fold change, mutant vs wild type tumors: 0.060.
  Mutant vs wild type adjusted p value: 4.08e-01.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.789. Selective essentiality: -0.272 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -1.344.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.857 (raw share of the gene's total tumor expression in malignant cells: 0.857; immune 0.081, stromal 0.039, other 0.022). Expressed in 25.3% of malignant cells, 3.5% of immune cells, 11.2% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.731 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.968. Cox hazard ratio per standard deviation of expression: 1.368 (log hazard ratio 0.313), p value 0.0161, n=132 mutant subgroup patients.

## Watch out

- No single weak point stands out in the tables; the usual caveats still apply (40 mutant lines behind the dependency claim, one single cell dataset, survival from 132 patients).

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.