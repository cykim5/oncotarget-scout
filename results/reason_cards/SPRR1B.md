# SPRR1B

Shortlist rank 5 of 20 (overall rank 331 of 28,328). Score 0.616. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by safety (0.963), dependency (0.582). Coverage is 4 of 5 layers (expression+dependency+safety+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.237. Log2 fold change, tumor vs healthy lung (shrunk): 6.023. Log2 fold change, mutant vs wild type tumors: -3.578.
  Mutant vs wild type adjusted p value: 5.92e-28.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.582. Selective essentiality: -0.201 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -0.819.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Not detected in the single cell dataset (zero counts in the LUAD cells used), so this layer does not contribute to the score and the gene's tumor location is unconfirmed.

### Safety (GTEx)

- Safety window layer value: 0.963 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.968. Cox hazard ratio per standard deviation of expression: 1.371 (log hazard ratio 0.315), p value 0.0142, n=132 mutant subgroup patients.

## Watch out

- Not detected in the single cell dataset, so its tumor location is unconfirmed; the rank rests on the other four layers.
- The time split check found no dependency signal before 2021 (selective essentiality 0.052 in DepMap 20Q2), so the dependency case rests on the newer data alone.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.