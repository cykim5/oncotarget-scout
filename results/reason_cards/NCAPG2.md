# NCAPG2

Shortlist rank 4 of 20 (overall rank 229 of 28,328). Score 0.633. It stays in the top 20 under all three weightings, so the rank does not depend on one weighting choice.

## What the evidence says

The rank is driven most by dependency (0.725), safety (0.861). Coverage is 5 of 5 layers (expression+dependency+safety+singlecell+outcome).

## Layer by layer

### Expression (TCGA)

- Tumor specificity layer value: 0.096. Log2 fold change, tumor vs healthy lung (shrunk): 1.768. Log2 fold change, mutant vs wild type tumors: 0.280.
  Mutant vs wild type adjusted p value: 4.07e-03.

### Dependency (DepMap CRISPR)

- Dependency layer value: 0.725. Selective essentiality: -0.250 (mean Chronos in STK11/KEAP1-mutant lung lines minus mean Chronos in other lung lines; more negative means the mutant lines need the gene more). Mean Chronos in mutant lung lines: -1.208.
  Behind this: 40 mutant lung lines and 86 other lung lines with CRISPR data. Common essential gene: no.

### Single cell location (CELLxGENE)

- Malignant cell share layer value: 0.834 (raw share of the gene's total tumor expression in malignant cells: 0.834; immune 0.120, stromal 0.012, other 0.033). Expressed in 28.7% of malignant cells, 5.9% of immune cells, 4.6% of stromal cells.

### Safety (GTEx)

- Safety window layer value: 0.861 (1 means quiet in heart, brain, and liver; 0 means loud).

### Survival (TCGA, mutant subgroup)

- Outcome layer value: 0.939. Cox hazard ratio per standard deviation of expression: 1.086 (log hazard ratio 0.083), p value 0.521, n=132 mutant subgroup patients.

## Watch out

- No significant survival link in the mutant subgroup, so this is a dependency hypothesis, not a prognosis story.

## Sources in this repository

results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/cox_phase3.csv, results/tables/dependency_layer.csv, results/tables/timesplit_comparison.csv, results/tables/dependency_layer_20q2.csv, results/tables/de_mutant_vs_wt.csv, results/tables/de_tumor_vs_normal_shrunk.csv.