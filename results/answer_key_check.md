# Answer key check v1 (Phase 1, full LUAD cohort)

Pre-committed check from the README: the full cohort ranking should
recover established LUAD target biology (EGFR, KRAS, CD274) near the
top without being steered to it. Ranks are reported wherever they land.

Ranking size: 35273 genes (partial score, expression + survival
layers, 2 of 5; see results/tables/target_scores_phase1.csv).

| Gene | Rank | Partial score | log2FoldChange (tumor vs normal) | padj |
| --- | --- | --- | --- | --- |
| EGFR | 282 | 0.148 | 0.413 | 0.0262 |
| KRAS | 155 | 0.221 | 0.650 | 4.14e-09 |
| CD274 | 457 | 0.105 | -0.744 | 0.000138 |


# Answer key check v2 (Phase 2, STK11/KEAP1 subgroup ranking)

Same pre-committed check, now on the Phase 2 ranking, which
scores genes for the STK11/KEAP1-mutant subgroup across four
layers (expression incl. mutant vs wild type, DepMap
dependency, GTEx safety, subgroup survival). Ranks are
reported wherever they land, untuned.

Ranking size: 28328 genes (4 of 5 layers; the single
cell layer arrives in Phase 3; see results/tables/
target_scores_phase2.csv).

| Gene | Rank | Score | log2FC mutant vs WT | Selective essentiality | Safety window |
| --- | --- | --- | --- | --- | --- |
| EGFR | 6166 | 0.289 | -0.819 | -0.007 | 0.737 |
| KRAS | 4993 | 0.351 | 0.123 | -0.062 | 0.715 |
| CD274 | 5646 | 0.309 | -1.001 | -0.005 | 0.845 |

Phase 1 ranks for comparison (full cohort, expression +
survival only): KRAS #155, EGFR #282, CD274 #457 of 35,273.


# Answer key check v3 (Phase 3, five-layer subgroup ranking)

Same pre-committed check, now on the Phase 3 ranking, which
scores genes for the STK11/KEAP1-mutant subgroup across all
five layers (expression incl. mutant vs wild type, DepMap
dependency, GTEx safety, single cell malignant share,
subgroup survival). Ranks are reported wherever they land,
untuned.

Ranking size: 28328 genes (5 of 5 layers; see
results/tables/target_scores_phase3.csv).

| Gene | Rank | Score | log2FC mutant vs WT | Selective essentiality | Safety window | Malignant share (raw) |
| --- | --- | --- | --- | --- | --- | --- |
| EGFR | 5182 | 0.382 | -0.819 | -0.007 | 0.737 | 0.907 |
| KRAS | 4104 | 0.429 | 0.123 | -0.062 | 0.715 | 0.870 |
| CD274 | 9694 | 0.308 | -1.001 | -0.005 | 0.845 | 0.298 |

Phase 2 ranks for comparison (four layers): KRAS #4993, EGFR #6166, CD274 #5646 of 28,328.
