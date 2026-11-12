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
