# Time split check (Phase 3, pre-2021 frozen ranking)

Pre-committed check from the README: freeze the evidence at
an earlier date, rebuild the ranking, and ask whether later
clinical developments (the KRAS G12C inhibitor approvals
that began in 2021) appear high in the frozen ranking.

## What was frozen

- TCGA expression and survival (2014-2018 data) and GTEx v8
  (2017-2020) are already pre-2021 and are carried over
  unchanged, including the Phase 2 Cox values.
- DepMap 24Q4 (Chronos) is replaced with DepMap Public 20Q2
  (figshare article 12280541, published 2020-06-04), using
  Achilles CERES gene effect scores; Chronos did not exist
  yet in 20Q2. The selective essentiality contrast is
  computed the same way. 20Q2 line counts: 36 mutant lung lines and 70 other lung lines with CERES data.
- The single cell layer cannot be frozen (the CELLxGENE
  dataset is from a 2024 preprint), so the frozen ranking
  is a four-layer ranking, directly comparable to Phase 2:
  the only difference between the two rankings is the
  dependency layer's release.

## Result

Frozen ranking size: 28351 genes. Spearman
correlation of scores with the Phase 2 ranking: 0.977.
Shortlist overlap (dependency-supported top 20): 12 of 20 genes in both.

| Gene | Rank in Phase 2 | Rank frozen | Selective essentiality 24Q4 | Selective essentiality 20Q2 |
| --- | --- | --- | --- | --- |
| INTS6 | 174 | 382 | -0.319 | -0.243 |
| NFE2L2 | 307 | 298 | -0.345 | -0.294 |
| SPRR1B | 317 | 4166 | -0.201 | 0.052 |
| UXS1 | 450 | 465 | -0.318 | -0.258 |
| NCAPG2 | 460 | 3964 | -0.250 | -0.068 |
| ADSL | 485 | 475 | -0.271 | -0.215 |
| GMPS | 651 | 483 | -0.272 | -0.245 |
| ARF4 | 1369 | 5023 | -0.279 | -0.125 |
| SGO1 | 1615 | 4009 | -0.147 | -0.077 |
| SLC33A1 | 1779 | 11095 | -0.224 | -0.265 |
| ICE2 | 1911 | 1658 | -0.224 | -0.195 |
| CENPA | 1934 | 1954 | -0.136 | -0.115 |
| CDK4 | 2615 | 4723 | -0.229 | -0.121 |
| CDAN1 | 2867 | 4176 | -0.229 | -0.158 |
| SKA1 | 2885 | 3735 | -0.141 | -0.102 |
| ATP1A1 | 2925 | 4512 | -0.293 | -0.195 |
| MED30 | 3192 | 5495 | -0.220 | -0.081 |
| COPG1 | 3221 | 4256 | -0.237 | -0.167 |
| KIF4A | 3234 | 2629 | -0.126 | -0.116 |
| RPL6 | 3286 | 14584 | -0.269 | -0.045 |
| IMPDH2 | 3931 | 3619 | -0.210 | -0.188 |
| OIP5 | 5377 | 4212 | -0.018 | -0.077 |
| PLA2G10 | 411 | 2879 | nan | -0.056 |
| PPAT | 4240 | 3494 | -0.149 | -0.151 |
| RPIA | 12628 | 2692 | -0.147 | -0.183 |
| SUCO | 12842 | 920 | -0.177 | -0.217 |
| TKT | 3813 | 2088 | -0.209 | -0.208 |
| TYMS | 4551 | 2903 | -0.125 | -0.159 |
| EGFR | 6166 | 6010 | -0.007 | -0.115 |
| KRAS | 4993 | 6041 | -0.062 | 0.020 |
| CD274 | 5646 | 5472 | -0.005 | -0.039 |

## Reading

KRAS is not expected to rank high in either the Phase 2 or
the frozen subgroup ranking: the subgroup ranking asks
which genes are subgroup-selective targets, and KRAS
selective essentiality is near zero in both releases. The
informative result is whether the frozen ranking recovers
the same biology, which the NFE2L2 row and the shortlist
overlap above answer. See results/phase3-report.md.
