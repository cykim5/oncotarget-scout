# Answer key check, corrected ranking

The three genes below are approved drug targets in lung
adenocarcinoma. They are never excluded and never up-weighted.
Where they land is a check on the method, not an input to it.

The rule being checked was chosen because the renormalized
ranking rewarded missing data, not because of where it put
these three. Both positions are shown so the change is visible.

| Gene | Before | After | Percentile after |
| --- | --- | --- | --- |
| EGFR | 5182 of 28328 | 588 of 16711 | 96.5 |
| KRAS | 4104 of 28328 | 157 of 16711 | 99.1 |
| CD274 | 9694 of 28328 | 4656 of 16711 | 72.1 |

KRAS is the gene this subgroup is defined against, and it now
sits in the top one percent without being told to. EGFR is in
the top four percent. CD274, which encodes PD-L1, stays far
down the list, which is the correct result: it is a biomarker
for who responds to immunotherapy, not a gene this subgroup
depends on, and a method that ranked it high would be wrong.
