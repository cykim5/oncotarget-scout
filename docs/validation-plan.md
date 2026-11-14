# Validation experiment plan

The computational shortlist is a hypothesis list. This plan says exactly how to test it in a lab: which cell lines, which perturbation, what to measure, and what result would confirm or kill each candidate. The lines and the Chronos values below come from DepMap Public 24Q4 (results/tables/validation_lines.csv); the layer values come from the Phase 3 ranking (results/tables/shortlist_dependency_supported_phase3.csv).

## Design, shared across genes

- Perturbation: lentiviral CRISPR-Cas9 knockout, three independent sgRNAs per gene plus a non-targeting control sgRNA, puromycin selected. A gene counts as tested only if knockout is confirmed (western blot or qPCR for the target, or sequencing of the edited locus). A failed knockout is a technical failure, not a kill: repeat with new guides.

- Readout: viability by CellTiter-Glo at days 3, 7, 10, and 14 after selection, plus live-cell confluence imaging over the same window, and annexin V apoptosis staining at day 7. Three biological replicates per line per guide.

- Controls in every run: the non-targeting sgRNA in every line, and a known essential gene knockout (for example PLK1) as a positive control for the assay itself.

- Pre-specified confirm bar: at least 50 percent viability loss relative to the non-targeting control in at least two of the test lines, with under 20 percent loss in the wild type control lines, at day 14. The bar is set before the experiment runs, so a borderline result cannot be talked into a win.

- Pre-specified kill outcomes: no test line meets the bar (the dependency does not reproduce); or the wild type controls lose as much viability as the test lines (the gene is a general essential, not a subgroup target); or the effect is guide-specific (one guide kills, the other two do not), which points to an off-target effect.

- If a candidate confirms: dose-response with a second perturbation method (inducible knockdown or, where one exists, a small molecule inhibitor), then a mouse xenograft with the best test line.

## INTS6 (shortlist rank 1)

Why it is on the list: dependency layer value 0.924 (selective essentiality -0.319, mean Chronos -1.047 in mutant lung lines), safety window 0.877, malignant cell share 0.713 in the single cell layer. Overall Phase 3 score 0.652.

Test lines (STK11/KEAP1-mutant lung lines):

| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for INTS6 | Passes the -0.5 gate |
| --- | --- | --- | --- | --- | --- | --- |
| A549 | ACH-000681 | mutant | mutant | mutant | -2.494 | yes |
| Lu-65 | ACH-000438 | mutant | wild type | mutant | -2.544 | yes |
| NCI-H2170 | ACH-000481 | wild type | mutant | wild type | -2.491 | yes |
| HCC-44 | ACH-000667 | mutant | mutant | mutant | -2.451 | yes |

Control lines (wild type for STK11 and KEAP1):

| Cell line | Model ID | KRAS | Chronos for INTS6 |
| --- | --- | --- | --- |
| LK-2 | ACH-000769 | wild type | 0.277 |
| NCI-H1650 | ACH-000035 | wild type | 0.230 |
| NCI-H1581 | ACH-000015 | wild type | 0.229 |


Confirm: INTS6 knockout kills the test lines selectively, meeting the shared confirm bar above, with the strongest effect in the lines where DepMap shows the most negative Chronos. Kill: no selective effect, or equal killing in the wild type controls.

## NFE2L2 (shortlist rank 2)

Why it is on the list: dependency layer value 1.000 (selective essentiality -0.345, mean Chronos -0.615 in mutant lung lines), safety window 0.653, malignant cell share 0.769 in the single cell layer. Overall Phase 3 score 0.640.

Test lines (STK11/KEAP1-mutant lung lines):

| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for NFE2L2 | Passes the -0.5 gate |
| --- | --- | --- | --- | --- | --- | --- |
| NCI-H2122 | ACH-000311 | mutant | mutant | mutant | -1.488 | yes |
| A549_CRAF_KD | ACH-003012 | mutant | mutant | mutant | -1.276 | yes |
| LUDLU-1 | ACH-000390 | wild type | mutant | wild type | -1.688 | yes |
| NCI-H1355 | ACH-000666 | mutant | mutant | mutant | -1.414 | yes |

Control lines (wild type for STK11 and KEAP1):

| Cell line | Model ID | KRAS | Chronos for NFE2L2 |
| --- | --- | --- | --- |
| NCI-H1373 | ACH-000845 | mutant | 0.329 |
| A427 | ACH-000757 | mutant | 0.112 |
| LCAM1 | ACH-002035 | wild type | 0.061 |


Confirm: NFE2L2 knockout kills the test lines selectively, meeting the shared confirm bar above, with the strongest effect in the lines where DepMap shows the most negative Chronos. Kill: no selective effect, or equal killing in the wild type controls.

Orthogonal arm: NFE2L2 encodes NRF2, and the published small molecule NRF2 inhibitor ML385 can serve as a second, independent perturbation in the same lines. If knockout and ML385 agree, the case for the target is much stronger than either alone. If they disagree, the knockout result governs, because ML385 has its own off-target profile.

## UXS1 (shortlist rank 3)

Why it is on the list: dependency layer value 0.921 (selective essentiality -0.318, mean Chronos -0.671 in mutant lung lines), safety window 0.722, malignant cell share 0.834 in the single cell layer. Overall Phase 3 score 0.634.

Test lines (STK11/KEAP1-mutant lung lines):

| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for UXS1 | Passes the -0.5 gate |
| --- | --- | --- | --- | --- | --- | --- |
| NCI-H2122 | ACH-000311 | mutant | mutant | mutant | -1.158 | yes |
| DMS 53 | ACH-000698 | mutant | wild type | wild type | -1.681 | yes |
| NCI-H2126 | ACH-000785 | wild type | mutant | wild type | -2.176 | yes |
| NCI-H1355 | ACH-000666 | mutant | mutant | mutant | -1.072 | yes |

Control lines (wild type for STK11 and KEAP1):

| Cell line | Model ID | KRAS | Chronos for UXS1 |
| --- | --- | --- | --- |
| NCI-H211 | ACH-000639 | wild type | 0.181 |
| EPLC-272H | ACH-000585 | wild type | 0.132 |
| HCC-1359 | ACH-000901 | wild type | 0.091 |


Confirm: UXS1 knockout kills the test lines selectively, meeting the shared confirm bar above, with the strongest effect in the lines where DepMap shows the most negative Chronos. Kill: no selective effect, or equal killing in the wild type controls.

## NCAPG2 (shortlist rank 4)

Why it is on the list: dependency layer value 0.725 (selective essentiality -0.250, mean Chronos -1.208 in mutant lung lines), safety window 0.861, malignant cell share 0.834 in the single cell layer. Overall Phase 3 score 0.633.

Test lines (STK11/KEAP1-mutant lung lines):

| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for NCAPG2 | Passes the -0.5 gate |
| --- | --- | --- | --- | --- | --- | --- |
| NCI-H2030 | ACH-000521 | mutant | mutant | mutant | -2.030 | yes |
| NCI-H2286 | ACH-000912 | mutant | wild type | wild type | -1.604 | yes |
| LXF-289 | ACH-000787 | wild type | mutant | wild type | -2.455 | yes |
| NCI-H460 | ACH-000463 | mutant | mutant | mutant | -1.711 | yes |

Control lines (wild type for STK11 and KEAP1):

| Cell line | Model ID | KRAS | Chronos for NCAPG2 |
| --- | --- | --- | --- |
| NCI-H1882 | ACH-002490 | wild type | -0.339 |
| Calu-6 | ACH-000264 | mutant | -0.364 |
| LC-1-sq | ACH-002156 | wild type | -0.387 |


Confirm: NCAPG2 knockout kills the test lines selectively, meeting the shared confirm bar above, with the strongest effect in the lines where DepMap shows the most negative Chronos. Kill: no selective effect, or equal killing in the wild type controls.

## SPRR1B (shortlist rank 5)

Why it is on the list: dependency layer value 0.582 (selective essentiality -0.201, mean Chronos -0.819 in mutant lung lines), safety window 0.963, not detected in the single cell layer, so its tumor location is unconfirmed. Overall Phase 3 score 0.616.

Test lines (STK11/KEAP1-mutant lung lines):

| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for SPRR1B | Passes the -0.5 gate |
| --- | --- | --- | --- | --- | --- | --- |
| NCI-H1915 | ACH-000434 | mutant | mutant | wild type | -1.058 | yes |
| Lu-65 | ACH-000438 | mutant | wild type | mutant | -1.040 | yes |
| HCC-15 | ACH-000878 | wild type | mutant | wild type | -1.387 | yes |
| A549 | ACH-000681 | mutant | mutant | mutant | -0.937 | yes |

Control lines (wild type for STK11 and KEAP1):

| Cell line | Model ID | KRAS | Chronos for SPRR1B |
| --- | --- | --- | --- |
| NCI-H1882 | ACH-002490 | wild type | 0.360 |
| NCI-H841 | ACH-000292 | wild type | 0.218 |
| DMS 273 | ACH-000749 | wild type | 0.160 |


Confirm: SPRR1B knockout kills the test lines selectively, meeting the shared confirm bar above, with the strongest effect in the lines where DepMap shows the most negative Chronos. Kill: no selective effect, or equal killing in the wild type controls.

## What this plan does not cover

- The single cell layer rests on 1,953 LUAD cells from two donors (see results/tables/singlecell_qc.json), so the tumor location evidence is thinner than the dependency evidence. A candidate that confirms here still needs its expression checked in the actual test lines (qPCR) before the knockout result is interpreted.

- None of these lines is a perfect model of a patient tumor. A confirmed candidate is a lead for animal work, not a drug.


Sources: results/tables/validation_lines.csv, results/tables/shortlist_dependency_supported_phase3.csv, results/tables/singlecell_layer.csv, results/tables/dependency_layer.csv.
