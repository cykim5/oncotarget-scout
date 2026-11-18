# What is at the top of the 2B69 screen

2738 compounds docked, 680 duplicate salt forms removed, 2058 distinct structures ranked.

## 1. Is the ranking still a size measurement

| Set | n | Spearman rho, score against heavy atoms | p |
| --- | --- | --- | --- |
| whole library | 2058 | -0.494 | 2.5e-127 |
| 31 to 43 heavy atoms, the control's window | 446 | +0.098 | 3.8e-02 |

The control reported +0.11 against size matched decoys. A strong
negative number across the whole library is expected and is not the
same failure as 4GLL: Vina sums over atom pairs, so more atoms earn
a more negative score almost mechanically. What matters is whether
the correlation survives inside a narrow size window, because that
is where a real binder has to win.

## 2. How many beat the substrate

97 of 2058 structures (4.7%) score better than UDP-glucuronate at -10.21.

In the broken 4GLL run about two thirds of the library beat the
substrate, which was the signal that the comparison was empty.

## 3. The three rankings

Raw score favours big molecules. Ligand efficiency, the score
divided by heavy atom count, favours small ones. The size matched
percentile asks the only fair question: how does this molecule do
against molecules its own size. Where the three disagree is where
the interesting compound usually is.

### Ranked by raw score

| # | Compound | ChEMBL | Heavy | Score | LE | Size matched pct |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Paliperidone | CHEMBL1621 | 31 | -12.15 | -0.392 | 100.0 |
| 2 | Pimozide | CHEMBL1423 | 34 | -12.14 | -0.357 | 99.8 |
| 3 | Rimonabant | CHEMBL111 | 30 | -11.99 | -0.400 | 99.7 |
| 4 | Darifenacin Hydrobromide | CHEMBL1200935 | 32 | -11.68 | -0.365 | 99.5 |
| 5 | Ezetimibe | CHEMBL1138 | 30 | -11.52 | -0.384 | 99.4 |
| 6 | Dabrafenib | CHEMBL2028663 | 35 | -11.51 | -0.329 | 99.2 |
| 7 | Ketanserin | CHEMBL51 | 29 | -11.45 | -0.395 | 99.4 |
| 8 | Ziprasidone Hydrochloride | CHEMBL1375743 | 28 | -11.41 | -0.407 | 99.3 |
| 9 | Mizolastine | CHEMBL94454 | 32 | -11.38 | -0.356 | 98.6 |
| 10 | Bictegravir Sodium | CHEMBL3989867 | 32 | -11.24 | -0.351 | 98.4 |
| 11 | Betamethasone Phosphoric Acid | CHEMBL1201207 | 32 | -11.15 | -0.349 | 98.2 |
| 12 | Avacopan | CHEMBL3989871 | 42 | -11.14 | -0.265 | 100.0 |
| 13 | Meclizine | CHEMBL1623 | 28 | -11.14 | -0.398 | 98.8 |
| 14 | Doxycycline | CHEMBL1200699 | 32 | -11.13 | -0.348 | 97.9 |
| 15 | Carbenicillin Indanyl Sodium | CHEMBL1200991 | 35 | -11.10 | -0.317 | 97.9 |
| 16 | Trifluperidol | CHEMBL15023 | 29 | -11.09 | -0.382 | 98.4 |
| 17 | Alvimopan | CHEMBL270190 | 31 | -11.09 | -0.358 | 97.7 |
| 18 | Almitrine | CHEMBL1183717 | 35 | -11.05 | -0.316 | 97.4 |
| 19 | Dicumarol | CHEMBL1466 | 25 | -11.02 | -0.441 | 99.5 |
| 20 | Atovaquone | CHEMBL1450 | 26 | -11.02 | -0.424 | 99.1 |
| 21 | Opicapone | CHEMBL1089318 | 27 | -10.97 | -0.406 | 98.7 |
| 22 | Lurasidone | CHEMBL1237021 | 35 | -10.96 | -0.313 | 97.2 |
| 23 | Flibanserin | CHEMBL231068 | 28 | -10.95 | -0.391 | 97.8 |
| 24 | Ticagrelor | CHEMBL398435 | 36 | -10.94 | -0.304 | 97.1 |
| 25 | Oxyphenisatine | CHEMBL245807 | 24 | -10.94 | -0.456 | 99.3 |

### Ranked by ligand efficiency

| # | Compound | ChEMBL | Heavy | Score | LE | Size matched pct |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Magnesium Chloride | CHEMBL3185229 | 1 | -1.68 | -1.676 | n/a |
| 2 | Strontium Chloride | CHEMBL2219640 | 1 | -1.20 | -1.195 | n/a |
| 3 | Potassium Iodide | CHEMBL1141 | 1 | -1.11 | -1.107 | n/a |
| 4 | Sodium Iodide I 131 | CHEMBL1201084 | 1 | -1.11 | -1.107 | n/a |
| 5 | Silver Nitrate | CHEMBL177367 | 4 | -3.59 | -0.896 | 30.8 |
| 6 | Calcium Acetate | CHEMBL1200800 | 4 | -3.33 | -0.832 | 23.1 |
| 7 | Talc | CHEMBL3989756 | 3 | -2.48 | -0.828 | n/a |
| 8 | Hydrotalcite Anhydrous | CHEMBL3833409 | 4 | -3.28 | -0.820 | 19.2 |
| 9 | Potassium Phosphate, Dibasic | CHEMBL1200459 | 5 | -4.05 | -0.811 | 41.4 |
| 10 | Zinc Sulfate | CHEMBL1200929 | 5 | -3.92 | -0.784 | 37.9 |
| 11 | Sodium Thiosulfate | CHEMBL2096650 | 5 | -3.63 | -0.725 | 31.0 |
| 12 | Foscarnet Sodium | CHEMBL754 | 7 | -5.05 | -0.721 | 31.9 |
| 13 | Methiodal | CHEMBL2110978 | 6 | -4.26 | -0.710 | 21.4 |
| 14 | Methiodal Sodium | CHEMBL2106907 | 6 | -4.17 | -0.695 | 17.9 |
| 15 | Chloral Hydrate | CHEMBL455917 | 7 | -4.75 | -0.679 | 26.6 |
| 16 | Calcium Oxybate | CHEMBL4594393 | 7 | -4.71 | -0.672 | 23.4 |
| 17 | Cysteine Hydrochloride | CHEMBL1200630 | 7 | -4.47 | -0.638 | 16.0 |
| 18 | Sulfanilamide | CHEMBL21 | 11 | -6.95 | -0.632 | 76.8 |
| 19 | Pyrithione Zinc | CHEMBL1200471 | 8 | -5.05 | -0.631 | 23.5 |
| 20 | Histamine Dihydrochloride | CHEMBL1533310 | 8 | -5.05 | -0.631 | 22.8 |
| 21 | Choline Salicylate | CHEMBL2104095 | 10 | -6.30 | -0.630 | 52.3 |
| 22 | Aminolevulinic Acid Hydrochloride | CHEMBL1200582 | 9 | -5.62 | -0.624 | 30.2 |
| 23 | Metformin Hydrochloride | CHEMBL1703 | 9 | -5.62 | -0.624 | 29.6 |
| 24 | Dextroamphetamine Adipate | CHEMBL1200782 | 10 | -6.21 | -0.621 | 46.8 |
| 25 | Amphetamine Aspartate | CHEMBL1200377 | 10 | -6.20 | -0.620 | 46.4 |

### Ranked by size matched percentile

| # | Compound | ChEMBL | Heavy | Score | LE | Size matched pct |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Paliperidone | CHEMBL1621 | 31 | -12.15 | -0.392 | 100.0 |
| 2 | Avacopan | CHEMBL3989871 | 42 | -11.14 | -0.265 | 100.0 |
| 3 | Carbamazepine | CHEMBL108 | 18 | -10.29 | -0.572 | 99.9 |
| 4 | Pimozide | CHEMBL1423 | 34 | -12.14 | -0.357 | 99.8 |
| 5 | Rimonabant | CHEMBL111 | 30 | -11.99 | -0.400 | 99.7 |
| 6 | Dicumarol | CHEMBL1466 | 25 | -11.02 | -0.441 | 99.5 |
| 7 | Darifenacin Hydrobromide | CHEMBL1200935 | 32 | -11.68 | -0.365 | 99.5 |
| 8 | Ketanserin | CHEMBL51 | 29 | -11.45 | -0.395 | 99.4 |
| 9 | Ezetimibe | CHEMBL1138 | 30 | -11.52 | -0.384 | 99.4 |
| 10 | Ziprasidone Hydrochloride | CHEMBL1375743 | 28 | -11.41 | -0.407 | 99.3 |
| 11 | Oxyphenisatine | CHEMBL245807 | 24 | -10.94 | -0.456 | 99.3 |
| 12 | Dabrafenib | CHEMBL2028663 | 35 | -11.51 | -0.329 | 99.2 |
| 13 | Atovaquone | CHEMBL1450 | 26 | -11.02 | -0.424 | 99.1 |
| 14 | Mianserin | CHEMBL6437 | 20 | -10.22 | -0.511 | 98.9 |
| 15 | Meclizine | CHEMBL1623 | 28 | -11.14 | -0.398 | 98.8 |
| 16 | Fluorescein Sodium | CHEMBL1628233 | 25 | -10.93 | -0.437 | 98.8 |
| 17 | Opicapone | CHEMBL1089318 | 27 | -10.97 | -0.406 | 98.7 |
| 18 | Paroxetine Mesylate | CHEMBL1200609 | 24 | -10.82 | -0.451 | 98.7 |
| 19 | Mizolastine | CHEMBL94454 | 32 | -11.38 | -0.356 | 98.6 |
| 20 | Niraparib | CHEMBL1094636 | 24 | -10.74 | -0.448 | 98.5 |
| 21 | Cyproheptadine Hydrochloride | CHEMBL1716 | 22 | -10.58 | -0.481 | 98.5 |
| 22 | Fezolinetant | CHEMBL3608680 | 25 | -10.91 | -0.436 | 98.5 |
| 23 | Pemirolast Potassium | CHEMBL1201050 | 17 | -9.60 | -0.565 | 98.5 |
| 24 | Indoramin Hydrochloride | CHEMBL2104941 | 26 | -10.93 | -0.420 | 98.5 |
| 25 | Trifluperidol | CHEMBL15023 | 29 | -11.09 | -0.382 | 98.4 |

### Agreement between the three top 25s

- raw and size matched percentile share 16 compounds
- raw and ligand efficiency share 0 compounds
- all three share 0 compounds

## 4. The natural ligands, for scale

| Ligand | Heavy | Score | Percentile against size matched decoys |
| --- | --- | --- | --- |
| UDP-xylose | 34 | -11.13 | 99.5 |
| UDP-galactose | 36 | -11.00 | 99.0 |
| UDP-glucuronate | 37 | -10.21 | 88.4 |
| UDP | 25 | -10.02 | 84.9 |
| UMP | 21 | -8.93 | 53.3 |

A drug worth looking at has to clear the substrate's standing
among molecules its own size, not just its raw number.

## 5. What clears the bar

13 structures sit at or above the 99th percentile among molecules their own size, which is where UDP-xylose and UDP-galactose sit.

| Compound | ChEMBL | Heavy | Score | LE | Pct | Peers |
| --- | --- | --- | --- | --- | --- | --- |
| Paliperidone | CHEMBL1621 | 31 | -12.15 | -0.392 | 100.0 | 612 |
| Pimozide | CHEMBL1423 | 34 | -12.14 | -0.357 | 99.8 | 454 |
| Rimonabant | CHEMBL111 | 30 | -11.99 | -0.400 | 99.7 | 648 |
| Darifenacin Hydrobromide | CHEMBL1200935 | 32 | -11.68 | -0.365 | 99.5 | 564 |
| Ezetimibe | CHEMBL1138 | 30 | -11.52 | -0.384 | 99.4 | 648 |
| Dabrafenib | CHEMBL2028663 | 35 | -11.51 | -0.329 | 99.2 | 386 |
| Ketanserin | CHEMBL51 | 29 | -11.45 | -0.395 | 99.4 | 697 |
| Ziprasidone Hydrochloride | CHEMBL1375743 | 28 | -11.41 | -0.407 | 99.3 | 733 |
| Avacopan | CHEMBL3989871 | 42 | -11.14 | -0.265 | 100.0 | 93 |
| Dicumarol | CHEMBL1466 | 25 | -11.02 | -0.441 | 99.5 | 812 |
| Atovaquone | CHEMBL1450 | 26 | -11.02 | -0.424 | 99.1 | 800 |
| Oxyphenisatine | CHEMBL245807 | 24 | -10.94 | -0.456 | 99.3 | 821 |
| Carbamazepine | CHEMBL108 | 18 | -10.29 | -0.572 | 99.9 | 719 |

## What this does not yet show

Nothing here looks at where the pose actually sits. A good score
with the molecule pressed against a box wall, or parked in some
other groove, is not a substrate site binder. That check runs
next, on the shortlist, and it is the one that killed the
4GLL shortlist.
