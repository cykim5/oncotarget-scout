# UXS1 virtual screen — approved-drug repurposing library

A small virtual screen of approved small-molecule drugs against the
UXS1 (UDP-xylose synthase / UDP-glucuronate decarboxylase 1) substrate
site, as a repurposing starting point after the docking run in
`../docking-uxs1/`. That run found no known UXS1 inhibitor (ChEMBL
target search returned no record) and established the holo setup used
here: the crystal NAD retained in the receptor, box centered on the
substrate site.

**This is a prediction, not a measurement.** Docking scores estimate
binding free energy from a scoring function. They do not prove a
compound inhibits UXS1, and they say nothing about cell permeability,
selectivity, or toxicity. Any hit here needs biochemical testing
(enzyme inhibition assay) before it means anything in the lab.

## Library

- **Source:** ChEMBL web services (molecule endpoint), queried Nov 16,
  2026: `max_phase=4` (approved), `molecule_type=Small molecule`.
- **Filters:** canonical SMILES present, preferred name present,
  150 <= full molecular weight <= 600 Da. Filtered pool: 2,757
  compounds.
- **Sample:** 100 compounds drawn with a fixed random seed (42), so the
  library is reproducible. Every compound is real, with its ChEMBL ID,
  preferred name, canonical SMILES, and molecular weight recorded in
  `library.csv`. No compound identities were invented.
- **Preparation:** RDKit ETKDGv3 3D conformer (one per compound, seed
  42), MMFF-optimized (UFF fallback), salts stripped to the largest
  fragment; Meeko `mk_prepare_ligand.py` (Gasteiger charges), the same
  preparation as the docking-uxs1 run. 100 of 100 compounds prepared
  successfully (see `logs/prep.log`).

## Method

- **Receptor:** `receptor/uxs1_holo.pdbqt`, copied from
  `../docking-uxs1/receptor/`: PDB 4GLL chain A with the crystal NAD
  kept as rigid receptor atoms. The NAD occupies its pocket, so ligands
  can only sample the substrate site.
- **Box:** center (60.0, 23.2, 39.0), size 24 x 24 x 24 A. Same as
  docking-uxs1 run 2: midpoint of the crystallographic UDP site (2B69
  UDP 900, transformed into the 4GLL frame) and the catalytic triad
  centroid (Thr202/Tyr231/Lys235 in 4GLL numbering).
- **Docking:** AutoDock Vina 1.2.7, 9 poses per ligand, energy range
  3 kcal/mol.
- **Pass 1:** all 100 ligands at exhaustiveness 4 (reduced, first
  pass). Driver: `run_screen_pass1.py`. Outputs in `outputs/pass1/`.
- **Pass 2:** the top 10 pass-1 hits re-docked at exhaustiveness 8.
  Driver: `run_redock_top10.py`. Outputs in `outputs/top10/`.
- **Reference points** from docking-uxs1 run 2 (same setup,
  exhaustiveness 16): the natural product UDP-xylose scored -8.796 and
  the natural substrate UDP-glucuronate -7.337 kcal/mol. A repurposing
  hit scoring near or better than the substrate is worth a second
  look, not a conclusion.

## Results — pass 1 (exhaustiveness 4)

100 of 100 compounds docked. Best Vina affinity per compound
(kcal/mol; more negative = stronger predicted binding). Top 10:

| Rank | Compound | ChEMBL ID | Best affinity (kcal/mol) | MW (Da) |
| --- | --- | --- | --- | --- |
| 1 | GLIMEPIRIDE | CHEMBL1481 | -9.735 | 490.63 |
| 2 | IRINOTECAN | CHEMBL481 | -9.402 | 586.69 |
| 3 | LAPATINIB | CHEMBL554 | -9.008 | 581.07 |
| 4 | AVAPRITINIB | CHEMBL4204794 | -8.530 | 498.57 |
| 5 | BALOXAVIR MARBOXIL | CHEMBL4297503 | -8.522 | 571.56 |
| 6 | PRALATREXATE | CHEMBL1201746 | -8.420 | 477.48 |
| 7 | TOVORAFENIB | CHEMBL3348923 | -8.362 | 506.30 |
| 8 | INDIGOTINDISULFONATE | CHEMBL1091250 | -8.355 | 422.40 |
| 9 | UMBRALISIB | CHEMBL3948730 | -8.279 | 571.56 |
| 10 | CICLESONIDE | CHEMBL2040682 | -8.248 | 540.70 |

The weakest scorers were all small molecules (bottom 5: -4.03 to
-4.94 kcal/mol, MW 186-258 Da), and all top 10 are MW 422-587 Da:
the size skew in Vina's scoring is visible in this screen, so part
of the top-10 ranking is a size effect, not binding specificity.

Full ranked list: `outputs/pass1/results_pass1.csv`.

## Results — pass 2 (top 10 re-docked at exhaustiveness 8)

| Rank | Compound | ChEMBL ID | Best affinity exh 8 (kcal/mol) | Pass 1 exh 4 (kcal/mol) | MW (Da) |
| --- | --- | --- | --- | --- | --- |
| 1 | IRINOTECAN | CHEMBL481 | -9.283 | -9.402 | 586.69 |
| 2 | LAPATINIB | CHEMBL554 | -9.196 | -9.008 | 581.07 |
| 3 | BALOXAVIR MARBOXIL | CHEMBL4297503 | -8.533 | -8.522 | 571.56 |
| 4 | AVAPRITINIB | CHEMBL4204794 | -8.515 | -8.530 | 498.57 |
| 5 | INDIGOTINDISULFONATE | CHEMBL1091250 | -8.363 | -8.355 | 422.40 |
| 6 | TOVORAFENIB | CHEMBL3348923 | -8.339 | -8.362 | 506.30 |
| 7 | PRALATREXATE | CHEMBL1201746 | -8.281 | -8.420 | 477.48 |
| 8 | UMBRALISIB | CHEMBL3948730 | -8.279 | -8.279 | 571.56 |
| 9 | CICLESONIDE | CHEMBL2040682 | -8.122 | -8.248 | 540.70 |
| 10 | GLIMEPIRIDE | CHEMBL1481 | -6.772 | -9.735 | 490.63 |

Full ranked list: `outputs/top10/results_top10.csv`.

### Reproducibility: 9 of 10 held, 1 did not

Nine of the ten re-docked within 0.19 kcal/mol of their pass-1 score.
The exception is glimepiride, the pass-1 leader: -9.735 at
exhaustiveness 4, but -6.772 at exhaustiveness 8, a 2.96 kcal/mol
drop that puts it last. Its pass-1 score was a lucky pose from the
coarse stochastic search, not a stable minimum. This is why the
re-dock step exists, and why the pass-1 ranking alone should not be
trusted for flexible ligands. Irinotecan and lapatinib are the
reproducible leaders (-9.283 and -9.196), both scoring better than
the natural product reference UDP-xylose (-8.796) — with the size
caveat above, since both are ~580 Da. Ranks 3-9 span only
0.41 kcal/mol, far inside Vina's error: do not over-read that order.

## What the numbers do and do not mean

- The pass-1 ranking at exhaustiveness 4 is a coarse first filter.
  Vina's search is stochastic and a low exhaustiveness can miss the
  best pose, so a compound ranked 11th is not meaningfully worse than
  one ranked 10th. The pass-2 re-dock at exhaustiveness 8 is the
  better estimate for the top 10.
- Compare hits against the substrate reference (-7.337) and product
  reference (-8.796) from the earlier run, not against each other in
  isolation. Do not over-read gaps under ~1 kcal/mol: that is well
  inside Vina's error.
- A good docking score is a reason to test a compound, not evidence
  it works. Repurposing candidates have known safety profiles, which
  is the practical advantage, but UXS1 inhibition still has to be
  measured.

## Limits

- Docking scores are predictions from Vina's empirical scoring
  function, not measured binding affinities. Expect errors of roughly
  2-3 kcal/mol.
- Vina scores skew more negative for larger ligands (more atoms make
  more contacts), so part of any ranking is a size effect. Molecular
  weights are in the results CSVs for that comparison.
- Rigid receptor: the protein does not move during docking. Real
  binding can involve side-chain and loop rearrangement.
- No explicit solvent, no entropy estimate, no dynamics (no molecular
  dynamics was run). Docking gives static poses, not a trajectory, so
  there is no video to make from this run.
- Ligands were docked in a single RDKit-generated 3D conformation
  with ChEMBL's as-drawn protonation state. At physiological pH,
  ionizable groups may carry different charges; that changes
  hydrogen-bond donor/acceptor assignment and can shift scores.
- Vina's search is stochastic (no fixed seed in the 1.2.7 Python API),
  so repeat runs can differ slightly in pose and score.
- Only chain A was used. UXS1 is a dimer in the crystal; the active
  site is within one subunit, but dimer-interface effects are not
  modeled.
- Only 100 of ~2,757 approved small molecules were screened, a random
  sample. A hit here is a starting point; the unscreened majority may
  hold better binders.

## Files

- `fetch_library.py` — ChEMBL query + sampling (reproducible, seed 42)
- `library.csv` — the 100-compound library (ChEMBL ID, name, SMILES, MW)
- `prepare_ligands.py` — RDKit 3D conformers + Meeko PDBQT preparation
- `ligands/<chembl_id>.sdf` — 3D structures
- `ligands/<chembl_id>.pdbqt` — prepared ligands for Vina
- `receptor/uxs1_holo.pdbqt` — holo receptor (copied from docking-uxs1)
- `run_screen_pass1.py` — pass-1 driver (exhaustiveness 4)
- `run_redock_top10.py` — pass-2 driver (top 10, exhaustiveness 8)
- `outputs/pass1/<chembl_id>_out.pdbqt` — docked poses, pass 1
- `outputs/pass1/<chembl_id>_energies.txt` — per-pose energies, pass 1
- `outputs/pass1/results_pass1.csv` — ranked pass-1 results
- `outputs/top10/<chembl_id>_out.pdbqt` — docked poses, pass 2
- `outputs/top10/<chembl_id>_energies.txt` — per-pose energies, pass 2
- `outputs/top10/results_top10.csv` — ranked pass-2 results
- `logs/prep.log` — ligand preparation log
- `logs/prep_failures.txt` — compounds that failed preparation (none)
- `logs/pass1.log` — pass-1 run log
- `logs/pass2.log` — pass-2 run log

## Environment

Linux VM, Python 3.12, AutoDock Vina 1.2.7, Meeko 0.8.0, RDKit
2026.09.2 (venv at /tmp/docking-venv/venv; ephemeral, the scripts and
configs above are enough to reproduce the run). Run Nov 16, 2026.

Runtime: ligand preparation ~2 min; pass 1 ~42 min (100 ligands at
exhaustiveness 4); pass 2 ~15 min (top 10 at exhaustiveness 8);
about 1 hour total. A runtime restart interrupted pass 1 at 95/100
and wiped /tmp, including the venv; the venv was recreated with
identical package versions and pass 1 resumed from its checkpointed
outputs (the driver skips ligands with existing outputs) and
completed all 100. No results were lost.
