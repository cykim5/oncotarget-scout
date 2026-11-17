# Full approved drug screen against UXS1: what it found, and what it did not

> **Correction, 2026-11-17.** The conclusion below rests on docking into
> PDB 4GLL, and 4GLL was solved with NAD bound but no UDP. With the
> substrate site empty, ARG272 relaxed into it: 14 of the 25 atoms of the
> crystallographic UDP position overlap that residue, the worst by 0.57 A,
> and 1.42 A of the overlap is backbone, which no flexible side chain
> setting can move. The superposition itself is sound, checked by the
> UDP to NAD geometry, 16.39 A in 2B69 against 17.33 A after transforming.
> So every ligand, including the natural ligands used as the positive
> control, was docked against a partly blocked pocket, and the corrected
> box was centred on space the receptor occupies.
>
> This does not overturn the finding, and it does mean the finding is not
> yet established. The re-test uses PDB 2B69, human UXS1 at 1.21 A with
> NAD and UDP both bound, so the site is already in the conformation that
> holds the ligand and nothing has to be superposed. Until that run is
> in, treat the conclusion here as pending rather than settled.


Run 2026-11-16 to 2026-11-17. 2,757 approved small molecules docked
against the UXS1 substrate site, scaling up the 100 compound pilot.

**The headline is not a hit list.** The screen ran as designed, and then
a control test showed the protocol cannot tell the enzyme's own
substrate from random drugs. That makes the ranked output unusable as a
source of candidates, and it is the most useful thing the exercise
produced.

## What was run

| Stage | Scope | Result |
| --- | --- | --- |
| Library | ChEMBL approved small molecules, 150 to 600 Da | 2,757 rows |
| Preparation | RDKit 3D conformer, Meeko | 2,733 prepared, 24 failed |
| Pass 1 | Vina, exhaustiveness 4, holo receptor | 2,738 scored in 6.7 h |
| Shortlist | top 40 after size matching and deduplication | re-docked at exhaustiveness 8 |
| Controls | 5 natural ligands in the same box | see below |

The 19 compounds that did not score are mostly boron containing drugs,
which AutoDock Vina has no atom type for, plus a few RDKit could not
generate a conformer for. Both are listed in `logs/`.

## Three problems found by checking rather than by reporting

### 1. Salt forms were inflating the counts

ChEMBL lists sodium, sulfate, besylate and similar forms as separate
molecules. Preparation strips them to the same parent, so the same
structure was docked repeatedly under different names. 2,757 library
rows collapse to **2,065 distinct structures**. Clopidogrel appears five
times.

After deduplication, 41 of 495 size matched structures beat the natural
product, **8.3 percent**, against 7.7 percent in the 100 compound pilot.
The pilot was representative, and scaling up changed the hit rate
hardly at all.

### 2. The search box was in the wrong place

The box was centred at (60.0, 23.2, 39.0). The crystallographic UDP in
4GLL sits at (62.3, 16.4, 34.2), **8.6 A away**. The crystal ligand fits
inside the box, but comes within 2.1 A of a wall.

The consequence shows up in the poses. Of the 40 re-docked compounds,
**20 come within 2 A of a box wall and 11 within 1 A**. A pose pressed
against a wall is a truncated search: the molecule would extend further
if it could, so its score reflects the box as much as the site.

A corrected run, centred on the crystal UDP position and widened to
28 A, is in `run_recentred.py`, with three random seeds per compound
because the first re-dock reused seed 42 from pass 1 and so could not
test pose stability independently.

### 3. The protocol fails its own positive controls

This is the finding that matters. The five natural ligands were docked
in the same box as the library, so they can be ranked against it
directly. A protocol with real discriminating power puts the native
ligands near the top of a random drug library.

| Natural ligand | Score | Size matched drugs | Beating it | Percentile |
| --- | --- | --- | --- | --- |
| UDP-xylose (product) | -8.796 | 649 | 44 | 93.2 |
| UDP | -7.480 | 1,152 | 263 | 77.2 |
| UDP-galactose | -7.391 | 523 | 324 | 38.0 |
| **UDP-glucuronate (substrate)** | **-7.337** | **439** | **299** | **31.9** |
| UMP | -6.997 | 1,128 | 304 | 73.0 |

**299 of 439 size matched approved drugs outscore the enzyme's actual
substrate.** On its own target, with its own natural ligand, the
protocol ranks the right answer in the bottom third.

So a drug scoring better than UDP-glucuronate has demonstrated nothing.
Roughly two thirds of random drugs do.

### Why this is not surprising in hindsight

UDP-glucuronate carries a diphosphate and a carboxylate and is strongly
anionic at physiological pH. Empirical scoring functions, Vina included,
are known to handle charged phosphate bearing ligands poorly, because
they do not model the cost of stripping water from a charged group well.
The site is built to pay that cost; the scoring function is not built to
account for it.

Meanwhile approved oral drugs are overwhelmingly neutral and lipophilic.
They are rewarded for burying hydrophobic surface, which this scoring
function counts generously. The screen was always going to favour them.

## Could the protocol be fixed? Three repairs tested, none passes

Before concluding that docking cannot work on this site, the three
cheapest explanations for the failed control were tested one at a time,
against the same size matched decoy set each time, so the percentiles
compare directly.

**The box.** Re-centring on the crystallographic ligand moved the
substrate from the 31.9th to the 58.3rd percentile. A real improvement,
and the largest single one, which is another way of saying the original
box was the biggest single error in the screen.

**The scoring function.** Vinardo weights polar contacts differently
from Vina, which is the obvious thing to try on a charged ligand. In the
corrected box it put the substrate at the 44.2nd percentile, worse than
Vina. 200 decoys were re-docked with each function to get this.

**The protonation state.** Every prepared natural ligand summed to a net
partial charge of about zero, so the phosphates and the carboxylate had
been modelled protonated. At physiological pH UDP-glucuronate carries
about -3. Deprotonating the acidic groups in place, keeping the docked
geometry, and re-running both functions gave the best result of the
whole exercise: 67.3rd percentile with Vina, 78.9th with Vinardo.

| Condition | UDP-glucuronate score | Percentile |
| --- | --- | --- |
| Vina, original box, neutral | -7.337 | 31.9 |
| Vinardo, corrected box, neutral | -5.258 | 44.2 |
| Vina, corrected box, neutral | -7.938 | 58.3 |
| Vina, corrected box, ionized | -8.138 | 67.3 |
| **Vinardo, corrected box, ionized** | **-6.008** | **78.9** |

So the repairs worked, in the sense that each one moved the number in a
direction that makes physical sense, and together they roughly halved
the error. The substrate still loses to 42 of 199 random approved drugs.
The bar set in advance was the top few percent. 78.9 is not that.

Two details make it worse rather than better. The order of the five
natural ligands reshuffles almost completely between conditions. UDP
goes from the 13.6th percentile to the 58.3rd on protonation state
alone, and UDP-galactose moves the wrong way in both functions. A
protocol whose ranking of five closely related molecules is that
unstable is not measuring the thing it is supposed to measure. And
UDP-xylose, the product, which looked like the one confident result at
93.2 in the original box, sits between the 37th and 54th percentile once
the box is right. That number was the misplaced box as well.

The conclusion is the one worth writing down once: this site is not
amenable to empirical scoring function docking at the level of
discrimination needed to pick compounds. Not because the software was
used badly, since the three obvious errors were found and corrected, but
because a deeply buried, strongly anionic, highly polar substrate site
is the case these functions are known to be worst at, and the correction
needed is larger than anything a box or a weighting change provides.

## The observation that did not survive

The sulfonylureas looked like the one thing worth keeping. Glyburide,
glipizide and glimepiride, three structurally related antidiabetics,
ranked together near the top, with poses reaching hydrogen bonding
distance of two of the three catalytic residues. A chemical series
recurring on its own is much harder to dismiss than a single compound.

Then the corrected box was run, with three random seeds per compound.

| Compound | Original box | Corrected box, mean of 3 seeds | Change | Seed spread |
| --- | --- | --- | --- | --- |
| Glyburide | -10.332 | -8.387 | **+1.945** | 0.916, unstable |
| Glipizide | -9.918 | -8.683 | **+1.235** | 0.504, unstable |
| Glimepiride | -8.981 | -8.918 | +0.063 | 0.222, stable |

Glyburide lost nearly two kilocalories and became the least stable
compound in the set. Both it and glipizide fall out of the top ten
entirely. The series was an artefact of a box centred 8.6 A from the
real site, and the apparent catalytic residue contacts went with it.

This is why the box check mattered. Without it, the writeup would have
led with a sulfonylurea series and a structure activity story that is
not there.

Across all 40 compounds, **14 were unstable across seeds**, with a
median spread of 0.36 and a maximum of 1.94 kcal/mol. Lapatinib, which
the 100 compound pilot singled out, has a spread of 1.00 and is among
the unstable ones.

What survives the corrected box and the seed test, with the enrichment
caveat below still applying to all of it: eravacycline, pralsetinib,
regorafenib and sorafenib, all stable within 0.5 kcal/mol, none of them
a story, all of them large kinase inhibitor shaped molecules of the kind
this scoring function reliably over rewards.

## A structural route that is now closed

A covalent inhibitor does not need the binding affinity a reversible one
does, so it is a standard fallback for a difficult site. It needs a
reactive residue lining the pocket. The nearest cysteine to the
crystallographic UDP position is **Cys230 at 12.7 A**, far outside the
site, and the next is at 21.6 A. There is no covalent handle at the
substrate site in this structure.

## What this changes about the project

Nothing about the target. UXS1 remains the most tractable of the top
five by structure, with a graded pocket, a solved complex, and
independent 2025 evidence that it is a selective dependency in
KEAP1-mutant lung cancer. What changed is the claim about chemistry:

- **Before:** a 100 compound screen produced a ranked list with
  lapatinib as the size matched leader.
- **After:** scaling to 2,757 moved lapatinib to rank 32 of 495 size
  matched structures, the corrected box showed its score varies by a
  full kilocalorie between random seeds, and the control test showed the
  ranking carries little information either way.

Reporting the negative control is a stronger result than reporting a top
ten, because anyone competent will ask for that control, and the answer
is already here.
