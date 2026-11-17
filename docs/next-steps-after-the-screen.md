# Where to take UXS1 next

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


Written 2026-11-17, after the full approved drug screen failed its
positive controls. This is an argument about what to do next and why,
not a task list.

## The thing to get straight first

The genetics and the chemistry are being asked to prove different
things, and only one of them is in trouble.

**The target hypothesis is in good shape.** UXS1 is a selective
dependency in KEAP1-mutant lung cancer, there is a solved structure, a
graded pocket, and an independent 2025 Cancer Research paper reaching
the same conclusion by a different route. None of that depended on the
screen.

**The chemistry is in trouble.** No compound is known to bind UXS1,
docking cannot rank its own substrate correctly, and there is no
cysteine near the site for a covalent approach.

The mistake would be to let the second drag down the first.

## The insight the biology is handing us

The published result is that **losing** UXS1 sensitizes KEAP1-mutant
tumors to checkpoint inhibitors. That came from knockout and knockdown,
which remove the protein.

An active site inhibitor does something different. It blocks catalysis
while the protein stays in the cell, still folded, still able to hold
whatever it holds and sit where it sits. For many enzymes that is the
same thing. For some it is not.

So the modality that matches the validated biology most closely is the
one that removes the protein:

1. **Degrader (PROTAC or molecular glue).** Removes the protein, which
   is what the genetics tested. It also sidesteps the hard part: a
   degrader needs a ligand that binds somewhere on the surface with
   modest affinity, not a potent active site inhibitor competing against
   a tight natural substrate.
2. **ASO or siRNA.** Removes the protein most directly of all. Lung is a
   plausible delivery target. This is how the genetic result would be
   reproduced pharmacologically with the least invention.
3. **Active site inhibitor.** The hardest of the three, competing with a
   charged substrate in a polar pocket, and the furthest from what was
   actually demonstrated.

The screen spent a day pursuing option 3, the hardest one, in the
chemical space least likely to contain it.

## Priorities

### 1. Run the genetic experiment. It needs no compound.

The validation plan already names the lines and the assays. CRISPR
knockout of UXS1 in STK11/KEAP1-mutant versus wild type lines, read out
on viability, then the checkpoint sensitization arm.

This is the experiment that makes or breaks the project's central claim,
it is a replication plus extension of a published 2025 result so the
expected outcome is known, and it is ordinary cell biology rather than
anything exotic. Everything else is downstream of it.

**Do this before any more chemistry.**

### 2. Publish the negative control, do not bury it

A screen that fails its own positive control is a finding. The
percentile table in
[screen-full-findings.md](screen-full-findings.md) is the kind of thing
that makes a reviewer trust the rest of the work, because the obvious
sceptical question has already been asked and answered.

The framing is: we built the screen, then built the test that showed the
screen could not work on this site, and we are reporting that rather
than a ranked list we do not believe.

### 3. If the chemistry continues, change the chemical space

Approved drugs were the wrong library for a nucleotide sugar site. Three
better options, in order of how well they fit the site:

- **Substrate analogues.** UDP-sugar chemistry is established. Analogues
  of UDP-glucuronate are the compounds most likely to bind, and they
  come with the obvious liabilities, poor cell permeability and poor
  selectivity against the many other UDP-sugar enzymes. Those liabilities
  are real but they are at least the normal problems of a normal
  programme.
- **NAD site fragments.** UXS1 is a short chain dehydrogenase with a
  Rossmann fold NAD site. That site is better precedented for small
  molecule binders than the sugar nucleotide site, and the holo
  structure already has NAD modelled.
- **Fragment screening, computationally or experimentally.** Fragments
  are the standard answer for a site with no known chemistry. They need
  a biophysical assay rather than docking, which brings us to the next
  point.

### 4. Do not try to fix the docking. It has been tried.

The rule was fixed in advance: a protocol is only used to rank compounds
once the native ligands sit in the top few percent of a size matched
decoy set. Three repairs have now been tested against that bar.

| Condition | Substrate percentile |
| --- | --- |
| Vina, original box, neutral ligands | 31.9 |
| Vinardo, corrected box, neutral | 44.2 |
| Vina, corrected box, neutral | 58.3 |
| Vina, corrected box, ionized | 67.3 |
| Vinardo, corrected box, ionized | 78.9 |

Every repair helped, and the best combination still leaves the enzyme's
own substrate losing to 42 of 199 random approved drugs. The five
natural ligands also reorder almost arbitrarily between conditions,
which says the remaining signal is not reliable even where the number
looks better.

What is left to try is not cheap: explicit solvent free energy methods,
or a physics based rescoring that models desolvation of a trianion
properly. Those are weeks of work and a different skill set, and they
would be spent making a prediction that a single binding assay would
settle in an afternoon. That is the wrong trade for this project.

So the decision is to stop here rather than keep tuning. If a compound
ever needs to be prioritized against this site, the route is an assay,
not a better scoring function.

### 5. What not to do

- **Do not chase the sulfonylureas.** They did not survive the corrected
  box: glyburide lost 1.945 kcal/mol and glipizide 1.235, and both became
  unstable across random seeds. Even had they survived, they have the
  physical profile of promiscuous binders, and glyburide and glipizide
  cause hypoglycaemia through KATP blockade, which rules them out as tool
  compounds and as clinical candidates alike.
- **Do not screen a bigger library with the same protocol.** Going from
  100 to 2,757 compounds moved the hit rate from 7.7 to 8.3 percent and
  changed no conclusion. Going to a million would produce a longer list
  with the same information content, which is close to none.
- **Do not present docking scores as evidence of binding** in any
  application or writeup. They are a prediction from a function that
  does not rank this site's own ligands correctly.

## The honest one line summary

The target survived the week and the chemistry did not. That is a normal
outcome in early discovery, and the project is better positioned than it
was, because it now knows which of its claims rest on evidence and which
rested on a scoring function.
