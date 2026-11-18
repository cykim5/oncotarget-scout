# The receptor was the problem, not the scoring function

Run 2026-11-17. This supersedes the conclusion in
[screen-full-findings.md](screen-full-findings.md) and
[next-steps-after-the-screen.md](next-steps-after-the-screen.md).

## What changed

Everything before this docked into PDB 4GLL, which was solved with NAD
bound and no UDP. With the substrate site empty, ARG272 had relaxed
into it: 14 of the 25 atoms of the crystallographic UDP position
overlap that residue, worst 0.57 A, and 1.42 A of the overlap is
backbone, so no flexible side chain setting could have fixed it.

This run uses PDB 2B69 instead: human UXS1 at 1.21 A with NAD and UDP
both bound. The site is already in the conformation that holds the
ligand, the box is centred on that ligand's own coordinates, and
nothing is superposed. The crystal ligand has zero clashes with the
prepared receptor, worst contact 2.69 A, against 14 clashing atoms in
4GLL.

Everything else is identical: same five natural ligands, same prepared
decoys from the same library, same size window of 31 to 43 heavy
atoms, same exhaustiveness, same seed.

## The result

| Ligand | Heavy atoms | 4GLL percentile | 2B69 percentile | Change |
| --- | --- | --- | --- | --- |
| UDP-xylose | 34 | 30.2 | **99.5** | +69.3 |
| UDP-galactose | 36 | 54.8 | **99.0** | +44.2 |
| UDP-glucuronate, the substrate | 37 | 58.3 | **88.4** | +30.1 |
| UDP | 25 | 13.6 | 84.9 | +71.3 |
| UMP | 21 | 15.1 | 53.3 | +38.2 |

UDP and UMP are outside the decoy size window, so their percentiles
compare molecules of different sizes and do not set the verdict. They
are shown because the movement is informative.

**Verdict: WEAK**, judged on the three size matched ligands. The worst
of them sits at the 88.4th percentile against a bar of 95 set before
the run. Two of the three clear the bar comfortably.

For comparison, the same three in 4GLL were at 30.2, 54.8 and 58.3.
The protocol has gone from unusable to borderline.

## The clearest evidence that the receptor was the cause

Score against heavy atom count, across the same decoys:

| Receptor | Spearman rho | p |
| --- | --- | --- |
| 4GLL, blocked site | **-0.67** | ~1e-30 |
| 2B69, open site | **+0.11** | 0.12 |

In 4GLL the ranking was substantially a measure of molecular size, and
that is exactly what a blocked pocket produces: nothing can form real
contacts, so the only thing separating molecules is how much bulk they
can pile up nearby. Open the pocket and shape complementarity takes
over. The size bias does not shrink, it disappears.

## What is left, and it is the thing we predicted

The substrate is the weakest of the three size matched ligands, and the
pattern tracks charge:

| Ligand | Formal charge at pH 7.4 | Percentile |
| --- | --- | --- |
| UDP-xylose | -2 | 99.5 |
| UDP-galactose | -2 | 99.0 |
| UDP-glucuronate | -3 | 88.4 |

Three points is not a trend, and it is the trend the scoring function's
missing desolvation term predicts: the more charge a ligand carries,
the more the function underestimates it. So the original diagnosis was
right about the mechanism and wrong about the magnitude. Charge costs
this protocol something. It was not what made the screen fail.

## What this means for the project

**The headline claim in the two earlier documents is withdrawn.** The
docking route is not closed. It was never properly tested, because the
pocket was occupied by the receptor's own backbone.

- The full screen of 2,757 approved drugs was run against a blocked
  site and its ranking should not be used, including the parts that
  looked like failures.
- The sulfonylurea result, already overturned once by the box
  correction, needs re-deriving from scratch, not rehabilitating.
- The enrichment failure that drove the whole negative conclusion was
  an artifact of the receptor.
- The screen is worth rerunning against 2B69.

**What survives unchanged:** the target. UXS1's selection came from the
five evidence layers and from independent work (Gebru et al., Cancer
Res 2025, PMID 40966319), neither of which touched docking.

**What also survives:** the discipline. Every one of these errors was
found by running a control rather than by reading the ranked list. The
screen never produced a compound worth believing at any point, and the
control said so every time, for a different reason each time.

## How to reproduce

```bash
cd control-2b69 && ./run_control.sh
```

Needs [dock-check](https://github.com/cykim5/dock-check) on PYTHONPATH,
and the prepared ligand library from `screening-full-uxs1/`.
One decoy, CHEMBL5315122, cannot be docked: it contains boron and
AutoDock has no atom type for it. 199 of 200 scored, 36.7 minutes.
