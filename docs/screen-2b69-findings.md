# The approved drug screen against 2B69

Run 2026-11-18. Supersedes [screen-full-findings.md](screen-full-findings.md),
which was run against 4GLL and whose conclusion was wrong for the
reason set out in [control-2b69-result.md](control-2b69-result.md).

2,738 compounds docked, 6 failures, all six boron compounds Meeko
cannot type. 6 hours 3 minutes. 2,058 distinct structures after salt
forms are deduplicated.

## The protocol works now

| Measurement | 4GLL | 2B69 |
| --- | --- | --- |
| Spearman rho, score against heavy atoms, inside the control's 31 to 43 window | -0.67 | **+0.098** |
| Fraction of the library beating the substrate | about two thirds | **4.7 percent** |

The control measured +0.11 against size matched decoys. The full
library gives +0.098 inside the same size window, on an entirely
different set of molecules. Two independent measurements landing
together is the evidence that the receptor was the problem.

Across the whole library rho is -0.494. That is the mechanical atom
count effect, not the 4GLL failure: Vina sums over atom pairs, so raw
score cannot rank a library spanning 14 to 100 heavy atoms. Rankings
are therefore reported three ways, in [analysis.md](../screening-2b69/outputs/analysis.md).

Raw score and ligand efficiency share **no compounds at all** in their
top 25. That is worth keeping as a demonstration of why the
normalisation choice is not cosmetic.

## The poses are in the right place

59 shortlisted structures re-docked with poses kept, against a
substrate site defined as the 18 residues within 4.5 A of the
crystallographic UDP. Defined by distance, not by residue number,
because 2B69 and 4GLL number this protein differently.

1 pose clipped by a box wall, 12 sitting in another groove, 38 clear on
both counts. See [pose_check.md](../screening-2b69/outputs/pose_check.md).
Geometry is no longer the failure mode.

## The chemistry is the answer

Filtering on an anionic group on the **parent** structure plus contact
with the four phosphate gripping arginines leaves 5 of 59:

| Compound | Score | Size matched pct | Anionic group | Arg contacts |
| --- | --- | --- | --- | --- |
| Betamethasone phosphate | -11.15 | 98.2 | phosphate ester | 4/4 |
| Carbenicillin indanyl | -11.10 | 97.9 | carboxylate | 4/4 |
| Alvimopan | -11.09 | 97.7 | carboxylate | 4/4 |
| Fluorescein | -10.93 | 98.8 | carboxylate, phenolate | 4/4 |
| Pemirolast | -9.60 | 98.5 | tetrazole | 4/4 |

Every survivor carries a phosphate, a carboxylate or a tetrazole, and
the scoring function has no electrostatics term. Recovering the site's
own pharmacophore from a library never designed for it, without being
told about charge, is what a working setup looks like. That coherence
is the result.

None of the five is a lead. Betamethasone phosphate and carbenicillin
indanyl are prodrugs whose anionic group exists to be cleaved. Alvimopan
is designed not to cross membranes and UXS1 is a Golgi enzyme.
Fluorescein is a dye and a known promiscuous binder, and notably the
PAINS A, B and C filters did not flag it. Pemirolast is the only one
worth a sentence: oral, not a prodrug, and a tetrazole is a standard
phosphate bioisostere. Its score sits below the substrate's, so it is
not a better binder than the natural ligand.

**Conclusion: the approved drug library contains no credible UXS1
substrate site binder, and that is unsurprising because oral drugs are
rarely polyanions. The finding is that the protocol can now tell the
difference. The 4GLL version demonstrably could not.**

## Two bugs in the plausibility filter, both the same mistake

Recorded because they are general and because any screen run on a
ChEMBL approved drug export has them.

First version passed zinc sulfate, dibasic potassium phosphate, calcium
acetate and hydrotalcite as chemically plausible. An anion is an anion,
and the filter never asked whether the molecule was a drug. They
reached the shortlist because ligand efficiency rewards being small.

Second and subtler: the filter **read counterions as pharmacophores**.
Paroxetine mesylate appeared to carry a sulfonate, which is the
mesylate that dissociates in solution. Amphetamine aspartate appeared
doubly anionic for the same reason. In the other direction, fluorescein
sodium and pemirolast potassium were thrown out as inorganic on account
of counterions, when both are good organic acids.

Judging the parent fragment rather than the formulated salt fixed both,
and took the surviving set from 19 partly nonsensical compounds to 5
coherent ones. Worth adding to dock-check as a check of its own.
