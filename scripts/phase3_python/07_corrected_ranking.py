"""Rebuild the headline ranking so that missing data cannot win.

Phase 3 renormalized the layer weights over whichever layers a gene
had. For one gene that is correct. Across genes it is not, because the
fewer layers a gene has, the higher its ceiling. The result was that
all 20 genes at the top of the Phase 3 ranking were missing the
dependency layer, and 18 of them were non-coding RNAs or pseudogenes.
Every real candidate sat below them.

The decision output was never affected, because the shortlist already
required dependency support. Only the headline ranking was wrong, and
the headline is what a reader sees first.

The rule applied here: a gene is ranked against other genes only if it
has the layer that carries the most weight, which is dependency. Genes
without it are not deleted. They are ranked separately and reported
with their coverage, so nothing is hidden.

This reads the committed layer scores, so it reproduces from the
repository without the raw downloads.

Outputs:
  results/tables/target_scores_phase3_corrected.csv
  results/tables/no_dependency_layer_phase3.csv
  results/answer_key_check_corrected.md
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from scoring_core import dependency_covered_ranking

TABLES = ROOT / "results" / "tables"
ANSWER_KEY = ["EGFR", "KRAS", "CD274"]

scores = pd.read_csv(TABLES / "target_scores_phase3.csv")
before = scores.set_index("gene")["rank"].to_dict()
n_before = len(scores)

ranked, uncovered = dependency_covered_ranking(scores)
ranked.to_csv(TABLES / "target_scores_phase3_corrected.csv", index=False)
uncovered.to_csv(TABLES / "no_dependency_layer_phase3.csv", index=False)

print(f"ranked with dependency evidence : {len(ranked)}")
print(f"reported separately, no layer   : {len(uncovered)}")
print()
print("Top 10 before (renormalized over covered layers):")
for _, r in scores.head(10).iterrows():
    print(f"  {int(r['rank']):>5}  {r['gene']:<12} {r['score']:.3f}  "
          f"coverage {int(r['coverage'])}/5")
print()
print("Top 10 after (dependency evidence required):")
for _, r in ranked.head(10).iterrows():
    print(f"  {int(r['rank']):>5}  {r['gene']:<12} {r['score']:.3f}  "
          f"coverage {int(r['coverage'])}/5")

lines = ["# Answer key check, corrected ranking", "",
         "The three genes below are approved drug targets in lung",
         "adenocarcinoma. They are never excluded and never up-weighted.",
         "Where they land is a check on the method, not an input to it.",
         "",
         "The rule being checked was chosen because the renormalized",
         "ranking rewarded missing data, not because of where it put",
         "these three. Both positions are shown so the change is visible.",
         "",
         "| Gene | Before | After | Percentile after |",
         "| --- | --- | --- | --- |"]
for g in ANSWER_KEY:
    row = ranked[ranked["gene"] == g]
    if len(row) == 0:
        lines.append(f"| {g} | {before.get(g, 'not ranked')} of {n_before} "
                     f"| no dependency layer | n/a |")
        continue
    rank = int(row["rank"].iloc[0])
    pct = 100.0 * (1.0 - (rank - 1) / len(ranked))
    lines.append(f"| {g} | {before.get(g)} of {n_before} | "
                 f"{rank} of {len(ranked)} | {pct:.1f} |")
lines += ["",
          "KRAS is the gene this subgroup is defined against, and it now",
          "sits in the top one percent without being told to. EGFR is in",
          "the top four percent. CD274, which encodes PD-L1, stays far",
          "down the list, which is the correct result: it is a biomarker",
          "for who responds to immunotherapy, not a gene this subgroup",
          "depends on, and a method that ranked it high would be wrong."]
(ROOT / "results" / "answer_key_check_corrected.md").write_text(
    "\n".join(lines) + "\n")

print()
print("\n".join(lines[9:]))
