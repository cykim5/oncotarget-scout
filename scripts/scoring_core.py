"""The scoring rules, in one place so they can be tested.

Phase 3 scored every gene by renormalizing the weights over whichever
layers that gene actually had. That is the right way to combine layers
for a single gene, and the wrong way to rank genes against each other,
because a gene with one layer can reach a score of 1.0 while a gene
with all five almost never can. Missing data was being rewarded.

This module holds the layer maths and the ranking rule that fixes it,
so that the pipeline scripts and the tests use exactly the same code.
"""
import numpy as np
import pandas as pd

LAYERS = ["tumor_specificity", "dependency", "safety_window",
          "malignant_share", "outcome_link"]

WEIGHTS = {"dependency": 0.30, "tumor_specificity": 0.25,
           "safety_window": 0.20, "malignant_share": 0.15,
           "outcome_link": 0.10}

DEP_ESSENTIALITY_GATE = -0.5   # mean Chronos in mutant lung lines


def scale01(s):
    """Scale a series to 0 to 1. All-equal or empty input gives NaN."""
    lo, hi = s.min(), s.max()
    if pd.isna(lo) == True or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


def weighted_score(df, weights):
    """Weighted mean over the layers a gene actually has.

    Missing layers are never imputed as zero. The weights are
    renormalized over the layers present, which is correct for
    describing one gene and not sufficient for ranking genes against
    each other. Use dependency_covered_ranking for that.
    """
    num = sum(weights[k] * df[k].fillna(0) for k in weights)
    den = sum(weights[k] * df[k].notna() for k in weights)
    return num / den.replace(0, np.nan)


def coverage(df):
    """How many of the five layers support this gene's score."""
    present = [c for c in LAYERS if c in df.columns]
    return df[present].notna().sum(axis=1)


def passes_essentiality_gate(mean_chronos, gate=DEP_ESSENTIALITY_GATE):
    """True where the gene is essential enough in mutant lung lines.

    A gene with no DepMap coverage is not a pass and not a fail. It is
    unknown, so it returns False here and is handled by the coverage
    rule instead.
    """
    return mean_chronos < gate


def dependency_covered_ranking(scores):
    """Rank only the genes that have the dependency layer.

    The question this project asks is which gene the subgroup depends
    on. Dependency carries the highest weight, and a gene with no
    DepMap coverage has no evidence on the layer that matters most, so
    it cannot be compared against one that does. Genes without the
    layer are not deleted, they are ranked separately and reported
    with their coverage.

    Returns (ranked, uncovered), both sorted by score, with a rank
    column starting at 1.
    """
    covered = scores[scores["dependency"].notna() == True].copy()
    uncovered = scores[scores["dependency"].isna() == True].copy()
    out = []
    for part in (covered, uncovered):
        part = part.sort_values("score", ascending=False).reset_index(drop=True)
        part["rank"] = np.arange(1, len(part) + 1)
        cols = ["rank"] + [c for c in part.columns if c != "rank"]
        out.append(part[cols])
    return out[0], out[1]
