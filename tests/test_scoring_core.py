"""Tests for the scoring rules.

The most important test here is test_missing_layers_do_not_win. It is
the regression test for the bug that put non-coding RNAs with two
layers above every real candidate with five.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from scoring_core import (LAYERS, WEIGHTS, coverage, dependency_covered_ranking,
                          passes_essentiality_gate, scale01, weighted_score)


def frame(rows):
    """Build a score table from dicts, filling absent layers with NaN."""
    df = pd.DataFrame(rows)
    for col in LAYERS:
        if col not in df.columns:
            df[col] = np.nan
    return df


# ------------------------------------------------------------- scale01
def test_scale01_maps_to_zero_and_one():
    out = scale01(pd.Series([2.0, 4.0, 6.0]))
    assert out.iloc[0] == 0.0
    assert out.iloc[2] == 1.0
    assert out.iloc[1] == 0.5


def test_scale01_all_equal_gives_nan():
    out = scale01(pd.Series([3.0, 3.0, 3.0]))
    assert out.isna().all() == True


def test_scale01_keeps_nan_as_nan():
    out = scale01(pd.Series([1.0, np.nan, 5.0]))
    assert pd.isna(out.iloc[1]) == True
    assert out.iloc[0] == 0.0
    assert out.iloc[2] == 1.0


# ------------------------------------------------------- weighted_score
def test_weighted_score_matches_hand_calculation():
    df = frame([{"gene": "A", "tumor_specificity": 1.0, "dependency": 0.0,
                 "safety_window": 1.0, "malignant_share": 0.0,
                 "outcome_link": 1.0}])
    got = weighted_score(df, WEIGHTS).iloc[0]
    want = 0.25 * 1.0 + 0.30 * 0.0 + 0.20 * 1.0 + 0.15 * 0.0 + 0.10 * 1.0
    assert got == pytest.approx(want)


def test_weighted_score_never_imputes_missing_as_zero():
    """A gene missing a layer must not be scored as if it had a 0 there."""
    full = frame([{"gene": "A", "tumor_specificity": 0.8, "dependency": 0.8,
                   "safety_window": 0.8, "malignant_share": 0.8,
                   "outcome_link": 0.8}])
    partial = frame([{"gene": "B", "tumor_specificity": 0.8,
                      "dependency": 0.8, "safety_window": 0.8}])
    assert weighted_score(full, WEIGHTS).iloc[0] == pytest.approx(0.8)
    assert weighted_score(partial, WEIGHTS).iloc[0] == pytest.approx(0.8)


def test_weighted_score_is_nan_when_no_layer_is_present():
    df = frame([{"gene": "A"}])
    assert pd.isna(weighted_score(df, WEIGHTS).iloc[0]) == True


# ------------------------------------------------------------ coverage
def test_coverage_counts_present_layers():
    df = frame([{"gene": "A", "tumor_specificity": 0.5, "dependency": 0.5},
                {"gene": "B", "tumor_specificity": 0.5, "dependency": 0.5,
                 "safety_window": 0.5, "malignant_share": 0.5,
                 "outcome_link": 0.5}])
    got = coverage(df)
    assert got.iloc[0] == 2
    assert got.iloc[1] == 5


# ---------------------------------------------------- essentiality gate
def test_gate_passes_strongly_essential_genes():
    s = pd.Series([-0.9, -0.5, -0.1])
    got = passes_essentiality_gate(s)
    assert got.iloc[0] == True
    assert got.iloc[1] == False   # the gate is strictly less than
    assert got.iloc[2] == False


def test_gate_treats_missing_coverage_as_not_a_pass():
    got = passes_essentiality_gate(pd.Series([np.nan]))
    assert got.iloc[0] == False


# ------------------------------------------- the ranking bug, regressed
def test_missing_layers_do_not_win():
    """A two-layer gene must not outrank a five-layer gene.

    This is the Phase 3 bug. Under renormalization alone, LINC00676
    scored 0.84 on expression, safety and outcome, and beat every
    candidate that had all five layers. The ranked output must now put
    the gene with dependency evidence first.
    """
    df = frame([
        {"gene": "LNC_LIKE", "tumor_specificity": 0.95,
         "safety_window": 0.95, "outcome_link": 0.95},
        {"gene": "REAL_TARGET", "tumor_specificity": 0.65,
         "dependency": 0.70, "safety_window": 0.60,
         "malignant_share": 0.60, "outcome_link": 0.55},
    ])
    df["score"] = weighted_score(df, WEIGHTS)
    assert df.loc[0, "score"] > df.loc[1, "score"]   # the bug is real

    ranked, uncovered = dependency_covered_ranking(df)
    assert ranked.loc[0, "gene"] == "REAL_TARGET"
    assert list(uncovered["gene"]) == ["LNC_LIKE"]


def test_ranking_keeps_every_gene_somewhere():
    df = frame([{"gene": "A", "dependency": 0.5, "score": 0.5},
                {"gene": "B", "score": 0.9},
                {"gene": "C", "dependency": 0.1, "score": 0.1}])
    ranked, uncovered = dependency_covered_ranking(df)
    assert len(ranked) + len(uncovered) == 3
    assert list(ranked["rank"]) == [1, 2]
    assert list(uncovered["rank"]) == [1]


def test_ranking_orders_by_score_within_each_group():
    df = frame([{"gene": "A", "dependency": 0.5, "score": 0.2},
                {"gene": "B", "dependency": 0.5, "score": 0.8}])
    ranked, _ = dependency_covered_ranking(df)
    assert list(ranked["gene"]) == ["B", "A"]
