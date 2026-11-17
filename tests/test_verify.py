"""Tests for the verification layer.

The headline test is test_catches_the_real_defect, which runs the
agents against the ranking as it stood before 2026-11-17 and asserts
that they flag it. The ranking after the fix is the negative control:
a checker that fires on everything is not a checker.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify import (Context, Finding, RuleBackend, run_decomposed,
                    run_single_pass)
from verify.agents import CrossCheckAgent, LayerBalanceAgent
from verify.context import LAYERS
from verify.tools import answer_key_position, layer_depletion

CARDS = "results/reason_cards"
SHORTLIST = "results/tables/shortlist_dependency_supported_phase3.csv"
BEFORE = "results/tables/target_scores_phase3.csv"
AFTER = "results/tables/target_scores_phase3_corrected.csv"


def fake_context(tmp_path, rows):
    df = pd.DataFrame(rows)
    for col in LAYERS:
        if col not in df.columns:
            df[col] = np.nan
    if "rank" not in df.columns:
        df.insert(0, "rank", np.arange(1, len(df) + 1))
    (tmp_path / "r.csv").write_text(df.to_csv(index=False))
    return Context(tmp_path, "r.csv")


# -------------------------------------------------------------- Finding
def test_finding_rejects_an_unknown_severity():
    with pytest.raises(ValueError):
        Finding(agent="a", severity="catastrophic", claim="c",
                evidence="e", where="w")


# -------------------------------------------------------- the invariant
def test_depletion_fires_when_the_top_lacks_a_layer(tmp_path):
    """Five genes at the top with no dependency, forty with it."""
    rows = [{"gene": f"TOP{i}", "tumor_specificity": 0.9} for i in range(5)]
    rows += [{"gene": f"REST{i}", "tumor_specificity": 0.5,
              "dependency": 0.5} for i in range(40)]
    ctx = fake_context(tmp_path, rows)
    got = {r["layer"]: r for r in layer_depletion(ctx, k=5)}
    assert got["dependency"]["present_in_top"] == 0
    assert got["dependency"]["p_depleted"] < 0.01


def test_depletion_stays_quiet_when_the_top_is_well_evidenced(tmp_path):
    rows = [{"gene": f"G{i}", "tumor_specificity": 0.5, "dependency": 0.5}
            for i in range(45)]
    ctx = fake_context(tmp_path, rows)
    for record in layer_depletion(ctx, k=5):
        assert record["p_depleted"] > 0.01


def test_depletion_does_not_fire_on_a_layer_nobody_has(tmp_path):
    """A layer absent from the whole list is not evidence of a defect."""
    rows = [{"gene": f"G{i}", "tumor_specificity": 0.5} for i in range(45)]
    ctx = fake_context(tmp_path, rows)
    got = {r["layer"]: r for r in layer_depletion(ctx, k=5)}
    assert got["dependency"]["population_rate"] == 0.0
    assert got["dependency"]["p_depleted"] == 1.0


# ------------------------------------------------------- answer key math
def test_answer_key_percentile(tmp_path):
    rows = [{"gene": "KRAS" if i == 0 else f"G{i}",
             "tumor_specificity": 0.5} for i in range(100)]
    ctx = fake_context(tmp_path, rows)
    got = {r["gene"]: r for r in answer_key_position(ctx)}
    assert got["KRAS"]["rank"] == 1
    assert got["KRAS"]["percentile"] == pytest.approx(100.0)
    assert got["EGFR"]["rank"] is None


# ------------------------------------------------------- the cross check
def test_cross_check_needs_two_independent_failures(tmp_path):
    ctx = fake_context(tmp_path, [{"gene": "A", "tumor_specificity": 0.5}])
    agent = CrossCheckAgent(RuleBackend())
    only_one = [Finding(agent="layer-balance", severity="critical",
                        claim="x", evidence="carry dependency, y", where="w")]
    assert agent.run(ctx, only_one) == []
    both = only_one + [Finding(agent="controls", severity="critical",
                               claim="KRAS is low", evidence="e", where="w")]
    assert len(agent.run(ctx, both)) == 1


# ------------------------------------------ the real defect, end to end
@pytest.mark.skipif((ROOT / BEFORE).exists() == False,
                    reason="result tables are not present")
def test_catches_the_real_defect():
    ctx = Context(ROOT, BEFORE, CARDS, SHORTLIST)
    result = run_decomposed(ctx, RuleBackend())
    claims = " ".join(f.claim for f in result["findings"])
    assert "dependency" in claims
    assert any(f.agent == "cross-check" for f in result["findings"])
    critical = [f for f in result["findings"] if f.severity == "critical"]
    assert len(critical) >= 4


@pytest.mark.skipif((ROOT / AFTER).exists() == False,
                    reason="result tables are not present")
def test_quiet_on_the_fixed_ranking():
    """The negative control.

    One finding survives, CD274 sitting at the 72nd percentile. That is
    a false positive and it is left in deliberately: the rule inherited
    the project's original wording, which expected every approved target
    near the top, and CD274 is a response biomarker rather than a gene
    this subgroup depends on. Silencing it by naming CD274 would be
    tuning the checker to the answer.
    """
    ctx = Context(ROOT, AFTER, CARDS, SHORTLIST)
    result = run_decomposed(ctx, RuleBackend())
    claims = [f.claim for f in result["findings"]]
    assert any("CD274" in c for c in claims)
    assert not any("dependency layer far more often" in c for c in claims)
    assert not any(f.agent == "cross-check" for f in result["findings"])


@pytest.mark.skipif((ROOT / BEFORE).exists() == False,
                    reason="result tables are not present")
def test_decomposition_adds_only_the_cross_check():
    ctx = Context(ROOT, BEFORE, CARDS, SHORTLIST)
    backend = RuleBackend()
    dec = {f.claim for f in run_decomposed(ctx, backend)["findings"]}
    sing = {f.claim for f in run_single_pass(ctx, backend)["findings"]}
    assert sing.issubset(dec)
    extra = dec - sing
    assert len(extra) == 1


# ------------------------------------------------- the three shapes agree
@pytest.mark.skipif((ROOT / BEFORE).exists() == False,
                    reason="result tables are not present")
def test_all_three_shapes_find_the_same_problems():
    """Cutting the work differently must not change what is found.

    If it did, the speedup comparison would be meaningless, because the
    shapes would be doing different work.
    """
    from verify.orchestrator import run_flat_parallel
    ctx = Context(ROOT, BEFORE, CARDS, SHORTLIST)
    backend = RuleBackend()
    by_agent = {f.claim for f in run_decomposed(ctx, backend)["findings"]}
    by_judgment = {f.claim for f in run_flat_parallel(ctx, backend)["findings"]}
    single = {f.claim for f in run_single_pass(ctx, backend)["findings"]}
    assert by_agent == by_judgment
    assert single.issubset(by_agent)


def test_collect_and_judge_split_is_faithful(tmp_path):
    """run() must equal collect() then judge(), or the pools diverge."""
    rows = [{"gene": f"TOP{i}", "tumor_specificity": 0.9} for i in range(5)]
    rows += [{"gene": f"REST{i}", "tumor_specificity": 0.5,
              "dependency": 0.5} for i in range(40)]
    ctx = fake_context(tmp_path, rows)
    agent = LayerBalanceAgent(RuleBackend(), k=5)
    direct = [f.claim for f in agent.run(ctx)]
    manual = []
    for kind, record in agent.collect(ctx):
        finding = agent.judge(ctx, kind, record)
        if finding is not None:
            manual.append(finding.claim)
    assert direct == manual
