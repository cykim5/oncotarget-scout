"""Would this have caught the bug, and is splitting the work up worth it?

The defect is real and dated. Until 2026-11-17 the ranking renormalized
the layer weights over whichever layers a gene had, so a gene with two
layers could outrank a gene with five. It survived three phases of my
own review. Both rankings are committed, so the same agents can run
against the broken one and the fixed one and the difference measured.

Three questions, which are not the same question.

1. Detection. Does the layer flag the known defect when it is present,
   and stay quiet when it is gone? A checker that fires on everything
   is not a checker.

2. What decomposition finds. Does splitting the work into specialists
   surface anything one worker running the same checks does not?

3. What decomposition costs, and at what level. Agent level
   parallelism is the intuitive design. Judgment level parallelism is
   the same work cut differently. They do not perform the same.

Writes results/verification-experiment.json
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify import Context, RuleBackend
from verify.agents import specialists
from verify.latency import CountingBackend, SlowBackend
from verify.orchestrator import (run_decomposed, run_flat_parallel,
                                 run_single_pass)

CARDS = "results/reason_cards"
SHORTLIST = "results/tables/shortlist_dependency_supported_phase3.csv"
RANKINGS = [
    ("before the fix", "results/tables/target_scores_phase3.csv"),
    ("after the fix", "results/tables/target_scores_phase3_corrected.csv"),
]
MODES = [run_single_pass, run_decomposed, run_flat_parallel]
REPEATS = 20

report = {"backend": "rules", "repeats": REPEATS, "runs": []}
backend = RuleBackend()

# ----------------------------------------------------------- 1 detection
for label, path in RANKINGS:
    ctx = Context(ROOT, path, CARDS, SHORTLIST)
    print(f"\n=== {label} ===")
    for runner in MODES:
        best = None
        for _ in range(REPEATS):
            result = runner(ctx, backend)
            if best is None or result["seconds"] < best["seconds"]:
                best = result
        sev = {}
        for f in best["findings"]:
            sev[f.severity] = sev.get(f.severity, 0) + 1
        print(f"  {best['mode']:<14} findings {len(best['findings']):>2}  "
              f"critical {sev.get('critical', 0):>2}  "
              f"{best['seconds'] * 1000:>7.1f} ms")
        report["runs"].append({
            "ranking": label, "file": path, "mode": best["mode"],
            "workers": best["workers"], "seconds": best["seconds"],
            "findings": [f.to_dict() for f in best["findings"]]})


def run_of(label, mode):
    return [r for r in report["runs"]
            if r["ranking"] == label and r["mode"] == mode][0]


before = run_of("before the fix", "flat-parallel")
after = run_of("after the fix", "flat-parallel")
report["detection"] = {
    "critical_before": len([f for f in before["findings"]
                            if f["severity"] == "critical"]),
    "critical_after": len([f for f in after["findings"]
                           if f["severity"] == "critical"]),
    "surviving_after": [f["claim"] for f in after["findings"]]}
print(f"\n--- detection ---")
print(f"  critical before the fix: {report['detection']['critical_before']}")
print(f"  critical after the fix : {report['detection']['critical_after']}")
for claim in report["detection"]["surviving_after"]:
    print(f"    still flagged: {claim}")

# ------------------------------------------- 2 what decomposition finds
dec = {f["claim"] for f in before["findings"]}
sing = {f["claim"] for f in run_of("before the fix", "single-pass")["findings"]}
report["only_in_decomposed"] = sorted(dec - sing)
print("\n--- what decomposition added ---")
if len(report["only_in_decomposed"]) == 0:
    print("  nothing")
for claim in report["only_in_decomposed"]:
    print(f"  + {claim}")

# -------------------------------------------------- 3 work balance, cost
ctx = Context(ROOT, RANKINGS[0][1], CARDS, SHORTLIST)
per = {}
for agent in specialists(RuleBackend()):
    counter = CountingBackend(RuleBackend())
    agent.backend = counter
    agent.run(ctx)
    per[agent.name] = counter.calls
total, slowest = sum(per.values()), max(per.values())
print("\n--- work per specialist ---")
for name, n in per.items():
    print(f"  {name:<16} {n} judgments")
print(f"  sequential {total}, busiest specialist {slowest}")
print(f"  so agent level parallelism caps at {total / slowest:.2f}x, "
      f"and judgment level at {total}x")
report["work_balance"] = {"per_agent": per, "sequential": total,
                          "busiest_specialist": slowest,
                          "agent_level_ceiling": total / slowest,
                          "judgment_level_ceiling": float(total)}

print("\n--- cost per check against wall clock ---")
print(f"  {'per check':>10} {'single':>10} {'by agent':>10} "
      f"{'by judgment':>12} {'best':>8}")
sweep = []
for ms in (0, 1, 10, 50, 100, 500):
    b = RuleBackend() if ms == 0 else SlowBackend(RuleBackend(), ms / 1000)
    n = 3 if ms < 100 else 1
    s = min(run_single_pass(ctx, b)["seconds"] for _ in range(n))
    d = min(run_decomposed(ctx, b)["seconds"] for _ in range(n))
    f = min(run_flat_parallel(ctx, b)["seconds"] for _ in range(n))
    sweep.append({"per_check_ms": ms, "single_pass_s": s,
                  "decomposed_s": d, "flat_parallel_s": f,
                  "agent_level_speedup": s / d, "judgment_level_speedup": s / f})
    print(f"  {ms:>8}ms {s * 1000:>9.1f}ms {d * 1000:>9.1f}ms "
          f"{f * 1000:>11.1f}ms {s / f:>7.2f}x")
report["latency_sweep"] = sweep

out = ROOT / "results" / "verification-experiment.json"
out.write_text(json.dumps(report, indent=2) + "\n")
print(f"\nwritten to {out.relative_to(ROOT)}")
