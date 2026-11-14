"""Phase 2, step 5: the dependency-supported shortlist.

The main Phase 2 ranking (script 05) renormalizes over covered layers,
exactly as docs/target-score.md specifies. A side effect, visible in
the results, is that the overall top 20 is dominated by non-coding
genes with no DepMap coverage at all: with the highest-weighted layer
missing, their scores are renormalized over expression, safety, and
outcome alone. The spec anticipates this ("genes are compared within
the same coverage where a decision depends on it"), and the README's
dependency check requires top candidates to be supported by the
DepMap CRISPR layer, not by expression alone.

This script therefore builds the decision shortlist: the top 20 genes
by main score among genes with dependency layer support
(dependency > 0, i.e. selectively essential in STK11/KEAP1-mutant
lung lines and actually essential there). Robustness is recomputed
within this subset under the same three weightings as script 05:
robust means top 20 in the subset under all three, fragile means
top 20 in the subset under exactly one.

It also draws the Kaplan-Meier figure for the top shortlist gene, in
the mutant subgroup survival cohort (same cohort definition as
script 05).

Outputs:
  results/tables/shortlist_dependency_supported.csv
  results/figures/km_shortlist_top.png

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/06_dependency_shortlist.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "xena"
CBIO = ROOT / "data" / "raw" / "cbioportal"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"

print("[1/2] Building the dependency-supported shortlist ...")
scores = pd.read_csv(TABLES / "target_scores_phase2.csv")
sub = scores[scores["dependency"] > 0].copy()
print(f"      {len(sub)} genes have dependency support (dependency > 0)")
for col, name in [("score", "main"), ("score_equal", "equal"),
                  ("score_expr_first", "expr_first")]:
    sub[f"subset_rank_{name}"] = sub[col].rank(ascending=False, method="min").astype(int)
top = {name: set(sub.nsmallest(20, f"subset_rank_{name}")["gene"])
       for name in ("main", "equal", "expr_first")}
robust = top["main"] & top["equal"] & top["expr_first"]
fragile = {g for g in (top["main"] | top["equal"] | top["expr_first"])
           if sum(g in s for s in top.values()) == 1}
sub["shortlist_robust_fragile"] = ""
sub.loc[sub["gene"].isin(robust), "shortlist_robust_fragile"] = "robust"
sub.loc[sub["gene"].isin(fragile), "shortlist_robust_fragile"] = "fragile"
shortlist = sub.sort_values("subset_rank_main").head(20)
shortlist.to_csv(TABLES / "shortlist_dependency_supported.csv", index=False)
print(f"      shortlist robust: {len(robust)}, fragile: {len(fragile)}")
for _, r in shortlist.head(5).iterrows():
    print(f"        #{int(r['subset_rank_main'])} {r['gene']} (overall rank "
          f"{int(r['rank'])}) score {r['score']:.3f} dep {r['dependency']:.3f} "
          f"{r['shortlist_robust_fragile'] or '-'}")

print("[2/2] Kaplan-Meier for the top shortlist gene ...")
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

km_gene = shortlist.iloc[0]["gene"]
clinical = pd.read_csv(RAW / "TCGA-LUAD.clinical.tsv.gz", sep="\t")
survival = pd.read_csv(RAW / "TCGA-LUAD.survival.tsv.gz", sep="\t")
expr = pd.read_csv(RAW / "TCGA-LUAD.star_counts.tsv.gz", sep="\t", index_col=0)
sample_type = expr.columns.str[13:15]
tumor_samples = expr.columns[sample_type == "01"].tolist()
first_by_patient = {}
for s in tumor_samples:
    first_by_patient.setdefault(s.replace(".", "-")[:12], s)

mut_patients = set()
for symbol in ("STK11", "KEAP1"):
    records = json.loads((CBIO / f"mutations_{symbol}.json").read_text())
    mut_patients.update(r["patientId"] for r in records
                        if r.get("mutationType") != "Silent")

clin = clinical.set_index("sample")
stage_map = {"Stage I": 1, "Stage IA": 1, "Stage IB": 1,
             "Stage II": 2, "Stage IIA": 2, "Stage IIB": 2,
             "Stage III": 3, "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
             "Stage IV": 4}
surv = survival.drop_duplicates("sample").set_index("sample")
rows = []
for s in first_by_patient.values():
    patient = s.replace(".", "-")[:12]
    if s not in surv.index or patient not in mut_patients:
        continue
    crow = clin.loc[clin.index.str[:12] == patient]
    if len(crow) == 0:
        continue
    crow = crow.iloc[0]
    stage = stage_map.get(str(crow.get("ajcc_pathologic_stage.diagnoses", "")), np.nan)
    age = crow.get("age_at_index.demographic", np.nan)
    srow = surv.loc[s]
    rows.append({"sample": s, "os_time": srow["OS.time"], "os_event": srow["OS"],
                 "stage_ord": stage, "age": age})
cohort = pd.DataFrame(rows).dropna(subset=["stage_ord", "age", "os_time"])
cohort = cohort[cohort["os_time"] > 0]

counts = np.rint(np.expm1(expr.to_numpy(dtype=np.float64) * np.log(2))).astype(np.int64)
counts = pd.DataFrame(counts, index=expr.index, columns=expr.columns)
probemap = pd.read_csv(RAW / "gencode.v36.annotation.gtf.gene.probemap", sep="\t")
probemap["ensembl_noversion"] = probemap["id"].str.split(".").str[0]
counts.index = counts.index.str.split(".").str[0]
symbol_map = probemap.drop_duplicates("ensembl_noversion").set_index("ensembl_noversion")["gene"]
counts.index = counts.index.map(symbol_map)
counts = counts[counts.index.notna()]
counts = counts.groupby(level=0).sum()
libsize = counts.sum(axis=0)
logexpr = np.log2(counts.div(libsize, axis=1) * 1e6 + 1)

x = logexpr.loc[km_gene, cohort["sample"]].astype(float)
median = x.median()
grp = np.where(x.values > median, f"{km_gene} high", f"{km_gene} low")
fig, ax = plt.subplots(figsize=(8, 6))
kmf = KaplanMeierFitter()
for label in [f"{km_gene} high", f"{km_gene} low"]:
    mask = grp == label
    kmf.fit(cohort["os_time"].values[mask], cohort["os_event"].values[mask], label=label)
    kmf.plot_survival_function(ax=ax)
res = logrank_test(
    cohort["os_time"].values[grp == f"{km_gene} high"],
    cohort["os_time"].values[grp == f"{km_gene} low"],
    cohort["os_event"].values[grp == f"{km_gene} high"],
    cohort["os_event"].values[grp == f"{km_gene} low"],
)
ax.set_xlabel("Days from diagnosis")
ax.set_ylabel("Overall survival probability")
ax.set_title(f"Overall survival by {km_gene} expression (median split,\n"
             f"STK11/KEAP1-mutant subgroup, n={len(cohort)})\n"
             f"log-rank p = {res.p_value:.3g}")
fig.tight_layout()
fig.savefig(FIGS / "km_shortlist_top.png", dpi=200)
plt.close(fig)
print(f"      KM gene: {km_gene} (log-rank p = {res.p_value:.3g}, "
      f"n={len(cohort)} mutant subgroup patients)")
print("Step 5 (dependency-supported shortlist) complete.")
