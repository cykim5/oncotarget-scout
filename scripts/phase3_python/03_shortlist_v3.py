"""Phase 3, step 3: the dependency-supported shortlist, v3.

Same rule as Phase 2 (scripts/phase2_python/06_dependency_shortlist.py):
the decision shortlist is the top 20 genes by main score among genes
with dependency layer support (dependency > 0, i.e. selectively
essential in STK11/KEAP1-mutant lung lines and actually essential
there). Robustness is recomputed within this subset under the same
three weightings as script 02: robust means top 20 in the subset
under all three, fragile means top 20 in the subset under exactly
one. The v2 shortlist file is kept as the Phase 2 record; this
script writes the v3 shortlist under its own name.

It also draws two figures: the Kaplan-Meier figure for the top
shortlist gene (mutant subgroup survival cohort, same cohort
definition as script 02), and a stacked compartment bar chart of
the single cell layer for the v3 shortlist genes plus the three
answer key genes.

Outputs:
  results/tables/shortlist_dependency_supported_phase3.csv
  results/figures/km_shortlist_top_phase3.png
  results/figures/singlecell_shortlist_shares.png

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/03_shortlist_v3.py
"""

import gc
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

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]

print("[1/3] Building the dependency-supported shortlist (v3) ...")
scores = pd.read_csv(TABLES / "target_scores_phase3.csv")
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
shortlist.to_csv(TABLES / "shortlist_dependency_supported_phase3.csv", index=False)
print(f"      shortlist robust: {len(robust)}, fragile: {len(fragile)}")
for _, r in shortlist.iterrows():
    print(f"        #{int(r['subset_rank_main'])} {r['gene']} (overall rank "
          f"{int(r['rank'])}) score {r['score']:.3f} dep {r['dependency']:.3f} "
          f"mal {r['malignant_share']:.3f} "
          f"{r['shortlist_robust_fragile'] or '-'}")

v2 = pd.read_csv(TABLES / "shortlist_dependency_supported.csv")
v2_genes = set(v2["gene"])
v3_genes = set(shortlist["gene"])
print(f"      vs Phase 2 shortlist: {len(v3_genes & v2_genes)} genes stay, "
      f"{len(v3_genes - v2_genes)} enter ({sorted(v3_genes - v2_genes)}), "
      f"{len(v2_genes - v3_genes)} leave ({sorted(v2_genes - v3_genes)})")

print("[2/3] Kaplan-Meier for the top shortlist gene ...")
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
del expr
gc.collect()
probemap = pd.read_csv(RAW / "gencode.v36.annotation.gtf.gene.probemap", sep="\t")
probemap["ensembl_noversion"] = probemap["id"].str.split(".").str[0]
counts.index = counts.index.str.split(".").str[0]
symbol_map = probemap.drop_duplicates("ensembl_noversion").set_index("ensembl_noversion")["gene"]
counts.index = counts.index.map(symbol_map)
counts = counts[counts.index.notna()]
counts = counts.groupby(level=0).sum()
libsize = counts.sum(axis=0)
logexpr = np.log2(counts.div(libsize, axis=1) * 1e6 + 1)
del counts
gc.collect()
logexpr = logexpr.loc[[km_gene], cohort["sample"]]

x = logexpr.loc[km_gene].astype(float)
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
fig.savefig(FIGS / "km_shortlist_top_phase3.png", dpi=200)
plt.close(fig)
print(f"      KM gene: {km_gene} (log-rank p = {res.p_value:.3g}, "
      f"n={len(cohort)} mutant subgroup patients)")

print("[3/3] Single cell compartment shares for the shortlist ...")
scl = pd.read_csv(TABLES / "singlecell_layer.csv").set_index("gene")
genes = list(shortlist["gene"]) + [g for g in ANSWER_KEY if g not in set(shortlist["gene"])]
rows = []
for g in genes:
    if g in scl.index:
        r = scl.loc[g]
        rows.append({"gene": g,
                     "malignant": r["malignant_share"],
                     "immune": r["immune_share"],
                     "stromal": r["stromal_share"],
                     "other": r["other_share"]})
    else:
        rows.append({"gene": g, "malignant": np.nan, "immune": np.nan,
                     "stromal": np.nan, "other": np.nan})
plot = pd.DataFrame(rows).fillna(0).sort_values("malignant")
fig, ax = plt.subplots(figsize=(8, 7))
left = np.zeros(len(plot))
colors = {"malignant": "#c0392b", "immune": "#2471a3",
          "stromal": "#1e8449", "other": "#b0b0b0"}
for comp in ("malignant", "immune", "stromal", "other"):
    ax.barh(plot["gene"], plot[comp], left=left, color=colors[comp], label=comp)
    left += plot[comp].values
ax.set_xlabel("Share of the gene's total tumor expression, by cell compartment")
ax.set_title("Where each shortlist gene is expressed (single cell layer)\n"
             "red = malignant cells; genes not detected score zero here, see table")
ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.08),
          ncol=4, frameon=False)
fig.tight_layout()
fig.savefig(FIGS / "singlecell_shortlist_shares.png", dpi=200)
plt.close(fig)
print("Step 3 (v3 shortlist, KM, compartment figure) complete.")
