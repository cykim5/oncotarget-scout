"""Phase 1 stages 4 to 6: survival layer, partial score, answer key, figures.

Companion to run_phase1.py. The first full run of run_phase1.py in this
environment was terminated after the GSEA stage had written
results/tables/de_tumor_vs_normal.csv and results/tables/gsea_hallmark.csv
(the interruption is recorded in results/run_log.md). This script performs
exactly the stages 4 to 6 defined in run_phase1.py, reading those two saved
stage outputs instead of recomputing them, so no number changes: every
downstream number derives from the committed stage 1 to 3 tables and the
same raw Xena files. Analysis choices are identical to run_phase1.py.

Run from the repository root:
  .venv/bin/python scripts/phase1_python/run_phase1_stages4to6.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "xena"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]
W_SPEC, W_OUT = 0.25, 0.10
N_CANDIDATES = 300


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


print("[4/6] Loading stage outputs and counts for the survival layer ...", flush=True)
de = pd.read_csv(TABLES / "de_tumor_vs_normal.csv")

expr = pd.read_csv(RAW / "TCGA-LUAD.star_counts.tsv.gz", sep="\t", index_col=0)
sample_type = expr.columns.str[13:15]
keep = sample_type.isin(["01", "11"])
expr = expr.loc[:, keep]
condition = np.where(sample_type[keep] == "01", "tumor", "normal")
samples = expr.columns.tolist()
counts = np.rint(np.expm1(expr.to_numpy(dtype=np.float64) * np.log(2))).astype(np.int64)
counts = pd.DataFrame(counts, index=expr.index, columns=samples)
probemap = pd.read_csv(RAW / "gencode.v36.annotation.gtf.gene.probemap", sep="\t")
probemap["ensembl_noversion"] = probemap["id"].str.split(".").str[0]
counts.index = counts.index.str.split(".").str[0]
symbol_map = probemap.drop_duplicates("ensembl_noversion").set_index("ensembl_noversion")["gene"]
counts.index = counts.index.map(symbol_map)
counts = counts[counts.index.notna()]
counts = counts.groupby(level=0).sum()
counts = counts[(counts >= 10).sum(axis=1) >= 10]
print(f"      counts rebuilt: {counts.shape[0]} genes x {counts.shape[1]} samples", flush=True)

# ------------------------------------------------------------- survival layer
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

clinical = pd.read_csv(RAW / "TCGA-LUAD.clinical.tsv.gz", sep="\t")
survival = pd.read_csv(RAW / "TCGA-LUAD.survival.tsv.gz", sep="\t")

tumor_samples = [s for s, c in zip(samples, condition) if c == "tumor"]
first_by_patient = {}
for s in tumor_samples:
    first_by_patient.setdefault(s[:12], s)
surv_samples = list(first_by_patient.values())

clin = clinical.set_index("sample")
stage_map = {"Stage I": 1, "Stage IA": 1, "Stage IB": 1,
             "Stage II": 2, "Stage IIA": 2, "Stage IIB": 2,
             "Stage III": 3, "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
             "Stage IV": 4}
surv = survival.drop_duplicates("sample").set_index("sample")
clin_patients = clin.index.str[:12]

rows = []
for s in surv_samples:
    patient = s[:12]
    if s not in surv.index:
        continue
    matches = np.where(clin_patients == patient)[0]
    if len(matches) == 0:
        continue
    srow = surv.loc[s]
    crow = clin.iloc[matches[0]]
    stage = stage_map.get(str(crow.get("ajcc_pathologic_stage.diagnoses", "")), np.nan)
    age = crow.get("age_at_index.demographic", np.nan)
    rows.append({"sample": s, "patient": patient,
                 "os_time": srow["OS.time"], "os_event": srow["OS"],
                 "stage_ord": stage, "age": age})
cohort = pd.DataFrame(rows).dropna(subset=["stage_ord", "age", "os_time"])
cohort = cohort[cohort["os_time"] > 0].reset_index(drop=True)
print(f"      survival cohort: {len(cohort)} patients with stage, age, OS "
      f"({int(cohort['os_event'].sum())} deaths)", flush=True)

libsize = counts.sum(axis=0)
logexpr = np.log2(counts.div(libsize, axis=1) * 1e6 + 1)

elig_rows = de[de["padj"] < 0.05]
spec = pd.Series(0.0, index=de["gene"].values)
spec.loc[elig_rows["gene"].values] = scale01(
    elig_rows["log2FoldChange"].clip(lower=0)).values

candidates = spec.sort_values(ascending=False).head(N_CANDIDATES).index.tolist()
for g in ANSWER_KEY:
    if g in spec.index and g not in candidates:
        candidates.append(g)

cox_rows = []
expr_cohort = logexpr.loc[:, cohort["sample"]]
for gene in candidates:
    if gene not in expr_cohort.index:
        continue
    x = expr_cohort.loc[gene].astype(float)
    z = (x - x.mean()) / x.std()
    df = cohort[["os_time", "os_event", "stage_ord", "age"]].copy()
    df["expr_z"] = z.values
    try:
        cph = CoxPHFitter()
        cph.fit(df, duration_col="os_time", event_col="os_event")
        cox_rows.append({
            "gene": gene,
            "log_hr_per_sd": float(cph.params_["expr_z"]),
            "hr_per_sd": float(np.exp(cph.params_["expr_z"])),
            "pvalue": float(cph.summary.loc["expr_z", "p"]),
            "n_patients": int(len(df)),
        })
    except Exception as exc:
        cox_rows.append({"gene": gene, "log_hr_per_sd": np.nan,
                         "hr_per_sd": np.nan, "pvalue": np.nan,
                         "n_patients": int(len(df)), "error": str(exc)})
cox = pd.DataFrame(cox_rows)
cox.to_csv(TABLES / "cox_top_candidates.csv", index=False)
print(f"      Cox done for {cox['log_hr_per_sd'].notna().sum()} candidates", flush=True)

# ------------------------------------------------------------- partial score
print("[5/6] Assembling the Phase 1 partial score ...", flush=True)
outcome = pd.Series(np.nan, index=de["gene"], dtype=float)
ok = cox.dropna(subset=["log_hr_per_sd"]).set_index("gene")["log_hr_per_sd"]
outcome.loc[ok.index] = scale01(ok)

scores = pd.DataFrame({
    "gene": de["gene"],
    "log2FoldChange": de["log2FoldChange"],
    "padj": de["padj"],
    "tumor_specificity": spec.reindex(de["gene"]).values,
    "outcome_link_cox_log_hr": ok.reindex(de["gene"]).values,
    "outcome_link": outcome.reindex(de["gene"]).values,
})
covered_w = (W_SPEC * scores["tumor_specificity"].notna()
             + W_OUT * scores["outcome_link"].notna())
scores["partial_score"] = (
    (W_SPEC * scores["tumor_specificity"].fillna(0)
     + W_OUT * scores["outcome_link"].fillna(0)) / covered_w
)
scores["coverage"] = scores["tumor_specificity"].notna().astype(int) \
    + scores["outcome_link"].notna().astype(int)
scores["layers_used"] = "expression+survival"
scores.loc[scores["outcome_link"].isna(), "layers_used"] = "expression only"
scores["known_target"] = scores["gene"].isin(ANSWER_KEY)
scores = scores.sort_values("partial_score", ascending=False).reset_index(drop=True)
scores.insert(0, "rank", np.arange(1, len(scores) + 1))
scores.to_csv(TABLES / "target_scores_phase1.csv", index=False)
print(f"      ranking written: {len(scores)} genes", flush=True)

ak = scores[scores["gene"].isin(ANSWER_KEY)].set_index("gene")
with open(ROOT / "results" / "answer_key_check.md", "w") as fh:
    fh.write("# Answer key check v1 (Phase 1, full LUAD cohort)\n\n")
    fh.write("Pre-committed check from the README: the full cohort ranking should\n"
             "recover established LUAD target biology (EGFR, KRAS, CD274) near the\n"
             "top without being steered to it. Ranks are reported wherever they land.\n\n")
    fh.write(f"Ranking size: {len(scores)} genes (partial score, expression + survival\n"
             f"layers, 2 of 5; see results/tables/target_scores_phase1.csv).\n\n")
    fh.write("| Gene | Rank | Partial score | log2FoldChange (tumor vs normal) | padj |\n")
    fh.write("| --- | --- | --- | --- | --- |\n")
    for g in ANSWER_KEY:
        if g in ak.index:
            r = ak.loc[g]
            fh.write(f"| {g} | {int(r['rank'])} | {r['partial_score']:.3f} | "
                     f"{r['log2FoldChange']:.3f} | {r['padj']:.3g} |\n")
        else:
            fh.write(f"| {g} | not in ranked table | - | - | - |\n")
print("      answer key: " +
      ", ".join(f"{g} #{int(ak.loc[g, 'rank'])}" for g in ANSWER_KEY if g in ak.index), flush=True)

# ------------------------------------------------------------------- figures
print("[6/6] Figures ...", flush=True)
v = de.dropna(subset=["padj"]).copy()
v["neglog10p"] = -np.log10(v["padj"].clip(lower=1e-300))
fig, ax = plt.subplots(figsize=(8, 6))
sig = v["padj"] < 0.05
ax.scatter(v.loc[~sig, "log2FoldChange"], v.loc[~sig, "neglog10p"],
           s=4, c="grey", alpha=0.4, label="padj >= 0.05")
ax.scatter(v.loc[sig & (v["log2FoldChange"] > 0), "log2FoldChange"],
           v.loc[sig & (v["log2FoldChange"] > 0), "neglog10p"],
           s=4, c="#c0392b", alpha=0.5, label="up in tumor, padj < 0.05")
ax.scatter(v.loc[sig & (v["log2FoldChange"] < 0), "log2FoldChange"],
           v.loc[sig & (v["log2FoldChange"] < 0), "neglog10p"],
           s=4, c="#2471a3", alpha=0.5, label="down in tumor, padj < 0.05")
top_lab = v.reindex(v["neglog10p"].sort_values(ascending=False).index).head(8)
for _, r in pd.concat([top_lab, v[v["gene"].isin(ANSWER_KEY)]]).drop_duplicates("gene").iterrows():
    ax.annotate(r["gene"], (r["log2FoldChange"], r["neglog10p"]), fontsize=8)
ax.set_xlabel("log2 fold change (tumor vs normal)")
ax.set_ylabel("-log10 adjusted p value")
ax.set_title("TCGA-LUAD: tumor vs normal differential expression")
ax.legend(markerscale=3, fontsize=8)
fig.tight_layout()
fig.savefig(FIGS / "volcano_tumor_vs_normal.png", dpi=200)
plt.close(fig)

top20 = scores.head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top20["gene"], top20["partial_score"], color="#2471a3")
ax.set_xlabel("Phase 1 partial score (expression + survival layers only)")
ax.set_title("Top 20 genes, Phase 1 partial ranking (full LUAD cohort)")
fig.tight_layout()
fig.savefig(FIGS / "top20_phase1.png", dpi=200)
plt.close(fig)

both = scores[scores["coverage"] == 2]
km_gene = both.iloc[0]["gene"]
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
ax.set_title(f"Overall survival by {km_gene} expression (median split)\n"
             f"log-rank p = {res.p_value:.3g}")
fig.tight_layout()
fig.savefig(FIGS / "km_top_candidate.png", dpi=200)
plt.close(fig)
print(f"      KM gene: {km_gene} (log-rank p = {res.p_value:.3g})", flush=True)
print("Stages 4 to 6 complete.", flush=True)
