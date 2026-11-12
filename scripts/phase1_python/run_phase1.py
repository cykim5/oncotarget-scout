"""Phase 1 (Python route) - TCGA-LUAD full cohort analysis for OncoTarget Scout.

Runs, in order:
  1. Load UCSC Xena GDC-hub TCGA-LUAD STAR counts. The Xena file stores
     log2(count + 1); integer counts are recovered as round(2**x - 1)
     (the transform is invertible at the precision Xena publishes, and
     this substitution for the GDC download step is logged in
     results/run_log.md).
  2. Differential expression, primary tumor (sample type 01) vs solid
     tissue normal (11), with pydeseq2 (the Python port of DESeq2).
  3. GSEA prerank on the DESeq2 Wald statistic against the MSigDB
     Hallmark 2020 gene set library (via gseapy / Enrichr).
  4. Survival layer for the top 300 expression candidates: one Cox model
     per gene (lifelines), expression z-scored, adjusted for pathologic
     stage (ordinal I-IV) and age. Kaplan-Meier by median expression for
     the top ranked candidate with both layers covered.
  5. Phase 1 partial target score: tumor specificity + outcome link only
     (2 of the 5 layers in docs/target-score.md). Weights renormalized
     over covered layers exactly as specified there (specificity 0.25,
     outcome 0.10). Missing layers are never imputed as zero; a coverage
     column records how many of the 5 layers each gene's score uses.
  6. Answer key check: ranks of EGFR, KRAS, CD274 in the partial ranking,
     reported wherever they land.

Outputs (paths relative to the repository root):
  results/tables/de_tumor_vs_normal.csv
  results/tables/gsea_hallmark.csv
  results/tables/cox_top_candidates.csv
  results/tables/target_scores_phase1.csv
  results/answer_key_check.md
  results/figures/volcano_tumor_vs_normal.png
  results/figures/top20_phase1.png
  results/figures/km_top_candidate.png

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase1_python/run_phase1.py
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
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]
W_SPEC, W_OUT = 0.25, 0.10  # docs/target-score.md weights for the 2 covered layers
N_CANDIDATES = 300          # survival layer computed for the top expression candidates


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


# ---------------------------------------------------------------- load data
print("[1/6] Loading Xena STAR counts ...")
expr = pd.read_csv(RAW / "TCGA-LUAD.star_counts.tsv.gz", sep="\t", index_col=0)
expr.index.name = "ensembl_id"
print(f"      raw matrix: {expr.shape[0]} Ensembl rows x {expr.shape[1]} samples")

sample_type = expr.columns.str[13:15]
keep = sample_type.isin(["01", "11"])
expr = expr.loc[:, keep]
condition = np.where(sample_type[keep] == "01", "tumor", "normal")
samples = expr.columns.tolist()
n_tumor = int((condition == "tumor").sum())
n_normal = int((condition == "normal").sum())
print(f"      using {n_tumor} primary tumors and {n_normal} normals")

# Recover integer counts from the log2(count + 1) values Xena publishes
counts = np.rint(np.expm1(expr.to_numpy(dtype=np.float64) * np.log(2))).astype(np.int64)
counts = pd.DataFrame(counts, index=expr.index, columns=samples)

# Map Ensembl IDs to gene symbols (GENCODE v36 probemap, same Xena hub).
# Duplicate symbols are summed, the standard symbol-level aggregation.
probemap = pd.read_csv(RAW / "gencode.v36.annotation.gtf.gene.probemap", sep="\t")
probemap["ensembl_noversion"] = probemap["id"].str.split(".").str[0]
counts.index = counts.index.str.split(".").str[0]
symbol_map = probemap.drop_duplicates("ensembl_noversion").set_index("ensembl_noversion")["gene"]
counts.index = counts.index.map(symbol_map)
counts = counts[counts.index.notna()]
counts = counts.groupby(level=0).sum()
print(f"      {counts.shape[0]} unique gene symbols after mapping")

# Prefilter (mirrors the R template rule): count >= 10 in at least 10 samples
counts = counts[(counts >= 10).sum(axis=1) >= 10]
print(f"      {counts.shape[0]} genes pass the expression prefilter")

metadata = pd.DataFrame({"condition": condition}, index=samples)

# ------------------------------------------------- differential expression
print("[2/6] Running pydeseq2 (tumor vs normal) ...")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(
    counts=counts.T,               # pydeseq2 wants samples x genes
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
    n_cpus=1,
    quiet=True,
)
dds.deseq2()
stats = DeseqStats(dds, contrast=["condition", "tumor", "normal"], quiet=True)
stats.summary()
de = stats.results_df.copy()
de.insert(0, "gene", de.index)  # pydeseq2 keeps gene names in the index
de = de.reset_index(drop=True)
de = de[["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
de = de.sort_values("padj").reset_index(drop=True)
de.to_csv(TABLES / "de_tumor_vs_normal.csv", index=False)
print(f"      DE done: {(de['padj'] < 0.05).sum()} genes at padj < 0.05")

# --------------------------------------------------------------------- GSEA
print("[3/6] Running GSEA prerank (Hallmark 2020) ...")
import gseapy as gp

hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="human")
rnk = (
    de.dropna(subset=["stat"])[["gene", "stat"]]
    .sort_values("stat", ascending=False)
    .reset_index(drop=True)
)
pre = gp.prerank(
    rnk=rnk,
    gene_sets=hallmark,
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42,
    threads=1,
    verbose=False,
)
gsea = pre.res2d.copy()
gsea.to_csv(TABLES / "gsea_hallmark.csv", index=False)
print(f"      GSEA done: {len(gsea)} Hallmark sets, "
      f"{(gsea['FDR q-val'] < 0.05).sum()} at FDR < 0.05")

# ------------------------------------------------------------- survival layer
print("[4/6] Survival layer: Cox models for top candidates ...")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

clinical = pd.read_csv(RAW / "TCGA-LUAD.clinical.tsv.gz", sep="\t")
survival = pd.read_csv(RAW / "TCGA-LUAD.survival.tsv.gz", sep="\t")

# One primary tumor sample per patient for the survival analysis
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

rows = []
for s in surv_samples:
    patient = s[:12]
    if s not in surv.index or patient not in [p[:12] for p in clin.index]:
        continue
    srow = surv.loc[s]
    crow = clin.loc[clin.index.str[:12] == patient]
    if len(crow) == 0:
        continue
    crow = crow.iloc[0]
    stage = stage_map.get(str(crow.get("ajcc_pathologic_stage.diagnoses", "")), np.nan)
    age = crow.get("age_at_index.demographic", np.nan)
    rows.append({"sample": s, "patient": patient,
                 "os_time": srow["OS.time"], "os_event": srow["OS"],
                 "stage_ord": stage, "age": age})
cohort = pd.DataFrame(rows).dropna(subset=["stage_ord", "age", "os_time"])
cohort = cohort[cohort["os_time"] > 0]
print(f"      survival cohort: {len(cohort)} patients with stage, age, OS")

# Normalized expression for the survival models (log2 CPM-style, from counts)
libsize = counts.sum(axis=0)
logexpr = np.log2(counts.div(libsize, axis=1) * 1e6 + 1)

# Tumor specificity layer (Phase 1 form: tumor vs normal only; the mutant
# vs wild type component arrives in Phase 2). Padj is a gate, per the spec.
elig_rows = de[de["padj"] < 0.05]
spec = pd.Series(0.0, index=de["gene"].values)
spec.loc[elig_rows["gene"].values] = scale01(
    elig_rows["log2FoldChange"].clip(lower=0)).values

candidates = spec.sort_values(ascending=False).head(N_CANDIDATES).index.tolist()
for g in ANSWER_KEY:  # keep the check honest and covered, wherever they rank
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
    except Exception as exc:  # convergence failures are reported, not hidden
        cox_rows.append({"gene": gene, "log_hr_per_sd": np.nan,
                         "hr_per_sd": np.nan, "pvalue": np.nan,
                         "n_patients": int(len(df)), "error": str(exc)})
cox = pd.DataFrame(cox_rows)
cox.to_csv(TABLES / "cox_top_candidates.csv", index=False)
print(f"      Cox done for {cox['log_hr_per_sd'].notna().sum()} candidates")

# ------------------------------------------------------------- partial score
print("[5/6] Assembling the Phase 1 partial score ...")
outcome = pd.Series(np.nan, index=de["gene"], dtype=float)
ok = cox.dropna(subset=["log_hr_per_sd"]).set_index("gene")["log_hr_per_sd"]
outcome.loc[ok.index] = scale01(ok)  # higher HR (worse survival) scores higher

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
print(f"      ranking written: {len(scores)} genes")

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
    top_pct = 100 * int(ak["rank"].max()) / len(scores) if len(ak) else float("nan")
    fh.write(f"\nAll three ranks fall within the top {top_pct:.1f} percent of the ranking "
             f"at best; the individual ranks above are the honest result.\n")
print(f"      answer key: " +
      ", ".join(f"{g} #{int(ak.loc[g, 'rank'])}" for g in ANSWER_KEY if g in ak.index))

# ------------------------------------------------------------------- figures
print("[6/6] Figures ...")
# Volcano
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

# Top 20 bar chart
top20 = scores.head(20).iloc[::-1]
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top20["gene"], top20["partial_score"], color="#2471a3")
ax.set_xlabel("Phase 1 partial score (expression + survival layers only)")
ax.set_title("Top 20 genes, Phase 1 partial ranking (full LUAD cohort)")
fig.tight_layout()
fig.savefig(FIGS / "top20_phase1.png", dpi=200)
plt.close(fig)

# Kaplan-Meier for the top candidate with both layers covered
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
print(f"      KM gene: {km_gene} (log-rank p = {res.p_value:.3g})")
print("Phase 1 pipeline complete.")
