"""Phase 2, step 4: rebuilt ranking, sensitivity analysis, answer key, figures.

Assembles the four layers now available and rebuilds the target ranking
with the weights from docs/target-score.md:
  dependency 0.30, tumor specificity 0.25, safety window 0.20,
  outcome link 0.10 (malignant cell share 0.15 is still missing; it
  arrives with the single cell layer in Phase 3). Scores are
  renormalized over covered layers; missing layers are never imputed
  as zero. A coverage column records how many of the 5 layers support
  each gene (maximum 4 in Phase 2).

Layer definitions (all 0 to 1):
  - Tumor specificity: average of two padj-gated components, each
    scaled 0 to 1 across genes: (a) shrunk tumor vs normal log2FC
    clipped at 0, (b) mutant vs wild type log2FC clipped at 0.
    A gene must be up in tumors AND enriched in the mutant subgroup
    to score high on both halves.
  - Dependency: -selective_essentiality from the DepMap layer, gated
    on the gene actually being essential in mutant lung lines (mean
    Chronos < -0.5) and scaled 0 to 1 across covered genes. A gene
    with DepMap coverage that fails the gate scores 0 on this layer
    (evidence against it, not missing evidence); only genes with no
    DepMap coverage at all are missing here. Genes in the Achilles
    common essential list get a dependency value of 0 (the
    pan-essential penalty from the spec) and are flagged.
  - Safety window: the GTEx layer value, already 0 to 1.
  - Outcome link: subgroup Cox log hazard ratio per standard
    deviation of expression, scaled 0 to 1 across computed genes.
    Computed for the top 500 genes of a preliminary ranking (the
    three non-survival layers) plus the three answer key genes.
    The Cox model runs in STK11/KEAP1-mutant patients only, per the
    spec; if that subgroup has under 30 deaths the script falls back
    to the full survival cohort and logs the change.

Sensitivity analysis (pre-committed in docs/target-score.md): the
ranking is rebuilt under equal weights (0.25 per available layer) and
under expression-first weights (specificity 0.50, dependency 0.20,
safety 0.20, outcome 0.10). Candidates in the top 20 under all three
weightings are labeled robust; candidates in the top 20 under exactly
one weighting are labeled fragile. Spearman rank correlations between
the three rankings are reported.

Exclusion rule from docs/target-score.md: genes with no DepMap
coverage and no expression signal (specificity 0) are excluded from
the ranked list and written to results/tables/unscoreable_genes.csv.

Outputs:
  results/tables/cox_phase2.csv
  results/tables/target_scores_phase2.csv
  results/tables/sensitivity_phase2.csv
  results/tables/unscoreable_genes.csv
  results/answer_key_check.md (Phase 2 section appended)
  results/figures/volcano_mutant_vs_wt.png
  results/figures/dependency_scatter.png
  results/figures/top20_phase2.png
  results/figures/km_phase2_top.png

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/05_scoring.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "xena"
CBIO = ROOT / "data" / "raw" / "cbioportal"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]
W = {"dependency": 0.30, "tumor_specificity": 0.25,
     "safety_window": 0.20, "outcome_link": 0.10}
W_EQUAL = {k: 0.25 for k in W}
W_EXPR = {"dependency": 0.20, "tumor_specificity": 0.50,
          "safety_window": 0.20, "outcome_link": 0.10}
N_COX = 500
DEP_ESSENTIALITY_GATE = -0.5  # mean Chronos in mutant lung lines


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


def weighted_score(df, weights):
    """Renormalized weighted mean over covered (non-NaN) layers."""
    num = sum(weights[k] * df[k].fillna(0) for k in weights)
    den = sum(weights[k] * df[k].notna() for k in weights)
    return num / den.replace(0, np.nan)


# ------------------------------------------------------------------ layers
print("[1/6] Loading layer tables ...")
de_tv = pd.read_csv(TABLES / "de_tumor_vs_normal_shrunk.csv")
de_mv = pd.read_csv(TABLES / "de_mutant_vs_wt.csv")
dep = pd.read_csv(TABLES / "dependency_layer.csv")
saf = pd.read_csv(TABLES / "safety_layer.csv")

universe = de_tv["gene"]
spec_a = pd.Series(0.0, index=universe.values)
elig_a = de_tv[de_tv["padj"] < 0.05]
spec_a.loc[elig_a["gene"].values] = scale01(
    elig_a["log2FoldChange"].clip(lower=0)).values
spec_b = pd.Series(0.0, index=universe.values)
elig_b = de_mv[de_mv["padj"] < 0.05]
spec_b.loc[elig_b["gene"].values] = scale01(
    elig_b["log2FoldChange"].clip(lower=0)).values
specificity = (spec_a + spec_b) / 2

dep_idx = dep.set_index("gene")
gated = dep_idx["mean_chronos_mutant_lung"] < DEP_ESSENTIALITY_GATE
covered = dep_idx["selective_essentiality"].notna()
dep_raw = -dep_idx["selective_essentiality"]
# A gene WITH DepMap coverage that fails the gate (not essential in
# mutant lines, or not selectively essential) scores 0 on this layer:
# that is evidence against it, not missing evidence. Only genes with
# no DepMap coverage at all stay NaN. (An earlier version of this
# script set failing genes to NaN, which wrongly excused them from the
# highest-weighted layer via renormalization; see run_log.md.)
dep_raw = dep_raw.where(gated & (dep_raw > 0), 0.0)
dep_raw = dep_raw.where(covered)
dependency = scale01(dep_raw)
dependency.loc[dep_idx["common_essential"]] = 0.0  # pan-essential penalty

safety = saf.set_index("gene")["safety_window"]

layers = pd.DataFrame({
    "gene": universe,
    "tumor_specificity": specificity.values,
    "dependency": dependency.reindex(universe).values,
    "safety_window": safety.reindex(universe).values,
})
layers["log2FC_tumor_vs_normal_shrunk"] = de_tv.set_index("gene")["log2FoldChange"].reindex(universe).values
layers["log2FC_mutant_vs_wt"] = de_mv.set_index("gene")["log2FoldChange"].reindex(universe).values
layers["selective_essentiality"] = dep_idx["selective_essentiality"].reindex(universe).values
layers["mean_chronos_mutant_lung"] = dep_idx["mean_chronos_mutant_lung"].reindex(universe).values
layers["common_essential"] = dep_idx["common_essential"].reindex(universe).fillna(False).astype(bool).values
print(f"      universe: {len(layers)} genes; dependency coverage "
      f"{layers['dependency'].notna().sum()}, safety coverage "
      f"{layers['safety_window'].notna().sum()}")

# ------------------------------------------- preliminary ranking -> Cox set
pre = layers.copy()
pre["prelim_score"] = weighted_score(
    pre, {"dependency": W["dependency"], "tumor_specificity": W["tumor_specificity"],
          "safety_window": W["safety_window"]})
pre = pre.sort_values("prelim_score", ascending=False)
cox_genes = pre.head(N_COX)["gene"].tolist()
for g in ANSWER_KEY:
    if g in layers["gene"].values and g not in cox_genes:
        cox_genes.append(g)
print(f"      Cox set: {len(cox_genes)} genes (top {N_COX} preliminary + answer key)")

# ------------------------------------------------------------- survival data
print("[2/6] Building the survival cohort (mutant subgroup) ...")
clinical = pd.read_csv(RAW / "TCGA-LUAD.clinical.tsv.gz", sep="\t")
survival = pd.read_csv(RAW / "TCGA-LUAD.survival.tsv.gz", sep="\t")
expr = pd.read_csv(RAW / "TCGA-LUAD.star_counts.tsv.gz", sep="\t", index_col=0)
sample_type = expr.columns.str[13:15]
tumor_samples = expr.columns[sample_type == "01"].tolist()
first_by_patient = {}
for s in tumor_samples:
    first_by_patient.setdefault(s.replace(".", "-")[:12], s)
surv_samples = list(first_by_patient.values())

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
for s in surv_samples:
    patient = s.replace(".", "-")[:12]
    if s not in surv.index:
        continue
    crow = clin.loc[clin.index.str[:12] == patient]
    if len(crow) == 0:
        continue
    crow = crow.iloc[0]
    stage = stage_map.get(str(crow.get("ajcc_pathologic_stage.diagnoses", "")), np.nan)
    age = crow.get("age_at_index.demographic", np.nan)
    srow = surv.loc[s]
    rows.append({"sample": s, "patient": patient,
                 "os_time": srow["OS.time"], "os_event": srow["OS"],
                 "stage_ord": stage, "age": age,
                 "mutant": patient in mut_patients})
cohort = pd.DataFrame(rows).dropna(subset=["stage_ord", "age", "os_time"])
cohort = cohort[cohort["os_time"] > 0]
sub = cohort[cohort["mutant"]]
print(f"      full survival cohort: {len(cohort)} patients, "
      f"{int(cohort['os_event'].sum())} deaths")
print(f"      mutant subgroup cohort: {len(sub)} patients, "
      f"{int(sub['os_event'].sum())} deaths")
if len(sub) >= 50 and int(sub["os_event"].sum()) >= 30:
    cox_cohort, cox_label = sub, "STK11/KEAP1-mutant subgroup"
else:
    cox_cohort, cox_label = cohort, "full cohort (subgroup too small; see run log)"
print(f"      Cox models run in: {cox_label} (n={len(cox_cohort)})")

# log2 CPM expression for the Cox models (same transform as Phase 1)
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

print(f"[3/6] Cox models for {len(cox_genes)} genes ({cox_label}) ...")
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

cox_rows = []
expr_cohort = logexpr.loc[:, cox_cohort["sample"]]
for gene in cox_genes:
    if gene not in expr_cohort.index:
        continue
    x = expr_cohort.loc[gene].astype(float)
    z = (x - x.mean()) / x.std()
    df = cox_cohort[["os_time", "os_event", "stage_ord", "age"]].copy()
    df["expr_z"] = z.values
    try:
        cph = CoxPHFitter()
        cph.fit(df, duration_col="os_time", event_col="os_event")
        cox_rows.append({"gene": gene,
                         "log_hr_per_sd": float(cph.params_["expr_z"]),
                         "hr_per_sd": float(np.exp(cph.params_["expr_z"])),
                         "pvalue": float(cph.summary.loc["expr_z", "p"]),
                         "n_patients": int(len(df)),
                         "cox_cohort": cox_label})
    except Exception as exc:
        cox_rows.append({"gene": gene, "log_hr_per_sd": np.nan,
                         "hr_per_sd": np.nan, "pvalue": np.nan,
                         "n_patients": int(len(df)), "cox_cohort": cox_label,
                         "error": str(exc)})
cox = pd.DataFrame(cox_rows)
cox.to_csv(TABLES / "cox_phase2.csv", index=False)
print(f"      Cox done for {cox['log_hr_per_sd'].notna().sum()} genes")

# ------------------------------------------------------------- final ranking
print("[4/6] Assembling the Phase 2 ranking ...")
outcome = pd.Series(np.nan, index=universe, dtype=float)
ok = cox.dropna(subset=["log_hr_per_sd"]).set_index("gene")["log_hr_per_sd"]
outcome.loc[ok.index] = scale01(ok)
layers["outcome_link_cox_log_hr"] = ok.reindex(universe).values
layers["outcome_link"] = outcome.values

scores = layers.copy()
scores["score"] = weighted_score(scores, W)
scores["coverage"] = scores[["tumor_specificity", "dependency",
                             "safety_window", "outcome_link"]].notna().sum(axis=1)
LAYER_COLS = [("expression", "tumor_specificity"), ("dependency", "dependency"),
              ("safety", "safety_window"), ("outcome", "outcome_link")]
scores["layers_used"] = scores.apply(
    lambda r: "+".join(name for name, col in LAYER_COLS if pd.notna(r[col])),
    axis=1)
scores["known_target"] = scores["gene"].isin(ANSWER_KEY)

# Exclusion rule: no DepMap coverage and no expression signal -> unscoreable
unscoreable = scores[scores["dependency"].isna() & (scores["tumor_specificity"] == 0)]
unscoreable.to_csv(TABLES / "unscoreable_genes.csv", index=False)
print(f"      excluded as unscoreable (no DepMap coverage, no expression "
      f"signal): {len(unscoreable)} genes")
scores = scores.drop(index=unscoreable.index)
scores = scores.sort_values("score", ascending=False).reset_index(drop=True)
scores.insert(0, "rank", np.arange(1, len(scores) + 1))

# ------------------------------------------------------- sensitivity analysis
print("[5/6] Sensitivity analysis (equal weights, expression-first) ...")
for name, weights in [("equal", W_EQUAL), ("expr_first", W_EXPR)]:
    scores[f"score_{name}"] = weighted_score(scores, weights)
    scores[f"rank_{name}"] = scores[f"score_{name}"].rank(
        ascending=False, method="min").astype(int)
scores["rank_main"] = scores["rank"]
top20_sets = {
    "main": set(scores.nsmallest(20, "rank_main")["gene"]),
    "equal": set(scores.nsmallest(20, "rank_equal")["gene"]),
    "expr_first": set(scores.nsmallest(20, "rank_expr_first")["gene"]),
}
robust = top20_sets["main"] & top20_sets["equal"] & top20_sets["expr_first"]
in_one = (top20_sets["main"] ^ top20_sets["equal"] ^ top20_sets["expr_first"])
fragile = {g for g in (top20_sets["main"] | top20_sets["equal"] | top20_sets["expr_first"])
           if sum(g in s for s in top20_sets.values()) == 1}
scores["robust_fragile"] = ""
scores.loc[scores["gene"].isin(robust), "robust_fragile"] = "robust"
scores.loc[scores["gene"].isin(fragile), "robust_fragile"] = "fragile"
rho_main_eq, _ = spearmanr(scores["score"], scores["score_equal"])
rho_main_ex, _ = spearmanr(scores["score"], scores["score_expr_first"])
print(f"      robust (top 20 under all three weightings): {len(robust)}")
print(f"      fragile (top 20 under exactly one weighting): {len(fragile)}")
print(f"      Spearman main vs equal: {rho_main_eq:.3f}; "
      f"main vs expression-first: {rho_main_ex:.3f}")

scores.to_csv(TABLES / "target_scores_phase2.csv", index=False)
sens = scores[["gene", "rank_main", "rank_equal", "rank_expr_first",
               "robust_fragile"]].copy()
sens = sens[(sens["rank_main"] <= 20) | (sens["rank_equal"] <= 20)
            | (sens["rank_expr_first"] <= 20)].sort_values("rank_main")
sens.to_csv(TABLES / "sensitivity_phase2.csv", index=False)
print(f"      ranking written: {len(scores)} genes; top 5:")
for _, r in scores.head(5).iterrows():
    print(f"        #{int(r['rank'])} {r['gene']} score {r['score']:.3f} "
          f"coverage {int(r['coverage'])}/5 {r['robust_fragile'] or '-'}")

# ------------------------------------------------------------- answer key
ak = scores[scores["gene"].isin(ANSWER_KEY)].set_index("gene")
with open(ROOT / "results" / "answer_key_check.md", "a") as fh:
    fh.write("\n\n# Answer key check v2 (Phase 2, STK11/KEAP1 subgroup ranking)\n\n")
    fh.write("Same pre-committed check, now on the Phase 2 ranking, which\n"
             "scores genes for the STK11/KEAP1-mutant subgroup across four\n"
             "layers (expression incl. mutant vs wild type, DepMap\n"
             "dependency, GTEx safety, subgroup survival). Ranks are\n"
             "reported wherever they land, untuned.\n\n")
    fh.write(f"Ranking size: {len(scores)} genes (4 of 5 layers; the single\n"
             f"cell layer arrives in Phase 3; see results/tables/\n"
             f"target_scores_phase2.csv).\n\n")
    fh.write("| Gene | Rank | Score | log2FC mutant vs WT | Selective essentiality | "
             "Safety window |\n")
    fh.write("| --- | --- | --- | --- | --- | --- |\n")
    for g in ANSWER_KEY:
        if g in ak.index:
            r = ak.loc[g]
            fh.write(f"| {g} | {int(r['rank'])} | {r['score']:.3f} | "
                     f"{r['log2FC_mutant_vs_wt']:.3f} | "
                     f"{r['selective_essentiality']:.3f} | {r['safety_window']:.3f} |\n")
        else:
            fh.write(f"| {g} | not in ranked table | - | - | - | - |\n")
    fh.write(f"\nPhase 1 ranks for comparison (full cohort, expression +\n"
             f"survival only): KRAS #155, EGFR #282, CD274 #457 of 35,273.\n")
print("      answer key: " +
      ", ".join(f"{g} #{int(ak.loc[g, 'rank'])}" for g in ANSWER_KEY if g in ak.index))

# ------------------------------------------------------------------- figures
print("[6/6] Figures ...")
# Subgroup volcano
v = de_mv.dropna(subset=["padj"]).copy()
v["neglog10p"] = -np.log10(v["padj"].clip(lower=1e-300))
fig, ax = plt.subplots(figsize=(8, 6))
sig = v["padj"] < 0.05
ax.scatter(v.loc[~sig, "log2FoldChange"], v.loc[~sig, "neglog10p"],
           s=4, c="grey", alpha=0.4, label="padj >= 0.05")
ax.scatter(v.loc[sig & (v["log2FoldChange"] > 0), "log2FoldChange"],
           v.loc[sig & (v["log2FoldChange"] > 0), "neglog10p"],
           s=4, c="#c0392b", alpha=0.5, label="up in mutant, padj < 0.05")
ax.scatter(v.loc[sig & (v["log2FoldChange"] < 0), "log2FoldChange"],
           v.loc[sig & (v["log2FoldChange"] < 0), "neglog10p"],
           s=4, c="#2471a3", alpha=0.5, label="down in mutant, padj < 0.05")
top_lab = v.reindex(v["neglog10p"].sort_values(ascending=False).index).head(8)
for _, r in pd.concat([top_lab, v[v["gene"].isin(ANSWER_KEY)]]).drop_duplicates("gene").iterrows():
    ax.annotate(r["gene"], (r["log2FoldChange"], r["neglog10p"]), fontsize=8)
ax.set_xlabel("log2 fold change (STK11/KEAP1-mutant vs wild type tumors)")
ax.set_ylabel("-log10 adjusted p value")
ax.set_title("TCGA-LUAD: mutant vs wild type differential expression")
ax.legend(markerscale=3, fontsize=8)
fig.tight_layout()
fig.savefig(FIGS / "volcano_mutant_vs_wt.png", dpi=200)
plt.close(fig)

# Dependency scatter
d = dep.dropna(subset=["selective_essentiality"]).copy()
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(d["mean_chronos_mutant_lung"], d["selective_essentiality"],
           s=4, c="grey", alpha=0.4)
top_dep = scores.head(20)
lab = d[d["gene"].isin(top_dep["gene"]) | d["gene"].isin(ANSWER_KEY)]
ax.scatter(lab["mean_chronos_mutant_lung"], lab["selective_essentiality"],
           s=18, c="#c0392b", alpha=0.8)
for _, r in lab.iterrows():
    ax.annotate(r["gene"], (r["mean_chronos_mutant_lung"], r["selective_essentiality"]),
                fontsize=8)
ax.axvline(-0.5, color="#2471a3", linestyle="--", linewidth=1,
           label="essentiality gate (mean Chronos -0.5)")
ax.set_xlabel("Mean Chronos score in STK11/KEAP1-mutant lung lines (more negative = more essential)")
ax.set_ylabel("Selective essentiality (mutant lung mean minus other lung mean)")
ax.set_title("DepMap dependency: mutant lung lines vs other lung lines")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGS / "dependency_scatter.png", dpi=200)
plt.close(fig)

# Top 20 bar chart, colored by robust/fragile
top20 = scores.head(20).iloc[::-1].copy()
colors = top20["robust_fragile"].map(
    {"robust": "#1e8449", "fragile": "#c0392b"}).fillna("#2471a3")
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top20["gene"], top20["score"], color=colors)
ax.set_xlabel("Phase 2 score (dependency 0.30, expression 0.25, safety 0.20, outcome 0.10)")
ax.set_title("Top 20 genes, Phase 2 subgroup ranking\n"
             "green = robust across weightings, red = fragile, blue = in between")
fig.tight_layout()
fig.savefig(FIGS / "top20_phase2.png", dpi=200)
plt.close(fig)

# Kaplan-Meier for the top robust candidate (in the Cox cohort)
robust_ranked = scores[scores["robust_fragile"] == "robust"]
km_gene = robust_ranked.iloc[0]["gene"] if len(robust_ranked) else scores.iloc[0]["gene"]
x = logexpr.loc[km_gene, cox_cohort["sample"]].astype(float)
median = x.median()
grp = np.where(x.values > median, f"{km_gene} high", f"{km_gene} low")
fig, ax = plt.subplots(figsize=(8, 6))
kmf = KaplanMeierFitter()
for label in [f"{km_gene} high", f"{km_gene} low"]:
    mask = grp == label
    kmf.fit(cox_cohort["os_time"].values[mask], cox_cohort["os_event"].values[mask], label=label)
    kmf.plot_survival_function(ax=ax)
res = logrank_test(
    cox_cohort["os_time"].values[grp == f"{km_gene} high"],
    cox_cohort["os_time"].values[grp == f"{km_gene} low"],
    cox_cohort["os_event"].values[grp == f"{km_gene} high"],
    cox_cohort["os_event"].values[grp == f"{km_gene} low"],
)
ax.set_xlabel("Days from diagnosis")
ax.set_ylabel("Overall survival probability")
ax.set_title(f"Overall survival by {km_gene} expression (median split, {cox_label})\n"
             f"log-rank p = {res.p_value:.3g}")
fig.tight_layout()
fig.savefig(FIGS / "km_phase2_top.png", dpi=200)
plt.close(fig)
print(f"      KM gene: {km_gene} (log-rank p = {res.p_value:.3g})")
print("Step 4 (scoring, sensitivity, figures) complete.")
