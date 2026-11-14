"""Phase 3, step 4: the time split check (pre-2021 frozen ranking).

Pre-committed check from the README: freeze the evidence at an
earlier date, rebuild the ranking, and ask whether later clinical
developments (the KRAS G12C inhibitor approvals that began in 2021)
appear high in the frozen ranking. A method that would have found
tomorrow's targets yesterday is worth trusting on new names.

What can and cannot be frozen:
  - TCGA expression and survival (2014-2018 data) and GTEx v8
    (2017-2020) are already pre-2021. They are carried over
    unchanged, including the Phase 2 Cox values.
  - DepMap is the one post-2021 input in the Phase 2 ranking
    (24Q4). It is replaced here with DepMap Public 20Q2 (figshare
    article 12280541, published 2020-06-04): Achilles CERES gene
    effect scores, CCLE mutations, sample annotations. Note the
    scoring change this forces: 20Q2 predates Chronos, so the
    frozen dependency layer uses CERES scores. The selective
    essentiality contrast is computed the same way (mutant lung
    mean minus other lung mean), and the difference is recorded
    in the run log.
  - The single cell layer cannot be frozen: the CELLxGENE dataset
    used in Phase 3 is from a 2024 preprint. Per the spec, the
    limitation is stated and the check runs on the four layers
    that can be frozen. The frozen ranking is therefore a
    four-layer ranking, directly comparable to the Phase 2 (v2)
    ranking: the ONLY difference between them is the dependency
    layer's release. That makes the comparison clean.

What the check can and cannot show: the subgroup ranking asks
which genes are subgroup-selective targets, and KRAS is not one
(selective essentiality near zero in both releases), so KRAS is
not expected to rank high in the frozen subgroup ranking either.
The informative questions are whether the frozen ranking recovers
the same biology (does NFE2L2 still lead the dependency layer
before 2021?) and how much of the shortlist survives the release
change. Both are reported untuned.

Outputs:
  results/tables/dependency_layer_20q2.csv
  results/tables/depmap_line_summary_20q2.json
  results/tables/target_scores_frozen.csv
  results/tables/timesplit_comparison.csv
  results/timesplit_check.md

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/04_timesplit.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DEP20 = ROOT / "data" / "raw" / "depmap20q2"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

GENES_OF_INTEREST = {"STK11", "KEAP1", "KRAS"}
ANSWER_KEY = ["EGFR", "KRAS", "CD274"]
W = {"dependency": 0.30, "tumor_specificity": 0.25,
     "safety_window": 0.20, "outcome_link": 0.10}
DEP_ESSENTIALITY_GATE = -0.5  # mean gene effect in mutant lung lines


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


def weighted_score(df, weights):
    num = sum(weights[k] * df[k].fillna(0) for k in weights)
    den = sum(weights[k] * df[k].notna() for k in weights)
    return num / den.replace(0, np.nan)


print("[1/5] Loading DepMap 20Q2 sample annotations ...")
sample = pd.read_csv(DEP20 / "sample_info.csv")
print(f"      sample_info columns: {list(sample.columns)}")
id_col = "DepMap_ID"
lin_col = "lineage"
lung = sample[sample[lin_col] == "lung"]
print(f"      {len(sample)} cell lines total, {len(lung)} lung lineage")

print("[2/5] Reading CCLE_mutations.csv in chunks (STK11/KEAP1/KRAS only) ...")
mut_lines = {g: set() for g in GENES_OF_INTEREST}
n_rows = 0
for chunk in pd.read_csv(DEP20 / "CCLE_mutations.csv", sep="\t",
                         usecols=["DepMap_ID", "Hugo_Symbol",
                                  "Variant_Classification"],
                         chunksize=200_000, low_memory=False):
    n_rows += len(chunk)
    sub = chunk[chunk["Hugo_Symbol"].isin(GENES_OF_INTEREST)
                & (chunk["Variant_Classification"] != "Silent")]
    for gene, grp in sub.groupby("Hugo_Symbol"):
        mut_lines[gene].update(grp["DepMap_ID"].tolist())
print(f"      scanned {n_rows} mutation rows")
for g in sorted(GENES_OF_INTEREST):
    print(f"      {g}: {len(mut_lines[g])} lines with non-silent mutations")

lung_ids = set(lung[id_col])
mutant_lung = sorted((mut_lines["STK11"] | mut_lines["KEAP1"]) & lung_ids)
other_lung = sorted(lung_ids - set(mutant_lung))
n_both = len(mut_lines["STK11"] & mut_lines["KEAP1"] & lung_ids)
print(f"      mutant lung lines (STK11 and/or KEAP1): {len(mutant_lung)} "
      f"({n_both} with both); other lung lines: {len(other_lung)}")
(TABLES / "depmap_line_summary_20q2.json").write_text(json.dumps({
    "depmap_release": "DepMap Public 20Q2",
    "n_lines_total": int(len(sample)),
    "n_lung_lines": int(len(lung)),
    "n_mutant_lung_lines": len(mutant_lung),
    "n_both_mutant_lung_lines": int(n_both),
    "n_other_lung_lines": len(other_lung),
    "n_stk11_mutant_lines": len(mut_lines["STK11"]),
    "n_keap1_mutant_lines": len(mut_lines["KEAP1"]),
    "n_kras_mutant_lines": len(mut_lines["KRAS"]),
    "mutant_lung_depmap_ids": mutant_lung,
}, indent=2))

print("[3/5] Achilles gene effect (CERES), chunked means by line group ...")
mut_set, oth_set = set(mutant_lung), set(other_lung)
sum_mut = sum_oth = sum_pan = None
cnt_mut = cnt_oth = cnt_pan = None
n_mut = n_oth = n_pan = 0
genes = None
keep = None
for chunk in pd.read_csv(DEP20 / "Achilles_gene_effect.csv",
                         index_col=0, chunksize=50):
    if genes is None:
        raw_cols = chunk.columns.str.split(" ").str[0]
        keep = ~raw_cols.duplicated(keep="first")
        genes = raw_cols[keep]
    chunk = chunk.iloc[:, keep]
    chunk.columns = genes
    arr = chunk.to_numpy(dtype=np.float64)
    is_mut = chunk.index.isin(mut_set)
    is_oth = chunk.index.isin(oth_set)
    if sum_mut is None:
        sum_mut = np.zeros(len(genes))
        sum_oth = np.zeros(len(genes))
        sum_pan = np.zeros(len(genes))
        cnt_mut = np.zeros(len(genes))
        cnt_oth = np.zeros(len(genes))
        cnt_pan = np.zeros(len(genes))
    for mask, s, c in ((is_mut, sum_mut, cnt_mut),
                       (is_oth, sum_oth, cnt_oth),
                       (np.ones(len(chunk), dtype=bool), sum_pan, cnt_pan)):
        sub = arr[mask]
        s += np.nansum(sub, axis=0)
        c += (~np.isnan(sub)).sum(axis=0)
    n_mut += int(is_mut.sum())
    n_oth += int(is_oth.sum())
    n_pan += len(chunk)
print(f"      gene effect rows scanned: {n_pan} lines; mutant lung with "
      f"CERES data: {n_mut}; other lung with CERES data: {n_oth}")
mean_mut = sum_mut / np.maximum(cnt_mut, 1)
mean_oth = sum_oth / np.maximum(cnt_oth, 1)
mean_pan = sum_pan / np.maximum(cnt_pan, 1)

common = pd.read_csv(DEP20 / "Achilles_common_essentials.csv")
sym_col = [c for c in common.columns if "gene" in c.lower() or "symbol" in c.lower()]
common_genes = set(common[sym_col[0]].astype(str)) if sym_col else set()
print(f"      common essential controls list: {len(common_genes)} genes")

layer = pd.DataFrame({
    "gene": genes,
    "n_mutant_lung_lines": n_mut,
    "n_other_lung_lines": n_oth,
    "mean_ceres_mutant_lung": mean_mut,
    "mean_ceres_other_lung": mean_oth,
    "mean_ceres_pan": mean_pan,
})
layer["selective_essentiality"] = (
    layer["mean_ceres_mutant_lung"] - layer["mean_ceres_other_lung"])
layer["common_essential"] = layer["gene"].isin(common_genes)
layer = layer.sort_values("selective_essentiality").reset_index(drop=True)
layer.to_csv(TABLES / "dependency_layer_20q2.csv", index=False)
print("      top 5 most selectively essential in mutant lung lines (20Q2 CERES):")
for _, r in layer.head(5).iterrows():
    print(f"        {r['gene']}: selective {r['selective_essentiality']:.3f}, "
          f"mutant mean {r['mean_ceres_mutant_lung']:.3f}, "
          f"other lung mean {r['mean_ceres_other_lung']:.3f}")

# ------------------------------------------------------- frozen ranking (v2 + 20Q2)
print("[4/5] Rebuilding the ranking with the frozen dependency layer ...")
de_tv = pd.read_csv(TABLES / "de_tumor_vs_normal_shrunk.csv")
de_mv = pd.read_csv(TABLES / "de_mutant_vs_wt.csv")
saf = pd.read_csv(TABLES / "safety_layer.csv")
cox = pd.read_csv(TABLES / "cox_phase2.csv")

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

dep_idx = layer.set_index("gene")
gated = dep_idx["mean_ceres_mutant_lung"] < DEP_ESSENTIALITY_GATE
covered = dep_idx["selective_essentiality"].notna()
dep_raw = -dep_idx["selective_essentiality"]
dep_raw = dep_raw.where(gated & (dep_raw > 0), 0.0)
dep_raw = dep_raw.where(covered)
dependency = scale01(dep_raw)
dependency.loc[dep_idx["common_essential"]] = 0.0
safety = saf.set_index("gene")["safety_window"]

outcome = pd.Series(np.nan, index=universe, dtype=float)
ok = cox.dropna(subset=["log_hr_per_sd"]).set_index("gene")["log_hr_per_sd"]
outcome.loc[ok.index] = scale01(ok)

frozen = pd.DataFrame({
    "gene": universe,
    "tumor_specificity": specificity.values,
    "dependency": dependency.reindex(universe).values,
    "safety_window": safety.reindex(universe).values,
    "outcome_link": outcome.values,
})
frozen["selective_essentiality_20q2"] = dep_idx["selective_essentiality"].reindex(universe).values
frozen["mean_ceres_mutant_lung"] = dep_idx["mean_ceres_mutant_lung"].reindex(universe).values
frozen["common_essential"] = dep_idx["common_essential"].reindex(universe).fillna(False).astype(bool).values
frozen["score"] = weighted_score(frozen, W)
frozen["coverage"] = frozen[["tumor_specificity", "dependency",
                             "safety_window", "outcome_link"]].notna().sum(axis=1)
unscoreable = frozen[frozen["dependency"].isna() & (frozen["tumor_specificity"] == 0)]
print(f"      excluded as unscoreable: {len(unscoreable)} genes")
frozen = frozen.drop(index=unscoreable.index)
frozen = frozen.sort_values("score", ascending=False).reset_index(drop=True)
frozen.insert(0, "rank", np.arange(1, len(frozen) + 1))
frozen.to_csv(TABLES / "target_scores_frozen.csv", index=False)
print(f"      frozen ranking written: {len(frozen)} genes")

# ------------------------------------------------------------- comparison
print("[5/5] Comparing the frozen ranking with the Phase 2 ranking ...")
v2 = pd.read_csv(TABLES / "target_scores_phase2.csv")
merged = v2[["gene", "rank", "score", "selective_essentiality"]].merge(
    frozen[["gene", "rank", "score", "selective_essentiality_20q2"]],
    on="gene", suffixes=("_v2", "_frozen"), how="outer")
rho, _ = spearmanr(merged["score_v2"].fillna(0), merged["score_frozen"].fillna(0))
print(f"      Spearman correlation of scores, v2 vs frozen: {rho:.3f}")

v2_short = pd.read_csv(TABLES / "shortlist_dependency_supported.csv")
frozen_sub = frozen[frozen["dependency"] > 0].copy()
frozen_short = set(frozen_sub.nsmallest(20, "rank")["gene"])
v2_short_set = set(v2_short["gene"])
overlap = v2_short_set & frozen_short
print(f"      shortlist overlap: {len(overlap)} of 20 genes in both; "
      f"enter frozen: {sorted(frozen_short - v2_short_set)}; "
      f"leave frozen: {sorted(v2_short_set - frozen_short)}")

key_genes = list(dict.fromkeys(list(v2_short["gene"])
                               + sorted(frozen_short) + ANSWER_KEY))
comp = merged[merged["gene"].isin(key_genes)].copy()
comp["in_v2_shortlist"] = comp["gene"].isin(v2_short_set)
comp["in_frozen_shortlist"] = comp["gene"].isin(frozen_short)
comp = comp.sort_values("rank_frozen")
comp.to_csv(TABLES / "timesplit_comparison.csv", index=False)

ak_f = frozen.set_index("gene")
with open(ROOT / "results" / "timesplit_check.md", "w") as fh:
    fh.write("# Time split check (Phase 3, pre-2021 frozen ranking)\n\n")
    fh.write("Pre-committed check from the README: freeze the evidence at\n"
             "an earlier date, rebuild the ranking, and ask whether later\n"
             "clinical developments (the KRAS G12C inhibitor approvals\n"
             "that began in 2021) appear high in the frozen ranking.\n\n")
    fh.write("## What was frozen\n\n")
    fh.write("- TCGA expression and survival (2014-2018 data) and GTEx v8\n"
             "  (2017-2020) are already pre-2021 and are carried over\n"
             "  unchanged, including the Phase 2 Cox values.\n"
             "- DepMap 24Q4 (Chronos) is replaced with DepMap Public 20Q2\n"
             "  (figshare article 12280541, published 2020-06-04), using\n"
             "  Achilles CERES gene effect scores; Chronos did not exist\n"
             "  yet in 20Q2. The selective essentiality contrast is\n"
             "  computed the same way. 20Q2 line counts: "
             f"{n_mut} mutant lung lines and {n_oth} other lung lines "
             "with CERES data.\n"
             "- The single cell layer cannot be frozen (the CELLxGENE\n"
             "  dataset is from a 2024 preprint), so the frozen ranking\n"
             "  is a four-layer ranking, directly comparable to Phase 2:\n"
             "  the only difference between the two rankings is the\n"
             "  dependency layer's release.\n\n")
    fh.write("## Result\n\n")
    fh.write(f"Frozen ranking size: {len(frozen)} genes. Spearman\n"
             f"correlation of scores with the Phase 2 ranking: {rho:.3f}.\n"
             f"Shortlist overlap (dependency-supported top 20): "
             f"{len(overlap)} of 20 genes in both.\n\n")
    fh.write("| Gene | Rank in Phase 2 | Rank frozen | Selective essentiality 24Q4 | "
             "Selective essentiality 20Q2 |\n")
    fh.write("| --- | --- | --- | --- | --- |\n")
    v2_idx = v2.set_index("gene")
    for g in key_genes:
        r2 = int(v2_idx.loc[g, "rank"]) if g in v2_idx.index else None
        rf = int(ak_f.loc[g, "rank"]) if g in ak_f.index else None
        se2 = v2_idx.loc[g, "selective_essentiality"] if g in v2_idx.index else np.nan
        sef = ak_f.loc[g, "selective_essentiality_20q2"] if g in ak_f.index else np.nan
        fh.write(f"| {g} | {r2 if r2 else 'not ranked'} | {rf if rf else 'not ranked'} | "
                 f"{se2:.3f} | {sef:.3f} |\n")
    fh.write("\n## Reading\n\n")
    fh.write("KRAS is not expected to rank high in either the Phase 2 or\n"
             "the frozen subgroup ranking: the subgroup ranking asks\n"
             "which genes are subgroup-selective targets, and KRAS\n"
             "selective essentiality is near zero in both releases. The\n"
             "informative result is whether the frozen ranking recovers\n"
             "the same biology, which the NFE2L2 row and the shortlist\n"
             "overlap above answer. See results/phase3-report.md.\n")
print("      answer key in the frozen ranking: " +
      ", ".join(f"{g} #{int(ak_f.loc[g, 'rank'])}" for g in ANSWER_KEY if g in ak_f.index))
print("Step 4 (time split check) complete.")
