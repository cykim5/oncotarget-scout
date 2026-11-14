"""Phase 2, step 1: shrinkage fix and STK11/KEAP1 subgroup differential expression.

Part A (shrinkage fix). Re-runs the Phase 1 tumor vs normal DESeq2 fit and
applies apeglm-style log2 fold change shrinkage (pydeseq2 DeseqStats
lfc_shrink, Cauchy prior, p values unchanged). The unshrunken Phase 1
fold changes let MAGEC1 reach log2FC 43.3 off a near-zero normal baseline;
shrinkage pulls low-information fold changes toward zero. The Phase 1
partial ranking is then rebuilt with the shrunk statistic
(results/tables/target_scores_phase1_shrunk.csv), reusing the Phase 1 Cox
table, and MAGEC1's new rank is reported. The original Phase 1 table and
report are untouched; the change is recorded as a dated correction note
in results/phase1-report.md.

Part B (subgroup DE). Splits the 528 primary tumors by STK11/KEAP1
mutation status from the cBioPortal mutation records (non-silent
mutations only; silent mutations do not count). Runs, within primary
tumors only:
  - combined mutant (STK11 and/or KEAP1) vs wild type (neither mutation)
  - STK11-only vs wild type, KEAP1-only vs wild type, both vs wild type
    (each only if the group has at least 20 tumors; the rule is stated
    in docs/project-spec.md)
Mutant vs wild type fold changes are left unshrunk: both groups are
tumors, so the near-zero baseline artifact that motivated shrinkage in
the tumor vs normal comparison does not apply the same way, and the
script reports the maximum |log2FC| of each subgroup table so any
residual artifact is visible. This choice is logged in run_log.md.

Outputs:
  results/tables/de_tumor_vs_normal_shrunk.csv
  results/tables/de_mutant_vs_wt.csv
  results/tables/de_stk11_vs_wt.csv / de_keap1_vs_wt.csv / de_both_vs_wt.csv
  results/tables/target_scores_phase1_shrunk.csv
  results/tables/mutation_groups.json

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/02_de_shrinkage_subgroup.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "xena"
CBIO = ROOT / "data" / "raw" / "cbioportal"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]
W_SPEC, W_OUT = 0.25, 0.10
MIN_GROUP = 20  # docs/project-spec.md: groups under 20 tumors are not analyzed alone


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


def load_counts():
    """Xena STAR counts, integer recovery, symbol mapping, prefilter.
    Identical to Phase 1 so the gene universe is comparable."""
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
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples "
          f"({int((condition == 'tumor').sum())} tumors, {int((condition == 'normal').sum())} normals)")
    return counts, samples, condition


def run_deseq(counts_sub, metadata, factor, contrast, shrink=False):
    """One DESeq2 fit + contrast; optionally with apeglm LFC shrinkage."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(counts=counts_sub.T, metadata=metadata,
                       design_factors=factor, refit_cooks=True,
                       n_cpus=1, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast, quiet=True)
    stats.summary()
    if shrink:
        coeff = [c for c in stats.LFC.columns if "[T." in c]
        assert len(coeff) == 1, f"expected one non-intercept coeff, got {coeff}"
        print(f"      shrinking LFCs for coeff {coeff[0]} ...")
        stats.lfc_shrink(coeff=coeff[0])
        stats.summary()
    de = stats.results_df.copy()
    de.insert(0, "gene", de.index)
    de = de.reset_index(drop=True)
    de = de[["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
    return de.sort_values("padj").reset_index(drop=True)


# ------------------------------------------------------------------ load
print("[1/4] Loading counts ...")
counts, samples, condition = load_counts()
tumor_samples = [s for s, c in zip(samples, condition) if c == "tumor"]

# ------------------------------------------------------- mutation strata
print("[2/4] Building STK11/KEAP1 mutation strata from cBioPortal records ...")
SILENT = {"Silent"}
mut_samples = {}
for symbol in ("STK11", "KEAP1"):
    records = json.loads((CBIO / f"mutations_{symbol}.json").read_text())
    pats = {r["patientId"] for r in records if r.get("mutationType") not in SILENT}
    mut_samples[symbol] = pats
    print(f"      {symbol}: {len(records)} records, {len(pats)} patients with "
          f"non-silent mutations")

def patient_of(sample_barcode):
    return sample_barcode.replace(".", "-")[:12]

groups = {}
for s in tumor_samples:
    p = patient_of(s)
    stk = p in mut_samples["STK11"]
    kea = p in mut_samples["KEAP1"]
    if stk and kea:
        groups[s] = "both"
    elif stk:
        groups[s] = "stk11_only"
    elif kea:
        groups[s] = "keap1_only"
    else:
        groups[s] = "wt"

sample_list = set(json.loads((CBIO / "sample_list.json").read_text()))
covered = {s for s in tumor_samples
           if patient_of(s) in {x[:12] for x in sample_list}}
unknown = [s for s in tumor_samples if s not in covered]
print(f"      tumors with cBioPortal mutation coverage: {len(covered)} of {len(tumor_samples)}; "
      f"no mutation data: {len(unknown)} (excluded from subgroup DE)")
from collections import Counter
sizes = Counter(groups[s] for s in covered)
print(f"      group sizes (covered tumors): {dict(sizes)}")
(TABLES / "mutation_groups.json").write_text(json.dumps({
    "n_primary_tumors": len(tumor_samples),
    "n_with_mutation_data": len(covered),
    "n_without_mutation_data": len(unknown),
    "group_sizes": dict(sizes),
    "rule": "non-silent mutations only; groups under 20 tumors are not analyzed alone",
}, indent=2))

# ------------------------------------------- Part A: shrinkage tumor vs normal
print("[3/4] Tumor vs normal DESeq2 with LFC shrinkage ...")
metadata = pd.DataFrame({"condition": condition}, index=samples)
de_shrunk = run_deseq(counts, metadata, "condition",
                      ["condition", "tumor", "normal"], shrink=True)
de_shrunk.to_csv(TABLES / "de_tumor_vs_normal_shrunk.csv", index=False)
mag = de_shrunk[de_shrunk["gene"] == "MAGEC1"]
print(f"      shrunk DE done: {(de_shrunk['padj'] < 0.05).sum()} genes at padj < 0.05; "
      f"max |log2FC| now {de_shrunk['log2FoldChange'].abs().max():.2f}")
if len(mag):
    print(f"      MAGEC1 shrunk log2FC: {mag.iloc[0]['log2FoldChange']:.2f} "
          f"(was 43.28 unshrunk), padj {mag.iloc[0]['padj']:.3g}")

# Rebuild the Phase 1 partial ranking with the shrunk statistic
cox = pd.read_csv(TABLES / "cox_top_candidates.csv")
elig = de_shrunk[de_shrunk["padj"] < 0.05]
spec = pd.Series(0.0, index=de_shrunk["gene"].values)
spec.loc[elig["gene"].values] = scale01(elig["log2FoldChange"].clip(lower=0)).values
outcome = pd.Series(np.nan, index=de_shrunk["gene"], dtype=float)
ok = cox.dropna(subset=["log_hr_per_sd"]).set_index("gene")["log_hr_per_sd"]
outcome.loc[ok.index] = scale01(ok)
scores = pd.DataFrame({
    "gene": de_shrunk["gene"],
    "log2FoldChange_shrunk": de_shrunk["log2FoldChange"],
    "padj": de_shrunk["padj"],
    "tumor_specificity": spec.reindex(de_shrunk["gene"]).values,
    "outcome_link_cox_log_hr": ok.reindex(de_shrunk["gene"]).values,
    "outcome_link": outcome.reindex(de_shrunk["gene"]).values,
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
scores.to_csv(TABLES / "target_scores_phase1_shrunk.csv", index=False)
mag_rank = scores[scores["gene"] == "MAGEC1"]
print(f"      Phase 1 ranking rebuilt with shrunk LFCs: new rank 1 is "
      f"{scores.iloc[0]['gene']} (score {scores.iloc[0]['partial_score']:.3f})")
if len(mag_rank):
    print(f"      MAGEC1 new rank: {int(mag_rank.iloc[0]['rank'])} of {len(scores)} "
          f"(was rank 1)")
for g in ANSWER_KEY:
    row = scores[scores["gene"] == g]
    if len(row):
        print(f"      {g} shrunk-rank: {int(row.iloc[0]['rank'])}")

# ------------------------------------------- Part B: subgroup DE (tumors only)
print("[4/4] Subgroup DE within primary tumors ...")
tumor_counts = counts.loc[:, sorted(covered)]
tumor_counts = tumor_counts[(tumor_counts >= 10).sum(axis=1) >= 10]
print(f"      tumor-only matrix: {tumor_counts.shape[0]} genes x {tumor_counts.shape[1]} samples")

mut2 = pd.Series([groups[s] for s in tumor_counts.columns], index=tumor_counts.columns)
mut2 = np.where(mut2.isin(["stk11_only", "keap1_only", "both"]), "mutant", "wt")
meta2 = pd.DataFrame({"mutgroup": mut2}, index=tumor_counts.columns)
n_mut = int((mut2 == "mutant").sum())
n_wt = int((mut2 == "wt").sum())
print(f"      combined mutant vs WT: {n_mut} mutant, {n_wt} WT")
de_mut = run_deseq(tumor_counts, meta2, "mutgroup", ["mutgroup", "mutant", "wt"])
de_mut.to_csv(TABLES / "de_mutant_vs_wt.csv", index=False)
print(f"      mutant vs WT done: {(de_mut['padj'] < 0.05).sum()} genes at padj < 0.05; "
      f"max |log2FC| {de_mut['log2FoldChange'].abs().max():.2f}")
top_mut = de_mut[de_mut["padj"] < 0.05].nlargest(5, "log2FoldChange")
print("      top 5 up in mutant tumors: " +
      ", ".join(f"{r['gene']} ({r['log2FoldChange']:.2f}, padj {r['padj']:.2g})"
                for _, r in top_mut.iterrows()))

# Per-gene contrasts from one 4-level fit (only for groups with >= 20 tumors)
meta4 = pd.DataFrame({"mutgroup4": [groups[s] for s in tumor_counts.columns]},
                     index=tumor_counts.columns)
levels = {"stk11_only": "de_stk11_vs_wt.csv",
          "keap1_only": "de_keap1_vs_wt.csv",
          "both": "de_both_vs_wt.csv"}
run_levels = {k: v for k, v in levels.items() if sizes.get(k, 0) >= MIN_GROUP}
for k in levels:
    if k not in run_levels:
        print(f"      {k}: n={sizes.get(k, 0)} < {MIN_GROUP}, skipped (per spec rule)")
if run_levels:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds4 = DeseqDataSet(counts=tumor_counts.T, metadata=meta4,
                        design_factors="mutgroup4", refit_cooks=True,
                        n_cpus=1, quiet=True)
    dds4.deseq2()
    for level, fname in run_levels.items():
        stats = DeseqStats(dds4, contrast=["mutgroup4", level, "wt"], quiet=True)
        stats.summary()
        de = stats.results_df.copy()
        de.insert(0, "gene", de.index)
        de = de.reset_index(drop=True)
        de = de[["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
        de = de.sort_values("padj").reset_index(drop=True)
        de.to_csv(TABLES / fname, index=False)
        print(f"      {level} vs WT done (n={sizes[level]}): "
              f"{(de['padj'] < 0.05).sum()} genes at padj < 0.05; "
              f"max |log2FC| {de['log2FoldChange'].abs().max():.2f}")

print("Step 1 (shrinkage + subgroup DE) complete.")
