"""Phase 2, step 1b: STK11/KEAP1 subgroup differential expression (tumors only).

Companion to 02_de_shrinkage_subgroup.py. The first run of script 02
completed Part A (shrinkage) and then crashed at the start of Part B
on a pandas set-indexer bug (fixed in script 02 since). Rather than
re-running the 16 minute tumor vs normal fit and shrinkage, this
script runs Part B alone with identical analysis code: it reloads the
counts, rebuilds the mutation strata from the same cBioPortal records,
and runs the mutant vs wild type DESeq2 fits within primary tumors.

  - combined mutant (STK11 and/or KEAP1) vs wild type (neither mutation)
  - STK11-only vs wild type, KEAP1-only vs wild type, both vs wild type
    (each only if the group has at least 20 tumors)

Fold changes are unshrunk here (see script 02's docstring for why),
and each table's maximum |log2FC| is printed so any residual artifact
is visible.

Outputs:
  results/tables/de_mutant_vs_wt.csv
  results/tables/de_stk11_vs_wt.csv / de_keap1_vs_wt.csv / de_both_vs_wt.csv

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/02b_subgroup_de.py
"""

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "xena"
CBIO = ROOT / "data" / "raw" / "cbioportal"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

MIN_GROUP = 20  # docs/project-spec.md: groups under 20 tumors are not analyzed alone


def load_counts():
    """Xena STAR counts, integer recovery, symbol mapping, prefilter.
    Identical to Phase 1 and to script 02."""
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


def run_deseq(counts_sub, metadata, factor, contrast):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(counts=counts_sub.T, metadata=metadata,
                       design_factors=factor, refit_cooks=True,
                       n_cpus=1, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast, quiet=True)
    stats.summary()
    de = stats.results_df.copy()
    de.insert(0, "gene", de.index)
    de = de.reset_index(drop=True)
    de = de[["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]]
    return de.sort_values("padj").reset_index(drop=True)


print("[1/3] Loading counts ...")
counts, samples, condition = load_counts()
tumor_samples = [s for s, c in zip(samples, condition) if c == "tumor"]

print("[2/3] Building STK11/KEAP1 mutation strata from cBioPortal records ...")
SILENT = {"Silent"}
mut_samples = {}
for symbol in ("STK11", "KEAP1"):
    records = json.loads((CBIO / f"mutations_{symbol}.json").read_text())
    mut_samples[symbol] = {r["patientId"] for r in records
                           if r.get("mutationType") not in SILENT}

def patient_of(barcode):
    return barcode.replace(".", "-")[:12]

groups = {}
for s in tumor_samples:
    p = patient_of(s)
    stk = p in mut_samples["STK11"]
    kea = p in mut_samples["KEAP1"]
    groups[s] = ("both" if stk and kea else
                 "stk11_only" if stk else
                 "keap1_only" if kea else "wt")

sample_list = set(json.loads((CBIO / "sample_list.json").read_text()))
covered = sorted(s for s in tumor_samples
                 if patient_of(s) in {x[:12] for x in sample_list})
sizes = Counter(groups[s] for s in covered)
print(f"      group sizes (covered tumors): {dict(sizes)}")

print("[3/3] Subgroup DE within primary tumors ...")
tumor_counts = counts.loc[:, covered]
tumor_counts = tumor_counts[(tumor_counts >= 10).sum(axis=1) >= 10]
print(f"      tumor-only matrix: {tumor_counts.shape[0]} genes x {tumor_counts.shape[1]} samples")

mut2 = pd.Series([groups[s] for s in tumor_counts.columns], index=tumor_counts.columns)
mut2 = np.where(mut2.isin(["stk11_only", "keap1_only", "both"]), "mutant", "wt")
meta2 = pd.DataFrame({"mutgroup": mut2}, index=tumor_counts.columns)
print(f"      combined mutant vs WT: {int((mut2 == 'mutant').sum())} mutant, "
      f"{int((mut2 == 'wt').sum())} WT")
de_mut = run_deseq(tumor_counts, meta2, "mutgroup", ["mutgroup", "mutant", "wt"])
de_mut.to_csv(TABLES / "de_mutant_vs_wt.csv", index=False)
print(f"      mutant vs WT done: {(de_mut['padj'] < 0.05).sum()} genes at padj < 0.05; "
      f"max |log2FC| {de_mut['log2FoldChange'].abs().max():.2f}")
top_mut = de_mut[de_mut["padj"] < 0.05].nlargest(5, "log2FoldChange")
print("      top 5 up in mutant tumors: " +
      ", ".join(f"{r['gene']} ({r['log2FoldChange']:.2f}, padj {r['padj']:.2g})"
                for _, r in top_mut.iterrows()))

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

print("Step 1b (subgroup DE) complete.")
