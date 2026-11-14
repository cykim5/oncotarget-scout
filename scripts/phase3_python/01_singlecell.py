"""Phase 3, step 1: single cell layer from CELLxGENE (malignant cell share).

Dataset: CELLxGENE 4ca422ae-e4ad-4884-9073-7c3799cd9143, 35,954 cells
from lung tumor samples (lung adenocarcinoma and squamous cell lung
carcinoma, 4 donors), Chromium 10x, from the collection
"Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx
DSP, and Chromium Methods Reveals Inter- and Intra-Patient
Heterogeneity" (bioRxiv 2024.11.01.621259). Chosen because its author
cell type labels include an explicit malignant cell compartment; the
other LUAD datasets found in the CELLxGENE scan either label tumor
cells only as lung epithelial cells or are multi-million cell atlases
too large for this environment. See results/run_log.md.

What this script does, in plain language: for every gene in the
Phase 2 ranking universe that the dataset measures, add up how much
of the gene's total expression sits in malignant cells, in immune
cells, in stromal cells, and elsewhere. The malignant share (the
fraction of the gene's total tumor expression found in malignant
cells) becomes the fifth evidence layer in script 02.

Design choices, all recorded here and in the run log:
  - Author labels are used as given; no re-clustering. The spec's
    risk mitigation says exactly this: choose a dataset with author
    cell type labels and treat the layer as a location check, not a
    discovery step. QC is still run and reported.
  - Counts, not log values, are summed: a share of log expression
    would be meaningless. The script checks whether X holds raw
    counts and, if it holds log1p values, recovers counts with
    expm1 before summing. Which case applied is printed and logged.
  - Cells are processed in chunks so the full matrix is never held
    in memory (this environment has about 6 GB).
  - LUAD cells are used where per-cell disease labels allow it;
    squamous cell carcinoma cells are excluded from the layer and
    the counts are reported. If no per-cell disease labels exist,
    all tumor cells are used and the mix is stated.
  - Compartment mapping is by keyword rules over the author labels
    and the full label to compartment table is printed, so the
    mapping is checkable. Non-malignant lung epithelial cells are
    their own "other" compartment, not stromal and not malignant.

Outputs:
  results/tables/singlecell_layer.csv
  results/tables/singlecell_qc.json
  results/tables/singlecell_compartments.csv (label to compartment map)

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/01_singlecell.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "data" / "raw" / "cellxgene" / "luad_sc_4ca422ae.h5ad"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

CHUNK = 4000  # cells per chunk (memory guard)
MIN_COUNTS = 500  # light QC: drop cells with fewer total counts
MIN_GENES = 200  # light QC: drop cells with fewer detected genes
MIN_GENE_COUNTS = 20  # a gene needs at least this many total counts in
# the dataset for its compartment shares to be defined; below that, a
# share of 1.0 could rest on a single read. Gated genes carry NaN.

# Keyword rules mapping author cell type labels to compartments.
# Order matters: malignant is checked first, then immune, then stromal.
IMMUNE_KEYS = ["t cell", "t-cell", "t-regulatory", "regulatory t",
               "regulatory", "b cell", "plasma cell", "natural killer",
               "nk cell", "monocyte", "macrophage", "dendritic", "mast cell",
               "granulocyte", "neutrophil", "lymphocyte", "myeloid",
               "cytotoxic t"]
STROMAL_KEYS = ["fibroblast", "endothelial", "pericyte", "smooth muscle"]


def compartment(label: str) -> str:
    low = label.lower()
    if "malignant" in low:
        return "malignant"
    if any(k in low for k in IMMUNE_KEYS):
        return "immune"
    if any(k in low for k in STROMAL_KEYS):
        return "stromal"
    return "other"


print("[1/5] Opening the h5ad in backed mode ...")
import anndata as ad

adata = ad.read_h5ad(H5AD, backed="r")
print(f"      shape: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"      obs columns: {list(adata.obs.columns)}")
print(f"      var columns: {list(adata.var.columns)}")
print(f"      X type: {type(adata.X).__name__}")

# Cell type column: CELLxGENE schema uses 'cell_type'; fall back to any
# column whose name contains 'cell_type' or 'author'.
obs = adata.obs
ct_col = None
for cand in ("cell_type", "author_cell_type", "cell_type_original"):
    if cand in obs.columns:
        ct_col = cand
        break
if ct_col is None:
    ct_col = [c for c in obs.columns if "cell_type" in c.lower()][0]
print(f"      cell type column: {ct_col}")
labels = obs[ct_col].astype(str)
print("      author labels:")
for lab, n in labels.value_counts().items():
    print(f"        {lab}: {n} cells -> {compartment(lab)}")

# Disease column, for the LUAD subset.
dis_col = None
for cand in ("disease", "disease_ontology_term_id", "donor_disease"):
    if cand in obs.columns:
        dis_col = cand
        break
if dis_col is not None:
    print(f"      disease column: {dis_col}")
    print("      disease values:")
    for val, n in obs[dis_col].astype(str).value_counts().items():
        print(f"        {val}: {n} cells")
else:
    print("      no per-cell disease column found")

# Gene symbols: CELLxGENE var index is Ensembl; symbols in feature_name.
var = adata.var
if "feature_name" in var.columns:
    symbols = var["feature_name"].astype(str)
else:
    symbols = pd.Series(var.index.astype(str), index=var.index)
print(f"      {symbols.nunique()} unique symbols across {len(symbols)} var rows")

# Is X raw counts or log1p? Sample the first chunk and inspect.
print("[2/5] Checking whether X holds raw counts ...")
sample = adata.X[:2000]
if sp.issparse(sample):
    sample = sample.toarray()
sample = np.asarray(sample, dtype=np.float64)
is_integer = np.all(sample == np.rint(sample))
print(f"      sample max {sample.max():.2f}, mean {sample.mean():.4f}, "
      f"all integer values: {is_integer}")
if is_integer and sample.max() > 50:
    counts_are_raw = True
    print("      X holds raw counts; summing directly.")
else:
    counts_are_raw = False
    print("      X does not look like raw counts; counts will be recovered "
          "with expm1 (log1p assumption) before summing.")
n_cells_total = int(adata.n_obs)

# ------------------------------------------------------------- QC and subset
print("[3/5] QC in chunks (total counts and genes detected per cell) ...")
tot_counts = np.zeros(n_cells_total, dtype=np.float64)
n_genes = np.zeros(n_cells_total, dtype=np.int64)
for start in range(0, n_cells_total, CHUNK):
    stop = min(start + CHUNK, n_cells_total)
    block = adata.X[start:stop]
    if sp.issparse(block):
        block = block.toarray()
    block = np.asarray(block, dtype=np.float64)
    if not counts_are_raw:
        block = np.expm1(block)
        block[block < 0] = 0
    tot_counts[start:stop] = block.sum(axis=1)
    n_genes[start:stop] = (block > 0).sum(axis=1)
    print(f"      ... {stop}/{n_cells_total} cells", flush=True)
keep_qc = (tot_counts >= MIN_COUNTS) & (n_genes >= MIN_GENES)
print(f"      QC: median total counts {np.median(tot_counts):.0f}, "
      f"median genes detected {int(np.median(n_genes))}")
print(f"      QC drops {(~keep_qc).sum()} of {n_cells_total} cells "
      f"(counts < {MIN_COUNTS} or genes < {MIN_GENES})")

# LUAD subset where per-cell disease labels allow it.
keep_disease = np.ones(n_cells_total, dtype=bool)
disease_note = "no per-cell disease labels; all tumor cells used (LUAD + LUSC mix)"
if dis_col is not None:
    dis_vals = obs[dis_col].astype(str).str.lower()
    luad_mask = dis_vals.str.contains("adenocarcinoma") & ~dis_vals.str.contains("squamous")
    if luad_mask.sum() > 0:
        keep_disease = luad_mask.to_numpy()
        disease_note = (f"restricted to LUAD cells by per-cell disease label: "
                        f"{int(luad_mask.sum())} of {n_cells_total} cells")
print(f"      {disease_note}")

keep = keep_qc & keep_disease
comp = np.array([compartment(lab) for lab in labels])
comp = comp[keep]
print(f"      cells kept: {int(keep.sum())} of {n_cells_total}")
for c in ("malignant", "immune", "stromal", "other"):
    print(f"        {c}: {int((comp == c).sum())} cells")
donor_note = "no donor column in obs"
if "donor_id" in obs.columns:
    donors_kept = obs.loc[keep, "donor_id"].astype(str).value_counts()
    donor_note = (f"LUAD cells used come from {len(donors_kept)} donor(s): "
                  + ", ".join(f"{d} ({n})" for d, n in donors_kept.items()))
    print(f"      {donor_note}")

# ------------------------------------------------- compartment sums per gene
print("[4/5] Summing expression per compartment, in chunks ...")
kept_idx = np.where(keep)[0]
compartments = ["malignant", "immune", "stromal", "other"]
sums = {c: np.zeros(adata.n_vars, dtype=np.float64) for c in compartments}
n_expr = {c: np.zeros(adata.n_vars, dtype=np.int64) for c in compartments}
n_cells_comp = {c: int((comp == c).sum()) for c in compartments}
for s in range(0, len(kept_idx), CHUNK):
    rows = kept_idx[s:s + CHUNK]
    block = adata.X[rows]
    if sp.issparse(block):
        block = block.toarray()
    block = np.asarray(block, dtype=np.float64)
    if not counts_are_raw:
        block = np.expm1(block)
        block[block < 0] = 0
    comp_block = comp[s:s + CHUNK]
    for c in compartments:
        m = comp_block == c
        if m.sum() == 0:
            continue
        sums[c] += block[m].sum(axis=0)
        n_expr[c] += (block[m] > 0).sum(axis=0)
    print(f"      ... {min(s + CHUNK, len(kept_idx))}/{len(kept_idx)} kept cells",
          flush=True)

# Aggregate duplicate symbols (sum), then restrict to the ranking universe.
frame = pd.DataFrame({
    "gene": symbols.values,
    "counts_malignant": sums["malignant"],
    "counts_immune": sums["immune"],
    "counts_stromal": sums["stromal"],
    "counts_other": sums["other"],
    "n_expr_malignant": n_expr["malignant"],
    "n_expr_immune": n_expr["immune"],
    "n_expr_stromal": n_expr["stromal"],
    "n_expr_other": n_expr["other"],
})
frame = frame.groupby("gene", as_index=False).sum(numeric_only=True)
frame["counts_total"] = (frame["counts_malignant"] + frame["counts_immune"]
                         + frame["counts_stromal"] + frame["counts_other"])
gated = frame["counts_total"] < MIN_GENE_COUNTS
frame["malignant_share"] = frame["counts_malignant"] / frame["counts_total"].replace(0, np.nan)
frame["immune_share"] = frame["counts_immune"] / frame["counts_total"].replace(0, np.nan)
frame["stromal_share"] = frame["counts_stromal"] / frame["counts_total"].replace(0, np.nan)
frame["other_share"] = frame["counts_other"] / frame["counts_total"].replace(0, np.nan)
frame.loc[gated, ["malignant_share", "immune_share", "stromal_share",
                  "other_share"]] = np.nan
frame["pct_malignant_expressing"] = 100 * frame["n_expr_malignant"] / max(n_cells_comp["malignant"], 1)
frame["pct_immune_expressing"] = 100 * frame["n_expr_immune"] / max(n_cells_comp["immune"], 1)
frame["pct_stromal_expressing"] = 100 * frame["n_expr_stromal"] / max(n_cells_comp["stromal"], 1)

universe = pd.read_csv(TABLES / "target_scores_phase2.csv", usecols=["gene"])
layer = frame[frame["gene"].isin(set(universe["gene"]))].copy()
layer = layer.sort_values("malignant_share", ascending=False).reset_index(drop=True)
layer.to_csv(TABLES / "singlecell_layer.csv", index=False)
n_covered = int(layer["malignant_share"].notna().sum())
print(f"      layer written: {len(layer)} genes in the ranking universe, "
      f"{n_covered} with a defined malignant share "
      f"({len(layer) - n_covered} not detected in the dataset)")

qc = {
    "dataset_id": "4ca422ae-e4ad-4884-9073-7c3799cd9143",
    "dataset_cells": n_cells_total,
    "cells_after_qc": int(keep_qc.sum()),
    "cells_used": int(keep.sum()),
    "disease_note": disease_note,
    "donor_note": donor_note,
    "counts_are_raw": bool(counts_are_raw),
    "cell_type_column": ct_col,
    "qc_thresholds": {"min_total_counts": MIN_COUNTS, "min_genes_detected": MIN_GENES},
    "min_gene_total_counts_for_share": MIN_GENE_COUNTS,
    "genes_gated_low_counts": int(gated.sum()),
    "median_total_counts": float(np.median(tot_counts)),
    "median_genes_detected": int(np.median(n_genes)),
    "n_cells_per_compartment": n_cells_comp,
    "genes_in_universe": int(len(layer)),
    "genes_with_malignant_share": n_covered,
}
(TABLES / "singlecell_qc.json").write_text(json.dumps(qc, indent=2))

label_map = (pd.DataFrame({"author_label": labels})
             .drop_duplicates()
             .assign(compartment=lambda d: d["author_label"].map(compartment),
                     n_cells=lambda d: d["author_label"].map(labels.value_counts())))
label_map.to_csv(TABLES / "singlecell_compartments.csv", index=False)

print("[5/5] Sanity checks on the shortlist and answer key ...")
check_genes = ["INTS6", "NFE2L2", "SPRR1B", "UXS1", "NCAPG2", "EGFR", "KRAS", "CD274"]
sub = layer.set_index("gene")
for g in check_genes:
    if g in sub.index:
        r = sub.loc[g]
        print(f"      {g}: malignant share {r['malignant_share']:.3f}, "
              f"immune {r['immune_share']:.3f}, stromal {r['stromal_share']:.3f}, "
              f"pct malignant cells expressing {r['pct_malignant_expressing']:.1f}%")
    else:
        print(f"      {g}: not in the single cell layer")
print("Step 1 (single cell layer) complete.")
