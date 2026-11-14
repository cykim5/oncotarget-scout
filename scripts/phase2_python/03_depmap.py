"""Phase 2, step 2: DepMap CRISPR dependency layer (DepMap Public 24Q4).

For every gene with Chronos scores, computes selective essentiality in
STK11/KEAP1-mutant lung cancer cell lines:
  selective_essentiality = mean Chronos in mutant lung lines
                         - mean Chronos in other lung lines
More negative means the mutant lines depend on the gene more than other
lung lines do. The pan-line mean is reported as a control. Genes in the
Achilles common essential controls list are flagged; the scoring script
penalizes them (dependency layer set to zero) because a gene every cell
needs offers no therapeutic window.

Mutation status per cell line comes from OmicsSomaticMutations.csv,
read in chunks and filtered to non-silent STK11 / KEAP1 / KRAS calls
only, so the 339 MB file is never loaded whole. KRAS status is recorded
for context (many STK11/KEAP1-mutant lines also carry KRAS) but does not
define the mutant group.

Line counts behind every claim are reported in the output table and in
stdout; a dependency claim with only a handful of mutant lines behind
it is a weak claim, and the report says so.

Outputs:
  results/tables/dependency_layer.csv
  results/tables/depmap_line_summary.json

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/03_depmap.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEPMAP = ROOT / "data" / "raw" / "depmap24q4"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

GENES_OF_INTEREST = {"STK11", "KEAP1", "KRAS"}
SILENT_TYPES = {"Silent", "Intron", "IGR", "5'UTR", "3'UTR", "5'Flank", "3'Flank"}

print("[1/4] Loading Model.csv (cell line annotations) ...")
model = pd.read_csv(DEPMAP / "Model.csv", usecols=["ModelID", "CellLineName",
                                                   "OncotreeLineage",
                                                   "OncotreePrimaryDisease"])
lung = model[model["OncotreeLineage"] == "Lung"]
print(f"      {len(model)} models total, {len(lung)} lung lineage")

print("[2/4] Reading OmicsSomaticMutations.csv in chunks (STK11/KEAP1/KRAS only) ...")
mut_lines = {g: set() for g in GENES_OF_INTEREST}
n_rows = 0
for chunk in pd.read_csv(DEPMAP / "OmicsSomaticMutations.csv",
                         usecols=["ModelID", "HugoSymbol", "VariantType"],
                         chunksize=200_000, low_memory=False):
    n_rows += len(chunk)
    sub = chunk[chunk["HugoSymbol"].isin(GENES_OF_INTEREST)
                & ~chunk["VariantType"].isin(SILENT_TYPES)]
    for gene, grp in sub.groupby("HugoSymbol"):
        mut_lines[gene].update(grp["ModelID"].tolist())
print(f"      scanned {n_rows} mutation rows")
for g in sorted(GENES_OF_INTEREST):
    print(f"      {g}: {len(mut_lines[g])} models with non-silent mutations")

lung_ids = set(lung["ModelID"])
mutant_lung = sorted((mut_lines["STK11"] | mut_lines["KEAP1"]) & lung_ids)
other_lung = sorted(lung_ids - set(mutant_lung))
n_both = len(mut_lines["STK11"] & mut_lines["KEAP1"] & lung_ids)
print(f"      mutant lung lines (STK11 and/or KEAP1): {len(mutant_lung)} "
      f"({n_both} with both); other lung lines: {len(other_lung)}")
(TABLES / "depmap_line_summary.json").write_text(json.dumps({
    "depmap_release": "DepMap Public 24Q4",
    "n_models_total": int(len(model)),
    "n_lung_models": int(len(lung)),
    "n_mutant_lung_lines": len(mutant_lung),
    "n_both_mutant_lung_lines": int(n_both),
    "n_other_lung_lines": len(other_lung),
    "n_stk11_mutant_models": len(mut_lines["STK11"]),
    "n_keap1_mutant_models": len(mut_lines["KEAP1"]),
    "n_kras_mutant_models": len(mut_lines["KRAS"]),
    "mutant_lung_model_ids": mutant_lung,
}, indent=2))

print("[3/4] Loading CRISPRGeneEffect.csv (Chronos scores) ...")
effect = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", index_col=0)
print(f"      gene effect matrix: {effect.shape[0]} models x {effect.shape[1]} genes")
# Columns look like "TP53 (7157)"; strip to the symbol.
effect.columns = effect.columns.str.split(" ").str[0]
effect = effect.loc[:, ~effect.columns.duplicated(keep="first")]

mut_idx = [m for m in mutant_lung if m in effect.index]
oth_idx = [m for m in other_lung if m in effect.index]
print(f"      mutant lung lines with CRISPR data: {len(mut_idx)}; "
      f"other lung lines with CRISPR data: {len(oth_idx)}")
mean_mut = effect.loc[mut_idx].mean(axis=0)
mean_oth = effect.loc[oth_idx].mean(axis=0)
mean_pan = effect.mean(axis=0)

print("[4/4] Assembling the dependency layer ...")
common = pd.read_csv(DEPMAP / "AchillesCommonEssentialControls.csv")
sym_col = [c for c in common.columns if "gene" in c.lower() or "symbol" in c.lower()]
common_genes = set(common[sym_col[0]].astype(str)) if sym_col else set()
print(f"      common essential controls list: {len(common_genes)} genes")

layer = pd.DataFrame({
    "gene": effect.columns,
    "n_mutant_lung_lines": len(mut_idx),
    "n_other_lung_lines": len(oth_idx),
    "mean_chronos_mutant_lung": mean_mut.values,
    "mean_chronos_other_lung": mean_oth.values,
    "mean_chronos_pan": mean_pan.values,
})
layer["selective_essentiality"] = (
    layer["mean_chronos_mutant_lung"] - layer["mean_chronos_other_lung"])
layer["common_essential"] = layer["gene"].isin(common_genes)
layer = layer.sort_values("selective_essentiality").reset_index(drop=True)
layer.to_csv(TABLES / "dependency_layer.csv", index=False)
print("      top 5 most selectively essential in mutant lung lines:")
for _, r in layer.head(5).iterrows():
    print(f"        {r['gene']}: selective {r['selective_essentiality']:.3f}, "
          f"mutant mean {r['mean_chronos_mutant_lung']:.3f}, "
          f"other lung mean {r['mean_chronos_other_lung']:.3f}")
print("Step 2 (DepMap dependency layer) complete.")
