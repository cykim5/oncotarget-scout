"""Phase 2, step 3: GTEx safety layer (GTEx v8 median TPM).

For every gene, computes a safety window from median expression in the
vital tissues named in docs/target-score.md: heart, brain, liver.
  safety_window = 1 - max(scaled heart, scaled brain, scaled liver)
where each organ value is the maximum log2(TPM + 1) across that organ's
tissues, scaled 0 to 1 across all genes. A gene quiet in all three
organs scores near 1; a gene loud in any one of them scores near 0.
Using the maximum across tissues within an organ (rather than the mean)
is the conservative choice for safety: one loud tissue is enough to
flag a risk.

Outputs:
  results/tables/safety_layer.csv

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase2_python/04_gtex.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GTEX = ROOT / "data" / "raw" / "gtex"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

ORGANS = {
    "heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle"],
    "brain": None,  # filled below: all "Brain - *" tissues
    "liver": ["Liver"],
}


def scale01(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or hi == lo:
        return s * np.nan
    return (s - lo) / (hi - lo)


print("Loading GTEx v8 median TPM ...")
gct = pd.read_csv(
    GTEX / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    sep="\t", skiprows=2, index_col=0)
symbols = gct["Description"]
tpm = gct.drop(columns=["Description"])
print(f"      {tpm.shape[0]} genes x {tpm.shape[1]} tissues")

ORGANS["brain"] = [c for c in tpm.columns if c.startswith("Brain - ")]
print(f"      brain tissues used ({len(ORGANS['brain'])}): {ORGANS['brain']}")

tpm.index = symbols
tpm = tpm.groupby(level=0).max()  # duplicate symbols: keep the loudest (conservative)
print(f"      {tpm.shape[0]} unique symbols after aggregation")

layer = pd.DataFrame({"gene": tpm.index})
for organ, tissues in ORGANS.items():
    raw = tpm[tissues].max(axis=1)
    layer[f"{organ}_tpm_max"] = raw.values
    layer[f"{organ}_scaled"] = scale01(np.log2(raw + 1)).values
layer["safety_window"] = 1 - layer[["heart_scaled", "brain_scaled",
                                    "liver_scaled"]].max(axis=1)
layer = layer[["gene", "heart_tpm_max", "brain_tpm_max", "liver_tpm_max",
               "safety_window"]]
layer = layer.sort_values("safety_window", ascending=False).reset_index(drop=True)
layer.to_csv(TABLES / "safety_layer.csv", index=False)
print(f"      safety layer written: {len(layer)} genes; "
      f"median safety window {layer['safety_window'].median():.3f}")
print("Step 3 (GTEx safety layer) complete.")
