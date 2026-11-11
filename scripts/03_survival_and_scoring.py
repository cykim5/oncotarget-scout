"""03_survival_and_scoring.py

Survival layer + score assembly for OncoTarget Scout.

Inputs (produced by the R steps and the documented manual downloads):
  results/tables/de_tumor_vs_normal.csv
  results/tables/de_stk11_keap1_mut_vs_wt.csv
  data/processed/tcga_luad_clinical.rds  -> export to CSV from R first
      (the R step that does this will be added on the first full run)
  data/raw/depmap/CRISPRGeneEffect.csv   (manual download, see data/README.md)
  data/raw/depmap/Model.csv
  data/raw/gtex_median_tpm.csv           (manual download, see data/README.md)

Output:
  results/tables/target_scores.csv  (one row per gene, layer values + score
                                     + coverage + known-target flag)
  results/answer_key_check.md       (rank positions of EGFR, KRAS, CD274)

Scoring follows docs/target-score.md exactly. Weights live in WEIGHTS below
and nowhere else, so the sensitivity analysis can rebuild the ranking under
alternate weights by changing one dictionary.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"

WEIGHTS = {
    "dependency": 0.30,
    "tumor_specificity": 0.25,
    "safety_window": 0.20,
    "malignant_share": 0.15,
    "outcome_link": 0.10,
}

ANSWER_KEY = ["EGFR", "KRAS", "CD274"]  # established LUAD target biology
VITAL_TISSUES = ["Heart", "Brain", "Liver"]  # GTEx safety window tissues


def scale01(series: pd.Series) -> pd.Series:
    """Min-max scale to 0..1, NaN preserved (missing layers are never zero)."""
    lo, hi = series.min(), series.max()
    if pd.isna(lo) or hi == lo:
        return series * np.nan
    return (series - lo) / (hi - lo)


def main() -> None:
    de_tumor = pd.read_csv(TABLES / "de_tumor_vs_normal.csv")
    de_sub = pd.read_csv(TABLES / "de_stk11_keap1_mut_vs_wt.csv")

    df = de_tumor[["gene", "log2FoldChange", "padj"]].rename(
        columns={"log2FoldChange": "lfc_tumor", "padj": "padj_tumor"}
    ).merge(
        de_sub[["gene", "log2FoldChange", "padj"]].rename(
            columns={"log2FoldChange": "lfc_subgroup", "padj": "padj_subgroup"}
        ),
        on="gene",
        how="outer",
    )

    # Tumor specificity layer: subgroup enrichment counts double, by design
    # (see docs/target-score.md). P values gate, they do not rank.
    df["tumor_specificity"] = scale01(
        df["lfc_tumor"].clip(lower=0) + 2 * df["lfc_subgroup"].clip(lower=0)
    )

    # TODO on first data run: join DepMap dependency contrast
    # (mutant lung lines vs all lines, Chronos), GTEx safety window,
    # single cell malignant share, and the Cox outcome layer.
    # Each join adds one column named exactly as in WEIGHTS, already scaled
    # 0..1, with NaN where the layer has no coverage for that gene.
    for layer in WEIGHTS:
        if layer not in df.columns:
            df[layer] = np.nan

    layer_cols = list(WEIGHTS)
    df["coverage"] = df[layer_cols].notna().sum(axis=1)

    # Weighted score over available layers, renormalized by covered weight,
    # so a gene is never punished for a missing layer, only compared fairly.
    w = pd.Series(WEIGHTS)
    covered_weight = df[layer_cols].notna().mul(w, axis=1).sum(axis=1)
    df["score"] = (
        df[layer_cols].fillna(0).mul(w, axis=1).sum(axis=1) / covered_weight
    )
    df["known_target"] = df["gene"].isin(ANSWER_KEY)

    out = df.sort_values("score", ascending=False)
    out.to_csv(TABLES / "target_scores.csv", index=False)

    lines = ["# Answer key check", ""]
    for gene in ANSWER_KEY:
        row = out.reset_index(drop=True)
        pos = row.index[row["gene"] == gene]
        rank = int(pos[0]) + 1 if len(pos) else None
        lines.append(f"- {gene}: rank {rank} of {len(out)} scored genes")
    (ROOT / "results" / "answer_key_check.md").write_text("\n".join(lines))

    print(f"Scored {len(out)} genes. Top 10:")
    print(out[["gene", "score", "coverage", "known_target"]].head(10))


if __name__ == "__main__":
    main()
