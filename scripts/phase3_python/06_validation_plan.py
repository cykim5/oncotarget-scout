"""Phase 3, step 6: the validation experiment plan.

Picks the cell lines and writes docs/validation-plan.md for the top
5 genes of the v3 dependency-supported shortlist, plus INTS6, NFE2L2,
and UXS1 if any of them fell out of the top 5 (the union is the plan
set). Every line named here comes from the DepMap Public 24Q4 data
already in data/raw/depmap24q4: mutation status from
OmicsSomaticMutations.csv (non-silent calls only) and per-line
Chronos values from CRISPRGeneEffect.csv.

Line selection rules, fixed in advance:
  - Test lines: lung lines carrying STK11 and/or KEAP1 non-silent
    mutations, with CRISPR data. Ranked with both-mutant lines
    first, then by most negative Chronos for that gene. Up to 4
    lines per gene, preferring at least one both-mutant, one
    STK11-only, and one KEAP1-only line where available. Lines with
    Chronos below -0.5 for the gene are preferred (the same
    essentiality gate the ranking uses); if fewer than 4 pass, the
    most negative remaining lines fill the slots and are flagged.
  - Control lines: lung lines wild type for both STK11 and KEAP1,
    with CRISPR data, where the gene's Chronos is above -0.5 (not
    essential). Two to three per gene, least essential first.
  - KRAS status is recorded for context (many STK11/KEAP1-mutant
    lines also carry KRAS) but does not define any group.

The plan itself is a fixed prose template filled with the selected
lines and the genes' layer values from the result tables, so every
number in it traces to a committed table. The perturbation,
readouts, and confirm/kill criteria are pre-specified plan
decisions, stated as such, not results.

Outputs:
  results/tables/validation_lines.csv
  docs/validation-plan.md

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/06_validation_plan.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEPMAP = ROOT / "data" / "raw" / "depmap24q4"
TABLES = ROOT / "results" / "tables"
DOCS = ROOT / "docs"

GENES_OF_INTEREST = {"STK11", "KEAP1", "KRAS"}
SILENT_TYPES = {"Silent", "Intron", "IGR", "5'UTR", "3'UTR", "5'Flank", "3'Flank"}
N_TEST = 4
N_CONTROL = 3
GATE = -0.5

print("[1/4] Loading the v3 shortlist and the target set ...")
short = pd.read_csv(TABLES / "shortlist_dependency_supported_phase3.csv")
plan_genes = list(short.head(5)["gene"])
for g in ("INTS6", "NFE2L2", "UXS1"):
    if g not in plan_genes:
        plan_genes.append(g)
print(f"      plan genes: {plan_genes}")
short_idx = short.set_index("gene")

print("[2/4] Mutation status per lung line (24Q4) ...")
model = pd.read_csv(DEPMAP / "Model.csv", usecols=["ModelID", "CellLineName",
                                                   "OncotreeLineage"])
lung = model[model["OncotreeLineage"] == "Lung"]
lung_ids = set(lung["ModelID"])
mut_lines = {g: set() for g in GENES_OF_INTEREST}
for chunk in pd.read_csv(DEPMAP / "OmicsSomaticMutations.csv",
                         usecols=["ModelID", "HugoSymbol", "VariantType"],
                         chunksize=200_000, low_memory=False):
    sub = chunk[chunk["HugoSymbol"].isin(GENES_OF_INTEREST)
                & ~chunk["VariantType"].isin(SILENT_TYPES)]
    for gene, grp in sub.groupby("HugoSymbol"):
        mut_lines[gene].update(grp["ModelID"].tolist())
mutant_lung = sorted((mut_lines["STK11"] | mut_lines["KEAP1"]) & lung_ids)
other_lung = sorted(lung_ids - set(mutant_lung))
print(f"      {len(lung)} lung lines; {len(mutant_lung)} mutant, "
      f"{len(other_lung)} wild type")

print("[3/4] Per-line Chronos values for the plan genes ...")
header = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv", nrows=0)
cols = {"DepMap_ID": "DepMap_ID"}
first_col = header.columns[0]
gene_cols = {}
for g in plan_genes:
    matches = [c for c in header.columns if c.split(" ")[0] == g]
    if not matches:
        raise SystemExit(f"gene {g} not found in CRISPRGeneEffect.csv")
    gene_cols[g] = matches[0]
effect = pd.read_csv(DEPMAP / "CRISPRGeneEffect.csv",
                     usecols=[first_col] + list(gene_cols.values()))
effect = effect.rename(columns={first_col: "ModelID",
                                **{v: k for k, v in gene_cols.items()}})
effect = effect.set_index("ModelID")
print(f"      gene effect subset: {effect.shape[0]} lines x {effect.shape[1]} genes")

names = lung.set_index("ModelID")["CellLineName"].to_dict()
rows = []
for gene in plan_genes:
    vals = effect[gene].dropna()
    test_pool = [m for m in mutant_lung if m in vals.index]
    ctrl_pool = [m for m in other_lung if m in vals.index]
    both = [m for m in test_pool
            if m in mut_lines["STK11"] and m in mut_lines["KEAP1"]]
    stk_only = [m for m in test_pool
                if m in mut_lines["STK11"] and m not in mut_lines["KEAP1"]]
    kea_only = [m for m in test_pool
                if m in mut_lines["KEAP1"] and m not in mut_lines["STK11"]]
    # Rank each pool by most negative Chronos, then take in the order
    # both, STK11-only, KEAP1-only, round-robin, to cover the strata.
    for pool in (both, stk_only, kea_only):
        pool.sort(key=lambda m: vals[m])
    ordered = []
    pools = [both, stk_only, kea_only]
    while len(ordered) < N_TEST and any(pools):
        for pool in pools:
            if pool and len(ordered) < N_TEST:
                ordered.append(pool.pop(0))
    # Prefer lines passing the essentiality gate; if the ordered set
    # has lines above the gate while gate-passing lines remain
    # unpicked, swap them in (keeps the strata mix where possible).
    picked = set(ordered)
    gate_pass = [m for m in test_pool if vals[m] < GATE and m not in picked]
    gate_pass.sort(key=lambda m: vals[m])
    for i, m in enumerate(ordered):
        if vals[m] >= GATE and gate_pass:
            ordered[i] = gate_pass.pop(0)
    for m in ordered:
        rows.append({"gene": gene, "role": "test",
                     "cell_line": names.get(m, m), "model_id": m,
                     "STK11_mutant": m in mut_lines["STK11"],
                     "KEAP1_mutant": m in mut_lines["KEAP1"],
                     "KRAS_mutant": m in mut_lines["KRAS"],
                     "chronos": float(vals[m]),
                     "passes_gate": bool(vals[m] < GATE)})
    ctrl_sorted = sorted([m for m in ctrl_pool if vals[m] > GATE],
                         key=lambda m: -vals[m])
    for m in ctrl_sorted[:N_CONTROL]:
        rows.append({"gene": gene, "role": "control",
                     "cell_line": names.get(m, m), "model_id": m,
                     "STK11_mutant": False, "KEAP1_mutant": False,
                     "KRAS_mutant": m in mut_lines["KRAS"],
                     "chronos": float(vals[m]),
                     "passes_gate": False})
lines = pd.DataFrame(rows)
lines.to_csv(TABLES / "validation_lines.csv", index=False)
print(f"      {len(lines)} line assignments written "
      f"({lines['role'].value_counts().to_dict()})")
for gene in plan_genes:
    sub = lines[lines["gene"] == gene]
    print(f"      {gene}: test "
          + ", ".join(f"{r['cell_line']} ({r['chronos']:.2f})"
                      for _, r in sub[sub["role"] == "test"].iterrows())
          + " | control "
          + ", ".join(f"{r['cell_line']} ({r['chronos']:.2f})"
                      for _, r in sub[sub["role"] == "control"].iterrows()))

print("[4/4] Writing docs/validation-plan.md ...")
scl = pd.read_csv(TABLES / "singlecell_layer.csv").set_index("gene")
parts = []
parts.append("# Validation experiment plan\n")
parts.append("The computational shortlist is a hypothesis list. This plan "
             "says exactly how to test it in a lab: which cell lines, "
             "which perturbation, what to measure, and what result would "
             "confirm or kill each candidate. The lines and the Chronos "
             "values below come from DepMap Public 24Q4 "
             "(results/tables/validation_lines.csv); the layer values "
             "come from the Phase 3 ranking "
             "(results/tables/shortlist_dependency_supported_phase3.csv).\n")
parts.append("## Design, shared across genes\n")
parts.append("- Perturbation: lentiviral CRISPR-Cas9 knockout, three "
             "independent sgRNAs per gene plus a non-targeting control "
             "sgRNA, puromycin selected. A gene counts as tested only if "
             "knockout is confirmed (western blot or qPCR for the target, "
             "or sequencing of the edited locus). A failed knockout is a "
             "technical failure, not a kill: repeat with new guides.\n")
parts.append("- Readout: viability by CellTiter-Glo at days 3, 7, 10, and "
             "14 after selection, plus live-cell confluence imaging over "
             "the same window, and annexin V apoptosis staining at day 7. "
             "Three biological replicates per line per guide.\n")
parts.append("- Controls in every run: the non-targeting sgRNA in every "
             "line, and a known essential gene knockout (for example "
             "PLK1) as a positive control for the assay itself.\n")
parts.append("- Pre-specified confirm bar: at least 50 percent viability "
             "loss relative to the non-targeting control in at least two "
             "of the test lines, with under 20 percent loss in the wild "
             "type control lines, at day 14. The bar is set before the "
             "experiment runs, so a borderline result cannot be talked "
             "into a win.\n")
parts.append("- Pre-specified kill outcomes: no test line meets the bar "
             "(the dependency does not reproduce); or the wild type "
             "controls lose as much viability as the test lines (the gene "
             "is a general essential, not a subgroup target); or the "
             "effect is guide-specific (one guide kills, the other two do "
             "not), which points to an off-target effect.\n")
parts.append("- If a candidate confirms: dose-response with a second "
             "perturbation method (inducible knockdown or, where one "
             "exists, a small molecule inhibitor), then a mouse xenograft "
             "with the best test line.\n")

for gene in plan_genes:
    r = short_idx.loc[gene]
    sub = lines[lines["gene"] == gene]
    test = sub[sub["role"] == "test"]
    ctrl = sub[sub["role"] == "control"]
    ms = scl.loc[gene, "malignant_share"] if gene in scl.index else np.nan
    parts.append(f"## {gene} (shortlist rank "
                 f"{int(r['subset_rank_main'])})\n")
    parts.append(f"Why it is on the list: dependency layer value "
                 f"{r['dependency']:.3f} (selective essentiality "
                 f"{r['selective_essentiality']:.3f}, mean Chronos "
                 f"{r['mean_chronos_mutant_lung']:.3f} in mutant lung "
                 f"lines), safety window {r['safety_window']:.3f}"
                 + (f", malignant cell share {ms:.3f} in the single cell "
                    f"layer" if pd.notna(ms) else ", not detected in the "
                    f"single cell layer, so its tumor location is "
                    f"unconfirmed")
                 + f". Overall Phase 3 score {r['score']:.3f}.\n")
    parts.append("Test lines (STK11/KEAP1-mutant lung lines):\n")
    test_rows = ["| Cell line | Model ID | STK11 | KEAP1 | KRAS | Chronos for "
                 f"{gene} | Passes the -0.5 gate |",
                 "| --- | --- | --- | --- | --- | --- | --- |"]
    for _, t in test.iterrows():
        test_rows.append(f"| {t['cell_line']} | {t['model_id']} | "
                         f"{'mutant' if t['STK11_mutant'] else 'wild type'} | "
                         f"{'mutant' if t['KEAP1_mutant'] else 'wild type'} | "
                         f"{'mutant' if t['KRAS_mutant'] else 'wild type'} | "
                         f"{t['chronos']:.3f} | "
                         f"{'yes' if t['passes_gate'] else 'no'} |")
    parts.append("\n".join(test_rows) + "\n")
    parts.append("Control lines (wild type for STK11 and KEAP1):\n")
    ctrl_rows = ["| Cell line | Model ID | KRAS | Chronos for "
                 f"{gene} |",
                 "| --- | --- | --- | --- |"]
    for _, t in ctrl.iterrows():
        ctrl_rows.append(f"| {t['cell_line']} | {t['model_id']} | "
                         f"{'mutant' if t['KRAS_mutant'] else 'wild type'} | "
                         f"{t['chronos']:.3f} |")
    parts.append("\n".join(ctrl_rows) + "\n")
    parts.append(f"\nConfirm: {gene} knockout kills the test lines "
                 f"selectively, meeting the shared confirm bar above, with "
                 f"the strongest effect in the lines where DepMap shows "
                 f"the most negative Chronos. Kill: no selective effect, "
                 f"or equal killing in the wild type controls.\n")
    if gene == "NFE2L2":
        parts.append("Orthogonal arm: NFE2L2 encodes NRF2, and the "
                     "published small molecule NRF2 inhibitor ML385 can "
                     "serve as a second, independent perturbation in the "
                     "same lines. If knockout and ML385 agree, the case "
                     "for the target is much stronger than either alone. "
                     "If they disagree, the knockout result governs, "
                     "because ML385 has its own off-target profile.\n")

parts.append("## What this plan does not cover\n")
parts.append("- The single cell layer rests on 1,953 LUAD cells from two "
             "donors (see results/tables/singlecell_qc.json), so the "
             "tumor location evidence is thinner than the dependency "
             "evidence. A candidate that confirms here still needs its "
             "expression checked in the actual test lines (qPCR) before "
             "the knockout result is interpreted.\n")
parts.append("- None of these lines is a perfect model of a patient "
             "tumor. A confirmed candidate is a lead for animal work, "
             "not a drug.\n")
parts.append("\nSources: results/tables/validation_lines.csv, "
             "results/tables/shortlist_dependency_supported_phase3.csv, "
             "results/tables/singlecell_layer.csv, "
             "results/tables/dependency_layer.csv.\n")
(DOCS / "validation-plan.md").write_text("\n".join(parts))
print("Step 6 (validation plan) complete.")
