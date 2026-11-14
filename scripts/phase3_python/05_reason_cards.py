"""Phase 3, step 5: plain-language reason cards for the top candidates.

Generates one short markdown file per gene in the v3
dependency-supported shortlist (results/reason_cards/<GENE>.md),
built only from the project's own result tables. No functional
claims are imported from outside the project: the two gene
descriptions used here (NFE2L2 encodes NRF2; UXS1 is an enzyme in
sugar metabolism) are the ones already stated in
results/phase2-report.md. Every number in a card traces to a table
named in the card's footer.

Each card states the gene's rank and score, what the evidence says
in one paragraph, the five layers with their numbers, and the weak
points ("Watch out") derived from the same tables.

Outputs:
  results/reason_cards/<GENE>.md (20 files)

Run from the repository root with the project venv:
  .venv/bin/python scripts/phase3_python/05_reason_cards.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "results" / "reason_cards"
OUT.mkdir(parents=True, exist_ok=True)

# Descriptions already stated in results/phase2-report.md. Nothing else
# is asserted about what a gene does.
KNOWN_DESCRIPTIONS = {
    "NFE2L2": "the gene that encodes NRF2, the transcription factor KEAP1 normally holds in check",
    "UXS1": "an enzyme in sugar metabolism",
}

print("Loading tables ...")
short = pd.read_csv(TABLES / "shortlist_dependency_supported_phase3.csv")
scl = pd.read_csv(TABLES / "singlecell_layer.csv").set_index("gene")
cox = pd.read_csv(TABLES / "cox_phase3.csv").set_index("gene")
dep = pd.read_csv(TABLES / "dependency_layer.csv").set_index("gene")
de_tv = pd.read_csv(TABLES / "de_tumor_vs_normal_shrunk.csv").set_index("gene")
de_mv = pd.read_csv(TABLES / "de_mutant_vs_wt.csv").set_index("gene")
ts = pd.read_csv(TABLES / "timesplit_comparison.csv").set_index("gene")
dep20 = pd.read_csv(TABLES / "dependency_layer_20q2.csv").set_index("gene")
n_ranked = len(pd.read_csv(TABLES / "target_scores_phase3.csv", usecols=["gene"]))

WEIGHTS = {"dependency": 0.30, "tumor_specificity": 0.25,
           "safety_window": 0.20, "malignant_share": 0.15,
           "outcome_link": 0.10}
LAYER_NAMES = {"dependency": "dependency", "tumor_specificity": "expression",
               "safety_window": "safety", "malignant_share": "single cell location",
               "outcome_link": "survival"}


def fmt(x, digits=3):
    return f"{x:.{digits}f}" if pd.notna(x) else "not measured"


def contribution(row):
    """Layer contributions (value x weight) for the lead sentence."""
    parts = []
    for col, w in WEIGHTS.items():
        v = row.get(col)
        if pd.notna(v):
            parts.append((LAYER_NAMES[col], float(v) * w, float(v)))
    parts.sort(key=lambda t: -t[1])
    return parts


def write_card(r):
    gene = r["gene"]
    rank = int(r["subset_rank_main"])
    overall = int(r["rank"])
    score = float(r["score"])
    robust = r["shortlist_robust_fragile"]
    robust_txt = ("It stays in the top 20 under all three weightings, so the "
                  "rank does not depend on one weighting choice." if robust == "robust"
                  else "It reaches the top 20 under only one weighting, so the "
                  "rank is sensitive to the weighting choice." if robust == "fragile"
                  else "It is in the top 20 under the main weighting and some "
                  "but not all alternate weightings.")

    contrib = contribution(r)
    lead_layers = ", ".join(f"{name} ({fmt(v)})" for name, _, v in contrib[:2])

    lines = []
    lines.append(f"# {gene}\n")
    desc = KNOWN_DESCRIPTIONS.get(gene)
    if desc:
        lines.append(f"{gene} is {desc}.\n")
    lines.append(f"Shortlist rank {rank} of 20 (overall rank {overall} of "
                 f"{n_ranked:,}). Score {fmt(score)}. {robust_txt}\n")
    lines.append("## What the evidence says\n")
    lines.append(f"The rank is driven most by {lead_layers}. "
                 f"Coverage is {int(r['coverage'])} of 5 layers "
                 f"({r['layers_used']}).\n")
    lines.append("## Layer by layer\n")
    lines.append("### Expression (TCGA)\n")
    lines.append(f"- Tumor specificity layer value: {fmt(r['tumor_specificity'])}. "
                 f"Log2 fold change, tumor vs healthy lung (shrunk): "
                 f"{fmt(r['log2FC_tumor_vs_normal_shrunk'])}. "
                 f"Log2 fold change, mutant vs wild type tumors: "
                 f"{fmt(r['log2FC_mutant_vs_wt'])}.")
    if gene in de_mv.index:
        lines.append(f"  Mutant vs wild type adjusted p value: "
                     f"{de_mv.loc[gene, 'padj']:.2e}.")
    lines.append("")
    lines.append("### Dependency (DepMap CRISPR)\n")
    lines.append(f"- Dependency layer value: {fmt(r['dependency'])}. "
                 f"Selective essentiality: {fmt(r['selective_essentiality'])} "
                 f"(mean Chronos in STK11/KEAP1-mutant lung lines minus mean "
                 f"Chronos in other lung lines; more negative means the "
                 f"mutant lines need the gene more). Mean Chronos in mutant "
                 f"lung lines: {fmt(r['mean_chronos_mutant_lung'])}.")
    if gene in dep.index:
        d = dep.loc[gene]
        lines.append(f"  Behind this: {int(d['n_mutant_lung_lines'])} mutant "
                     f"lung lines and {int(d['n_other_lung_lines'])} other "
                     f"lung lines with CRISPR data. Common essential gene: "
                     f"{'yes' if bool(d['common_essential']) else 'no'}.")
    lines.append("")
    lines.append("### Single cell location (CELLxGENE)\n")
    if gene in scl.index and pd.notna(scl.loc[gene, "malignant_share"]):
        s = scl.loc[gene]
        lines.append(f"- Malignant cell share layer value: "
                     f"{fmt(r['malignant_share'])} (raw share of the gene's "
                     f"total tumor expression in malignant cells: "
                     f"{fmt(s['malignant_share'])}; immune "
                     f"{fmt(s['immune_share'])}, stromal "
                     f"{fmt(s['stromal_share'])}, other "
                     f"{fmt(s['other_share'])}). Expressed in "
                     f"{fmt(s['pct_malignant_expressing'], 1)}% of malignant "
                     f"cells, {fmt(s['pct_immune_expressing'], 1)}% of immune "
                     f"cells, {fmt(s['pct_stromal_expressing'], 1)}% of "
                     f"stromal cells.")
    else:
        lines.append("- Not detected in the single cell dataset (zero counts "
                     "in the LUAD cells used), so this layer does not "
                     "contribute to the score and the gene's tumor location "
                     "is unconfirmed.")
    lines.append("")
    lines.append("### Safety (GTEx)\n")
    lines.append(f"- Safety window layer value: {fmt(r['safety_window'])} "
                 f"(1 means quiet in heart, brain, and liver; 0 means loud).")
    lines.append("")
    lines.append("### Survival (TCGA, mutant subgroup)\n")
    if gene in cox.index and pd.notna(cox.loc[gene, "log_hr_per_sd"]):
        c = cox.loc[gene]
        lines.append(f"- Outcome layer value: {fmt(r['outcome_link'])}. "
                     f"Cox hazard ratio per standard deviation of expression: "
                     f"{fmt(c['hr_per_sd'])} (log hazard ratio "
                     f"{fmt(c['log_hr_per_sd'])}), p value "
                     f"{c['pvalue']:.3g}, n={int(c['n_patients'])} mutant "
                     f"subgroup patients.")
    else:
        lines.append("- Survival was not modeled for this gene (it was "
                     "outside the Cox set), so this layer does not "
                     "contribute to the score.")
    lines.append("")
    lines.append("## Watch out\n")
    flags = []
    if pd.notna(r["safety_window"]) and r["safety_window"] < 0.6:
        flags.append(f"Safety is the weak point: the safety window is "
                     f"{fmt(r['safety_window'])}, so the gene is active in "
                     f"normal heart, brain, or liver.")
    if robust == "fragile":
        flags.append("The rank is fragile: it reaches the top 20 under only "
                     "one of the three weightings.")
    if gene in scl.index and pd.notna(scl.loc[gene, "malignant_share"]):
        ms = scl.loc[gene, "malignant_share"]
        if ms < 0.3:
            flags.append(f"Most of its expression in the tumor sits outside "
                         f"the malignant cells (malignant share "
                         f"{fmt(ms)}), so a drug aimed at it would act "
                         f"largely on non-tumor cells.")
    else:
        flags.append("Not detected in the single cell dataset, so its "
                     "tumor location is unconfirmed; the rank rests on the "
                     "other four layers.")
    if gene in ts.index:
        sel20 = ts.loc[gene, "selective_essentiality_20q2"]
        if pd.notna(sel20) and sel20 >= 0:
            flags.append(f"The time split check found no dependency signal "
                         f"before 2021 (selective essentiality "
                         f"{fmt(sel20)} in DepMap 20Q2), so the dependency "
                         f"case rests on the newer data alone.")
        elif not bool(ts.loc[gene, "in_frozen_shortlist"]):
            if (gene in dep20.index
                    and pd.notna(dep20.loc[gene, "mean_ceres_mutant_lung"])
                    and dep20.loc[gene, "mean_ceres_mutant_lung"] >= -0.5):
                flags.append(f"The time split check did not keep it in the "
                             f"frozen shortlist: its 20Q2 selective "
                             f"essentiality was {fmt(sel20)}, but its mean "
                             f"CERES in mutant lung lines "
                             f"({fmt(dep20.loc[gene, 'mean_ceres_mutant_lung'])}) "
                             f"fell short of the -0.5 essentiality gate, so "
                             f"it sits on the threshold.")
            else:
                flags.append("The time split check did not keep it in the "
                             "frozen pre-2021 shortlist, so its case is "
                             "weaker in the older data.")
    if r["tumor_specificity"] == 0:
        flags.append("Expression does not distinguish mutant tumors from "
                     "wild type tumors; the rank rests on the other layers.")
    if gene in cox.index and pd.notna(cox.loc[gene, "pvalue"]):
        if cox.loc[gene, "pvalue"] >= 0.05:
            flags.append("No significant survival link in the mutant "
                         "subgroup, so this is a dependency hypothesis, not "
                         "a prognosis story.")
    if not flags:
        flags.append("No single weak point stands out in the tables; the "
                     "usual caveats still apply (40 mutant lines behind the "
                     "dependency claim, one single cell dataset, survival "
                     "from 132 patients).")
    for f in flags:
        lines.append(f"- {f}")
    lines.append("")
    lines.append("## Sources in this repository\n")
    lines.append("results/tables/shortlist_dependency_supported_phase3.csv, "
                 "results/tables/singlecell_layer.csv, "
                 "results/tables/cox_phase3.csv, "
                 "results/tables/dependency_layer.csv, "
                 "results/tables/timesplit_comparison.csv, "
                 "results/tables/dependency_layer_20q2.csv, "
                 "results/tables/de_mutant_vs_wt.csv, "
                 "results/tables/de_tumor_vs_normal_shrunk.csv.")
    (OUT / f"{gene}.md").write_text("\n".join(lines))
    return gene


print(f"Writing reason cards for {len(short)} shortlist genes ...")
for _, row in short.iterrows():
    write_card(row)
    print(f"      {row['gene']}: rank {int(row['subset_rank_main'])}")
print(f"Step 5 (reason cards) complete: {len(short)} cards in "
      f"results/reason_cards/.")
