"""What is actually at the top of the 2B69 screen?

The raw Vina score grows with molecular size, so sorting the library by
raw score mostly sorts it by atom count. This script does the three
checks that have to pass before any compound name is worth saying out
loud:

  1. deduplicate salt forms, which score identically because ligand
     preparation strips them to the same parent
  2. measure how much of the ranking is still a size measurement
  3. rank three different ways and see whether the answers agree

Reference scores come from the 2B69 control, which docked the five
natural ligands under exactly this protocol.

Writes outputs/analysis.md and outputs/ranked.csv
"""
import csv
import os
import statistics

from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "outputs", "results.csv")
LIBRARY = os.path.join(os.path.dirname(BASE), "screening-full-uxs1", "library.csv")
OUT_MD = os.path.join(BASE, "outputs", "analysis.md")
OUT_CSV = os.path.join(BASE, "outputs", "ranked.csv")

# from control-2b69/report/dock-check.md, same receptor, same protocol
NATURAL = [
    ("UDP-xylose", 34, -11.126, 99.5),
    ("UDP-galactose", 36, -10.996, 99.0),
    ("UDP-glucuronate", 37, -10.213, 88.4),
    ("UDP", 25, -10.021, 84.9),
    ("UMP", 21, -8.927, 53.3),
]
SUBSTRATE_SCORE = -10.213
SIZE_WINDOW = 4          # heavy atoms either side, for the size matched percentile
MIN_PEERS = 25           # below this the local percentile is not trustworthy


def load():
    lib = {}
    for row in csv.DictReader(open(LIBRARY)):
        lib[row["chembl_id"]] = (row["pref_name"], row["canonical_smiles"])
    rows = []
    for row in csv.DictReader(open(RESULTS)):
        cid = row["id"]
        name, smiles = lib.get(cid, ("", ""))
        rows.append({"id": cid, "name": name, "smiles": smiles,
                     "score": float(row["score"]), "heavy": int(row["heavy"])})
    return rows


def parent_key(smiles):
    """InChIKey of the largest fragment, so salt forms collapse together."""
    if smiles == "":
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def deduplicate(rows):
    """One row per parent structure, keeping the best score."""
    seen = {}
    no_key = []
    for r in rows:
        k = parent_key(r["smiles"])
        if k is None:
            no_key.append(r)
            continue
        r["parent"] = k
        if k not in seen or r["score"] < seen[k]["score"]:
            if k in seen:
                seen[k]["dropped"] = seen[k].get("dropped", 0) + 1
            seen[k] = r
    kept = list(seen.values()) + no_key
    return kept, len(rows) - len(kept)


def local_percentile(row, rows):
    """Percentile among library compounds of a similar size.

    This is the same comparison the control made for the natural
    ligands: how does this molecule do against molecules its own size,
    rather than against the whole library.
    """
    lo = row["heavy"] - SIZE_WINDOW
    hi = row["heavy"] + SIZE_WINDOW
    peers = [r["score"] for r in rows
             if lo <= r["heavy"] <= hi and r["id"] != row["id"]]
    if len(peers) < MIN_PEERS:
        return None, len(peers)
    beaten = sum(1 for s in peers if s < row["score"])
    return 100.0 * (len(peers) - beaten) / len(peers), len(peers)


def table(rows, n, key, label, reverse=False):
    out = [f"### {label}", "",
           "| # | Compound | ChEMBL | Heavy | Score | LE | Size matched pct |",
           "| --- | --- | --- | --- | --- | --- | --- |"]
    picked = sorted(rows, key=key, reverse=reverse)[:n]
    for i, r in enumerate(picked, 1):
        pct = "n/a" if r["pct"] is None else f"{r['pct']:.1f}"
        name = r["name"] if r["name"] != "" else "(unnamed)"
        out.append(f"| {i} | {name.title()} | {r['id']} | {r['heavy']} | "
                   f"{r['score']:.2f} | {r['le']:.3f} | {pct} |")
    out.append("")
    return out, picked


def main():
    rows = load()
    n_raw = len(rows)
    rows, n_dropped = deduplicate(rows)
    print(f"{n_raw} docked rows, {n_dropped} duplicate salt forms removed, "
          f"{len(rows)} structures", flush=True)

    scores = [r["score"] for r in rows]
    heavies = [r["heavy"] for r in rows]

    rho_all, p_all = spearmanr(heavies, scores)
    window = [r for r in rows if 31 <= r["heavy"] <= 43]
    rho_win, p_win = spearmanr([r["heavy"] for r in window],
                               [r["score"] for r in window])

    for r in rows:
        r["le"] = r["score"] / r["heavy"]
        r["pct"], r["peers"] = local_percentile(r, rows)

    lines = ["# What is at the top of the 2B69 screen", "",
             f"{n_raw} compounds docked, {n_dropped} duplicate salt forms removed, "
             f"{len(rows)} distinct structures ranked.", "",
             "## 1. Is the ranking still a size measurement", "",
             "| Set | n | Spearman rho, score against heavy atoms | p |",
             "| --- | --- | --- | --- |",
             f"| whole library | {len(rows)} | {rho_all:+.3f} | {p_all:.1e} |",
             f"| 31 to 43 heavy atoms, the control's window | {len(window)} | "
             f"{rho_win:+.3f} | {p_win:.1e} |", "",
             "The control reported +0.11 against size matched decoys. A strong",
             "negative number across the whole library is expected and is not the",
             "same failure as 4GLL: Vina sums over atom pairs, so more atoms earn",
             "a more negative score almost mechanically. What matters is whether",
             "the correlation survives inside a narrow size window, because that",
             "is where a real binder has to win.", ""]

    beat = [r for r in rows if r["score"] < SUBSTRATE_SCORE]
    lines += ["## 2. How many beat the substrate", "",
              f"{len(beat)} of {len(rows)} structures ({100.0*len(beat)/len(rows):.1f}%) "
              f"score better than UDP-glucuronate at {SUBSTRATE_SCORE:.2f}.", "",
              "In the broken 4GLL run about two thirds of the library beat the",
              "substrate, which was the signal that the comparison was empty.", ""]

    lines += ["## 3. The three rankings", "",
              "Raw score favours big molecules. Ligand efficiency, the score",
              "divided by heavy atom count, favours small ones. The size matched",
              "percentile asks the only fair question: how does this molecule do",
              "against molecules its own size. Where the three disagree is where",
              "the interesting compound usually is.", ""]

    t1, raw_top = table(rows, 25, lambda r: r["score"], "Ranked by raw score")
    t2, le_top = table(rows, 25, lambda r: r["le"], "Ranked by ligand efficiency")
    scored = [r for r in rows if r["pct"] is not None]
    t3, pct_top = table(scored, 25, lambda r: -r["pct"],
                        "Ranked by size matched percentile")
    lines += t1 + t2 + t3

    raw_ids = set(r["id"] for r in raw_top)
    le_ids = set(r["id"] for r in le_top)
    pct_ids = set(r["id"] for r in pct_top)
    lines += ["### Agreement between the three top 25s", "",
              f"- raw and size matched percentile share {len(raw_ids & pct_ids)} compounds",
              f"- raw and ligand efficiency share {len(raw_ids & le_ids)} compounds",
              f"- all three share {len(raw_ids & le_ids & pct_ids)} compounds", ""]
    all_three = raw_ids & le_ids & pct_ids
    if len(all_three) > 0:
        names = [r["name"].title() for r in rows if r["id"] in all_three]
        lines += ["Present in all three: " + ", ".join(sorted(names)), ""]

    lines += ["## 4. The natural ligands, for scale", "",
              "| Ligand | Heavy | Score | Percentile against size matched decoys |",
              "| --- | --- | --- | --- |"]
    for name, heavy, score, pct in NATURAL:
        lines.append(f"| {name} | {heavy} | {score:.2f} | {pct} |")
    lines += ["",
              "A drug worth looking at has to clear the substrate's standing",
              "among molecules its own size, not just its raw number.", ""]

    strong = [r for r in scored if r["pct"] >= 99.0]
    lines += ["## 5. What clears the bar", "",
              f"{len(strong)} structures sit at or above the 99th percentile "
              f"among molecules their own size, which is where UDP-xylose and "
              f"UDP-galactose sit.", ""]
    if len(strong) > 0:
        lines += ["| Compound | ChEMBL | Heavy | Score | LE | Pct | Peers |",
                  "| --- | --- | --- | --- | --- | --- | --- |"]
        for r in sorted(strong, key=lambda r: r["score"])[:40]:
            name = r["name"] if r["name"] != "" else "(unnamed)"
            lines.append(f"| {name.title()} | {r['id']} | {r['heavy']} | "
                         f"{r['score']:.2f} | {r['le']:.3f} | {r['pct']:.1f} | "
                         f"{r['peers']} |")
        lines.append("")

    lines += ["## What this does not yet show", "",
              "Nothing here looks at where the pose actually sits. A good score",
              "with the molecule pressed against a box wall, or parked in some",
              "other groove, is not a substrate site binder. That check runs",
              "next, on the shortlist, and it is the one that killed the",
              "4GLL shortlist.", ""]

    open(OUT_MD, "w").write("\n".join(lines))

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "name", "heavy", "score", "le",
                                           "pct", "peers", "smiles"])
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["score"]):
            w.writerow({"id": r["id"], "name": r["name"], "heavy": r["heavy"],
                        "score": round(r["score"], 3), "le": round(r["le"], 4),
                        "pct": "" if r["pct"] is None else round(r["pct"], 1),
                        "peers": r["peers"], "smiles": r["smiles"]})

    print(f"rho whole library {rho_all:+.3f}, rho in window {rho_win:+.3f}")
    print(f"{len(beat)} beat the substrate, {len(strong)} at or above the 99th pct")
    print(f"written to {OUT_MD}")


if __name__ == "__main__":
    main()
