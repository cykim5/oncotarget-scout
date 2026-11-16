"""Size correct the pass 1 screen.

Vina's scoring function adds up favourable contacts, so a bigger molecule
makes more contacts and scores better whether or not it fits the site in
any meaningful way. This script measures how strong that effect is in our
own screen, then compares compounds against the natural product on a like
for like basis instead.

Reads:  outputs/pass1/results_pass1.csv, ligands/*.sdf
Writes: outputs/pass1/results_size_corrected.csv
"""

import csv
import os
import statistics

BASE = os.path.dirname(os.path.abspath(__file__))

# The natural product from the docking-uxs1 run, same receptor and box.
REF_NAME = "UDP-xylose"
REF_AFFINITY = -8.796
REF_HEAVY = 34

# A compound is only a fair comparison against the reference if it is
# roughly the same size. This band is the reference plus or minus 6 heavy
# atoms, which is wide enough to keep a useful number of compounds.
BAND_LOW = 30
BAND_HIGH = 40


def heavy_atom_count(sdf_path):
    """Count non-hydrogen atoms in a V2000 SDF atom block."""
    lines = open(sdf_path).read().splitlines()
    n_atoms = int(lines[3][0:3])
    count = 0
    for line in lines[4:4 + n_atoms]:
        element = line[31:34].strip()
        if element != "H":
            count = count + 1
    return count


def pearson(xs, ys):
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    top = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    bottom = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    return top / bottom


rows = []
results_csv = os.path.join(BASE, "outputs", "pass1", "results_pass1.csv")
with open(results_csv) as f:
    for record in csv.DictReader(f):
        sdf = os.path.join(BASE, "ligands", record["chembl_id"] + ".sdf")
        if os.path.exists(sdf) == False:
            continue
        heavy = heavy_atom_count(sdf)
        affinity = float(record["best_affinity_kcal_per_mol"])
        rows.append({
            "chembl_id": record["chembl_id"],
            "pref_name": record["pref_name"],
            "best_affinity_kcal_per_mol": affinity,
            "heavy_atoms": heavy,
            "ligand_efficiency": round(-affinity / heavy, 4),
            "in_size_band": BAND_LOW <= heavy <= BAND_HIGH,
            "beats_reference": affinity < REF_AFFINITY,
        })

r = pearson([x["heavy_atoms"] for x in rows],
            [x["best_affinity_kcal_per_mol"] for x in rows])

band = [x for x in rows if x["in_size_band"] == True]
band.sort(key=lambda x: x["best_affinity_kcal_per_mol"])
winners = [x for x in band if x["beats_reference"] == True]

out_csv = os.path.join(BASE, "outputs", "pass1", "results_size_corrected.csv")
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    for x in sorted(rows, key=lambda v: v["best_affinity_kcal_per_mol"]):
        writer.writerow(x)

print(f"compounds scored: {len(rows)}")
print(f"heavy atoms vs affinity: r = {r:.2f}  "
      f"(size alone explains about {r * r * 100:.0f}% of the ranking)")
print(f"size band {BAND_LOW} to {BAND_HIGH} heavy atoms: {len(band)} compounds")
print(f"beating {REF_NAME} ({REF_AFFINITY} kcal/mol, {REF_HEAVY} heavy atoms): "
      f"{len(winners)}")
for x in winners:
    print(f"  {x['pref_name']}  {x['best_affinity_kcal_per_mol']}  "
          f"{x['heavy_atoms']} heavy atoms")
print(f"written to {out_csv}")
