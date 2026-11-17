"""Re-dock the shortlist in a box centred on the real substrate site.

Two problems with the original box, both measured in analyze_poses.py:

  1. It is centred at (60.0, 23.2, 39.0), which is 8.6 A away from the
     crystallographic UDP position at (62.3, 16.4, 34.2). The crystal
     UDP fits inside, but only just, coming within 2.1 A of a wall.
  2. Half the re-docked poses come within 2 A of a wall and 11 of 40
     come within 1 A, meaning the search was cut off. A ligand that
     wants to extend further cannot, so its score is an artefact of the
     box rather than of the site.

This run centres the box on the crystal UDP centroid and widens it from
24 to 28 A, which clears the walls. It also runs three different random
seeds per compound, because the earlier re-dock reused seed 42 from
pass 1 and so could not test pose stability independently.

Writes outputs/recentred/<cid>_seed<N>_energies.txt and
outputs/results_recentred.csv
"""
import csv
import os
import statistics
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
OUT_DIR = os.path.join(BASE, "outputs", "recentred")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [62.3, 16.4, 34.2]      # crystal UDP centroid, from 4GLL
SIZE = [28, 28, 28]
EXHAUSTIVENESS = 8
N_POSES = 9
SEEDS = [11, 42, 97]

shortlist = list(csv.DictReader(
    open(os.path.join(BASE, "outputs", "shortlist_size_matched.csv"))))
print(f"Re-centred run: {len(shortlist)} compounds, {len(SEEDS)} seeds, "
      f"box {SIZE[0]} A at {CENTER}", flush=True)

rows = []
t_start = time.time()
for i, item in enumerate(shortlist, start=1):
    cid = item["chembl_id"]
    name = item["pref_name"]
    old = float(item["best_affinity_kcal_per_mol"])
    lig = os.path.join(BASE, "ligands", cid + ".pdbqt")
    if os.path.exists(lig) == False:
        continue

    scores = []
    for seed in SEEDS:
        ef = os.path.join(OUT_DIR, f"{cid}_seed{seed}_energies.txt")
        if os.path.exists(ef) == True:
            with open(ef) as f:
                scores.append(float(f.readline().split()[0]))
            continue
        v = Vina(sf_name="vina", verbosity=0, seed=seed)
        v.set_receptor(RECEPTOR)
        v.set_ligand_from_file(lig)
        v.compute_vina_maps(center=CENTER, box_size=SIZE)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        scores.append(energies[0][0])
        with open(ef, "w") as f:
            for e in energies:
                f.write(" ".join(f"{x:.3f}" for x in e) + "\n")
        if seed == SEEDS[0]:
            v.write_poses(os.path.join(OUT_DIR, cid + "_out.pdbqt"),
                          n_poses=1, overwrite=True)

    spread = max(scores) - min(scores)
    rows.append({
        "chembl_id": cid,
        "pref_name": name,
        "mean_affinity": round(statistics.mean(scores), 3),
        "best_affinity": round(min(scores), 3),
        "seed_spread": round(spread, 3),
        "stable": spread <= 0.5,
        "old_box_affinity": old,
        "change_vs_old_box": round(statistics.mean(scores) - old, 3),
        "heavy_atoms": item["heavy_atoms"],
    })
    print(f"  [{i}/{len(shortlist)}] {name[:24]}: mean {statistics.mean(scores):.3f}, "
          f"spread {spread:.3f} (old box {old:.3f})", flush=True)

rows.sort(key=lambda r: r["mean_affinity"])
out = os.path.join(BASE, "outputs", "results_recentred.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

stable = sum(1 for r in rows if r["stable"] == True)
print(f"\nstable across seeds (spread <= 0.5): {stable} of {len(rows)}", flush=True)
print(f"written to {out}", flush=True)
print(f"wall time {(time.time() - t_start)/60:.1f} min", flush=True)
