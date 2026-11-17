"""Re-dock the size matched shortlist at exhaustiveness 8.

The 100 compound screen showed why this step exists: glimepiride led
pass 1 at exhaustiveness 4 and then lost 2.96 kcal/mol when re-docked,
because its pass 1 pose was a lucky result of the coarse search rather
than a stable minimum. Anything that does not hold up here is noise.

Reads:  outputs/shortlist_size_matched.csv
Writes: outputs/redock/<cid>_energies.txt, outputs/results_redock.csv
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
OUT_DIR = os.path.join(BASE, "outputs", "redock")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 8
N_POSES = 9
SEED = 42

shortlist = list(csv.DictReader(
    open(os.path.join(BASE, "outputs", "shortlist_size_matched.csv"))))
print(f"Re-dock: {len(shortlist)} compounds at exhaustiveness {EXHAUSTIVENESS}",
      flush=True)

results = []
t_start = time.time()
for i, row in enumerate(shortlist, start=1):
    cid = row["chembl_id"]
    name = row["pref_name"]
    pass1 = float(row["best_affinity_kcal_per_mol"])
    lig = os.path.join(BASE, "ligands", cid + ".pdbqt")
    energy_file = os.path.join(OUT_DIR, cid + "_energies.txt")

    if os.path.exists(energy_file) == True:
        with open(energy_file) as f:
            best = float(f.readline().split()[0])
    else:
        t0 = time.time()
        v = Vina(sf_name="vina", verbosity=0, seed=SEED)
        v.set_receptor(RECEPTOR)
        v.set_ligand_from_file(lig)
        v.compute_vina_maps(center=CENTER, box_size=SIZE)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        best = energies[0][0]
        with open(energy_file, "w") as f:
            for e in energies:
                f.write(" ".join(f"{x:.3f}" for x in e) + "\n")
        v.write_poses(os.path.join(OUT_DIR, cid + "_out.pdbqt"),
                      n_poses=1, overwrite=True)
        print(f"  [{i}/{len(shortlist)}] {name}: {best:.3f} "
              f"(pass 1 {pass1:.3f}, {time.time() - t0:.0f}s)", flush=True)

    results.append({
        "chembl_id": cid,
        "pref_name": name,
        "redock_exh8": round(best, 3),
        "pass1_exh4": pass1,
        "shift": round(best - pass1, 3),
        "heavy_atoms": row["heavy_atoms"],
        "held": abs(best - pass1) <= 0.5,
    })

results.sort(key=lambda r: r["redock_exh8"])
out = os.path.join(BASE, "outputs", "results_redock.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    for r in results:
        w.writerow(r)

held = sum(1 for r in results if r["held"] == True)
print(f"\nheld within 0.5 kcal/mol: {held} of {len(results)}", flush=True)
print(f"written to {out}", flush=True)
print(f"wall time {(time.time() - t_start) / 60:.1f} min", flush=True)
