"""Pass 1 over the full approved drug library, same setup as the 100.

Receptor, box, scoring function and exhaustiveness are identical to
screening-uxs1/run_screen_pass1.py and to docking-uxs1 run 2, so scores
are comparable to the natural ligand references from that run.

Resumable: a compound with an existing energies file is skipped, so the
run can be stopped and restarted without losing work.
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
LIGAND_DIR = os.path.join(BASE, "ligands")
LIBRARY = os.path.join(BASE, "library.csv")
OUT_DIR = os.path.join(BASE, "outputs")
RESULTS = os.path.join(OUT_DIR, "results_full.csv")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 4
N_POSES = 9
SEED = 42

with open(LIBRARY, newline="") as f:
    meta = {r["chembl_id"]: r for r in csv.DictReader(f)}

ligands = sorted(meta)
print(f"Full screen: {len(ligands)} ligands, exhaustiveness {EXHAUSTIVENESS}",
      flush=True)

t_start = time.time()
done = 0
skipped = 0
failed = 0

for i, cid in enumerate(ligands, start=1):
    name = meta[cid]["pref_name"]
    lig = os.path.join(LIGAND_DIR, cid + ".pdbqt")
    energy_file = os.path.join(OUT_DIR, cid + "_energies.txt")

    if os.path.exists(energy_file) == True:
        skipped = skipped + 1
        continue
    if os.path.exists(lig) == False:
        failed = failed + 1
        continue

    t0 = time.time()
    try:
        v = Vina(sf_name="vina", verbosity=0, seed=SEED)
        v.set_receptor(RECEPTOR)
        v.set_ligand_from_file(lig)
        v.compute_vina_maps(center=CENTER, box_size=SIZE)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        best = energies[0][0]
        with open(energy_file, "w") as f:
            for row in energies:
                f.write(" ".join(f"{x:.3f}" for x in row) + "\n")
        v.write_poses(os.path.join(OUT_DIR, cid + "_out.pdbqt"),
                      n_poses=1, overwrite=True)
        done = done + 1
    except Exception as err:
        failed = failed + 1
        with open(os.path.join(BASE, "logs", "dock_failures.txt"), "a") as f:
            f.write(f"{cid}\t{name}\t{str(err)[:200]}\n")
        continue

    if done % 25 == 0:
        spent = time.time() - t_start
        rate = spent / done
        left = (len(ligands) - i) * rate / 3600
        print(f"  [{i}/{len(ligands)}] {done} docked, {failed} failed, "
              f"{rate:.1f}s each, about {left:.1f} h left", flush=True)

# Collect everything that has an energies file into one ranked table.
rows = []
for cid in ligands:
    energy_file = os.path.join(OUT_DIR, cid + "_energies.txt")
    if os.path.exists(energy_file) == False:
        continue
    with open(energy_file) as f:
        first = f.readline().split()
    if len(first) == 0:
        continue
    rows.append({
        "chembl_id": cid,
        "pref_name": meta[cid]["pref_name"],
        "best_affinity_kcal_per_mol": float(first[0]),
        "full_mwt": meta[cid]["full_mwt"],
    })

rows.sort(key=lambda r: r["best_affinity_kcal_per_mol"])
with open(RESULTS, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["rank", "chembl_id", "pref_name",
                                      "best_affinity_kcal_per_mol", "full_mwt"])
    w.writeheader()
    for n, r in enumerate(rows, start=1):
        r["rank"] = n
        w.writerow(r)

print(f"docked {done}, skipped {skipped}, failed {failed}", flush=True)
print(f"ranked {len(rows)} compounds into {RESULTS}", flush=True)
print(f"total wall time {(time.time() - t_start) / 3600:.2f} h", flush=True)
