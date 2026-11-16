"""Virtual screen, pass 2: re-dock the top 10 pass-1 hits at exhaustiveness 8.

Reads outputs/pass1/results_pass1.csv, takes the top N_TOP ligands, and
re-docks them into the same holo substrate-site setup with
exhaustiveness 8 (the standard setting from the docking-uxs1 run).
Outputs: outputs/top10/<chembl_id>_out.pdbqt, *_energies.txt,
results_top10.csv (ranked by best affinity at exhaustiveness 8).
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
LIGAND_DIR = os.path.join(BASE, "ligands")
PASS1 = os.path.join(BASE, "outputs", "pass1", "results_pass1.csv")
OUT_DIR = os.path.join(BASE, "outputs", "top10")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 8
N_POSES = 9
N_TOP = 10

with open(PASS1, newline="") as f:
    rows = list(csv.DictReader(f))
top = rows[:N_TOP]
print(f"Pass 2: re-docking top {len(top)} at exhaustiveness {EXHAUSTIVENESS}", flush=True)

results = []
t_start = time.time()
for i, row in enumerate(top, start=1):
    cid = row["chembl_id"]
    name = row["pref_name"]
    lig = os.path.join(LIGAND_DIR, f"{cid}.pdbqt")
    out_pdbqt = os.path.join(OUT_DIR, f"{cid}_out.pdbqt")
    energy_file = os.path.join(OUT_DIR, f"{cid}_energies.txt")
    print(f"[{i}/{len(top)}] {cid} ({name}) ...", flush=True)
    t0 = time.time()
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(RECEPTOR)
    v.set_ligand_from_file(lig)
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    energies = v.energies(n_poses=N_POSES)
    v.write_poses(out_pdbqt, n_poses=N_POSES, energy_range=3.0, overwrite=True)
    with open(energy_file, "w") as f:
        f.write(f"ligand: {cid} ({name})\n")
        f.write(f"receptor: {RECEPTOR}\n")
        f.write(f"center: {CENTER}\nbox_size: {SIZE}\n")
        f.write(f"exhaustiveness: {EXHAUSTIVENESS}\nn_poses: {N_POSES}\n")
        f.write("mode | affinity = total (kcal/mol) | intermolecular (kcal/mol) | intramolecular (kcal/mol) | torsional (kcal/mol) | intra best pose (kcal/mol)\n")
        for j, erow in enumerate(energies, start=1):
            f.write(f"{j} | {erow[0]:.3f} | {erow[1]:.3f} | {erow[2]:.3f} | {erow[3]:.3f} | {erow[4]:.3f}\n")
    best = energies[0][0]
    rt = time.time() - t0
    results.append((cid, name, best, float(row["best_affinity_kcal_per_mol"]), rt))
    print(f"  best affinity: {best:.3f} kcal/mol (pass 1: {float(row['best_affinity_kcal_per_mol']):.3f}, {rt:.0f}s)", flush=True)

results.sort(key=lambda r: r[2])
with open(os.path.join(OUT_DIR, "results_top10.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "chembl_id", "pref_name", "best_affinity_exh8_kcal_per_mol",
                "best_affinity_pass1_exh4_kcal_per_mol", "full_mwt", "runtime_s"])
    for i, (cid, name, best, p1, rt) in enumerate(results, start=1):
        mw = next(r["full_mwt"] for r in top if r["chembl_id"] == cid)
        w.writerow([i, cid, name, f"{best:.3f}", f"{p1:.3f}", mw, f"{rt:.0f}"])

print(f"\nPass 2 complete in {time.time()-t_start:.0f}s", flush=True)
print("Top 10 (exhaustiveness 8):", flush=True)
for i, (cid, name, best, p1, rt) in enumerate(results, start=1):
    print(f"  {i}. {name} ({cid}): {best:.3f}", flush=True)
