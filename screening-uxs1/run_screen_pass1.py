"""Virtual screen, pass 1: ~100 approved drugs vs the UXS1 substrate site.

Setup (identical to docking-uxs1 run 2, the holo run):
- Receptor: uxs1_holo.pdbqt (4GLL chain A + crystal NAD as rigid atoms)
- Box: center (60.0, 23.2, 39.0), size 24 x 24 x 24 A (substrate site)
- AutoDock Vina 1.2.7, exhaustiveness 4 (reduced for the first pass),
  9 poses per ligand, energy range 3 kcal/mol.

Ligands: ligands/<chembl_id>.pdbqt from prepare_ligands.py (ChEMBL
approved drugs, RDKit 3D conformers, Meeko/Gasteiger).
Outputs: outputs/pass1/<chembl_id>_out.pdbqt, *_energies.txt,
results_pass1.csv (ranked by best affinity).

Resumable: ligands with an existing output PDBQT are skipped.
"""
import csv
import os
import sys
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
LIGAND_DIR = os.path.join(BASE, "ligands")
LIBRARY = os.path.join(BASE, "library.csv")
OUT_DIR = os.path.join(BASE, "outputs", "pass1")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 4
N_POSES = 9

with open(LIBRARY, newline="") as f:
    meta = {r["chembl_id"]: r for r in csv.DictReader(f)}

ligands = sorted(meta)
print(f"Pass 1: {len(ligands)} ligands, exhaustiveness {EXHAUSTIVENESS}", flush=True)

results = []
t_start = time.time()
for i, cid in enumerate(ligands, start=1):
    name = meta[cid]["pref_name"]
    lig = os.path.join(LIGAND_DIR, f"{cid}.pdbqt")
    out_pdbqt = os.path.join(OUT_DIR, f"{cid}_out.pdbqt")
    energy_file = os.path.join(OUT_DIR, f"{cid}_energies.txt")
    if os.path.exists(out_pdbqt) and os.path.exists(energy_file):
        print(f"[{i}/{len(ligands)}] {cid} ({name}): skipped, output exists", flush=True)
        continue
    if not os.path.exists(lig):
        print(f"[{i}/{len(ligands)}] {cid} ({name}): no PDBQT, skipped", flush=True)
        continue
    print(f"[{i}/{len(ligands)}] {cid} ({name}) ...", flush=True)
    t0 = time.time()
    try:
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(RECEPTOR)
        v.set_ligand_from_file(lig)
        v.compute_vina_maps(center=CENTER, box_size=SIZE)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        v.write_poses(out_pdbqt, n_poses=N_POSES, energy_range=3.0, overwrite=True)
    except Exception as e:
        print(f"  FAILED: {e}", flush=True)
        with open(os.path.join(BASE, "logs", "dock_failures_pass1.txt"), "a") as f:
            f.write(f"{cid}\t{name}\t{e}\n")
        continue
    with open(energy_file, "w") as f:
        f.write(f"ligand: {cid} ({name})\n")
        f.write(f"receptor: {RECEPTOR}\n")
        f.write(f"center: {CENTER}\nbox_size: {SIZE}\n")
        f.write(f"exhaustiveness: {EXHAUSTIVENESS}\nn_poses: {N_POSES}\n")
        f.write("mode | affinity = total (kcal/mol) | intermolecular (kcal/mol) | intramolecular (kcal/mol) | torsional (kcal/mol) | intra best pose (kcal/mol)\n")
        for j, row in enumerate(energies, start=1):
            f.write(f"{j} | {row[0]:.3f} | {row[1]:.3f} | {row[2]:.3f} | {row[3]:.3f} | {row[4]:.3f}\n")
    best = energies[0][0]
    rt = time.time() - t0
    results.append((cid, name, best, rt))
    print(f"  best affinity: {best:.3f} kcal/mol ({rt:.0f}s, elapsed {time.time()-t_start:.0f}s)", flush=True)

# Re-read any previously completed results (resumed run) from energies files
for cid in ligands:
    if any(r[0] == cid for r in results):
        continue
    energy_file = os.path.join(OUT_DIR, f"{cid}_energies.txt")
    if not os.path.exists(energy_file):
        continue
    with open(energy_file) as f:
        for line in f:
            if line.startswith("1 |"):
                best = float(line.split("|")[1].strip())
                results.append((cid, meta[cid]["pref_name"], best, float("nan")))
                break

results.sort(key=lambda r: r[2])
with open(os.path.join(OUT_DIR, "results_pass1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "chembl_id", "pref_name", "best_affinity_kcal_per_mol", "full_mwt", "runtime_s"])
    for i, (cid, name, best, rt) in enumerate(results, start=1):
        w.writerow([i, cid, name, f"{best:.3f}", meta[cid]["full_mwt"],
                    f"{rt:.0f}" if rt == rt else ""])

print(f"\nPass 1 complete: {len(results)} ligands docked in {time.time()-t_start:.0f}s", flush=True)
print("Top 10 (pass 1, exhaustiveness 4):", flush=True)
for i, (cid, name, best, rt) in enumerate(results[:10], start=1):
    print(f"  {i}. {name} ({cid}): {best:.3f}", flush=True)
