"""Dock candidate ligands into the UXS1 (PDB 4GLL, chain A) active site with AutoDock Vina.

Box centered on the midpoint between the catalytic triad centroid
(Thr202/Tyr231/Lys235 in 4GLL numbering; Thr118/Tyr147/Lys151 in the
JBC paper's construct numbering, offset +84) and the co-crystallized NAD
centroid. Size 24 x 24 x 24 A covers the NAD site and the UDP-sugar site.

Run: /tmp/docking-venv/venv/bin/python run_docking.py
Outputs: outputs/<ligand>_out.pdbqt (poses), outputs/<ligand>_energies.txt,
         outputs/results.csv (best affinity per ligand).
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_chainA.pdbqt")
LIGAND_DIR = os.path.join(BASE, "ligands")
OUT_DIR = os.path.join(BASE, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [53.3, 27.7, 42.3]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 16
N_POSES = 9

LIGANDS = ["NAD", "UDP", "UDP-glucuronate", "UDP-xylose", "UDP-galactose", "UMP"]

results = []
for name in LIGANDS:
    lig = os.path.join(LIGAND_DIR, f"{name}.pdbqt")
    out_pdbqt = os.path.join(OUT_DIR, f"{name}_out.pdbqt")
    energy_file = os.path.join(OUT_DIR, f"{name}_energies.txt")
    print(f"=== {name} ===", flush=True)
    t0 = time.time()
    v = Vina(sf_name="vina", verbosity=1)
    v.set_receptor(RECEPTOR)
    v.set_ligand_from_file(lig)
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    energies = v.energies(n_poses=N_POSES)
    v.write_poses(out_pdbqt, n_poses=N_POSES, energy_range=3.0, overwrite=True)
    with open(energy_file, "w") as f:
        f.write(f"ligand: {name}\n")
        f.write(f"receptor: {RECEPTOR}\n")
        f.write(f"center: {CENTER}\nbox_size: {SIZE}\n")
        f.write(f"exhaustiveness: {EXHAUSTIVENESS}\nn_poses: {N_POSES}\n")
        f.write("mode | affinity = total (kcal/mol) | intermolecular (kcal/mol) | intramolecular (kcal/mol) | torsional (kcal/mol) | intra best pose (kcal/mol)\n")
        for i, row in enumerate(energies, start=1):
            f.write(f"{i} | {row[0]:.3f} | {row[1]:.3f} | {row[2]:.3f} | {row[3]:.3f} | {row[4]:.3f}\n")
    best = energies[0][0]
    results.append((name, best, time.time() - t0))
    print(f"  best affinity: {best:.3f} kcal/mol ({time.time()-t0:.0f}s)", flush=True)

results.sort(key=lambda r: r[1])
with open(os.path.join(OUT_DIR, "results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ligand", "best_affinity_kcal_per_mol", "runtime_s"])
    for name, best, rt in results:
        w.writerow([name, f"{best:.3f}", f"{rt:.0f}"])

print("\nRanked by best affinity (kcal/mol, lower is better):")
for i, (name, best, rt) in enumerate(results, start=1):
    print(f"  {i}. {name}: {best:.3f}")
