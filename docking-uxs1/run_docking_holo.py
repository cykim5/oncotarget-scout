"""Follow-up docking: UDP ligands into the NAD-bound (holo) UXS1 receptor.

Why: in the apo run, every UDP ligand fell into the deep NAD pocket instead
of the crystallographic UDP site (2B69 UDP 900). In the real enzyme NAD is
bound when the substrate binds, so the holo receptor (NAD kept as rigid
receptor atoms) is the correct model for the substrate site.

Receptor: uxs1_holo.pdbqt = Meeko apo receptor + crystal NAD atoms
(from OpenBabel conversion of 4GLL chain A + NAD).
Box: center (60.0, 23.2, 39.0) = midpoint of the 2B69 UDP 900 centroid
(62.27, 16.41, 34.18, transformed into the 4GLL frame) and the catalytic
triad centroid (57.66, 29.97, 43.73). Size 24 x 24 x 24 A.
Ligands: the 5 UDP ligands only (NAD's site is occupied by the receptor NAD).
Outputs: outputs/holo/<ligand>_out.pdbqt, *_energies.txt, results_holo.csv.
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
LIGAND_DIR = os.path.join(BASE, "ligands")
OUT_DIR = os.path.join(BASE, "outputs", "holo")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
EXHAUSTIVENESS = 16
N_POSES = 9

LIGANDS = ["UDP", "UDP-glucuronate", "UDP-xylose", "UDP-galactose", "UMP"]

results = []
for name in LIGANDS:
    lig = os.path.join(LIGAND_DIR, f"{name}.pdbqt")
    out_pdbqt = os.path.join(OUT_DIR, f"{name}_out.pdbqt")
    energy_file = os.path.join(OUT_DIR, f"{name}_energies.txt")
    print(f"=== {name} (holo) ===", flush=True)
    t0 = time.time()
    v = Vina(sf_name="vina", verbosity=1)
    v.set_receptor(RECEPTOR)
    v.set_ligand_from_file(lig)
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    energies = v.energies(n_poses=N_POSES)
    v.write_poses(out_pdbqt, n_poses=N_POSES, energy_range=3.0, overwrite=True)
    with open(energy_file, "w") as f:
        f.write(f"ligand: {name} (holo receptor, NAD retained)\n")
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
with open(os.path.join(OUT_DIR, "results_holo.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ligand", "best_affinity_kcal_per_mol", "runtime_s"])
    for name, best, rt in results:
        w.writerow([name, f"{best:.3f}", f"{rt:.0f}"])

print("\nRanked by best affinity (kcal/mol, lower is better):")
for i, (name, best, rt) in enumerate(results, start=1):
    print(f"  {i}. {name}: {best:.3f}")
