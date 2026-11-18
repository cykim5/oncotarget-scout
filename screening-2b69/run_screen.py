"""The full approved drug screen, rerun against a receptor that works.

The first screen docked into PDB 4GLL, which has NAD bound and no UDP.
With the substrate site empty, ARG272 had relaxed into it, backbone
included, so every compound was scored against a partly blocked
pocket. The control failed and the screen was abandoned.

Against PDB 2B69, human UXS1 at 1.21 A with NAD and UDP both bound,
the same control goes from FAIL to WEAK: the natural ligands move from
the 30th, 55th and 58th percentiles to the 99.5th, 99.0th and 88.4th,
and the correlation between score and molecular size collapses from
rho -0.67 to +0.11. See docs/control-2b69-result.md.

So the protocol is worth pointing at the library again.

Exhaustiveness is 8, not the 4 the first screen used, because 8 is
what the control validated. Screening with a protocol other than the
one that passed its control is how the first screen went wrong.

Resumable: a compound with an existing energies file is skipped.
Compounds Vina cannot parse, boron drugs in particular, are recorded
and skipped rather than ending the run.

Writes outputs/results.csv and outputs/failures.tsv
"""
import csv
import os
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "2b69_holo.pdbqt")
LIGAND_DIR = os.path.join(os.path.dirname(BASE), "screening-full-uxs1", "ligands")
LIBRARY = os.path.join(os.path.dirname(BASE), "screening-full-uxs1", "library.csv")
OUT_DIR = os.path.join(BASE, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [50.48, 21.14, 23.36]   # the UDP 900 centroid, in 2B69's own frame
SIZE = [28, 28, 28]
EXHAUSTIVENESS = 8               # what the control validated
SEED = 42
N_POSES = 9

names = {}
if os.path.exists(LIBRARY) == True:
    for row in csv.DictReader(open(LIBRARY)):
        names[row.get("id") or row.get("chembl_id")] = row.get("name", "")

ligands = sorted(f for f in os.listdir(LIGAND_DIR) if f.endswith(".pdbqt"))
print(f"screen: {len(ligands)} ligands, exhaustiveness {EXHAUSTIVENESS}, "
      f"box {SIZE[0]} A at {CENTER}", flush=True)

results_path = os.path.join(OUT_DIR, "results.csv")
fail_path = os.path.join(OUT_DIR, "failures.tsv")
done_ids = set()
if os.path.exists(results_path) == True:
    for row in csv.DictReader(open(results_path)):
        done_ids.add(row["id"])
    print(f"  resuming: {len(done_ids)} already scored", flush=True)

new_file = os.path.exists(results_path) == False
out = open(results_path, "a", newline="")
w = csv.DictWriter(out, fieldnames=["id", "name", "score", "heavy"])
if new_file == True:
    w.writeheader()
fails = open(fail_path, "a")

started = time.time()
done = len(done_ids)
failed = 0
for i, fn in enumerate(ligands):
    cid = fn[:-6]
    if cid in done_ids:
        continue
    path = os.path.join(LIGAND_DIR, fn)
    try:
        v = Vina(sf_name="vina", verbosity=0, seed=SEED)
        v.set_receptor(RECEPTOR)
        v.set_ligand_from_file(path)
        v.compute_vina_maps(center=CENTER, box_size=SIZE)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        score = float(v.energies(n_poses=1)[0][0])
    except Exception as err:
        failed = failed + 1
        fails.write(f"{cid}\t{str(err).strip().splitlines()[-1][:140]}\n")
        fails.flush()
        continue
    heavy = 0
    for l in open(path):
        if l[:6].strip() in ("ATOM", "HETATM"):
            t = l.split()[-1].upper()
            if t not in ("H", "HD", "HS"):
                heavy = heavy + 1
    w.writerow({"id": cid, "name": names.get(cid, ""),
                "score": round(score, 3), "heavy": heavy})
    out.flush()
    done = done + 1
    if done % 100 == 0:
        mins = (time.time() - started) / 60.0
        rate = (done - len(done_ids)) / max(mins, 0.01)
        left = (len(ligands) - done) / max(rate, 0.01)
        print(f"  {done}/{len(ligands)} scored, {failed} failed, "
              f"{mins:.0f} min elapsed, about {left:.0f} min left", flush=True)

out.close()
fails.close()
print(f"done: {done} scored, {failed} failed, "
      f"{(time.time() - started) / 60:.1f} minutes", flush=True)
