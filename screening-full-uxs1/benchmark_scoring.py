"""Can any scoring function tell the real substrate from random drugs?

The Vina enrichment test failed badly: UDP-glucuronate, the enzyme's
actual substrate, ranked below two thirds of size matched approved
drugs. Before concluding that docking cannot work on this site, the
obvious thing to rule out is the scoring function itself. Vina's
empirical function is known to handle charged, phosphate bearing
ligands poorly, and every natural ligand here is an anion.

This re-scores a sample in the re-centred box with both Vina and
Vinardo, which weights polar contacts differently, and reports where
the natural ligands land in each.

Writes outputs/scoring_benchmark.csv
"""
import csv
import os
import random
import statistics
import time

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
NAT_DIR = os.path.join(BASE, "..", "docking-uxs1", "ligands")
OUT_DIR = os.path.join(BASE, "outputs", "benchmark")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [62.3, 16.4, 34.2]      # the corrected box
SIZE = [28, 28, 28]
EXHAUSTIVENESS = 8
SEED = 42
N_DECOYS = 200                   # random drugs, size matched to the substrate
HEAVY_LOW, HEAVY_HIGH = 31, 43   # UDP-glucuronate has 37
NATURALS = ["UDP-glucuronate", "UDP-xylose", "UDP", "UDP-galactose", "UMP"]


def heavy_sdf(path):
    L = open(path).read().splitlines()
    n = int(L[3][0:3])
    return sum(1 for l in L[4:4+n] if l[31:34].strip() != "H")


def dock(ligand_path, sf):
    v = Vina(sf_name=sf, verbosity=0, seed=SEED)
    v.set_receptor(RECEPTOR)
    v.set_ligand_from_file(ligand_path)
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
    return v.energies(n_poses=1)[0][0]


# pick the decoy sample once, reproducibly
pool = []
for fn in sorted(os.listdir(os.path.join(BASE, "ligands"))):
    if fn.endswith(".sdf") == False:
        continue
    cid = fn[:-4]
    pdbqt = os.path.join(BASE, "ligands", cid + ".pdbqt")
    if os.path.exists(pdbqt) == False:
        continue
    h = heavy_sdf(os.path.join(BASE, "ligands", fn))
    if HEAVY_LOW <= h <= HEAVY_HIGH:
        pool.append((cid, h))
random.seed(42)
decoys = random.sample(pool, min(N_DECOYS, len(pool)))
print(f"decoy pool {len(pool)}, sampled {len(decoys)}", flush=True)

rows = []
t0 = time.time()
for sf in ("vina", "vinardo"):
    for name in NATURALS:
        p = os.path.join(NAT_DIR, name + ".pdbqt")
        if os.path.exists(p) == False:
            continue
        cache = os.path.join(OUT_DIR, f"nat_{name}_{sf}.txt")
        if os.path.exists(cache) == True:
            score = float(open(cache).read().strip())
        else:
            score = dock(p, sf)
            open(cache, "w").write(f"{score:.3f}")
        rows.append({"kind": "natural", "id": name, "sf": sf,
                     "score": round(score, 3),
                     "heavy": heavy_sdf(os.path.join(NAT_DIR, name + ".sdf"))})
        print(f"  {sf:<8} {name:<18} {score:.3f}", flush=True)

    for i, (cid, h) in enumerate(decoys, 1):
        cache = os.path.join(OUT_DIR, f"dec_{cid}_{sf}.txt")
        if os.path.exists(cache) == True:
            score = float(open(cache).read().strip())
        else:
            try:
                score = dock(os.path.join(BASE, "ligands", cid + ".pdbqt"), sf)
            except Exception:
                continue
            open(cache, "w").write(f"{score:.3f}")
        rows.append({"kind": "decoy", "id": cid, "sf": sf,
                     "score": round(score, 3), "heavy": h})
        if i % 50 == 0:
            print(f"  {sf}: {i}/{len(decoys)} decoys, "
                  f"{(time.time()-t0)/60:.0f} min elapsed", flush=True)

out = os.path.join(BASE, "outputs", "scoring_benchmark.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["kind", "id", "sf", "score", "heavy"])
    w.writeheader()
    for r in rows:
        w.writerow(r)

print("\nPERCENTILE OF EACH NATURAL LIGAND AMONG THE DECOYS")
for sf in ("vina", "vinardo"):
    dec = sorted(r["score"] for r in rows if r["kind"] == "decoy" and r["sf"] == sf)
    print(f"\n  {sf}  ({len(dec)} decoys, median {statistics.median(dec):.2f})")
    for name in NATURALS:
        hit = [r for r in rows if r["kind"] == "natural"
               and r["id"] == name and r["sf"] == sf]
        if len(hit) == 0:
            continue
        s = hit[0]["score"]
        better = sum(1 for d in dec if d < s)
        pct = 100.0 * (len(dec) - better) / len(dec)
        print(f"    {name:<18}{s:>8.3f}   beaten by {better:>3} of {len(dec)}"
              f"   percentile {pct:>5.1f}")
print(f"\nwritten to {out}")
print(f"wall time {(time.time()-t0)/60:.1f} min")
