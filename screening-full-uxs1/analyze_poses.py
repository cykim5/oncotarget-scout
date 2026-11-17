"""Where does each docked pose actually sit?

A good Vina score is not the same as binding in the right place. This
checks every re-docked pose against three physical references:

  1. the crystallographic UDP position, taken from 4GLL, which is where
     the real substrate sits
  2. the catalytic triad, Thr202, Tyr231 and Lys235 in 4GLL numbering
  3. the NAD cofactor held in the holo receptor

A compound that scores well but sits far from all three is occupying
some other groove on the protein and is not a substrate site binder,
whatever its score says.

Writes outputs/pose_analysis.csv
"""
import csv
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))
REDOCK = os.path.join(BASE, "outputs", "redock")
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
CRYSTAL_UDP = os.path.join(BASE, "..", "docking-uxs1", "receptor",
                           "crystal_UDP900_in4GLL.pdb")

CENTER = [60.0, 23.2, 39.0]
SIZE = [24, 24, 24]
TRIAD = {"THR": 202, "TYR": 231, "LYS": 235}
CONTACT = 4.0          # angstroms, a generous contact distance


def read_atoms(path, element_column=True):
    """(x, y, z, resname, resseq) for ATOM and HETATM lines."""
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) == False:
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        resname = line[17:20].strip()
        try:
            resseq = int(line[22:26])
        except ValueError:
            resseq = -1
        out.append((x, y, z, resname, resseq))
    return out


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms)/n,
            sum(a[1] for a in atoms)/n,
            sum(a[2] for a in atoms)/n)


def min_dist(pose, target):
    best = 1e9
    for p in pose:
        for t in target:
            d = dist(p, t)
            if d < best:
                best = d
    return best


udp = read_atoms(CRYSTAL_UDP)
udp_centre = centroid(udp)

receptor = read_atoms(RECEPTOR)
nad = [a for a in receptor if a[3] == "NAD"]
triad = {}
for name, seq in TRIAD.items():
    triad[f"{name}{seq}"] = [a for a in receptor if a[3] == name and a[4] == seq]

lo = [CENTER[i] - SIZE[i]/2 for i in range(3)]
hi = [CENTER[i] + SIZE[i]/2 for i in range(3)]

redock = {}
results_redock = os.path.join(BASE, "outputs", "results_redock.csv")
if os.path.exists(results_redock) == True:
    for r in csv.DictReader(open(results_redock)):
        redock[r["chembl_id"]] = r

rows = []
for fn in sorted(os.listdir(REDOCK)):
    if fn.endswith("_out.pdbqt") == False:
        continue
    cid = fn.replace("_out.pdbqt", "")
    pose = read_atoms(os.path.join(REDOCK, fn))
    if len(pose) == 0:
        continue
    # first MODEL only: Vina writes the best pose first
    pc = centroid(pose)
    near_udp = sum(1 for p in pose
                   if min(dist(p, u) for u in udp) <= CONTACT)
    wall = min(min(p[i] - lo[i], hi[i] - p[i]) for p in pose for i in range(3))
    row = {
        "chembl_id": cid,
        "pref_name": redock.get(cid, {}).get("pref_name", ""),
        "redock_exh8": redock.get(cid, {}).get("redock_exh8", ""),
        "centroid_to_crystal_udp_A": round(dist(pc, udp_centre), 2),
        "pct_atoms_overlapping_udp_site": round(100.0 * near_udp / len(pose), 1),
        "min_dist_to_NAD_A": round(min_dist(pose, nad), 2) if nad else "",
        "closest_box_wall_A": round(wall, 2),
    }
    for label, atoms in triad.items():
        row[f"min_dist_{label}_A"] = round(min_dist(pose, atoms), 2) if atoms else ""
    rows.append(row)

rows.sort(key=lambda r: float(r["redock_exh8"]) if r["redock_exh8"] else 0)
out = os.path.join(BASE, "outputs", "pose_analysis.csv")
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"crystal UDP centroid: {tuple(round(v,1) for v in udp_centre)}")
print(f"box centre {CENTER}, so the true substrate site is "
      f"{dist(udp_centre, CENTER):.1f} A off centre")
print(f"analysed {len(rows)} poses, written to {out}")
print()
hdr = f"  {'compound':<24}{'score':>8}{'d(UDP)':>9}{'%overlap':>10}{'d(NAD)':>8}{'wall':>7}"
print(hdr)
for r in rows[:15]:
    print(f"  {r['pref_name'][:23]:<24}{r['redock_exh8']:>8}"
          f"{r['centroid_to_crystal_udp_A']:>9}{r['pct_atoms_overlapping_udp_site']:>10}"
          f"{r['min_dist_to_NAD_A']:>8}{r['closest_box_wall_A']:>7}")
