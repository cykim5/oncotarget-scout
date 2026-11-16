"""Quantitative check of a docked pose vs the crystal NAD reference.

Reports: centroid distance between docked pose and crystal NAD, and the
fraction of docked heavy atoms within 2 A of any crystal NAD heavy atom.
Usage: venv python check_pose.py <ligand> [model_number, default 1]
"""
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


def read_atoms(path, model=1):
    atoms = []
    cur = 0
    in_model = model is None
    for ln in open(path):
        if ln.startswith("MODEL"):
            cur = int(ln.split()[1])
            in_model = (cur == model)
            continue
        if ln.startswith("ENDMDL"):
            in_model = False
            continue
        if in_model and ln.startswith(("ATOM", "HETATM")):
            el = ln[76:78].strip() or ln[12:16].strip()[0]
            if el == "H":
                continue
            atoms.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return atoms


lig = sys.argv[1]
model = int(sys.argv[2]) if len(sys.argv) > 2 else 1
docked = read_atoms(os.path.join(BASE, "outputs", f"{lig}_out.pdbqt"), model)
crystal = read_atoms(os.path.join(BASE, "receptor", "crystal_NAD_chainA.pdb"), None)

def centroid(pts):
    n = len(pts)
    return tuple(sum(p[i] for p in pts) / n for i in range(3))

dc, cc = centroid(docked), centroid(crystal)
print(f"{lig} model {model}: {len(docked)} heavy atoms")
print(f"  centroid distance to crystal NAD: {math.dist(dc, cc):.2f} A")
near = sum(1 for d in docked if min(math.dist(d, c) for c in crystal) < 2.0)
print(f"  heavy atoms within 2 A of a crystal NAD atom: {near}/{len(docked)} ({100*near/len(docked):.0f}%)")
# distance from docked centroid to catalytic triad centroid
triad = (57.66, 29.97, 43.73)
print(f"  centroid distance to catalytic triad: {math.dist(dc, triad):.2f} A")
