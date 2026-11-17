"""Were the natural ligands handicapped by being docked uncharged?

Every prepared natural ligand sums to a net partial charge of about
zero, which means the phosphates and the carboxylate were modelled
protonated. At physiological pH UDP-glucuronate carries about -3. The
drugs it was ranked against are genuinely neutral, so the control may
have been scored in a protonation state it never occupies.

This deprotonates the acidic groups in place, keeping the 3D geometry
that was already docked, re-prepares with Meeko, and re-docks the five
natural ligands in the same corrected box with both scoring functions.
The decoy distributions from benchmark_scoring.py are reused, so the
percentiles are directly comparable.

Writes outputs/ionization_check.csv
"""
import csv
import os
import statistics

from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "receptor", "uxs1_holo.pdbqt")
NAT_DIR = os.path.join(BASE, "..", "docking-uxs1", "ligands")
OUT_DIR = os.path.join(BASE, "outputs", "ionized")
BENCH = os.path.join(BASE, "outputs", "scoring_benchmark.csv")
os.makedirs(OUT_DIR, exist_ok=True)

CENTER = [62.3, 16.4, 34.2]
SIZE = [28, 28, 28]
EXHAUSTIVENESS = 8
SEED = 42
NATURALS = ["UDP-glucuronate", "UDP-xylose", "UDP", "UDP-galactose", "UMP"]

# acidic OH groups that are deprotonated at pH 7.4
ACIDIC = [
    Chem.MolFromSmarts("[OX2H1]-[PX4]"),        # phosphate / diphosphate OH
    Chem.MolFromSmarts("[OX2H1]-[CX3]=[OX1]"),  # carboxylic acid OH
]


def deprotonate(mol):
    """Set the acidic oxygens to -1 and drop their hydrogens, in place."""
    rw = Chem.RWMol(mol)
    targets = set()
    for patt in ACIDIC:
        for match in rw.GetSubstructMatches(patt):
            targets.add(match[0])
    drop = []
    for idx in sorted(targets):
        atom = rw.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "H":
                drop.append(nb.GetIdx())
        atom.SetFormalCharge(-1)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out, len(targets)


def prepare(name):
    sdf = os.path.join(NAT_DIR, name + ".sdf")
    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        return None, 0, "RDKit could not read the SDF"
    ion, n = deprotonate(mol)
    charge = Chem.GetFormalCharge(ion)
    try:
        setups = MoleculePreparation().prepare(ion)
    except Exception as err:
        return None, charge, f"Meeko prepare failed: {str(err)[:100]}"
    if len(setups) == 0:
        return None, charge, "Meeko produced no setup"
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if ok == False:
        return None, charge, f"Meeko write failed: {str(err)[:100]}"
    path = os.path.join(OUT_DIR, name + "_ion.pdbqt")
    with open(path, "w") as f:
        f.write(text)
    return path, charge, ""


def dock(path, sf):
    v = Vina(sf_name=sf, verbosity=0, seed=SEED)
    v.set_receptor(RECEPTOR)
    v.set_ligand_from_file(path)
    v.compute_vina_maps(center=CENTER, box_size=SIZE)
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
    return v.energies(n_poses=1)[0][0]


# reuse the decoy distributions already measured
decoys = {"vina": [], "vinardo": []}
neutral = {}
with open(BENCH) as f:
    for row in csv.DictReader(f):
        if row["kind"] == "decoy":
            decoys[row["sf"]].append(float(row["score"]))
        else:
            neutral[(row["id"], row["sf"])] = float(row["score"])

rows = []
for name in NATURALS:
    path, charge, err = prepare(name)
    if path is None:
        print(f"  {name}: {err}")
        continue
    print(f"  {name}: formal charge {charge:+d}")
    for sf in ("vina", "vinardo"):
        score = dock(path, sf)
        rows.append({
            "id": name,
            "sf": sf,
            "formal_charge": charge,
            "score_neutral": neutral.get((name, sf), ""),
            "score_ionized": round(score, 3),
        })

with open(os.path.join(BASE, "outputs", "ionization_check.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("\nNEUTRAL VERSUS IONIZED, SAME BOX, SAME DECOY SETS\n")
for sf in ("vina", "vinardo"):
    pool = sorted(decoys[sf])
    n = len(pool)
    print(f"  {sf}  ({n} decoys, median {statistics.median(pool):.2f})")
    print(f"    {'ligand':<20} {'neutral':>9} {'pctile':>7}   {'ionized':>9} {'pctile':>7}   {'shift':>7}")
    for r in rows:
        if r["sf"] != sf:
            continue
        pn = 100.0 * sum(1 for d in pool if d > r["score_neutral"]) / n
        pi = 100.0 * sum(1 for d in pool if d > r["score_ionized"]) / n
        shift = r["score_ionized"] - r["score_neutral"]
        print(f"    {r['id']:<20} {r['score_neutral']:>9.3f} {pn:>6.1f}%   "
              f"{r['score_ionized']:>9.3f} {pi:>6.1f}%   {shift:>+7.3f}")
    print()
