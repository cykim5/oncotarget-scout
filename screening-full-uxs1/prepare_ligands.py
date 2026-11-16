"""Build a 3D PDBQT for every compound in library.csv.

Same recipe as the 100 compound screen: largest fragment, explicit
hydrogens, one ETKDGv3 conformer at a fixed seed, force field clean up,
then Meeko. The difference is that Meeko is called through its Python
API rather than a command line script, so there is no hardcoded path to
anyone's virtual environment.

Resumable: a compound with an existing PDBQT is skipped.
"""
import csv
import os
import sys
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
LIBRARY = os.path.join(BASE, "library.csv")
LIGAND_DIR = os.path.join(BASE, "ligands")
LOG_DIR = os.path.join(BASE, "logs")
os.makedirs(LIGAND_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

SEED = 42


def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    return max(frags, key=lambda f: f.GetNumHeavyAtoms())


def build(row):
    """Return (ok, reason). Writes <cid>.sdf and <cid>.pdbqt."""
    cid = row["chembl_id"]
    sdf = os.path.join(LIGAND_DIR, cid + ".sdf")
    pdbqt = os.path.join(LIGAND_DIR, cid + ".pdbqt")
    if os.path.exists(pdbqt) == True:
        return True, "already prepared"

    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        return False, "RDKit SMILES parse failed"
    mol = largest_fragment(mol)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.numThreads = 1
    conf_id = AllChem.EmbedMolecule(mol, params)
    if conf_id == -1:
        params.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(mol, params)
    if conf_id == -1:
        return False, "ETKDG embedding failed"

    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass  # an unoptimised conformer is still usable

    mol.SetProp("_Name", cid)
    with Chem.SDWriter(sdf) as w:
        w.write(mol)

    try:
        setups = MoleculePreparation().prepare(mol)
    except Exception as err:
        return False, f"Meeko prepare failed: {str(err)[:120]}"
    if len(setups) == 0:
        return False, "Meeko produced no setup"
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if ok == False:
        return False, f"Meeko write failed: {str(err)[:120]}"
    with open(pdbqt, "w") as f:
        f.write(text)
    return True, "prepared"


with open(LIBRARY, newline="") as f:
    rows = list(csv.DictReader(f))

print(f"Preparing {len(rows)} ligands", flush=True)
failures = []
done = 0
t0 = time.time()
for i, row in enumerate(rows, start=1):
    ok, reason = build(row)
    if ok == True:
        done = done + 1
    else:
        failures.append((row["chembl_id"], row["pref_name"], reason))
    if i % 100 == 0:
        rate = (time.time() - t0) / i
        left = (len(rows) - i) * rate / 60
        print(f"  [{i}/{len(rows)}] ok {done}, failed {len(failures)}, "
              f"about {left:.0f} min left", flush=True)

with open(os.path.join(LOG_DIR, "prep_failures.txt"), "w") as f:
    for cid, name, reason in failures:
        f.write(f"{cid}\t{name}\t{reason}\n")

print(f"Prepared {done} / {len(rows)}. Failures: {len(failures)}", flush=True)
