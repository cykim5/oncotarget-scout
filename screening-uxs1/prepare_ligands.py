"""Prepare screening ligands: SMILES -> RDKit 3D conformer -> Meeko PDBQT.

Reads library.csv (ChEMBL approved drugs). For each compound:
- parse SMILES with RDKit, keep the largest fragment (desalt),
- add hydrogens, embed one 3D conformer with ETKDGv3 (fixed seed),
- MMFF-optimize (UFF fallback, then unoptimized fallback),
- write ligands/<chembl_id>.sdf, then run Meeko mk_prepare_ligand.py
  (Gasteiger charges, same as the docking-uxs1 run) to make the PDBQT.

Compounds that fail parsing, embedding, or Meeko are recorded in
logs/prep_failures.txt and skipped, never invented or substituted.
"""
import csv
import os
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

BASE = os.path.dirname(os.path.abspath(__file__))
LIBRARY = os.path.join(BASE, "library.csv")
LIGAND_DIR = os.path.join(BASE, "ligands")
LOG_DIR = os.path.join(BASE, "logs")
os.makedirs(LIGAND_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

VENV_BIN = "/tmp/docking-venv/venv/bin"
MK_PREPARE = os.path.join(VENV_BIN, "mk_prepare_ligand.py")

SEED = 42


def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    return max(frags, key=lambda f: f.GetNumHeavyAtoms())


def main():
    with open(LIBRARY, newline="") as f:
        rows = list(csv.DictReader(f))

    failures = []
    prepared = []
    for row in rows:
        cid = row["chembl_id"]
        name = row["pref_name"]
        sdf = os.path.join(LIGAND_DIR, f"{cid}.sdf")
        pdbqt = os.path.join(LIGAND_DIR, f"{cid}.pdbqt")
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            failures.append((cid, name, "RDKit SMILES parse failed"))
            continue
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
            failures.append((cid, name, "ETKDG embedding failed"))
            continue
        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is not None:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass  # keep the unoptimized conformer
        mol.SetProp("_Name", cid)
        with Chem.SDWriter(sdf) as w:
            w.write(mol)
        r = subprocess.run(
            [MK_PREPARE, "-i", sdf, "-o", pdbqt],
            capture_output=True, text=True, timeout=300,
        )
        if r.returncode != 0 or not os.path.exists(pdbqt):
            failures.append((cid, name, f"Meeko failed: {(r.stderr or '').strip()[:200]}"))
            continue
        prepared.append(cid)
        print(f"  prepared {cid} ({name})", flush=True)

    with open(os.path.join(LOG_DIR, "prep_failures.txt"), "w") as f:
        for cid, name, reason in failures:
            f.write(f"{cid}\t{name}\t{reason}\n")

    print(f"\nPrepared {len(prepared)} / {len(rows)} ligands.", flush=True)
    print(f"Failures: {len(failures)} (see logs/prep_failures.txt)", flush=True)
    return 0 if prepared else 1


if __name__ == "__main__":
    sys.exit(main())
