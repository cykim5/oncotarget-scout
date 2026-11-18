"""Could this molecule chemically bind this site at all?

The substrate site of UXS1 exists to grip a pyrophosphate. Four
arginines, two asparagines, a threonine, a tyrosine and a glutamine
line it, and the natural ligands all carry two or three negative
charges at pH 7.4. A site like that pays for charge pairing, not for
bulk hydrophobic surface.

Vina has no electrostatics term. It has a hydrogen bond term and a
hydrophobic term, and the hydrophobic term is the easier one to earn.
So a large greasy molecule can collect a good score in a polar site by
filling it, without forming a single interaction the real pocket is
built for. That artifact is well known and it is the main thing that
separates a docking hit list from a shortlist worth testing.

This script asks, for each shortlisted pose:

  - what is the molecule's formal charge at pH 7.4
  - does it carry any anionic group at all
  - how lipophilic is it
  - does the pose actually contact the phosphate gripping arginines

Writes outputs/plausibility.csv and outputs/plausibility.md
"""
import csv
import math
import os

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
APO = os.path.join(os.path.dirname(BASE), "control-2b69", "receptor", "2b69_chainA_NAD.pdb")
POSE_DIR = os.path.join(BASE, "outputs", "poses")
POSE_CHECK = os.path.join(BASE, "outputs", "pose_check.csv")
OUT_CSV = os.path.join(BASE, "outputs", "plausibility.csv")
OUT_MD = os.path.join(BASE, "outputs", "plausibility.md")

# the residues that hold the pyrophosphate, read off the crystal contacts
PHOSPHATE_GRIP = {("ARG", "144"), ("ARG", "188"), ("ARG", "213"), ("ARG", "277")}
CONTACT = 4.5

ANIONIC = [
    ("carboxylate", "[CX3](=O)[OX1H0-,OX2H1]"),
    ("phosphate or phosphonate", "[PX4](=O)([OX1H0-,OX2H1])[OX1H0-,OX2H1]"),
    ("sulfonate or sulfate", "[SX4](=O)(=O)[OX1H0-,OX2H1]"),
    ("tetrazole", "c1nnn[nH1,n-]1"),
    ("acyl sulfonamide", "[CX3](=O)[NX3H1][SX4](=O)(=O)"),
]


def read_atoms(path):
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) == False:
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        out.append((x, y, z, line[17:20].strip(), line[22:26].strip()))
    return out


def pose_atoms(path):
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            try:
                out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
        if line.startswith("ENDMDL"):
            break
    return out


def charge_at_ph74(mol):
    """A simple count: strong acids lose a proton, basic amines gain one.

    This is a rule of thumb, not a pKa calculation. It is here to
    separate polyanions from cations, which is a difference large
    enough that a rough rule gets it right.
    """
    neg = 0
    for _, smarts in ANIONIC:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            neg += len(mol.GetSubstructMatches(patt))
    basic = Chem.MolFromSmarts("[NX3;H0,H1,H2;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[#16](=O)=O)]")
    pos = 0
    if basic is not None:
        pos = len(mol.GetSubstructMatches(basic))
    pos += len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX4+]")))
    return -neg + pos, neg, pos


def pains_catalog():
    """Known frequent hitters, the PAINS A/B/C filter sets.

    An anionic group in the right place is necessary but not
    sufficient. Several of the compounds that pass the charge filter
    are dyes and quinones, which hit in every assay they are put in,
    and a screen that reports them has learned nothing.
    """
    params = FilterCatalog.FilterCatalogParams()
    for cat in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
                FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
                FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C):
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


def main():
    rows = list(csv.DictReader(open(POSE_CHECK)))
    pains = pains_catalog()
    rec = read_atoms(APO)
    grip = [a for a in rec if (a[3], a[4]) in PHOSPHATE_GRIP]
    print(f"{len(rows)} poses, {len(grip)} atoms in the phosphate grip", flush=True)

    out = []
    for r in rows:
        smiles = r.get("smiles", "")
        if smiles == "":
            smiles = None
        mol = None
        if smiles is not None:
            mol = Chem.MolFromSmiles(smiles)
        # Everything below has to be asked of the PARENT, not the salt.
        # Two bugs came from not doing this. Paroxetine mesylate looked
        # like it carried a sulfonate: that is the mesylate counterion,
        # which dissociates in solution and is not part of the drug.
        # And fluorescein sodium looked inorganic because of its sodium,
        # when the sodium is just the counterion of a fine organic acid.
        # The docked species is the parent, so the parent is what gets
        # judged.
        if mol is not None:
            try:
                mol = rdMolStandardize.FragmentParent(mol)
                # FragmentParent hands back a mol with no ring info, so
                # anything asking about rings throws until it is
                # sanitised again.
                Chem.SanitizeMol(mol)
            except Exception:
                pass
        if mol is None:
            out.append({**r, "charge": "", "anionic_groups": "", "basic_groups": "",
                        "clogp": "", "grip_contacts": "", "pains": "",
                        "verdict": "no structure"})
            continue

        net, neg, pos = charge_at_ph74(mol)
        clogp = round(Crippen.MolLogP(mol), 2)

        pose_path = os.path.join(POSE_DIR, r["id"] + ".pdbqt")
        n_grip = 0
        if os.path.exists(pose_path) == True:
            atoms = pose_atoms(pose_path)
            touched = set()
            for a in grip:
                if min((math.dist(a[:3], p) for p in atoms), default=99) <= CONTACT:
                    touched.add((a[3], a[4]))
            n_grip = len(touched)

        reasons = []

        # A filter that only asks "is it anionic" passes mineral salts.
        # Zinc sulfate and dibasic potassium phosphate both cleared the
        # charge test on the first run, because an anion is an anion.
        # They are excipients, not candidate binders, and they only
        # reached the shortlist because ligand efficiency divides by
        # heavy atom count and so rewards being tiny.
        hard = []
        metals = {"Na", "K", "Ca", "Mg", "Zn", "Fe", "Al", "Li", "Ba", "Sr"}
        has_metal = any(a.GetSymbol() in metals for a in mol.GetAtoms())
        n_carbon = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
        n_ring = rdMolDescriptors.CalcNumRings(mol)
        heavy = mol.GetNumHeavyAtoms()
        if has_metal == True:
            hard.append("inorganic salt or counterion")
        if n_carbon < 8 or n_ring == 0:
            hard.append(f"not drug like ({n_carbon} C, {n_ring} rings)")
        reasons += hard

        if neg == 0:
            reasons.append("no anionic group")
        if net > 0:
            reasons.append(f"net positive ({net:+d})")
        if clogp > 4.0:
            reasons.append(f"lipophilic (clogP {clogp})")
        if n_grip == 0:
            reasons.append("no contact with the phosphate arginines")

        pains_hit = ""
        entry = pains.GetFirstMatch(mol)
        if entry is not None:
            pains_hit = entry.GetDescription()
            reasons.append(f"known frequent hitter ({pains_hit})")

        if len(hard) > 0:
            verdict = "implausible: " + ", ".join(reasons)
        elif len(reasons) == 0:
            verdict = "plausible"
        elif neg >= 1 and n_grip >= 1:
            verdict = "worth a look"
        else:
            verdict = "implausible: " + ", ".join(reasons)

        if len(hard) == 0 and pains_hit != "" and verdict == "worth a look":
            verdict = f"anionic but a known frequent hitter ({pains_hit})"

        out.append({**r, "charge": net, "anionic_groups": neg, "basic_groups": pos,
                    "clogp": clogp, "grip_contacts": n_grip, "pains": pains_hit,
                    "verdict": verdict})

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    plaus = [r for r in out if str(r["verdict"]).startswith("implausible") == False
             and r["verdict"] != "no structure"]
    implaus = [r for r in out if str(r["verdict"]).startswith("implausible") == True]

    lines = ["# Could the shortlist chemically bind this site", "",
             "The substrate site grips a pyrophosphate. Four arginines line it,",
             "ARG144, ARG188, ARG213 and ARG277, and every natural ligand carries",
             "two or three negative charges at pH 7.4. Vina has no electrostatics",
             "term, so it cannot reward charge pairing and cannot punish its",
             "absence. It does have a hydrophobic term, which a large greasy",
             "molecule can earn by filling the site without forming a single",
             "interaction the pocket exists for.", "",
             "That is the artifact this page checks for.", "",
             "## Summary", "",
             f"- {len(out)} shortlisted structures",
             f"- {len(implaus)} are chemically implausible for this site",
             f"- {len(plaus)} carry an anionic group and reach the arginines", ""]

    def block(rs, title, note):
        if len(rs) == 0:
            return [f"### {title}", "", "None.", ""]
        b = [f"### {title}", "", note, "",
             "| Compound | ChEMBL | Score | Pct | Charge | Anionic | clogP | "
             "Arg contacts | Read |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for r in sorted(rs, key=lambda r: float(r["score"])):
            name = r["name"] if r["name"] != "" else "(unnamed)"
            b.append(f"| {name.title()} | {r['id']} | {r['score']} | "
                     f"{r['pct'] or 'n/a'} | {r['charge']:+} | {r['anionic_groups']} | "
                     f"{r['clogp']} | {r['grip_contacts']}/4 | {r['verdict']} |")
        b.append("")
        return b

    lines += block(plaus, "Chemically plausible for this site",
                   "An anionic group, and a pose that reaches the arginines.")
    lines += block(implaus[:30], "Implausible, highest scoring first",
                   "These are the scores the hydrophobic term paid for.")

    open(OUT_MD, "w").write("\n".join(lines))
    print(f"{len(plaus)} plausible, {len(implaus)} implausible")
    print(f"written to {OUT_MD}")


if __name__ == "__main__":
    main()
