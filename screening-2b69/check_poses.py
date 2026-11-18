"""Where does each shortlisted pose actually sit?

A good score is not the same as binding in the right place. Three ways
a high scoring pose can be worthless:

  1. it is pressed against a box wall, which means the search was cut
     off and the real optimum is outside the box
  2. it is parked in some other groove, nowhere near the substrate site
  3. it is sitting on the NAD cofactor rather than in the substrate site

The substrate site is defined empirically here, as the 18 receptor
residues within 4.5 A of the crystallographic UDP in 2B69, rather than
by residue number. 2B69 and 4GLL number the same protein differently,
and hard coding numbers from one into the other is how you end up
measuring the wrong residue.

Writes outputs/pose_check.csv and outputs/pose_check.md
"""
import csv
import math
import os
import sys

from vina import Vina

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.join(BASE, "2b69_holo.pdbqt")
UDP = os.path.join(os.path.dirname(BASE), "control-2b69", "receptor", "2b69_UDP900.pdb")
APO = os.path.join(os.path.dirname(BASE), "control-2b69", "receptor", "2b69_chainA_NAD.pdb")
LIGANDS = os.path.join(os.path.dirname(BASE), "screening-full-uxs1", "ligands")
RANKED = os.path.join(BASE, "outputs", "ranked.csv")
POSE_DIR = os.path.join(BASE, "outputs", "poses")
OUT_CSV = os.path.join(BASE, "outputs", "pose_check.csv")
OUT_MD = os.path.join(BASE, "outputs", "pose_check.md")

CENTER = [50.48, 21.14, 23.36]
SIZE = [28, 28, 28]
EXHAUSTIVENESS = 8
SEED = 42
N_POSES = 9
CONTACT = 4.5            # angstroms, the same cutoff that defined the site
WALL_WARN = 2.0          # a pose atom this close to a face is clipped


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


def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms) / n,
            sum(a[1] for a in atoms) / n,
            sum(a[2] for a in atoms) / n)


def site_residues(udp, rec):
    """The residues that touch the crystal ligand."""
    keys = set()
    for a in rec:
        if a[3] == "UDP":
            continue
        d = min(math.dist(a[:3], u[:3]) for u in udp)
        if d <= CONTACT:
            keys.add((a[3], a[4]))
    return keys


def pose_atoms(pdbqt_text):
    """Atoms of model 1 only."""
    out = []
    in_model = False
    for line in pdbqt_text.splitlines():
        if line.startswith("MODEL"):
            if in_model == True:
                break
            in_model = True
        if line.startswith(("ATOM", "HETATM")):
            try:
                out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
        if line.startswith("ENDMDL"):
            break
    return out


def wall_margin(atoms):
    """Smallest distance from any pose atom to any box face."""
    best = 1e9
    for x, y, z in atoms:
        for v, c, s in ((x, CENTER[0], SIZE[0]), (y, CENTER[1], SIZE[1]),
                        (z, CENTER[2], SIZE[2])):
            best = min(best, (c + s / 2.0) - v, v - (c - s / 2.0))
    return best


def main():
    shortlist = []
    rows = list(csv.DictReader(open(RANKED)))
    by_score = sorted(rows, key=lambda r: float(r["score"]))[:25]
    by_le = sorted(rows, key=lambda r: float(r["le"]))[:25]
    with_pct = [r for r in rows if r["pct"] != ""]
    by_pct = sorted(with_pct, key=lambda r: -float(r["pct"]))[:25]
    high_pct = [r for r in with_pct if float(r["pct"]) >= 99.0]
    seen = set()
    for group, tag in ((by_score, "raw"), (by_le, "le"),
                       (by_pct, "pct"), (high_pct, "pct99")):
        for r in group:
            if r["id"] in seen:
                continue
            seen.add(r["id"])
            shortlist.append(r)
    limit = int(os.environ.get("SHORTLIST_LIMIT", "0"))
    if limit > 0:
        shortlist = shortlist[:limit]
    print(f"shortlist of {len(shortlist)} structures", flush=True)

    udp = read_atoms(UDP)
    rec = read_atoms(APO)
    site = site_residues(udp, rec)
    udp_c = centroid(udp)
    nad = [a for a in rec if a[3] == "NAD"]
    print(f"site defined by {len(site)} residues, {len(nad)} NAD atoms", flush=True)

    os.makedirs(POSE_DIR, exist_ok=True)
    out = []
    for i, r in enumerate(shortlist, 1):
        lig = os.path.join(LIGANDS, r["id"] + ".pdbqt")
        if os.path.exists(lig) == False:
            print(f"  {r['id']}: no prepared ligand, skipped", flush=True)
            continue
        pose_path = os.path.join(POSE_DIR, r["id"] + ".pdbqt")
        if os.path.exists(pose_path) == False:
            v = Vina(sf_name="vina", seed=SEED, verbosity=0)
            v.set_receptor(RECEPTOR)
            v.set_ligand_from_file(lig)
            v.compute_vina_maps(center=CENTER, box_size=SIZE)
            v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
            v.write_poses(pose_path, n_poses=N_POSES, overwrite=True)
        atoms = pose_atoms(open(pose_path).read())
        if len(atoms) == 0:
            continue

        pc = centroid(atoms)
        touched = set()
        for a in rec:
            if a[3] == "UDP":
                continue
            key = (a[3], a[4])
            if key not in site:
                continue
            if min(math.dist(a[:3], p) for p in atoms) <= CONTACT:
                touched.add(key)
        nad_min = min((math.dist(a[:3], p) for a in nad for p in atoms),
                      default=99.0)
        out.append({
            "id": r["id"], "name": r["name"], "heavy": r["heavy"],
            "score": r["score"], "le": r["le"], "pct": r["pct"],
            "wall_margin": round(wall_margin(atoms), 2),
            "centroid_to_udp": round(math.dist(pc, udp_c), 2),
            "site_residues_touched": len(touched),
            "site_fraction": round(len(touched) / len(site), 2),
            "nad_closest": round(nad_min, 2),
            "smiles": r.get("smiles", ""),
        })
        print(f"  [{i}/{len(shortlist)}] {r['id']:<16} wall {out[-1]['wall_margin']:>5} "
              f"udp {out[-1]['centroid_to_udp']:>5} site {len(touched)}/{len(site)}",
              flush=True)

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    clipped = [r for r in out if r["wall_margin"] < WALL_WARN]
    elsewhere = [r for r in out if r["site_fraction"] < 0.3]
    good = [r for r in out if r["wall_margin"] >= WALL_WARN
            and r["site_fraction"] >= 0.5]

    lines = ["# Where the shortlisted poses actually sit", "",
             f"{len(out)} structures re-docked with poses kept, under the same",
             "protocol as the screen. The substrate site is the 18 receptor",
             "residues within 4.5 A of the crystallographic UDP in 2B69.", "",
             "## Summary", "",
             f"- {len(clipped)} poses come within {WALL_WARN} A of a box face, "
             f"so the search was cut off and the score is not trustworthy",
             f"- {len(elsewhere)} poses touch less than 30% of the site residues, "
             f"so they are sitting somewhere else on the protein",
             f"- {len(good)} poses are clear of the walls and touch at least half "
             f"the site residues", ""]

    def block(rs, title, note):
        if len(rs) == 0:
            return [f"### {title}", "", "None.", ""]
        b = [f"### {title}", "", note, "",
             "| Compound | ChEMBL | Heavy | Score | Pct | Wall margin | "
             "Centroid to UDP | Site residues | NAD closest |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for r in sorted(rs, key=lambda r: float(r["score"])):
            name = r["name"] if r["name"] != "" else "(unnamed)"
            b.append(f"| {name.title()} | {r['id']} | {r['heavy']} | {r['score']} | "
                     f"{r['pct'] or 'n/a'} | {r['wall_margin']} | "
                     f"{r['centroid_to_udp']} | {r['site_residues_touched']}/18 | "
                     f"{r['nad_closest']} |")
        b.append("")
        return b

    lines += block(good, "Poses that survive both checks",
                   "These are the only ones worth a sentence.")
    lines += block(clipped, "Poses clipped by a box wall",
                   "The search space ran out before the pose did.")
    lines += block(elsewhere, "Poses somewhere other than the substrate site",
                   "A score earned in a different groove.")

    open(OUT_MD, "w").write("\n".join(lines))
    print(f"\n{len(clipped)} clipped, {len(elsewhere)} off site, {len(good)} survive")
    print(f"written to {OUT_MD}")


if __name__ == "__main__":
    main()
