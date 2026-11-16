"""Render the best holo-run pose (UDP-xylose) in the UXS1 substrate site.

Shows: protein cartoon (gray), docked UDP-xylose mode 1 (cyan sticks),
crystal NAD retained in the receptor (yellow lines), crystal UDP 900 from
2B69 transformed into the 4GLL frame (magenta lines), catalytic triad
Thr202/Tyr231/Lys235 (red sticks).

Usage: LIGAND=UDP-xylose PYTHONPATH=/tmp/imp_shim pymol -c render_pose_holo.py
Run from the docking-uxs1 directory. Output: outputs/holo/<ligand>_best_pose_holo.png
"""
import os

BASE = os.getcwd()
LIGAND = os.environ.get("LIGAND", "UDP-xylose")

from pymol import cmd

cmd.set("opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.load(os.path.join(BASE, "receptor", "uxs1_chainA.pdb"), "uxs1")
cmd.load(os.path.join(BASE, "outputs", "holo", f"{LIGAND}_out.pdbqt"), "docked_all")
cmd.create("docked", "docked_all", 1)
cmd.delete("docked_all")
cmd.load(os.path.join(BASE, "receptor", "crystal_NAD_chainA.pdb"), "crystal_NAD")
cmd.load(os.path.join(BASE, "receptor", "crystal_UDP900_in4GLL.pdb"), "crystal_UDP900")

cmd.hide("everything", "uxs1")
cmd.show("cartoon", "uxs1")
cmd.color("gray80", "uxs1")
cmd.show("sticks", "resi 202+231+235 and uxs1")
cmd.color("red", "resi 202+231+235 and uxs1")
cmd.show("sticks", "docked")
cmd.color("cyan", "docked")
cmd.show("lines", "crystal_NAD")
cmd.color("yellow", "crystal_NAD")
cmd.show("lines", "crystal_UDP900")
cmd.color("magenta", "crystal_UDP900")
cmd.set("stick_radius", 0.25, "docked")
cmd.center("docked")
cmd.zoom("docked", 14)
cmd.bg_color("white")
cmd.ray(1200, 900)
cmd.png(os.path.join(BASE, "outputs", "holo", f"{LIGAND}_best_pose_holo.png"))
print(f"wrote {LIGAND}_best_pose_holo.png")
