"""Render the best-scoring docked pose in the UXS1 pocket with PyMOL.

Usage: LIGAND=NAD PYTHONPATH=/tmp/imp_shim pymol -c render_pose.py
The ligand's mode-1 pose is taken from outputs/<ligand>_out.pdbqt.
Outputs: outputs/<ligand>_best_pose.png
"""
import os

# NOTE: PyMOL's execfile does not set __file__ to the script path, so use
# the working directory. Run this script from the docking-uxs1 directory.
BASE = os.getcwd()
LIGAND = os.environ.get("LIGAND", "NAD")

# Load the multi-model PDBQT directly (PyMOL reads PDBQT natively and
# keeps the atom typing; extracting ATOM lines to PDB breaks bond
# perception). State 1 is the best-scoring pose.
src = os.path.join(BASE, "outputs", f"{LIGAND}_out.pdbqt")
holo_src = os.path.join(BASE, "outputs", "holo", f"{LIGAND}_out.pdbqt")
if not os.path.exists(src) and os.path.exists(holo_src):
    src = holo_src

from pymol import cmd

cmd.set("opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.load(os.path.join(BASE, "receptor", "uxs1_chainA.pdb"), "uxs1")
cmd.load(src, "docked_all")
cmd.create("docked", "docked_all", 1)
cmd.delete("docked_all")
cmd.load(os.path.join(BASE, "receptor", "crystal_NAD_chainA.pdb"), "crystal_NAD")

cmd.hide("everything", "uxs1")
cmd.show("cartoon", "uxs1")
cmd.color("gray80", "uxs1")
cmd.show("sticks", "resi 202+231+235 and uxs1")
cmd.color("red", "resi 202+231+235 and uxs1")
cmd.show("sticks", "docked")
cmd.color("cyan", "docked")
cmd.show("lines", "crystal_NAD")
cmd.color("yellow", "crystal_NAD")
cmd.set("stick_radius", 0.25, "docked")
cmd.center("docked")
cmd.zoom("docked", 12)
cmd.bg_color("white")
cmd.ray(1200, 900)
cmd.png(os.path.join(BASE, "outputs", f"{LIGAND}_best_pose.png"))
print(f"wrote {LIGAND}_best_pose.png")
