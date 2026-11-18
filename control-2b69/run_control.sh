#!/bin/bash
cd "$HOME/uxs1-2b69"
export PYTHONPATH="$HOME/dock-check"
L="$HOME/oncotarget-scout-work/docking-uxs1/ligands"
/opt/anaconda3/envs/docking/bin/python -m dockcheck run \
  --receptor receptor/2b69_holo.pdbqt \
  --native "$L/UDP-glucuronate.pdbqt" \
  --native "$L/UDP-xylose.pdbqt" \
  --native "$L/UDP-galactose.pdbqt" \
  --native "$L/UDP.pdbqt" \
  --native "$L/UMP.pdbqt" \
  --decoys decoys/ \
  --center 50.48 21.14 23.36 \
  --size 28 \
  --exhaustiveness 8 \
  --seed 42 \
  --out report/
