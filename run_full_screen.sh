#!/bin/bash
# Prepare every ligand, then dock the whole library. Both steps are
# resumable, so this can be stopped and started again safely.
set -u
PY=/opt/anaconda3/envs/docking/bin/python
cd "$(dirname "$0")/screening-full-uxs1"
mkdir -p logs
echo "=== prepare started $(date) ==="
"$PY" -I prepare_ligands.py
echo "=== screen started $(date) ==="
"$PY" -I run_screen.py
echo "=== finished $(date) ==="
