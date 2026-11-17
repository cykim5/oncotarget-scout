#!/bin/bash
# Wait for the re-centred run, then start the scoring benchmark.
cd "$(dirname "$0")"
while pgrep -f run_recentred.py > /dev/null; do sleep 20; done
/opt/anaconda3/envs/docking/bin/python -I benchmark_scoring.py
