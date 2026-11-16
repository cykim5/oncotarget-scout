#!/bin/bash
# Is the UXS1 screen still going, and how far along is it?
D="$HOME/oncotarget-scout-work/screening-full-uxs1"
TOTAL=$(( $(wc -l < "$D/library.csv") - 1 ))
DONE=$(ls "$D/outputs" 2>/dev/null | grep -c "_energies.txt")

if pgrep -f "run_screen.py" > /dev/null; then
  STATE="RUNNING"
else
  if [ -f "$D/outputs/results_full.csv" ]; then
    STATE="FINISHED"
  else
    STATE="STOPPED before finishing, rerun ./run_full_screen.sh to resume"
  fi
fi

echo "status : $STATE"
echo "docked : $DONE of $TOTAL  ($(( DONE * 100 / TOTAL ))%)"

FIRST=$(ls -tr "$D/outputs"/*_energies.txt 2>/dev/null | head -1)
if [ -n "$FIRST" ] && [ "$DONE" -gt 0 ]; then
  START=$(stat -f %m "$FIRST")
  NOW=$(date +%s)
  ELAPSED=$(( NOW - START ))
  echo "elapsed: $(( ELAPSED / 60 )) min"
  if [ "$STATE" = "RUNNING" ] && [ "$DONE" -gt 5 ]; then
    LEFT=$(( (TOTAL - DONE) * ELAPSED / DONE ))
    echo "left   : about $(( LEFT / 3600 ))h $(( (LEFT % 3600) / 60 ))m"
    echo "eta    : $(date -v +${LEFT}S '+%a %H:%M')"
  fi
fi

if [ "$STATE" = "FINISHED" ]; then
  echo
  echo "top 5 by raw affinity (remember this ranking is half molecular weight):"
  head -6 "$D/outputs/results_full.csv" | column -t -s,
fi
