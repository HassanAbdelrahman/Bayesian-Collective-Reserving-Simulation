#!/usr/bin/env bash
set -euo pipefail

N_REPS="${1:-3}"
SCENARIO="${2:-rich}"
MAX_PARALLEL="${3:-1}"
PROFILE="${4:-pilot}"

export SIM_ROOT="$(cd "$(dirname "$0")" && pwd)"
export SIM_PROFILE="$PROFILE"
export SIM_SCENARIO="$SCENARIO"

mkdir -p "$SIM_ROOT/results/$SCENARIO/logs"

run_one() {
  local r="$1"
  echo "Starting replication $r / $N_REPS ($SCENARIO, $PROFILE)"
  Rscript "$SIM_ROOT/R/run_replication.R" "$r" "$SCENARIO" \
    > "$SIM_ROOT/results/$SCENARIO/logs/rep_${r}.log" 2>&1
}
export -f run_one
export N_REPS SCENARIO PROFILE SIM_ROOT SIM_PROFILE

active=0
pids=()
for r in $(seq 1 "$N_REPS"); do
  run_one "$r" &
  pids+=("$!")
  active=$((active + 1))
  if (( active >= MAX_PARALLEL )); then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
    active=$((active - 1))
  fi
done
for pid in "${pids[@]}"; do wait "$pid"; done

Rscript "$SIM_ROOT/R/summarize_results.R"
