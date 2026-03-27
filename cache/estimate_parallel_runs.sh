#!/bin/bash
# estimate_parallel_runs.sh - Launch independent estimation runs
#
# Usage: bash code/scripts/estimate_parallel_runs.sh

set -e
cd "$(dirname "$0")/.."

N_RUNS=4
THREADS_PER_RUN=7

echo "Launching ${N_RUNS} runs × ${THREADS_PER_RUN} threads = $((N_RUNS * THREADS_PER_RUN)) cores"

PIDS=()
for i in $(seq 1 $N_RUNS); do
    mkdir -p "results/run_${i}"
    echo "Starting run ${i}..."
    julia -t ${THREADS_PER_RUN} --project=. scripts/estimate_parameters_run.jl ${i} > "results/run_${i}/log.txt" 2>&1 &
    PIDS+=($!)
done

echo "PIDs: ${PIDS[@]}"
echo ""
echo "Monitor with: tail -f code/results/run_*/log.txt"
echo "After completion: julia --project=code code/scripts/combine_ensembles.jl"
echo ""

# Wait and report
for i in "${!PIDS[@]}"; do
    RUN=$((i + 1))
    wait ${PIDS[$i]} && echo "Run ${RUN} finished OK" || echo "Run ${RUN} FAILED"
done

echo "All runs complete."
