#!/usr/bin/env bash
# Generic Slurm pipeline: sweep -> merge -> design -> merge -> analytic -> merge.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export PEP_REPO_ROOT="$root"
config=${1:-"$root/config/paper.toml"}
out=${2:-"$root/results/cluster"}
config=$(cd -- "$(dirname -- "$config")" && pwd)/$(basename -- "$config")
mkdir -p -- "$out/logs"
out=$(cd -- "$out" && pwd)
groups=$(julia --project="$root" "$root/scripts/run_experiments.jl" --config "$config" --stage count)
[[ "$groups" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid group count\n' >&2; exit 1; }
export PEP_ARRAY_TASK_COUNT=$groups
export PEP_SOLVER="${PEP_SOLVER:-auto}"
parallel=${PEP_CONCURRENCY:-2}
cpus=${PEP_CPUS:-1}
memory=${PEP_MEMORY:-64G}
walltime=${PEP_TIME:-24:00:00}
dependency=()
for stage in sweep design analytic; do
  job=$(sbatch --parsable --export=ALL --cpus-per-task="$cpus" --mem="$memory" \
    --time="$walltime" --array="0-$((groups-1))%$parallel" \
    --output="$out/logs/%x-%A_%a.out" "${dependency[@]}" \
    "$root/cluster/worker.slurm" "$config" "$out" "$stage")
  job=${job%%;*}
  # The merger is a single job; remove the array identity before running it.
  merge=$(sbatch --parsable --export=ALL --cpus-per-task=1 --mem=8G --time=00:30:00 \
    --dependency="afterok:$job" --output="$out/logs/merge-%j.out" \
    "$root/cluster/worker.slurm" "$config" "$out" "merge-$stage")
  merge=${merge%%;*}
  dependency=(--dependency="afterok:$merge")
  printf '%s array=%s merge=%s\n' "$stage" "$job" "$merge"
done
printf 'After the last merge, see docs/CLUSTER.md for plotting commands.\n'
