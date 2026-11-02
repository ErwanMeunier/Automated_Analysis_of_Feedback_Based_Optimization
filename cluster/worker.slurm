#!/usr/bin/env bash
#SBATCH --job-name=pep-review
# No account, partition, mail address, hostname, or license path is embedded.
set -euo pipefail
root=${PEP_REPO_ROOT:?Use cluster/submit.sh or export PEP_REPO_ROOT}
config=$1
out=$2
stage=$3
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export GKSwstype=100
export PEP_SOLVER="${PEP_SOLVER:-auto}"
# Load the Julia module appropriate to your cluster before submission.
julia --project="$root" --threads="$JULIA_NUM_THREADS" "$root/scripts/run_experiments.jl" \
  --config "$config" --out "$out" --stage "$stage" --solver "$PEP_SOLVER"
