# Generic Slurm execution

The launchers contain no cluster-specific account, partition, email address, hostname, scratch path, or license location. Configure those through your scheduler/environment. The repository must be on a filesystem shared by all assigned nodes.

## Setup

Load Julia 1.12.6 using your cluster's module system, clone/copy this repository, and instantiate once on a login/setup node:

```sh
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
export PEP_SOLVER=mosek
```

For MOSEK, arrange a license visible to compute nodes using its normal installation or `MOSEKLM_LICENSE_FILE` in your shell. License contents are never copied into this repository. Use `PEP_SOLVER=clarabel` for license-free execution or `auto` for availability-based selection. Set a shared `JULIA_DEPOT_PATH` yourself if your cluster requires one.

## Submit

From the repository root:

```sh
# A small scheduler/pipeline check.
bash cluster/submit.sh config/smoke.toml results/cluster_smoke

# Full paper plus supplementary grid.
PEP_CONCURRENCY=2 PEP_CPUS=1 PEP_MEMORY=64G PEP_TIME=24:00:00 bash cluster/submit.sh config/paper.toml results/cluster_full
```

Use the scheduler's `SBATCH_ACCOUNT` and `SBATCH_PARTITION` environment variables where required. The resource values above are configurable allocation defaults, not measured guarantees for every horizon. High horizons can require substantial memory and time; begin with the smoke configuration. Increasing Julia threads can multiply concurrent SDP memory. Each solver itself uses one thread by default.

The submitter reads the configuration to count (T,L_A,κ_f) groups and creates:

```text
sweep array -> merge-sweep -> design array -> merge-design -> analytic array -> merge-analytic
```

There is one array task per group, and dependencies require successful completion before the next stage. Logs go to `OUTPUT/logs/`. `PEP_REPO_ROOT` is exported by the submitter because Slurm executes a spooled copy of the worker script; the worker does not infer the repository from that copy's location.

## Resume and inspect

Each task writes crash-safe shards under its stage's `tasks/` directory. Resubmitting the same configuration into the same output directory reuses stored points. Keep solver selection and configuration unchanged when resuming. A nonzero analytical job exit prevents downstream merging and leaves successful shards available for diagnosis/retry. Failed/uncertified sweep and design rows are retained, and the completion summary identifies them; successful scheduling does not imply numerical acceptance of every point.

For targeted reruns, pass the original task index to `sbatch --array=INDEX ... cluster/worker.slurm CONFIG OUTPUT STAGE` with `PEP_ARRAY_TASK_COUNT` and `PEP_REPO_ROOT` set consistently. After all tasks succeed, a stage can be merged manually:

```sh
julia --project=. scripts/run_experiments.jl --config config/paper.toml --out results/cluster_full --stage merge-analytic
```

To plot a completed full run, use the direct plotter command in [REPRODUCIBILITY.md](REPRODUCIBILITY.md), replacing `results/full` with your cluster output directory. Copy the merged Arrow files to a workstation if preferred; plotting does not need a MOSEK license. Historical T=50 refinements and bounded reference repair are documented in the same guide.
