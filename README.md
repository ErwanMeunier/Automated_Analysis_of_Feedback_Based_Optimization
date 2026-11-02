# Automated Analysis of Feedback-Based Optimization

Julia code and reference data for reproducing the paper's numerical experiments with performance estimation problems (PEPs). The main experiment uses **κ_f = μ_f/L_f = 0**. Results for **κ_f = 0.1** are included as supplementary material. All supplied data use normalized **L_f = L_B = R_x = ρ = 1**; see [the normalization and numerical conventions](docs/REPRODUCIBILITY.md).

## Start here

Use **Julia 1.12.6**, the version recorded in `Manifest.toml`, and run these commands from this repository's root. The environment is pinned; no external checkout of the PEPDesign toolbox is needed.

```sh
julia --project=. -e "import Pkg; Pkg.instantiate()"
julia --project=. scripts/inspect_results.jl
julia --project=. scripts/reproduce_figures.jl
```

The last two commands **do not load an SDP solver or require a MOSEK license**. They read the supplied Arrow tables and generate CSV reports, PDF figures, and PNG previews under `results/`:

- `results/inspection/`: readable data tables, coverage, and unresolved numerical points.
- `results/figures/paper_kappa_0/`: the three paper experiments and auxiliary views.
- `results/figures/supplement_kappa_0p1/`: the same experiments for κ_f = 0.1.

Ready-made plots are also supplied in [data/figures](data/figures). The three main outputs are `paper_figure3.pdf` (step-size sweep at T=20), `best_PEP_and_alpha_vs_T_split.pdf` (horizon dependence), and `paper_figure2.pdf` (step size and performance versus L_A).

Only numerically accepted analytical results are plotted by default. To inspect the historical uncertified estimates as well, use the explicit option below; they receive red cross markers. Every plot directory contains the plotted-value tables and its settings.

```sh
julia --project=. scripts/reproduce_figures.jl --kappa 0 --out results/paper
julia --project=. scripts/reproduce_figures.jl --include-uncertified --out results/diagnostic_figures
```

## Recompute a small experiment

```sh
julia --project=. test/runtests.jl
julia --project=. --threads=1 scripts/run_experiments.jl --config config/smoke.toml --out results/smoke
julia --project=. scripts/check_smoke.jl results/smoke
```

The smoke configuration runs T=3 for both κ values, with four sampled step sizes and two FOM design iterations. It checks the entire sweep → design → analytical-comparison path. The test suite separately checks a T=2 SDP and the L_f scaling identity. To check only the supplied data, use `julia --project=. test/runtests.jl --data-only`. Results and remaining exceptions are recorded in [the validation report](docs/VALIDATION.md).

MOSEK is preferred and is the solver used for the archived experiments. In the default `auto` mode, a small availability/license probe selects MOSEK; if it is unavailable, the code uses the open-source Clarabel solver. The selected backend is logged and saved in each experiment's `config.txt`. To select a solver explicitly:

```sh
julia --project=. scripts/run_experiments.jl --solver clarabel --out results/smoke_clarabel
julia --project=. scripts/run_experiments.jl --solver mosek --out results/smoke_mosek
julia --project=. test/runtests.jl --solver clarabel
```

Explicit `mosek` mode reports an error if MOSEK is unavailable. An individual PEP's numerical failure is recorded rather than silently switching solvers. Clarabel uses the same residual and gap acceptance checks; equality with the archived MOSEK values is subject to numerical tolerance. See the official [MOSEK license setup](https://docs.mosek.com/latest/licensing/quickstart.html) and [Clarabel JuMP interface](https://clarabel.org/stable/julia/jump/).

## Full experiment and cluster execution

```sh
julia --project=. --threads=1 scripts/run_experiments.jl --config config/paper.toml --out results/full
```

This is a substantial computation, including horizons through T=50; it is not necessary to regenerate figures from the supplied data. The [Slurm instructions](docs/CLUSTER.md) provide a resumable array pipeline with configurable resources and no embedded cluster accounts, hosts, or license locations. The [reproduction guide](docs/REPRODUCIBILITY.md) describes configurations, historical T=50 refinements, and bounded repairs.

## Repository contents

| Location | Contents |
|---|---|
| `src/` | Control PEP, plotting/data helpers, and analytical formulas |
| `PEP_NumericalDesign/src/` | Vendored PEPDesign source needed by the control model |
| `scripts/` | Reproduction, inspection, experiment, repair, and merge entry points |
| `config/` | Small and full scientific grids |
| `data/reference/` | Original latest Arrow tables and SHA-256 provenance |
| `data/corrections/` | Separately recorded, bounded retry results and diagnostics |
| `data/figures/` | Generated reviewer figures and plotted-value CSVs |
| `docs/` | Scientific conventions, data schema, cluster use, and validation |
| `test/` | Data integrity and small numerical checks |

The paper package contains Julia and optional Bash/Slurm launchers. The original Python implementation, duplicate checkouts, personal remote-access scripts, manuscript drafts, and old experiments are excluded.

The code is distributed under the [MIT license](LICENSE). Vendored code retains its notices; see [third-party notices](THIRD_PARTY_NOTICES.md). Solver packages retain their respective licenses.
