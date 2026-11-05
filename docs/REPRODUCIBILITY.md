# Reproducing the experiments

## Scientific configuration

The primary case is κ_f=0 (smooth convex); κ_f=0.1 is supplementary. The supplied results fix L_f=L_B=R_x=1 and ρ=L_B R_u/R_x=1, and vary L_A from 0 to 0.9 in increments of 0.1. The full configuration is [config/paper.toml](../config/paper.toml).

The normalized step size is α=L_f L_B² η, where η is the physical step size. Convert a normalized objective value back by multiplying by L_f R_x². In particular, physical L_f=10 corresponds to η=α/10 when L_B=1. The archived calculations and the plotted step-size axes use normalized L_f=1.

The main metric is **f(G u_T) − f(x_*)**, with G=(I−A)⁻¹B. It is not the transient state gap f(x_T)−f(x_*). The implementation keeps T step-size coefficients according to its trajectory indexing. The scripts and figures use the stored `T` column without shifting it.

Best-point and horizon plots use W(α,T)/W(0,T), where W(0,T)=ρ²/[2(1−L_A)²] in normalized units. The α-sweep remains in raw regret units. The uniform analytical step size is

```text
α_analytic = 2(1 − L_A)^3 / (L_f L_B^2 (3 − L_A)).
```

Its performance is evaluated with a fresh PEP at this exact step size, rather than a nearest sampled step size. Main figures use T=20 for the α-sweep and T∈{5,10,20,30,40,50} for the norm-bound comparison. The horizon view uses all available design horizons.

## Stored data and figure regeneration

Run `scripts/reproduce_figures.jl` as shown in the root README. It loads the original tables, restricts them to 0≤L_A≤0.9, applies only accepted repair rows, and plots κ=0 and κ=0.1 into separate directories. The original `.arrow` files remain byte-for-byte unchanged. `--kappa 0` selects just the paper case.

Analytical rows marked `UNCERTIFIED` are excluded by default. With `--include-uncertified`, they are displayed with red cross markers. Missing finite designs are excluded from the default design-derived curves, and `unresolved.csv` reports them. Optional diagnostic views can use sampled substitutes; their `best_source` field identifies them explicitly. A finite sample minimum is not silently presented as an optimized design in the default view.

## Fresh computation

The entry point `scripts/run_experiments.jl` sets one shared configuration for all three stages:

1. `sweep`: solve the PEP on the constant-step grid (100 points from 0.001 to 2 for the full run).
2. `design`: use the lowest finite sampled value as the starting point for a positive constant-policy FOM design; the initializer is 0.25 only when no matching sample is available.
3. `analytic`: evaluate the uniform analytical step size on the same class PEP.

Use `--stage sweep`, `--stage design`, or `--stage analytic` to run an individual stage. Outputs live in `sweep/`, `design/`, and `analytic/`, each with shards, `merged.arrow`, and configuration records. Keep the same configuration, solver, and output directory when resuming; use a new output directory when changing any of them. Stored scientific configurations are checked before resuming.

To plot a newly computed full run:

```sh
julia --project=. scripts/plot_preliminary_sweep.jl results/full/sweep/merged.arrow --design results/full/design/merged.arrow --analytic results/full/analytic/merged.arrow --kappa 0 --curve-T 20 --horizons 5,10,20,30,40,50 --norm-a 0.0,0.1,0.3,0.5,0.7 --include-design-only --normalize-worst-cases --sweep-xmax 1 --split-horizon-panels --out results/full/figures
```

All fresh runs use explicit facial reduction and original-space residual checks. Defaults accept a relative gap ≤1e−5 and scalar/cone violations ≤3e−6. `certified` means passing these floating-point numerical checks; it is not an exact-arithmetic proof. Non-finite, failed, or rejected values are retained as diagnostic rows and are excluded from plotted finite results.

FOM produces a locally optimized positive constant step size. Its output depends on initialization and the iteration budget; this code does not certify a global minimum over all step sizes. Hardware, solver versions, and numerical conditioning can change which points pass tolerance checks.

## Historical T=50 refinement

The latest reference design incorporates six 60-iteration T=50 reruns, replacing earlier 30-iteration designs. Their starts came from T=40 designs. The affected (L_A,κ_f) pairs are (0,0), (0.1,0), (0.2,0), (0,0.1), (0.1,0.1), and (0.4,0.1). Their exact starts are stored as `alpha_init` in the supplied design table.

`scripts/replay_design.jl` replays selected stored starts, with 60 iterations for these six historical cases and 30 otherwise. This is the route for comparing a rerun against a particular archived design. Its default is a single small T=3 case; larger replay requests must be selected explicitly. It is separate from a fresh grid experiment, whose initializers come from its new sweep.

## Bounded repair attempts

```sh
julia --project=. scripts/retry_reference_points.jl --max-T 12 --limit 2 --design-iters 3 --seconds 60
```

This identifies missing designs and uncertified analytical points, attempts only the requested bounded subset, and writes `data/corrections/`. A retry is overlaid on a reference row only if it passes the numerical acceptance checks. The original tables are never overwritten. The runtime limit applies to each solver call; compilation and several design calls can add time.

`retry_attempts.csv` records computed and deferred points, objectives, solver, gaps, residuals, and status. `unresolved.csv` lists remaining groups. Increase the horizon, count, time, and design iteration limits explicitly for broader repairs, preferably on a cluster. This package's validation does not rerun the full T=50 experiment.
