# Validation of the reviewer package

Validated on 2026-11-02 with Julia 1.12.6 and the pinned root environment.

| Check | Result |
|---|---|
| Original Arrow SHA-256 digests, row counts, κ convention, and correction policy | 15 checks passed |
| T=2 SDP acceptance and L_f scaling, κ=0 and κ=0.1 | 10 numerical checks passed with MOSEK and with Clarabel |
| Automatic fallback with an unavailable MOSEK license | Selected Clarabel and passed the small SDP certificate checks |
| T=3 end-to-end smoke, MOSEK | 8 accepted sweep rows, 2 accepted designs, 2 accepted analytical evaluations |
| T=3 end-to-end smoke, Clarabel | Same accepted row counts; residual and gap checks passed |
| Archived design replay | A T=3, L_A=0.3, κ=0 run with one design iteration returned an accepted result |
| Scientific configuration reuse | A different configuration in an existing output directory was rejected |
| Paper and supplementary plots | Regenerated from supplied data; main PNG layouts inspected |
| Slurm launchers | Bash syntax checked; real Slurm submission was not executed on this Windows workstation |

Clarabel's equilibration and chordal decomposition are disabled for this model: default equilibration produced numerical failures in the reduced degenerate SDP. Acceptance tolerances were retained. The generic solver path computes a primal/dual objective gap if the solver does not expose the MOI relative-gap attribute.

The small solver comparison is saved in [smoke_solver_comparison.csv](../data/diagnostics/smoke_solver_comparison.csv). Its largest objective difference, divided by `max(1,abs(MOSEK value))`, is about 2.06e−5 for the sweep, 1.55e−7 for the designs, and 2.17e−7 for the analytical evaluations. Design comparisons are separate local optimization runs rather than fixed-step evaluations.

All 435 archived analytical rows marked `OK` satisfy the package's configured gap and residual thresholds. Bounded retries recovered one missing T=11 design; the T=12 analytical retry remained uncertified. Seven groups remain unresolved, listed in [the diagnostics](../data/diagnostics/unresolved.csv). Five analytical estimates remain excluded from the default figures. The full high-horizon experiment was not rerun as part of this local validation.

To repeat the small checks:

```sh
julia --project=. test/runtests.jl --solver mosek
julia --project=. test/runtests.jl --solver clarabel
julia --project=. scripts/run_experiments.jl --solver mosek --out results/smoke_mosek
julia --project=. scripts/check_smoke.jl results/smoke_mosek
julia --project=. scripts/run_experiments.jl --solver clarabel --out results/smoke_clarabel
julia --project=. scripts/check_smoke.jl results/smoke_clarabel
```

The data-only test works without loading either solver. The root README gives the solver-free inspection and figure-regeneration commands.
