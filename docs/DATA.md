# Reference data and parsing

The three original Arrow tables total approximately 6.4 MB. Their provenance and SHA-256 digests are recorded in [PROVENANCE.toml](../data/reference/PROVENANCE.toml). Run `julia --project=. test/runtests.jl --data-only` to verify them.

| Table | Archived rows | Coverage |
|---|---:|---|
| `merged_sweep.arrow` | 46,200 | T=3:20,25,30,40; κ=0,0.1; 100 α samples per group |
| `merged_design.arrow` | 484 | T=3:20,25,30,40,50; κ=0,0.1 |
| `merged_analytic_theorem1.arrow` | 440 | Same design horizons, stable L_A≤0.9 |

The raw sweep/design archives include L_A=1, outside the paper's stable domain. The loader preserves the originals but excludes these rows from reviewer figures. In the stable domain, the uncorrected design table has 436 finite `OK` rows and four failures. The analytical table has 435 `OK` rows and five `UNCERTIFIED` rows. The sweep contains 26,830 finite accepted rows and 15,170 rejected rows; the rejected values are not substituted into curves. There is no archived T=50 α-sweep; its design and analytical points are available independently.

## Read a table directly

```julia
using Arrow, DataFrames
df = DataFrame(Arrow.Table("data/reference/merged_design.arrow"))
paper = df[(df.kappa .== 0.0) .& (df.norm_A .<= 0.9), :]
show(first(paper, 10); allcols=true)
```

For domain filtering and accepted correction overlays:

```julia
include("src/ReviewData.jl")
using .ReviewData, DataFrames
tables = load_reference()
issues = audit_reference(tables)
issues[issues.needs_repair, :]
```

The command `julia --project=. scripts/inspect_results.jl` exports all three stable-domain tables to CSV and prints concise status counts. `coverage.csv` has one row per expected design/analytical group; `unresolved.csv` contains the remaining exceptions. It does not require a solver.

## Columns

All tables share `T`, `n_steps=T+1`, `kappa=mu_f/L_f`, `mu_f`, `L_f`, `norm_A`, `norm_B`, `R_x`, `R_u`, `status`, and `solve_time` in seconds.

| Table | Additional useful columns |
|---|---|
| Sweep | `alpha`, `pep` |
| Design | `alpha_init`, `optimizer`, `alpha_best`, `pep_best`, `n_solves` |
| Analytical | `alpha_analytic`, `pep_analytic`, `analytic_relative_gap`, `analytic_scalar_violation`, `analytic_cone_violation`, `analytic_solver_status` |

The analytical schema retains historical recovery-related fields, such as `gram_rank` and `recovered_norm_A`. They are unused/NaN in the uniform Theorem-1 experiment and should not be interpreted as recovered plant data. Empty `realization_file` entries are expected.

Original sweep status strings may contain the individual rejection diagnostics. Original design rows contain status and objective but do not carry every solver residual; reading the archive does not rerun numerical certification. Repair attempts have separate residual records. `OK` denotes a numerically accepted result, `UNCERTIFIED` a rejected estimate, and `ERROR` a failed calculation. `SLOW_PROGRESS` must be interpreted together with the producing solver/certificate path, rather than treated as a proof of optimality by itself.

## Generated plot tables

Each plot directory includes `best_points.csv` (including horizon values), `sampled_best.csv`, `group_coverage.csv`, `alpha_sweep.csv`, and `plot_summary.txt`. Tables are restricted to that directory's selected κ. The α-sweep CSV retains the full sampled range even when the displayed axis is clipped. Normalized columns retain matching `*_unnormalized` values and the denominator `W_zero`. `best_source`, `analytic_source`, and `analytic_certified` describe provenance and acceptance.

Original reference tables and separately generated correction/diagnostic artifacts are versioned. Large new experiment outputs go under the ignored `results/` directory.
