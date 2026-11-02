# Data-quality reports

- `coverage_original.csv`: coverage and numerical status before bounded repairs.
- `coverage_with_repairs.csv`: the same groups after accepted repair overlays.
- `unresolved.csv`: remaining incomplete/uncertified groups.
- `sweep_group_quality.csv`: sampled and finite-row counts for every stable-domain sweep group.
- `smoke_solver_comparison.csv`: MOSEK and Clarabel objectives and step sizes from the small end-to-end run; design values use each solver's own optimized step size.

Detailed repair attempts and rejected retry values are retained in
[`../corrections/`](../corrections/). The raw archives are unchanged. These
reports are generated from the supplied data, rather than inferred from the
appearance of curves. `scripts/inspect_results.jl` regenerates a current
coverage report under `results/inspection/`.
