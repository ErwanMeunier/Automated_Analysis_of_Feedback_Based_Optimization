# Bounded repairs of reference points

The original tables in `data/reference/` remain unchanged. The loader overlays
only repair rows with `status=OK` and matching physical parameters.

On 2026-11-02, a bounded MOSEK run attempted two points, using T≤12, three
design iterations, and a 45-second limit per solver call:

- **T=11, L_A=0.9, κ=0:** recovered a finite positive FOM design with a passing
  numerical certificate. Its stored step is about 0.001196805 and its objective
  is about 20.68564. This short repair is separate from the original 30-iteration
  experiment and is not a global step-size optimality certificate.
- **T=12, L_A=0.9, κ=0.1:** the analytical retry remained uncertified (cone
  violation about 3.51e−6, above the 3e−6 threshold). It does not replace the
  original row.

Seven groups remain unresolved: three lack an accepted design and five lack an
accepted analytical result, with one group counted in both categories. Larger
retries were deferred to respect the bounded local validation budget.
`retry_attempts.csv` records attempts and deferrals; `unresolved.csv` gives the
remaining group list. Accepted correction data are stored in
`merged_design.arrow`; `merged_analytic_theorem1.arrow` retains the rejected
retry for audit, and is ignored by the overlay.

Re-running the repair script appends diagnostics and attempts the still
unresolved subset. Use explicit CLI limits for any larger computation.
