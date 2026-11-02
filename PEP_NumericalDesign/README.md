# Vendored PEPDesign source

This directory contains the PEPDesign source required by `src/PEPControl.jl`. It is loaded using a relative `include`, so no external repository, network checkout, or path dependency is required. Use the root `Project.toml` and `Manifest.toml` for every reproduction command.

The mathematical modeling, compilation, facial reduction, and design implementation were copied from the supporting `PEP_Control` repository. The reviewer package adds optional MOSEK loading, automatic MOSEK/Clarabel backend selection, a generic primal/dual objective gap fallback, and optional solver time limits. Default backend selection is shared by direct PEP evaluation and FOM design.

The full source dependency closure is retained; unrelated toolbox experiment suites, development scripts, slides, and archived outputs are omitted. The local [MIT license](LICENSE) and root [third-party notices](../THIRD_PARTY_NOTICES.md) apply. Consult the root [reproduction guide](../docs/REPRODUCIBILITY.md) for the numerical acceptance conventions used by this paper package.
