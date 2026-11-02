# Check the numerical products of the configured small end-to-end run.
import Pkg
Pkg.activate(joinpath(@__DIR__,"..");io=devnull)
using Arrow, DataFrames, TOML
root=normpath(joinpath(@__DIR__,".."))
length(ARGS)==1 || error("usage: julia --project=. scripts/check_smoke.jl RESULTS_DIR")
cfg=TOML.parsefile(joinpath(root,"config","smoke.toml"))
groups=length(cfg["horizons"])*cfg["n_A"]*length(cfg["kappas"])
for (stage,value,expected) in (("sweep",:pep,groups*cfg["n_alpha"]),
                              ("design",:pep_best,groups),
                              ("analytic",:pep_analytic,groups))
    df=DataFrame(Arrow.Table(joinpath(ARGS[1],stage,"merged.arrow")))
    nrow(df)==expected || error("$stage has $(nrow(df)) rows, expected $expected")
    all(isfinite,df[!,value]) || error("$stage contains non-finite values")
    if stage=="sweep"
        all(in(("OPTIMAL","ALMOST_OPTIMAL","SLOW_PROGRESS","OK")),df.status) ||
            error("sweep contains rejected rows")
    else
        all(==("OK"),df.status) || error("$stage contains rejected rows")
    end
    if stage=="analytic"
        all(df.analytic_relative_gap .<= cfg["gap_tol"]) || error("analytical gap failure")
        all(df.analytic_scalar_violation .<= cfg["cert_tol"]) || error("analytical scalar residual failure")
        all(df.analytic_cone_violation .<= cfg["cert_tol"]) || error("analytical cone residual failure")
    end
    println("$stage: $expected finite accepted rows")
end
println("Small end-to-end experiment passed")
