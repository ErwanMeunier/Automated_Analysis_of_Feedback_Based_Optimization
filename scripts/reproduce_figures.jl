# Solver-free: reads reference Arrow tables and writes figures and their CSVs.
import Pkg
Pkg.activate(joinpath(@__DIR__,".."); io=devnull)
include(joinpath(@__DIR__,"..","src","ReviewData.jl"))
using .ReviewData, Arrow, CSV, DataFrames

function main()
    opts = Dict("data"=>joinpath(ROOT,"data","reference"),
        "out"=>joinpath(ROOT,"results","figures"),"kappa"=>"both",
        "include-uncertified"=>false)
    i=1
    while i<=length(ARGS)
        key=replace(ARGS[i],r"^--"=>"")
        haskey(opts,key) || error("unknown option $(ARGS[i])")
        if key=="include-uncertified"
            opts[key]=true; i+=1
        else
            i<length(ARGS) || error("$(ARGS[i]) needs a value")
            opts[key]=ARGS[i+1]; i+=2
        end
    end
    opts["kappa"] in ("both","0","0.0","0.1") || error("kappa must be 0, 0.1 or both")
    tables=load_reference(directory=opts["data"])
    out=abspath(opts["out"])
    prepared=joinpath(out,"tables")
    mkpath(prepared)
    # Curated copies stay separate from immutable, original reference tables.
    for (kind,file) in ((:sweep,"merged_sweep.arrow"),(:design,"merged_design.arrow"),
                        (:analytic,"merged_analytic_theorem1.arrow"))
        Arrow.write(joinpath(prepared,file),getproperty(tables,kind))
    end
    audit=audit_reference(tables)
    CSV.write(joinpath(out,"coverage.csv"),audit)
    CSV.write(joinpath(out,"unresolved.csv"),audit[audit.needs_repair,:])
    kappas=opts["kappa"]=="both" ? ("0","0.1") : (opts["kappa"],)
    for kappa in kappas
        folder=parse(Float64,kappa)==0 ? "paper_kappa_0" : "supplement_kappa_0p1"
        args=[joinpath(prepared,"merged_sweep.arrow"),
            "--design",joinpath(prepared,"merged_design.arrow"),
            "--analytic",joinpath(prepared,"merged_analytic_theorem1.arrow"),
            "--analytic-mode","theorem1","--include-design-only",
            "--allow-partial-design","--horizons","5,10,20,30,40,50",
            "--curve-T","20","--T-plot-horizons","all","--split-horizon-panels",
            "--norm-a","0.0,0.1,0.3,0.5,0.7","--colormap","viridis",
            "--normalize-worst-cases","--sweep-ymax","2","--sweep-xmax","1","--no-titles",
            "--axis-font-size","18","--legend-font-size","13",
            "--tick-font-size","13","--star-size","7","--line-width","4",
            "--kappa",kappa,"--out",joinpath(out,folder)]
        if opts["include-uncertified"]
            append!(args,["--include-uncertified-analytic","--mark-uncertified-analytic"])
        end
        cmd=`$(Base.julia_cmd()) --project=$ROOT $(joinpath(@__DIR__,"plot_preliminary_sweep.jl")) $args`
        process=run(ignorestatus(addenv(cmd,"GKSwstype"=>"100")))
        success(process) || error("plotting failed for kappa=$kappa; see diagnostics above")
    end
    println("Figures, plotted-value CSVs, and diagnostics: ",relpath(out,ROOT))
end
main()
