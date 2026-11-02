import Pkg
Pkg.activate(joinpath(@__DIR__,".."); io=devnull)
include(joinpath(@__DIR__,"..","src","ReviewData.jl"))
using .ReviewData, DataFrames, CSV, Arrow

function main()
    opts = Dict("data"=>joinpath(ROOT,"data","reference"),
                "out"=>joinpath(ROOT,"results","inspection"))
    for i in 1:2:length(ARGS)
        i < length(ARGS) || error("$(ARGS[i]) needs a value")
        key = replace(ARGS[i],r"^--"=>"")
        haskey(opts,key) || error("unknown option $(ARGS[i])")
        opts[key] = ARGS[i+1]
    end
    tables = load_reference(directory=opts["data"])
    mkpath(opts["out"])
    for kind in (:sweep,:design,:analytic)
        table = getproperty(tables,kind)
        println("$kind: $(nrow(table)) stable-domain rows; horizons=$(sort(unique(table.T))); kappa=$(sort(unique(table.kappa)))")
        classes = first.(split.(String.(table.status),":"))
        counts = combine(groupby(DataFrame(status=classes),:status),nrow=>:rows)
        show(stdout,counts; allrows=true,summary=false); println()
        CSV.write(joinpath(opts["out"],"$(kind).csv"),table)
    end
    audit = audit_reference(tables)
    CSV.write(joinpath(opts["out"],"coverage.csv"),audit)
    issues = audit[audit.needs_repair,:]
    CSV.write(joinpath(opts["out"],"unresolved.csv"),issues)
    println("Unresolved groups: $(nrow(issues)) / $(nrow(audit))")
    show(stdout,issues; allrows=true,summary=false); println()
    println("CSV tables and coverage report: ",relpath(opts["out"],ROOT))
end
main()
