import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPControl.jl"))
using .PEPControl
using DataFrames

function main()
    isempty(ARGS) && error("usage: merge_array_results.jl RUN_DIR [EXPECTED_TASKS] [COPY_DEST]")
    run_dir = abspath(ARGS[1])
    expected = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
    copy_dest = length(ARGS) >= 3 ? abspath(ARGS[3]) : nothing
    taskroot = joinpath(run_dir, "tasks")
    isdir(taskroot) || error("no array task directory found at $taskroot")

    taskdirs = sort(filter(isdir, readdir(taskroot; join = true)))
    complete = filter(d -> isfile(joinpath(d, "_SUCCESS")), taskdirs)
    expected === nothing || length(complete) == expected ||
        error("expected $expected successful tasks, found $(length(complete))")
    length(complete) == length(taskdirs) ||
        error("$(length(taskdirs) - length(complete)) task directories lack _SUCCESS")

    out = merge_array_parts(run_dir)
    out === nothing && error("no Arrow shards found below $taskroot")
    df = load_array_results(run_dir)
    println("merged $(nrow(df)) rows from $(length(complete)) tasks into $out")

    if copy_dest !== nothing
        mkpath(dirname(copy_dest))
        temporary = copy_dest * ".tmp"
        cp(out, temporary; force = true)
        mv(temporary, copy_dest; force = true)
        println("copied merged result to $copy_dest")
    end
end

main()
