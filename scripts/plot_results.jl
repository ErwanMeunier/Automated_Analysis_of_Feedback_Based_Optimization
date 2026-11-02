# Plot any stored quantity against any other, straight from the Arrow runs
# written by `scripts/sweep_alpha.jl` and `scripts/design_best_alpha.jl`.
# No solver is loaded, so this starts in a couple of seconds.
#
# Discover what is in a run:
#   julia --project=Julia Julia/scripts/plot_results.jl Julia/results/sweep_XXXXXXXX --list
#
# Worst case vs step size, one curve per normalized strong-convexity ratio:
#   ... plot_results.jl RUN --x alpha --y pep --series kappa --filter "T=2,norm_A=0.75"
#
# Best-over-α worst case vs horizon:
#   ... plot_results.jl RUN --x T --y pep --series kappa --agg min
#
# Worst case over the (α, ‖A‖) plane at fixed horizon:
#   ... plot_results.jl RUN --mode heatmap --x alpha --y norm_A --z pep --filter "T=2"
#
# Designed vs initial step size (design runs):
#   ... plot_results.jl RUN --mode scatter --x alpha_init --y alpha_best --series kappa
#
# Options: --mode xy|heatmap|scatter  --x/--y/--z COL  --series COL
#          --filter "COL=VAL[,COL>=VAL…]" (repeatable)  --agg min|max|mean|median|first
#          --xscale/--yscale log  --title STR  --out FILE  --csv FILE
#          --keep-all (keep failed solves)  --contour  --no-diagonal  --list

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPPlots.jl"))
using .PEPPlots
using DataFrames, CSV, Printf

function parse_args(args)
    isempty(args) && error("usage: plot_results.jl RUN_DIR [options] (see the header)")
    opts = Dict{String,Any}("filter" => String[])
    run = nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if !startswith(a, "--")
            run === nothing || error("unexpected positional argument \"$a\"")
            run = a
            i += 1
        elseif a in ("--list", "--keep-all", "--contour", "--no-diagonal")
            opts[a[3:end]] = true
            i += 1
        else
            i < length(args) || error("option $a expects a value")
            key, val = a[3:end], args[i + 1]
            key == "filter" ? append!(opts["filter"], split(val, ',')) :
                              (opts[key] = val)
            i += 2
        end
    end
    run === nothing && error("missing RUN_DIR (a results directory or .arrow file)")
    return run, opts
end

getsym(opts, key) = haskey(opts, key) ? Symbol(opts[key]) : nothing

function main()
    run, opts = parse_args(ARGS)
    df = load_run(run)

    if get(opts, "list", false)
        println("run: ", run)
        describe_run(df)
        return
    end

    is_design = hasproperty(df, :alpha_best)
    mode = Symbol(get(opts, "mode", is_design ? "scatter" : "xy"))
    x = something(getsym(opts, "x"), is_design ? :alpha_init : :alpha)
    y = something(getsym(opts, "y"), is_design ? :alpha_best : :pep)
    z = getsym(opts, "z")
    series = getsym(opts, "series")
    agg = Symbol(get(opts, "agg", "min"))
    xscale = Symbol(get(opts, "xscale", "identity"))
    yscale = Symbol(get(opts, "yscale", "identity"))
    title = get(opts, "title", "")

    for c in filter(!isnothing, [x, y, z, series])
        hasproperty(df, c) ||
            error("unknown column :$c — available: $(join(names(df), ", "))")
    end

    df = clean_rows(df, (mode === :heatmap ? [z] : [x, y]);
                    keep_all = get(opts, "keep-all", false))
    df = apply_filters(df, parse_filter.(opts["filter"]))
    nrow(df) == 0 && error("no rows left after filtering/cleaning — relax the " *
                           "filters, or pass --keep-all to keep failed solves")

    plt, tbl = if mode === :xy
        xy_plot(df; x, y, series, agg, xscale, yscale, title)
    elseif mode === :heatmap
        z === nothing && error("--mode heatmap requires --z COL")
        heatmap_plot(df; x, y, z, agg, xscale, yscale, title,
                     contour = get(opts, "contour", false))
    elseif mode === :scatter
        scatter_plot(df; x, y, series, xscale, yscale, title,
                     diagonal = !get(opts, "no-diagonal", false))
    else
        error("unknown --mode $mode (xy, heatmap, scatter)")
    end

    name = mode === :heatmap ? "$(mode)_$(z)_over_$(x)_$(y)" : "$(mode)_$(y)_vs_$(x)"
    series === nothing || (name *= "_by_$(series)")
    out = get(opts, "out", joinpath(isdir(run) ? run : dirname(run), "plots", "$name.png"))
    save_plot(plt, out)
    haskey(opts, "csv") && CSV.write(opts["csv"], tbl)

    @info "plot_results: done" rows = nrow(df) figure = out csv = get(opts, "csv", nothing)
end

main()
