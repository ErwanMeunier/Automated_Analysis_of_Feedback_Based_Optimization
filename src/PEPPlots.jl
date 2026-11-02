# PEPPlots — figures from the stored Arrow sweeps.
#
# Deliberately independent of `PEPControl.jl` (and therefore of PEPDesign and
# Mosek): plotting a finished run must not pay a solver-stack load. The only
# shared knowledge is the shard layout written by `PEPControl.ResultWriter`.

module PEPPlots

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Arrow, DataFrames, LaTeXStrings, Plots, Printf, Statistics

export load_run, clean_rows, apply_filters, parse_filter, aggregate_xy,
       xy_plot, heatmap_plot, scatter_plot, save_plot, describe_run

const OK_STATUS = ("OPTIMAL", "SLOW_PROGRESS", "OK", "ALMOST_OPTIMAL")

# ── loading ──────────────────────────────────────────────────────────────────

"""
    load_run(path) -> DataFrame

Load a run directory (`merged.arrow` if present, otherwise all `part-*.arrow`
shards, including Slurm task subdirectories) or a single `.arrow` file.
"""
function load_run(path::AbstractString)
    if isfile(path)
        return DataFrame(Arrow.Table(path))
    end
    isdir(path) || error("no such run directory or Arrow file: $path")
    merged = joinpath(path, "merged.arrow")
    isfile(merged) && return DataFrame(Arrow.Table(merged))
    files = String[]
    for (root, _, names) in walkdir(path)
        append!(files, (joinpath(root, f) for f in names
                        if occursin(r"^part-\d+\.arrow$", f)))
    end
    sort!(files)
    isempty(files) && error("no Arrow shards under $path")
    return reduce(vcat, DataFrame(Arrow.Table(f)) for f in files)
end

"""Drop failed solves: bad `status` and non-finite values in `valuecols`."""
function clean_rows(df::DataFrame, valuecols; keep_all::Bool = false)
    keep_all && return df
    out = hasproperty(df, :status) ?
          df[[String(s) in OK_STATUS for s in df.status], :] : df
    for c in valuecols
        (c === nothing || !hasproperty(out, c)) && continue
        eltype(out[!, c]) <: Real || continue
        out = out[isfinite.(out[!, c]), :]
    end
    return out
end

# ── filtering ────────────────────────────────────────────────────────────────

"""Parse `"T=2"`, `"norm_A>=0.5"`, `"kappa!=1"` into `(column, op, value)`."""
function parse_filter(s::AbstractString)
    m = match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(>=|<=|!=|=|>|<)\s*(.+?)\s*$", s)
    m === nothing && error("cannot parse filter \"$s\" (use e.g. T=2, norm_A>=0.5)")
    col, op, raw = Symbol(m.captures[1]), Symbol(m.captures[2]), m.captures[3]
    val = something(tryparse(Float64, raw), raw)
    return (col, op, val)
end

# Float grid values are sampled, so `=` matches the nearest sampled value.
function _snap(col::AbstractVector{<:Real}, v::Real, name::Symbol)
    u = sort(unique(col))
    any(==(v), u) && return v
    i = argmin(abs.(u .- v))
    span = maximum(u) - minimum(u)
    tol = max(1e-8, 0.01 * (span == 0 ? abs(v) + 1 : span))
    abs(u[i] - v) <= tol || error("no sampled $name near $v; available: " *
        join((@sprintf("%.4g", w) for w in first(u, 12)), ", ") *
        (length(u) > 12 ? ", …" : ""))
    return u[i]
end

"""Apply `(column, op, value)` triples; `=` on float columns snaps to the grid."""
function apply_filters(df::DataFrame, filters)
    out = df
    for (col, op, val) in filters
        hasproperty(out, col) ||
            error("unknown column :$col (available: $(join(names(out), ", ")))")
        v = out[!, col]
        if val isa AbstractString || !(eltype(v) <: Real)
            sv = string(val)
            out = out[[(op === :!= ? string(x) != sv : string(x) == sv) for x in v], :]
            continue
        end
        target = op === :(=) ? _snap(v, val, col) : val
        pred = op === :(=)  ? (x -> x == target) :
               op === :!=   ? (x -> x != target) :
               op === :<    ? (x -> x < target)  :
               op === :<=   ? (x -> x <= target) :
               op === :>    ? (x -> x > target)  :
                              (x -> x >= target)
        out = out[pred.(v), :]
    end
    return out
end

# ── aggregation ──────────────────────────────────────────────────────────────

_reducer(agg::Symbol) =
    agg === :min ? minimum : agg === :max ? maximum : agg === :mean ? mean :
    agg === :median ? median : agg === :first ? first :
    error("unknown aggregator :$agg (min, max, mean, median, first)")

"""
    aggregate_xy(df, x, y; agg) -> DataFrame

Collapse duplicate `x` values (i.e. marginalize the axes not shown) with `agg`,
sorted by `x`. This is what turns "best PEP over α versus T" into one call.
"""
function aggregate_xy(df::DataFrame, x::Symbol, y::Symbol; agg::Symbol = :min)
    g = combine(groupby(df, x), y => _reducer(agg) => y)
    return sort!(g, x)
end

_scale(s::Symbol) = s === :log ? :log10 : s

function _series_groups(df::DataFrame, series)
    series === nothing && return [(nothing, df)]
    hasproperty(df, series) || error("unknown series column :$series")
    gs = groupby(df, series; sort = true)
    return [(first(g[!, series]), DataFrame(g)) for g in gs]
end

const LATEX_AXIS_LABELS = Dict(
    :T => L"T",
    :alpha => L"\alpha",
    :alpha_init => L"\alpha_{init}",
    :alpha_best => L"\alpha_{\mathrm{best}}",
    :pep => L"W(\alpha,T)",
    :pep_best => L"W(\alpha_{\mathrm{best}},T)",
    :norm_A => L"\Vert A\Vert",
    :norm_B => L"\Vert B\Vert",
    :kappa => L"\kappa",
    :mu_f => L"\mu_f",
    :L_f => L"L_f",
)

_latex_identifier(name) = latexstring("\\mathrm{" * replace(string(name), "_" => "\\_") * "}")
_axis_label(name::Symbol) = get(LATEX_AXIS_LABELS, name, _latex_identifier(name))
_latex_title(title::AbstractString) = isempty(title) ? "" :
    latexstring("\\mathrm{" * replace(title, " " => "\\ ", "_" => "\\_") * "}")
_latex_tick(x::Real) = latexstring(@sprintf("%.3g", x))

_label(series, v) = v === nothing ? "" :
    v isa Real ? latexstring(string(_axis_label(series)), " = ", @sprintf("%.4g", v)) :
    latexstring(string(_axis_label(series)), " = ", string(v))

function _sampled_ticks(values; max_ticks = 11, all_integer = false)
    vals = sort(unique(Float64.(filter(isfinite, values))))
    isempty(vals) && return :auto
    integer_axis = all_integer && all(isinteger, vals)
    if integer_axis || length(vals) <= max_ticks
        selected = vals
    else
        indices = unique(round.(Int, range(1, length(vals); length = max_ticks)))
        selected = vals[indices]
    end
    labels = integer_axis ? latexstring.(string.(Int.(round.(selected)))) :
                            _latex_tick.(selected)
    return (selected, labels)
end

function _linear_ticks(values; count = 11)
    vals = Float64.(filter(isfinite, values))
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    lo == hi && return ([lo], [_latex_tick(lo)])
    raw_step = (hi - lo) / max(count - 1, 1)
    magnitude = 10.0^floor(log10(raw_step))
    factor = first(candidate for candidate in (1.0, 2.0, 2.5, 5.0, 10.0)
                   if candidate * magnitude >= raw_step)
    step = factor * magnitude
    tick_lo = lo >= 0 ? 0.0 : floor(lo / step) * step
    tick_hi = ceil(hi / step) * step
    selected = collect(tick_lo:step:tick_hi)
    return (selected, _latex_tick.(selected))
end

function _log_ticks(values)
    vals = Float64.(filter(x -> isfinite(x) && x > 0, values))
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    lo == hi && return ([lo], [_latex_tick(lo)])
    selected = Float64[]
    for exponent in floor(Int, log10(lo)):ceil(Int, log10(hi)), multiplier in (1.0, 2.0, 5.0)
        tick = multiplier * 10.0^exponent
        lo <= tick <= hi && push!(selected, tick)
    end
    isempty(selected) && (selected = collect(exp10.(range(log10(lo), log10(hi); length = 7))))
    return (selected, _latex_tick.(selected))
end

_axis_ticks(values, scale::Symbol; all_integer = false) =
    scale === :log ? _log_ticks(values) :
    all_integer ? _sampled_ticks(values; all_integer = true) : _linear_ticks(values)

_tick_limits(ticks, scale::Symbol) =
    scale === :log || !(ticks isa Tuple) ? :auto : extrema(first(ticks))

# ── plots ────────────────────────────────────────────────────────────────────

"""
    xy_plot(df; x, y, series, agg, ...) -> (plot, plotted_table)

One curve of `y` against `x` per value of `series`, with duplicate `x` values
reduced by `agg`.
"""
function xy_plot(df::DataFrame; x::Symbol, y::Symbol, series = nothing,
                 agg::Symbol = :min, xscale::Symbol = :identity,
                 yscale::Symbol = :identity, title::String = "",
                 marker::Bool = true)
    plt = plot(; xlabel = _axis_label(x), ylabel = _axis_label(y),
               title = _latex_title(title), fontfamily = "Computer Modern",
               xscale = _scale(xscale), yscale = _scale(yscale),
               legend = series === nothing ? :none : :best)
    tables = DataFrame[]
    for (v, sub) in _series_groups(df, series)
        g = aggregate_xy(sub, x, y; agg)
        nrow(g) == 0 && continue
        plot!(plt, g[!, x], g[!, y]; label = _label(series, v),
              marker = marker ? :circle : :none, markersize = 3)
        series === nothing || (g[!, series] = fill(v, nrow(g)))
        push!(tables, g)
    end
    isempty(tables) && error("nothing to plot: every series is empty")
    table = reduce(vcat, tables)
    integer_x = eltype(table[!, x]) <: Integer || all(isinteger, table[!, x])
    x_ticks = _axis_ticks(table[!, x], xscale; all_integer = integer_x)
    y_ticks = _axis_ticks(table[!, y], yscale)
    plot!(plt; xticks = x_ticks, xlims = _tick_limits(x_ticks, xscale),
          xrotation = integer_x && length(unique(table[!, x])) > 12 ? 45.0 : 0.0,
          yticks = y_ticks, ylims = _tick_limits(y_ticks, yscale))
    return plt, table
end

"""
    heatmap_plot(df; x, y, z, agg, ...) -> (plot, plotted_table)

`z` over the (`x`, `y`) grid; cells with no sample stay `NaN` (blank).
"""
function heatmap_plot(df::DataFrame; x::Symbol, y::Symbol, z::Symbol,
                      agg::Symbol = :min, xscale::Symbol = :identity,
                      yscale::Symbol = :identity, title::String = "",
                      contour::Bool = false)
    g = combine(groupby(df, [x, y]), z => _reducer(agg) => z)
    nrow(g) == 0 && error("nothing to plot: the filtered table is empty")
    xs, ys = sort(unique(g[!, x])), sort(unique(g[!, y]))
    xi = Dict(v => i for (i, v) in enumerate(xs))
    yi = Dict(v => i for (i, v) in enumerate(ys))
    Z = fill(NaN, length(ys), length(xs))
    for r in eachrow(g)
        Z[yi[r[y]], xi[r[x]]] = r[z]
    end
    f = contour ? contourf : Plots.heatmap
    integer_x = eltype(xs) <: Integer || all(isinteger, xs)
    plt = f(xs, ys, Z; xlabel = _axis_label(x), ylabel = _axis_label(y),
            colorbar_title = _axis_label(z), title = _latex_title(title),
            fontfamily = "Computer Modern",
            xticks = _axis_ticks(xs, xscale; all_integer = integer_x),
            xrotation = integer_x && length(xs) > 12 ? 45.0 : 0.0,
            yticks = _axis_ticks(ys, yscale),
            xscale = _scale(xscale), yscale = _scale(yscale))
    return plt, sort!(g, [x, y])
end

"""Raw scatter of `y` against `x` (no aggregation), optionally with `y = x`."""
function scatter_plot(df::DataFrame; x::Symbol, y::Symbol, series = nothing,
                      xscale::Symbol = :identity, yscale::Symbol = :identity,
                      title::String = "", diagonal::Bool = true)
    plt = plot(; xlabel = _axis_label(x), ylabel = _axis_label(y),
               title = _latex_title(title), fontfamily = "Computer Modern",
               xscale = _scale(xscale), yscale = _scale(yscale),
               legend = series === nothing ? :none : :best)
    for (v, sub) in _series_groups(df, series)
        nrow(sub) == 0 && continue
        scatter!(plt, sub[!, x], sub[!, y]; label = _label(series, v),
                 markersize = 4, markerstrokewidth = 0)
    end
    if diagonal && nrow(df) > 0
        lo = min(minimum(df[!, x]), minimum(df[!, y]))
        hi = max(maximum(df[!, x]), maximum(df[!, y]))
        plot!(plt, [lo, hi], [lo, hi]; linestyle = :dash, color = :black,
              label = L"y=x")
    end
    integer_x = eltype(df[!, x]) <: Integer || all(isinteger, df[!, x])
    x_ticks = _axis_ticks(df[!, x], xscale; all_integer = integer_x)
    y_ticks = _axis_ticks(df[!, y], yscale)
    plot!(plt; xticks = x_ticks, xlims = _tick_limits(x_ticks, xscale),
          xrotation = integer_x && length(unique(df[!, x])) > 12 ? 45.0 : 0.0,
          yticks = y_ticks, ylims = _tick_limits(y_ticks, yscale))
    return plt, select(df, unique(filter(!isnothing, [x, y, series])))
end

function save_plot(plt, out::AbstractString)
    mkpath(dirname(out))
    savefig(plt, out)
    return out
end

# ── discovery ────────────────────────────────────────────────────────────────

"""Print the schema, status histogram and the sampled values of each grid axis."""
function describe_run(df::DataFrame)
    println("rows: ", nrow(df))
    println("columns:")
    for n in names(df)
        @printf("  %-12s %s\n", n, eltype(df[!, n]))
    end
    if hasproperty(df, :status)
        println("status:")
        for r in eachrow(sort!(combine(groupby(df, :status), nrow), :nrow; rev = true))
            @printf("  %-40s %d\n", first(String(r.status), 40), r.nrow)
        end
    end
    println("grid axes:")
    for n in names(df)
        c = df[!, n]
        eltype(c) <: Real || continue
        u = sort(unique(c))
        (1 < length(u) <= 40) || continue
        @printf("  %-12s (%d) %s\n", n, length(u),
                join((@sprintf("%.4g", v) for v in u), ", "))
    end
    return nothing
end

end # module
