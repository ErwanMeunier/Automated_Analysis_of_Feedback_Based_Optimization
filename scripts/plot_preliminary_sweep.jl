# Paper-style figures from sweep_alpha Arrow shards, optionally refined with a
# complete FOM design_best_alpha run.
#
# This script never loads or modifies PEPControl/PEPDesign. It approximates
# alpha_best by the minimum over the sampled alpha grid when no design run is
# supplied, and excludes groups whose alpha sweep is incomplete.
#
# Usage:
#   julia --project=Julia Julia/scripts/plot_preliminary_sweep.jl RUN_DIR \
#       [--design DESIGN_DIR] [--out DIR] [--kappa VALUE] \
#       [--analytic ANALYTIC_RESULTS] \
#       [--analytic-mode former|recovered|theorem1] \
#       [--include-uncertified-analytic] \
#       [--mark-uncertified-analytic] \
#       [--include-design-only] \
#       [--colormap MATPLOTLIB_NAME] \
#       [--normalize-worst-cases] \
#       [--sweep-ymax VALUE|auto] \
#       [--no-titles] \
#       [--title-font-size N] [--axis-font-size N] \
#       [--legend-font-size N] [--tick-font-size N] \
#       [--star-size N] [--line-width N] \
#       [--horizons 10,15,20] [--curve-T 10] \
#       [--T-plot-horizons all|3,4,5,15,20,30,40,50] \
#       [--split-horizon-panels] \
#       [--norm-a 0.2,0.4,0.6,0.8] [--allow-partial] [--strict-design]
#
# When RUN_DIR is a merged_sweep.arrow file, a sibling merged_design.arrow is
# used automatically. By default, unmatched design groups are excluded from
# design-derived figures; pass --strict-design to require complete matching.
# In `theorem1` mode (the default), the sharp uniform Theorem-1 step size and
# its fresh PEP evaluation are detected automatically. `recovered` loads the
# plant-specific result obtained from worst-case matrices. `former` uses the
# same norm-bound step size as Theorem 1 but retains the legacy nearest sampled
# sweep evaluation.
# `--include-design-only` includes finite design groups that have no complete
# alpha sweep in the best-point and horizon figures. It does not invent sweep
# samples or add them to the alpha-sweep figure.
# `--include-uncertified-analytic` includes finite but uncertified analytical
# PEP estimates. Included uncertified points always receive red cross markers.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPPlots.jl"))
using .PEPPlots
using CSV, DataFrames, LaTeXStrings, Plots, Printf, Statistics

const GROUP_KEY = [:T, :norm_A, :kappa]

function parse_args(args)
    isempty(args) && error("usage: plot_preliminary_sweep.jl RUN_DIR [options]")
    run = args[1]
    opts = Dict{String,Any}(
        "horizons" => "5,10,20,30,40,50",
        "T-plot-horizons" => "all",
        "split-horizon-panels" => false,
        "curve-T" => "20",
        "kappa" => "0",
        "norm-a" => "0.2,0.4,0.6,0.8",
        "analytic-mode" => "theorem1",
        "colormap" => "viridis",
        "normalize-worst-cases" => false,
        "include-uncertified-analytic" => false,
        "mark-uncertified-analytic" => false,
        "include-design-only" => false,
        "no-titles" => false,
        "sweep-ymax" => "default",
        "title-font-size" => "16",
        "axis-font-size" => "16",
        "legend-font-size" => "11",
        "tick-font-size" => "12",
        "star-size" => "7",
        "line-width" => "2",
        "allow-partial-design" => true,
    )
    i = 2
    while i <= length(args)
        arg = args[i]
        if arg in ("--allow-partial", "--allow-partial-design", "--sweep-only",
                   "--normalize-worst-cases", "--no-titles",
                   "--include-uncertified-analytic", "--include-design-only",
                   "--mark-uncertified-analytic", "--split-horizon-panels")
            opts[arg[3:end]] = true
            i += 1
        elseif arg == "--strict-design"
            opts["allow-partial-design"] = false
            i += 1
        elseif startswith(arg, "--")
            i < length(args) || error("$arg expects a value")
            opts[arg[3:end]] = args[i + 1]
            i += 2
        else
            error("unexpected argument $arg")
        end
    end
    return abspath(run), opts
end

function normalize_worst_case_table(df, columns)
    required = vcat([:norm_A, :R_u], collect(columns))
    missingcols = filter(c -> !hasproperty(df, c), required)
    isempty(missingcols) ||
        error("cannot normalize worst cases; table lacks columns $missingcols")

    result = copy(df)
    norm_A = Float64.(result.norm_A)
    rho = Float64.(result.R_u)
    valid = isfinite.(norm_A) .& isfinite.(rho) .& (rho .> 0)
    all(valid) || error(
        "cannot normalize worst cases: norm_A must be finite and rho=R_u " *
        "must be finite and positive (invalid rows: $(count(.!valid)))")

    # W(0,T) = rho^2 / (2(1-L_A)^2). Multiplying by its reciprocal
    # remains well defined at L_A=1, where the baseline diverges and the
    # normalized value has limiting value zero.
    result.W_zero = rho .^ 2 ./ (2 .* (1 .- norm_A) .^ 2)
    multiplier = 2 .* (1 .- norm_A) .^ 2 ./ rho .^ 2
    for column in columns
        raw_column = Symbol(column, "_unnormalized")
        result[!, raw_column] = copy(result[!, column])
        result[!, column] = Float64.(result[!, column]) .* multiplier
    end
    return result
end

worst_case_label(normalized) = normalized ?
    L"W(\alpha,T)/W(0,T)" : L"W(\alpha,T)"

csv_floats(s) = parse.(Float64, filter(!isempty, strip.(split(s, ','))))
csv_ints(s) = parse.(Int, filter(!isempty, strip.(split(s, ','))))

function parse_sweep_ymax(value; normalized = false)
    text = lowercase(strip(String(value)))
    text == "default" && return normalized ? nothing : 2.5
    text in ("auto", "none") && return nothing
    ymax = tryparse(Float64, text)
    ymax === nothing && error(
        "--sweep-ymax must be a positive number or 'auto', got '$value'")
    isfinite(ymax) && ymax > 0 || error(
        "--sweep-ymax must be finite and positive, got '$value'")
    return ymax
end

function parse_positive_plot_value(opts, name)
    value = tryparse(Float64, String(opts[name]))
    value !== nothing && isfinite(value) && value > 0 ||
        error("--$name must be finite and positive, got '$(opts[name])'")
    return value
end

function parse_plot_style(opts)
    return (
        show_titles = !opts["no-titles"],
        title_font_size = parse_positive_plot_value(opts, "title-font-size"),
        axis_font_size = parse_positive_plot_value(opts, "axis-font-size"),
        legend_font_size = parse_positive_plot_value(opts, "legend-font-size"),
        tick_font_size = parse_positive_plot_value(opts, "tick-font-size"),
        star_size = parse_positive_plot_value(opts, "star-size"),
        line_width = parse_positive_plot_value(opts, "line-width"),
    )
end

# ColorSchemes.jl renames a few Matplotlib maps to avoid collisions with
# existing Julia schemes. Keep the command-line interface aligned with the
# names users see in Matplotlib's colormap catalogue.
const MATPLOTLIB_COLORMAP_ALIASES = Dict(
    "autumn" => :autumn1,
    "gray" => :gray1,
    "jet" => :jet1,
    "rainbow" => :rainbow1,
)

function matplotlib_palette(name, n::Integer)
    n > 0 || error("a colormap must be sampled at least once")
    requested = strip(String(name))
    isempty(requested) && error("--colormap cannot be empty")
    reversed = endswith(requested, "_r")
    base = reversed ? chop(requested; tail = 2) : requested

    # Stop in Viridis' green region: the top 20% is already yellowish in
    # rendered paper figures. Reversal uses the same cropped colors.
    if base == "viridis"
        gradient = cgrad(:viridis)
        positions = range(0.0, 0.80; length = n)
        colors = [gradient[position] for position in positions]
        reversed && reverse!(colors)
        return colors
    end

    # ColorSchemes ships Matplotlib's twilight map but not its half-cycle
    # shifted alias, so construct that standard Matplotlib variant directly.
    if base == "twilight_shifted"
        gradient = cgrad(:twilight)
        positions = range(0.0, 1.0; length = n)
        colors = [gradient[mod(position + 0.5, 1.0)] for position in positions]
        reversed && reverse!(colors)
        return colors
    end

    scheme = get(MATPLOTLIB_COLORMAP_ALIASES, base, Symbol(base))
    try
        return palette(scheme, n; rev = reversed)
    catch err
        error("unknown Matplotlib colormap '$requested'; names are case-sensitive " *
              "and may use the '_r' suffix for reversal (for example, plasma_r). " *
              "Underlying Plots.jl error: $(sprint(showerror, err))")
    end
end

function config_int(run, name)
    isdir(run) || return nothing
    file = joinpath(run, "config.txt")
    isfile(file) || return nothing
    m = match(Regex("\\b$name\\s*=\\s*(\\d+)"), read(file, String))
    return m === nothing ? nothing : parse(Int, m.captures[1])
end

function config_int_list(run, name)
    isdir(run) || return nothing
    file = joinpath(run, "config.txt")
    isfile(file) || return nothing
    m = match(Regex("\\b$name\\s*=\\s*\\[([^\\]]*)\\]"), read(file, String))
    m === nothing && return nothing
    values = filter(!isempty, strip.(split(m.captures[1], ',')))
    return parse.(Int, values)
end

function expected_shape(run, df)
    tlo, thi = config_int(run, "T_lo"), config_int(run, "T_hi")
    tvalues = config_int_list(run, "T_values")
    nA = config_int(run, "n_A")
    nalpha = config_int(run, "n_alpha")
    nkappa = config_int(run, "n_kappa")
    observed_counts = combine(groupby(df, GROUP_KEY),
        :alpha => (x -> length(unique(x))) => :n_alpha)
    inferred_nalpha = isempty(observed_counts.n_alpha) ? 0 : maximum(observed_counts.n_alpha)
    nalpha = something(nalpha, inferred_nalpha)
    nt = tvalues === nothing ?
         (any(isnothing, (tlo, thi)) ? nothing : thi - tlo + 1) :
         length(tvalues)
    expected_groups = any(isnothing, (nt, nA, nkappa)) ? nothing :
                      nt * nA * nkappa
    expected_rows = expected_groups === nothing ? nothing : expected_groups * nalpha
    return (; nalpha, expected_groups, expected_rows, observed_counts)
end

keytuple(row) = (Int(row.T), Float64(row.norm_A), Float64(row.kappa))

function complete_sweep_rows(raw, nalpha)
    counts = combine(groupby(raw, GROUP_KEY),
        :alpha => (x -> length(unique(x))) => :n_alpha)
    complete = Set(keytuple(row) for row in eachrow(counts) if row.n_alpha == nalpha)
    mask = [keytuple(row) in complete for row in eachrow(raw)]
    return raw[mask, :], counts, complete
end

function analytic_alpha(norm_A, L_f, norm_B = 1.0)
    a = Float64(norm_A)
    L = Float64(L_f)
    b = Float64(norm_B)
    (!isfinite(a) || !isfinite(L) || !isfinite(b) ||
     L <= 0 || a >= 1 || b <= 0) && return NaN
    # The legacy expression simplifies exactly to the sharp uniform threshold
    # in Theorem 1. Legacy mode differs only in evaluating it at the nearest
    # sampled sweep alpha rather than by a fresh PEP solve.
    return 2 * (1 - a)^3 / (L * b^2 * (3 - a))
end

function sampled_best_table(df)
    rows = NamedTuple[]
    for group in groupby(df, GROUP_KEY; sort = true)
        isempty(group.pep) && continue
        best_i = argmin(group.pep)
        ref_alpha = analytic_alpha(first(group.norm_A), first(group.L_f),
                                   first(group.norm_B))
        ref_i = isfinite(ref_alpha) ? argmin(abs.(group.alpha .- ref_alpha)) : nothing
        push!(rows, (
            T = Int(first(group.T)),
            norm_A = Float64(first(group.norm_A)),
            kappa = Float64(first(group.kappa)),
            mu_f = Float64(first(group.mu_f)),
            L_f = Float64(first(group.L_f)),
            norm_B = Float64(first(group.norm_B)),
            R_x = Float64(first(group.R_x)),
            R_u = Float64(first(group.R_u)),
            alpha_best_sampled = Float64(group.alpha[best_i]),
            pep_best_sampled = Float64(group.pep[best_i]),
            alpha_analytic_norm_bound = ref_alpha,
            alpha_analytic_sampled = ref_i === nothing ? NaN : Float64(group.alpha[ref_i]),
            pep_analytic_sampled = ref_i === nothing ? NaN : Float64(group.pep[ref_i]),
        ))
    end
    return DataFrame(rows)
end

function pep_analytic_table(df; include_uncertified = false)
    required = [:T, :norm_A, :kappa, :alpha_analytic, :pep_analytic, :status]
    missingcols = filter(c -> !hasproperty(df, c), required)
    isempty(missingcols) || error("analytical-from-PEP data lack columns $missingcols")
    statuses = String.(df.status)
    valid = ((statuses .== "OK") .|
             (include_uncertified .& (statuses .== "UNCERTIFIED"))) .&
            isfinite.(df.alpha_analytic) .& (df.alpha_analytic .> 0) .&
            isfinite.(df.pep_analytic)
    clean = df[valid, :]
    rows = NamedTuple[]
    for group in groupby(clean, GROUP_KEY; sort = true)
        row = last(group)
        parameter(field) = hasproperty(row, field) ? Float64(row[field]) : NaN
        push!(rows, (
            T = Int(row.T), norm_A = Float64(row.norm_A),
            kappa = Float64(row.kappa),
            mu_f = parameter(:mu_f), L_f = parameter(:L_f),
            norm_B = parameter(:norm_B), R_x = parameter(:R_x),
            R_u = parameter(:R_u),
            alpha_analytic_pep = Float64(row.alpha_analytic),
            pep_analytic_pep = Float64(row.pep_analytic),
            analytic_certified_pep = String(row.status) == "OK",
        ))
    end
    return DataFrame(rows)
end

function attach_pep_analytical(best, analytical, source)
    result = copy(best)
    result.alpha_analytic = fill(NaN, nrow(result))
    result.pep_analytic = fill(NaN, nrow(result))
    result.analytic_certified = fill(false, nrow(result))
    result.analytic_source = fill("unavailable", nrow(result))
    analytical === nothing && return result, 0

    lookup = Dict(keytuple(row) => row for row in eachrow(analytical))
    matched = 0
    for row in eachrow(result)
        evaluated = get(lookup, keytuple(row), nothing)
        evaluated === nothing && continue
        row.alpha_analytic = evaluated.alpha_analytic_pep
        row.pep_analytic = evaluated.pep_analytic_pep
        row.analytic_certified = evaluated.analytic_certified_pep
        row.analytic_source = evaluated.analytic_certified_pep ?
            source : source * " (uncertified)"
        matched += 1
    end
    @info "matched analytical PEP results" source matched total = nrow(result)
    return result, matched
end

function attach_former_analytical(best)
    result = copy(best)
    result.alpha_analytic = copy(result.alpha_analytic_norm_bound)
    result.pep_analytic = copy(result.pep_analytic_sampled)
    result.analytic_certified = fill(true, nrow(result))
    result.analytic_source = fill(
        "former norm-bound formula + nearest sampled sweep", nrow(result))
    matched = count(isfinite.(result.alpha_analytic) .&
                    isfinite.(result.pep_analytic))
    @info "using former norm-bound analytical step size" matched total = nrow(result)
    return result, matched
end

function parse_analytic_mode(value)
    mode = Symbol(lowercase(strip(String(value))))
    mode in (:former, :recovered, :theorem1) ||
        error("--analytic-mode must be former, recovered, or theorem1, got '$value'")
    return mode
end

function fom_design_table(df; include_parameters = false)
    required = [:T, :norm_A, :kappa, :optimizer, :alpha_best, :pep_best]
    include_parameters && append!(required, [:mu_f, :L_f, :norm_B, :R_x, :R_u])
    missingcols = filter(c -> !hasproperty(df, c), required)
    isempty(missingcols) || error("design data lack columns $missingcols")
    fom = df[lowercase.(String.(df.optimizer)) .== "fom", :]
    fom = clean_rows(fom, [:alpha_best, :pep_best])
    rows = NamedTuple[]
    for group in groupby(fom, GROUP_KEY; sort = true)
        best_i = argmin(group.pep_best)
        parameter(field) = hasproperty(group, field) ? Float64(first(group[!, field])) : NaN
        push!(rows, (
            T = Int(first(group.T)),
            norm_A = Float64(first(group.norm_A)),
            kappa = Float64(first(group.kappa)),
            mu_f = parameter(:mu_f),
            L_f = parameter(:L_f),
            norm_B = parameter(:norm_B),
            R_x = parameter(:R_x),
            R_u = parameter(:R_u),
            alpha_best_design = Float64(group.alpha_best[best_i]),
            pep_best_design = Float64(group.pep_best[best_i]),
        ))
    end
    return DataFrame(rows)
end

function select_best_points(sampled, design;
                            allow_partial_design = false,
                            include_design_only = false)
    best = copy(sampled)
    best.alpha_best = copy(best.alpha_best_sampled)
    best.pep_best = copy(best.pep_best_sampled)
    best.best_source = fill("sampled sweep", nrow(best))
    design === nothing && return best, 0

    lookup = Dict(keytuple(row) => row for row in eachrow(design))
    matched = 0
    keep = falses(nrow(best))
    for (i, row) in enumerate(eachrow(best))
        drow = get(lookup, keytuple(row), nothing)
        drow === nothing && continue
        row.alpha_best = drow.alpha_best_design
        row.pep_best = drow.pep_best_design
        row.best_source = "FOM design"
        keep[i] = true
        matched += 1
    end
    unmatched = nrow(best) - matched
    if unmatched > 0 && !allow_partial_design
        error("FOM design covers $matched/$(nrow(best)) valid sweep groups; " *
              "refusing partial design data (omit --design or pass --allow-partial-design)")
    end
    unmatched > 0 && @warn "excluding sweep groups without a matched finite FOM design result" matched unmatched
    matched > 0 || error("no sweep groups have a matched finite FOM design result")
    result = best[keep, :]
    if include_design_only
        sampled_keys = Set(keytuple(row) for row in eachrow(sampled))
        for drow in eachrow(design)
            keytuple(drow) in sampled_keys && continue
            push!(result, (
                T = drow.T, norm_A = drow.norm_A, kappa = drow.kappa,
                mu_f = drow.mu_f, L_f = drow.L_f,
                norm_B = drow.norm_B, R_x = drow.R_x, R_u = drow.R_u,
                alpha_best_sampled = NaN, pep_best_sampled = NaN,
                alpha_analytic_norm_bound = analytic_alpha(drow.norm_A,
                                                           drow.L_f,
                                                           drow.norm_B),
                alpha_analytic_sampled = NaN, pep_analytic_sampled = NaN,
                alpha_best = drow.alpha_best_design,
                pep_best = drow.pep_best_design,
                best_source = "FOM design",
            ))
        end
        added = nrow(result) - matched
        added > 0 && @info "included FOM design groups without complete alpha sweeps" added
        sort!(result, GROUP_KEY)
    end
    return result, matched
end

function add_analytic_only_sampled(best, sampled, analytical)
    analytical === nothing && return best, 0
    analytic_keys = Set(keytuple(row) for row in eachrow(analytical))
    best_keys = Set(keytuple(row) for row in eachrow(best))
    missing = [keytuple(row) in analytic_keys && !(keytuple(row) in best_keys)
               for row in eachrow(sampled)]
    added = copy(sampled[missing, :])
    result = copy(best)
    if nrow(added) > 0
        added.alpha_best = copy(added.alpha_best_sampled)
        added.pep_best = copy(added.pep_best_sampled)
        added.best_source = fill("sampled sweep (design unavailable)", nrow(added))
        result = vcat(result, added; cols = :setequal)
        @warn "including analytical PEP groups without finite FOM design; sampled sweep used for comparison" added = nrow(added)
    end
    result_keys = Set(keytuple(row) for row in eachrow(result))
    analytic_only = 0
    for row in eachrow(analytical)
        keytuple(row) in result_keys && continue
        all(isfinite, (row.mu_f, row.L_f, row.norm_B, row.R_x, row.R_u)) ||
            error("analytical-only group lacks parameters needed for plotting: $(keytuple(row))")
        push!(result, (
            T = row.T, norm_A = row.norm_A, kappa = row.kappa,
            mu_f = row.mu_f, L_f = row.L_f,
            norm_B = row.norm_B, R_x = row.R_x, R_u = row.R_u,
            alpha_best_sampled = NaN, pep_best_sampled = NaN,
            alpha_analytic_norm_bound = analytic_alpha(row.norm_A, row.L_f,
                                                       row.norm_B),
            alpha_analytic_sampled = NaN, pep_analytic_sampled = NaN,
            alpha_best = NaN, pep_best = NaN,
            best_source = "analytic PEP only (no finite best)",
        ))
        push!(result_keys, keytuple(row))
        analytic_only += 1
    end
    analytic_only > 0 &&
        @warn "including analytical PEP groups without a finite design or sampled comparison" added = analytic_only
    sort!(result, GROUP_KEY)
    return result, nrow(added) + analytic_only
end

function nearest_value(values, target; label = "value", tolerance = Inf)
    vals = sort(unique(Float64.(values)))
    isempty(vals) && return nothing
    value = vals[argmin(abs.(vals .- target))]
    abs(value - target) <= tolerance || return nothing
    value == target || @warn "snapped requested $label to sampled grid" requested = target selected = value
    return value
end

function style_kwargs(; xlabel = L"", ylabel = L"", title = L"",
                      plot_style, size = (800, 500))
    # An empty LaTeXString renders as literal "$$" in the GR backend.
    title = plot_style.show_titles ? title : ""
    return (; xlabel, ylabel, title, grid = true, framestyle = :box,
            linewidth = plot_style.line_width, markersize = 7,
            markerstrokewidth = 1,
            markerstrokecolor = :black,
            guidefontsize = plot_style.axis_font_size,
            tickfontsize = plot_style.tick_font_size,
            legendfontsize = plot_style.legend_font_size,
            titlefontsize = plot_style.title_font_size,
            fontfamily = "Computer Modern",
            left_margin = 8 * Plots.mm, right_margin = 4 * Plots.mm,
            bottom_margin = 13 * Plots.mm, top_margin = 3 * Plots.mm,
            size, dpi = 180)
end

latex_tick(x::Real) = latexstring(@sprintf("%.3g", x))

function sampled_ticks(values; max_ticks = 11, all_integer = false)
    vals = sort(unique(Float64.(filter(isfinite, values))))
    isempty(vals) && return :auto
    integer_axis = all_integer && all(isinteger, vals)
    if length(vals) <= max_ticks
        selected = vals
    else
        indices = unique(round.(Int, range(1, length(vals); length = max_ticks)))
        selected = vals[indices]
    end
    labels = integer_axis ? latexstring.(string.(Int.(round.(selected)))) :
                            latex_tick.(selected)
    return (selected, labels)
end

function linear_ticks(values; count = 11)
    vals = Float64.(filter(isfinite, values))
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    lo == hi && return ([lo], [latex_tick(lo)])
    raw_step = (hi - lo) / max(count - 1, 1)
    magnitude = 10.0^floor(log10(raw_step))
    factor = first(candidate for candidate in (1.0, 2.0, 2.5, 5.0, 10.0)
                   if candidate * magnitude >= raw_step)
    step = factor * magnitude
    tick_lo = lo >= 0 ? 0.0 : floor(lo / step) * step
    tick_hi = ceil(hi / step) * step
    selected = collect(tick_lo:step:tick_hi)
    return (selected, latex_tick.(selected))
end

function log_ticks(values)
    vals = Float64.(filter(x -> isfinite(x) && x > 0, values))
    isempty(vals) && return :auto
    lo, hi = extrema(vals)
    lo == hi && return ([lo], [latex_tick(lo)])
    candidates = Float64[]
    for exponent in floor(Int, log10(lo)):ceil(Int, log10(hi)), multiplier in (1.0, 2.0, 5.0)
        tick = multiplier * 10.0^exponent
        lo <= tick <= hi && push!(candidates, tick)
    end
    isempty(candidates) && (candidates = collect(exp10.(range(log10(lo), log10(hi); length = 7))))
    return (candidates, latex_tick.(candidates))
end

tick_limits(ticks) = ticks isa Tuple ? extrema(first(ticks)) : :auto

function save_pdf_png(plt, outdir, stem)
    png = joinpath(outdir, stem * ".png")
    pdf = joinpath(outdir, stem * ".pdf")
    savefig(plt, png)
    try
        savefig(plt, pdf)
    catch err
        @warn "could not overwrite PDF; close any PDF viewer and rerun (PNG was written)" pdf exception = (err, catch_backtrace())
    end
end

function paper_figure2(best, horizons, kappa, outdir;
                       partial = false, colormap = "viridis",
                       normalized = false, mark_uncertified = false,
                       plot_style)
    sub = best[(best.kappa .== kappa) .& in.(best.T, Ref(horizons)), :]
    nrow(sub) > 0 || error("no complete groups for Figure 2 selections")
    colors = matplotlib_palette(colormap, max(length(horizons), 2))
    p_alpha = plot(; style_kwargs(xlabel = L"L_A", ylabel = L"\alpha",
                                   title = L"\mathrm{Best\ step\ size\ versus\ system\ dynamics}",
                                   plot_style = plot_style)...,
                   yscale = :log10, legend = :bottomleft)
    p_pep = plot(; style_kwargs(xlabel = L"L_A",
                                  ylabel = worst_case_label(normalized),
                                  title = L"\mathrm{Best\ worst\ case\ versus\ system\ dynamics}",
                                  plot_style = plot_style,
                                  size = (800, 575))...,
                 yscale = :log10, legend = :top, legend_column = 2)
    analytic_pep_plotted = false
    uncertified_pep_plotted = false
    for (j, T) in enumerate(horizons)
        curve = sort(sub[sub.T .== T, :], :norm_A)
        nrow(curve) == 0 && (@warn "no complete Figure 2 groups" T; continue)
        plot!(p_alpha, curve.norm_A, curve.alpha_best;
              color = colors[j], marker = :star5,
              markersize = plot_style.star_size,
              label = latexstring(raw"\alpha_{\mathrm{best}};\ T=", T))
        plot!(p_pep, curve.norm_A, curve.pep_best;
              color = colors[j], marker = :star5,
              markersize = plot_style.star_size,
              label = latexstring(raw"W(\alpha_{\mathrm{best}},T=", T, ")"))
        analytic_pep_curve = curve[isfinite.(curve.pep_analytic) .&
                                   (curve.pep_analytic .> 0), :]
        if nrow(analytic_pep_curve) > 0
            plot!(p_pep, analytic_pep_curve.norm_A,
                  analytic_pep_curve.pep_analytic;
                  color = colors[j], linestyle = :dash, marker = :diamond,
                  markersize = 5, label = "")
            analytic_pep_plotted = true
            uncertified_curve = analytic_pep_curve[.!analytic_pep_curve.analytic_certified, :]
            if mark_uncertified && nrow(uncertified_curve) > 0
                scatter!(p_pep, uncertified_curve.norm_A,
                         uncertified_curve.pep_analytic;
                         color = :red, marker = :xcross,
                         markersize = plot_style.star_size + 2,
                         markerstrokewidth = 2, label = "")
                uncertified_pep_plotted = true
            end
        end
    end
    analytic_pep_plotted &&
        plot!(p_pep, [NaN], [NaN]; color = :black, linestyle = :dash,
              marker = :diamond, markersize = 5,
              label = L"W(\alpha_{\mathrm{analytic}},T=\cdot)")
    uncertified_pep_plotted &&
        scatter!(p_pep, [NaN], [NaN]; color = :red, marker = :xcross,
                 markersize = plot_style.star_size + 2,
                 label = "uncertified PEP")
    # The analytical step size is horizon-invariant, so show one reference
    # curve rather than repeating the same line and legend entry for every T.
    analytic_alpha_curve = sub[isfinite.(sub.alpha_analytic) .&
                               (sub.alpha_analytic .> 0),
                               [:T, :norm_A, :alpha_analytic]]
    sort!(analytic_alpha_curve, [:norm_A, :T])
    unique!(analytic_alpha_curve, :norm_A)
    nrow(analytic_alpha_curve) > 0 &&
        plot!(p_alpha, analytic_alpha_curve.norm_A,
              analytic_alpha_curve.alpha_analytic;
              color = :red, linestyle = :dash,
              label = L"\alpha_{\mathrm{analytic}}")
    norm_A_ticks = linear_ticks(sub.norm_A)
    plot!(p_alpha; xticks = norm_A_ticks,
          xlims = (minimum(sub.norm_A) - 0.015, maximum(sub.norm_A) + 0.015))
    positive_pep = filter(x -> isfinite(x) && x > 0,
                          vcat(sub.pep_best, sub.pep_analytic))
    isempty(positive_pep) && error("no positive worst-case values for Figure 2")
    # Leave room on the log axis for the centered legend without covering the
    # performance curves in either the standalone or combined figure.
    legend_entries = length(horizons) + (analytic_pep_plotted ? 1 : 0) +
                     (uncertified_pep_plotted ? 1 : 0)
    upper_padding = max(3.0, 1.1 * ceil(Int, legend_entries / 2))
    plot!(p_pep; xticks = norm_A_ticks,
          xlims = (minimum(sub.norm_A) - 0.02,
                   maximum(sub.norm_A) + 0.04),
          yticks = log_ticks(positive_pep),
          ylims = (minimum(positive_pep) * 0.8,
                   maximum(positive_pep) * upper_padding))
    combined = plot(p_alpha, p_pep; layout = (1, 2), size = (1700, 625),
                    left_margin = 10 * Plots.mm, right_margin = 6 * Plots.mm,
                    bottom_margin = 18 * Plots.mm, top_margin = 5 * Plots.mm,
                    dpi = 180)
    save_pdf_png(p_alpha, outdir, "best_alpha_vs_normA")
    save_pdf_png(p_pep, outdir, "PEPbest_alpha_vs_normA")
    save_pdf_png(combined, outdir, "paper_figure2")
end

function relative_improvement_figure(best, horizons, kappa, outdir;
                                     colormap = "viridis", normalized = false,
                                     mark_uncertified = false, plot_style)
    selected = best[(best.kappa .== kappa) .& in.(best.T, Ref(horizons)), :]
    valid = isfinite.(selected.pep_best) .&
            isfinite.(selected.pep_analytic) .&
            (selected.pep_analytic .> 0) .&
            (selected.best_source .== "FOM design")
    improvement = copy(selected[valid, :])
    if nrow(improvement) == 0
        @warn "no analytical rows for the relative-improvement figure"
        return nothing
    end
    improvement.relative_improvement_percent = 100 .* (
        improvement.pep_analytic .- improvement.pep_best) ./
        improvement.pep_analytic

    excluded = nrow(selected) - nrow(improvement)
    excluded > 0 && @warn "excluded rows without matched finite FOM and analytic values from relative-improvement figure" excluded

    colors = matplotlib_palette(colormap, max(length(horizons), 2))
    improvement_style = style_kwargs(
        xlabel = L"\Vert A\Vert",
        ylabel = L"\mathrm{Relative\ improvement}\;[\%]",
        title = L"\mathrm{Relative\ improvement\ of\ the\ PEP\ optimal\ step\ size}",
        plot_style = plot_style)
    p = plot(; improvement_style..., legend = :topleft, legend_column = 2)
    for (j, T) in enumerate(horizons)
        curve = sort(improvement[improvement.T .== T, :], :norm_A)
        nrow(curve) == 0 && (@warn "no relative-improvement rows for horizon" T; continue)
        plot!(p, curve.norm_A, curve.relative_improvement_percent;
              color = colors[j], marker = :circle,
              label = latexstring(raw"T=", T))
    end
    uncertified = improvement[.!improvement.analytic_certified, :]
    mark_uncertified && nrow(uncertified) > 0 &&
        scatter!(p, uncertified.norm_A,
                 uncertified.relative_improvement_percent;
                 color = :red, marker = :xcross,
                 markersize = plot_style.star_size + 2,
                 markerstrokewidth = 2, label = "uncertified PEP")
    hline!(p, [0.0]; color = :black, linestyle = :dash, label = "")
    x_ticks = linear_ticks(improvement.norm_A)
    y_ticks = linear_ticks(improvement.relative_improvement_percent)
    y_upper = 1.1 * maximum(improvement.relative_improvement_percent)
    y_lower = minimum(tick_limits(y_ticks))
    plot!(p; xticks = x_ticks,
          xlims = (minimum(improvement.norm_A) - 0.02,
                   maximum(improvement.norm_A) + 0.04),
          yticks = y_ticks,
          ylims = (y_lower, y_upper))

    columns = [:T, :norm_A, :kappa, :alpha_analytic, :alpha_best,
               :pep_analytic, :pep_best, :best_source, :analytic_source,
               :analytic_certified,
               :relative_improvement_percent]
    normalized && append!(columns,
        [:W_zero, :pep_analytic_unnormalized, :pep_best_unnormalized])
    CSV.write(joinpath(outdir, "relative_improvement.csv"),
              sort(select(improvement, columns), [:T, :norm_A]))
    save_pdf_png(p, outdir, "relative_improvement_vs_normA")
end

function paper_figure3(clean, best, curve_T, requested_A, kappa, outdir;
                       partial = false, colormap = "viridis",
                       sweep_ymax = nothing, sweep_xmax = nothing, plot_style)
    available = best[(best.T .== curve_T) .& (best.kappa .== kappa), :norm_A]
    selected = Float64[]
    for target in requested_A
        value = nearest_value(available, target; label = "norm_A", tolerance = 0.06)
        value === nothing ? @warn("no complete group near requested norm_A", requested = target) : push!(selected, value)
    end
    unique!(selected)
    isempty(selected) && error("no complete groups for Figure 3 selections")

    p = plot(; style_kwargs(xlabel = L"\alpha",
                             ylabel = worst_case_label(false),
                             title = L"\mathrm{Worst\ case\ versus\ step\ size}",
                             plot_style = plot_style)...,
             legend = :topright)
    colors = matplotlib_palette(colormap, max(length(selected), 2))
    analytic_reference_plotted = false
    for (j, normA) in enumerate(selected)
        curve = sort(clean[(clean.T .== curve_T) .& (clean.kappa .== kappa) .&
                           (clean.norm_A .== normA), :], :alpha)
        nrow(curve) == 0 && continue
        plot!(p, curve.alpha, curve.pep; color = colors[j],
              label = latexstring(raw"L_A = ", @sprintf("%.2f", normA)))
        optimum = best[(best.T .== curve_T) .& (best.kappa .== kappa) .&
                       (best.norm_A .== normA), :]
        nrow(optimum) > 0 && scatter!(p, [first(optimum.alpha_best)],
            [first(optimum.pep_best)]; color = colors[j], marker = :star5,
            markersize = plot_style.star_size, label = "")
        if nrow(optimum) > 0
            ref = first(optimum.alpha_analytic)
            if isfinite(ref)
                vline!(p, [ref]; color = colors[j], linestyle = :dash, label = "")
                analytic_reference_plotted = true
            end
        end
    end
    analytic_reference_plotted &&
        plot!(p, [NaN], [NaN]; color = :black, linestyle = :dash,
              label = L"\alpha_{analytic}")
    # A=0 is the delayed-GD limiting case used in the paper, when available.
    zeroA = nearest_value(available, 0.0; label = "delayed-GD norm_A", tolerance = 1e-10)
    if zeroA !== nothing && zeroA ∉ selected
        curve = sort(clean[(clean.T .== curve_T) .& (clean.kappa .== kappa) .&
                           (clean.norm_A .== zeroA), :], :alpha)
        nrow(curve) > 0 && plot!(p, curve.alpha, curve.pep;
                                color = :red, linestyle = :dashdot,
                                label = L"\mathrm{Delayed\ GD}")
    end
    plotted = clean[(clean.T .== curve_T) .& (clean.kappa .== kappa) .&
                    in.(clean.norm_A, Ref([selected; zeroA])), :]
    CSV.write(joinpath(outdir, "alpha_sweep.csv"), sort(plotted, [:norm_A, :alpha]))
    alpha_ticks = sweep_xmax === nothing ? linear_ticks(plotted.alpha) :
                  linear_ticks([0.0, sweep_xmax])
    if sweep_ymax === nothing
        plot!(p; xticks = alpha_ticks, xlims = tick_limits(alpha_ticks))
    else
        plot!(p; xticks = alpha_ticks, xlims = tick_limits(alpha_ticks),
              yticks = linear_ticks([0.0, sweep_ymax]),
              ylims = (0.0, sweep_ymax))
    end
    save_pdf_png(p, outdir, "paper_figure3")
end

function horizon_figures(best, requested_A, kappa, outdir;
                         partial = false, colormap = "viridis",
                         normalized = false, mark_uncertified = false,
                         plot_horizons = nothing,
                         split_horizon_panels = false,
                         return_plots = false,
                         plot_style)
    if split_horizon_panels
        # Preserve the historical single-axis filenames and add the split
        # alternatives alongside them for direct comparison.
        horizon_figures(best, requested_A, kappa, outdir;
                        partial, colormap, normalized, mark_uncertified,
                        plot_horizons, plot_style)
        return split_horizon_figures(best, requested_A, kappa, outdir;
                                     partial, colormap, normalized,
                                     mark_uncertified, plot_horizons, plot_style)
    end
    horizon_mask = plot_horizons === nothing ? trues(nrow(best)) :
                   in.(best.T, Ref(plot_horizons))
    available = best[best.kappa .== kappa, :norm_A]
    selected = Float64[]
    for target in requested_A
        value = nearest_value(available, target; label = "norm_A", tolerance = 0.06)
        value === nothing || push!(selected, value)
    end
    unique!(selected)
    isempty(selected) && return
    colors = matplotlib_palette(colormap, max(length(selected), 2))
    pep_style = style_kwargs(xlabel = L"T",
                             ylabel = worst_case_label(normalized),
                             title = L"\mathrm{Best\ worst\ case\ versus\ horizon}",
                             plot_style = plot_style)
    p_pep_linear = plot(; pep_style..., yscale = :identity,
                        legend = :topright)
    p_pep_log10 = plot(; pep_style..., yscale = :log10,
                       legend = :topright)
    # Keep the alpha legend outside the axes so it cannot cover a curve.
    alpha_style = merge(style_kwargs(xlabel = L"T",
                                      ylabel = L"\alpha",
                                      title = L"\mathrm{Best\ step\ size\ versus\ horizon}",
                                      plot_style = plot_style),
                        (; right_margin = 20 * Plots.mm))
    p_alpha = plot(; alpha_style..., legend = :outerright)
    pep_linear_values = Float64[]
    analytic_linear_points = 0
    analytic_log10_points = 0
    uncertified_points = 0
    sampled_without_design = 0
    has_design = any(best.best_source .== "FOM design")
    for (j, normA) in enumerate(selected)
        curve = sort(best[(best.kappa .== kappa) .& (best.norm_A .== normA) .&
                          horizon_mask, :], :T)
        nrow(curve) == 0 && continue
        label = latexstring(raw"L_A = ", @sprintf("%.2f", normA))

        best_linear = curve[isfinite.(curve.pep_best) .&
                            (.!has_design .| (curve.best_source .== "FOM design")), :]
        if nrow(best_linear) > 0
            plot!(p_pep_linear, best_linear.T, best_linear.pep_best;
                  color = colors[j], marker = :circle, label)
            append!(pep_linear_values, best_linear.pep_best)
        end
        best_log10 = best_linear[best_linear.pep_best .> 0, :]
        if nrow(best_log10) > 0
            plot!(p_pep_log10, best_log10.T, best_log10.pep_best;
                  color = colors[j], marker = :circle, label)
        end
        sampled_only = has_design ?
            curve[curve.best_source .== "sampled sweep (design unavailable)", :] :
            curve[1:0, :]
        if nrow(sampled_only) > 0
            for p in (p_pep_linear, p_pep_log10)
                scatter!(p, sampled_only.T, sampled_only.pep_best;
                         color = :gray, marker = :square,
                         markersize = plot_style.star_size, label = "")
            end
            sampled_without_design += nrow(sampled_only)
        end

        analytic_linear = curve[isfinite.(curve.pep_analytic), :]
        if nrow(analytic_linear) > 0
            plot!(p_pep_linear, analytic_linear.T, analytic_linear.pep_analytic;
                  color = colors[j], linestyle = :dash, marker = :diamond,
                  markersize = 5, label = "")
            append!(pep_linear_values, analytic_linear.pep_analytic)
            analytic_linear_points += nrow(analytic_linear)
        end
        analytic_log10 = analytic_linear[analytic_linear.pep_analytic .> 0, :]
        if nrow(analytic_log10) > 0
            plot!(p_pep_log10, analytic_log10.T, analytic_log10.pep_analytic;
                  color = colors[j], linestyle = :dash, marker = :diamond,
                  markersize = 5, label = "")
            analytic_log10_points += nrow(analytic_log10)
        end
        uncertified = analytic_log10[.!analytic_log10.analytic_certified, :]
        if mark_uncertified && nrow(uncertified) > 0
            for p in (p_pep_linear, p_pep_log10)
                scatter!(p, uncertified.T, uncertified.pep_analytic;
                         color = :red, marker = :xcross,
                         markersize = plot_style.star_size + 2,
                         markerstrokewidth = 2, label = "")
            end
            uncertified_points += nrow(uncertified)
        end

        alpha_best_curve = has_design ?
            curve[curve.best_source .== "FOM design", :] : curve
        plot!(p_alpha, alpha_best_curve.T, alpha_best_curve.alpha_best;
              color = colors[j], marker = :circle, label)
        nrow(sampled_only) > 0 &&
            scatter!(p_alpha, sampled_only.T, sampled_only.alpha_best;
                     color = :gray, marker = :square,
                     markersize = plot_style.star_size, label = "")
        analytic_alpha_curve = curve[isfinite.(curve.alpha_analytic), :]
        nrow(analytic_alpha_curve) > 0 &&
            plot!(p_alpha, analytic_alpha_curve.T,
                  analytic_alpha_curve.alpha_analytic;
                  color = colors[j], linestyle = :dash, marker = :none,
                  label = "")
    end

    analytic_linear_points == 0 &&
        @warn "no finite analytical-step-size performance is available for the horizon plot"
    analytic_log10_points < analytic_linear_points &&
        @warn "omitting non-positive analytical performance values from log10 horizon plot" omitted = analytic_linear_points - analytic_log10_points

    # Colours identify ||A||; these two legend entries identify the step-size method.
    for p in (p_pep_linear, p_pep_log10)
        # In the split, two-column legend, a spacer keeps the two method
        # entries together on the last row instead of straddling columns.
        if return_plots && analytic_linear_points > 0 && isodd(length(selected))
            plot!(p, [NaN], [NaN]; color = :white, linealpha = 0,
                  marker = :none, label = " ")
        end
        plot!(p, [NaN], [NaN]; color = :black, marker = :circle,
              label = L"\alpha_{\mathrm{best}}")
        analytic_linear_points > 0 &&
            plot!(p, [NaN], [NaN]; color = :black, linestyle = :dash,
                  marker = :diamond, markersize = 5,
                  label = L"\alpha_{\mathrm{analytic}}")
        uncertified_points > 0 &&
            scatter!(p, [NaN], [NaN]; color = :red, marker = :xcross,
                     markersize = plot_style.star_size + 2,
                     label = "uncertified PEP")
        sampled_without_design > 0 &&
            scatter!(p, [NaN], [NaN]; color = :gray, marker = :square,
                     markersize = plot_style.star_size,
                     label = "sampled (no design)")
    end
    # In the split two-column legend, pad an odd number of L_A entries so
    # the two step-size methods share the final row (best left, analytic right).
    if return_plots && analytic_linear_points > 0 && isodd(length(selected))
        plot!(p_alpha, [NaN], [NaN]; color = :white, linealpha = 0,
              marker = :none, label = " ")
    end
    plot!(p_alpha, [NaN], [NaN]; color = :black, marker = :circle,
          label = L"\alpha_{\mathrm{best}}")
    analytic_linear_points > 0 &&
        plot!(p_alpha, [NaN], [NaN]; color = :black, linestyle = :dash,
              marker = :none,
              label = L"\alpha_{\mathrm{analytic}}")
    sampled_without_design > 0 &&
        scatter!(p_alpha, [NaN], [NaN]; color = :gray, marker = :square,
                 markersize = plot_style.star_size,
                 label = "sampled (no design)")
    selected_rows = best[(best.kappa .== kappa) .&
                         in.(best.norm_A, Ref(selected)) .&
                         horizon_mask, :]
    horizons = sort(unique(selected_rows.T))
    if plot_horizons === nothing
        horizon_ticks = sampled_ticks(horizons; max_ticks = 8,
                                      all_integer = true)
        if horizon_ticks !== :auto
            positions = sort(unique(vcat(horizon_ticks[1],
                                         Float64.(filter(in(horizons), [30, 40, 50])))))
            horizon_ticks = (positions,
                             latexstring.(string.(Int.(round.(positions)))))
        end
    else
        # Explicit horizon selection for paper plots: avoid crowding the labels
        # between 10 and 15 while showing every requested long-horizon anchor.
        anchors = filter(in(horizons), [3, 5, 10, 15, 20, 30, 40, 50])
        if isempty(anchors)
            horizon_ticks = sampled_ticks(horizons; max_ticks = 8,
                                          all_integer = true)
        else
            positions = Float64.(sort(unique(vcat(anchors,
                                                   [first(horizons), last(horizons)]))))
            horizon_ticks = (positions,
                             latexstring.(string.(Int.(round.(positions)))))
        end
    end
    rotation = 0.0
    pep_linear_ticks = linear_ticks(pep_linear_values)
    plot!(p_pep_linear; xticks = horizon_ticks, xrotation = rotation,
          yticks = pep_linear_ticks, ylims = tick_limits(pep_linear_ticks))
    # Keep Plots.jl's native log10 ticks/labels so the logarithmic scale is
    # explicit on the rendered axis instead of looking like a linear axis
    # with manually chosen decimal labels.
    plot!(p_pep_log10; xticks = horizon_ticks, xrotation = rotation)
    alpha_values = copy(selected_rows.alpha_best)
    append!(alpha_values, filter(isfinite, selected_rows.alpha_analytic))
    alpha_ticks = linear_ticks(alpha_values)
    plot!(p_alpha; xticks = horizon_ticks, xrotation = rotation,
          yticks = alpha_ticks, ylims = tick_limits(alpha_ticks))
    return_plots && return (p_pep_linear, p_pep_log10, p_alpha)
    # Keep the historical name as a linear-scale alias for existing workflows.
    save_pdf_png(p_pep_linear, outdir, "best_PEP_vs_T")
    save_pdf_png(p_pep_linear, outdir, "best_PEP_vs_T_linear")
    save_pdf_png(p_pep_log10, outdir, "best_PEP_vs_T_log10")
    save_pdf_png(p_alpha, outdir, "best_alpha_vs_T")
end

function split_horizon_figures(best, requested_A, kappa, outdir;
                               partial, colormap, normalized,
                               mark_uncertified, plot_horizons, plot_style)
    chosen = plot_horizons === nothing ? sort(unique(best[best.kappa .== kappa, :T])) :
             plot_horizons
    early = filter(T -> T < 20, chosen)
    late = filter(T -> T >= 20, chosen)
    isempty(early) && error("--split-horizon-panels needs a plotted T <= 15")
    isempty(late) && error("--split-horizon-panels needs a plotted T >= 20")

    left = horizon_figures(best, requested_A, kappa, outdir;
                           partial, colormap, normalized, mark_uncertified,
                           plot_horizons = early, return_plots = true, plot_style)
    right = horizon_figures(best, requested_A, kappa, outdir;
                            partial, colormap, normalized, mark_uncertified,
                            plot_horizons = late, return_plots = true, plot_style)
    selected_A = [nearest_value(best[best.kappa .== kappa, :norm_A], target;
                                label = "norm_A", tolerance = 0.06)
                  for target in requested_A]
    filter!(x -> x !== nothing, selected_A)
    rows = best[(best.kappa .== kappa) .& in.(best.norm_A, Ref(selected_A)) .&
                in.(best.T, Ref(vcat(early, late))), :]
    pep_values = filter(isfinite, vcat(rows.pep_best, rows.pep_analytic))
    positive_pep = filter(>(0), pep_values)
    alpha_values = filter(isfinite, vcat(rows.alpha_best, rows.alpha_analytic))
    isempty(positive_pep) && error("split horizon plot has no positive PEP values")
    linear_yticks = linear_ticks(pep_values)
    alpha_yticks = linear_ticks(alpha_values)
    y_settings = ((tick_limits(linear_yticks), linear_yticks),
                  ((minimum(positive_pep) * 0.8, maximum(positive_pep) * 2.5), :auto),
                  (tick_limits(alpha_yticks), alpha_yticks))
    stems = ("best_PEP_vs_T_linear_split", "best_PEP_vs_T_log10_split",
             "best_alpha_vs_T_split")
    for i in 1:3
        ylims, yticks = y_settings[i]
        left_plot, right_plot = left[i], right[i]
        if i == 3
            # The shorter alpha legend needs space above the curves when it
            # moves inside the wider early-horizon panel.
            ylo, yhi = ylims
            ylims = (ylo, yhi + 0.25 * (yhi - ylo))
        end
        # GR scales legend samples to the subplot width. The late-horizon
        # panel is deliberately narrow, so its native samples are tiny.
        plot!(left_plot; legend = :topright, legend_column = 2,
              xlims = (minimum(early) - 0.5, maximum(early) + 0.5),
              ylims, yticks, right_margin = 2 * Plots.mm)
        plot!(right_plot; legend = false,
              xlims = (19.0, maximum(late) + 2.0), ylims,
              yticks = (Float64[], String[]), ylabel = "",
              left_margin = 2 * Plots.mm, right_margin = 4 * Plots.mm)
        combined = plot(left_plot, right_plot;
                        layout = grid(1, 2; widths = [0.62, 0.38]),
                        size = (800, 500), dpi = 180)
        save_pdf_png(combined, outdir, stems[i])
        i == 1 && save_pdf_png(combined, outdir, "best_PEP_vs_T_split")
    end

    # Additional paper layout: performance above step size, with exactly the
    # same horizontal plotting area for both quantities in each T range.
    early_ticks = Float64.(sort(unique(vcat(first(early),
        filter(in(early), [3, 5, 10, 15]), last(early)))))
    late_ticks = Float64.(sort(unique(vcat(first(late),
        filter(in(late), [20, 30, 40, 50]), last(late)))))
    top_left, top_right = left[2], right[2]
    bottom_left, bottom_right = left[3], right[3]
    plot!(top_left; xlabel = "",
          xticks = (early_ticks, fill("", length(early_ticks))),
          left_margin = 12 * Plots.mm, bottom_margin = -1.5 * Plots.mm)
    plot!(top_right; xlabel = "",
          xticks = (late_ticks, fill("", length(late_ticks))),
          bottom_margin = -1.5 * Plots.mm)
    plot!(bottom_left; legend = false, left_margin = 12 * Plots.mm,
          top_margin = -1.5 * Plots.mm, ylims = (0.0, 1.0))
    plot!(bottom_right; legend = false, top_margin = -1.5 * Plots.mm,
          ylims = (0.0, 1.0))
    # Plots scales marker strokes when four plots are composed on one canvas;
    # thin them only for this layout so their fill colors remain visible.
    for panel in (top_left, top_right, bottom_left, bottom_right)
        for series in panel.series_list
            series[:markersize] *= 0.8
            series[:markerstrokewidth] = 0.25
        end
    end
    stacked = plot(top_left, top_right, bottom_left, bottom_right;
                   layout = grid(2, 2; widths = [0.62, 0.38]),
                   size = (900, 900), dpi = 180)

    # Magnify the L_A=0.7 endpoint comparison without removing the full
    # T=20:50 curves. The inset has no T labels, as the bottom row owns them.
    zoom_normA = nearest_value(selected_A, 0.7;
                               label = "inset L_A", tolerance = 0.06)
    if zoom_normA !== nothing && 50 in late
        zoom = sort(rows[(rows.norm_A .== zoom_normA) .&
                         in.(rows.T, Ref([40, 50])), :], :T)
        if nrow(zoom) == 2 && all(isfinite.(zoom.pep_best)) &&
           all(isfinite.(zoom.pep_analytic))
            zoom_color = matplotlib_palette(colormap,
                max(length(selected_A), 2))[findfirst(==(zoom_normA), selected_A)]
            zoom_lo, zoom_hi = extrema(vcat(zoom.pep_best, zoom.pep_analytic))
            zoom_pad = 0.12 * (zoom_hi - zoom_lo)
            zoom_ylims = (max(0.0, zoom_lo - zoom_pad), zoom_hi + zoom_pad)
            zoom_yticks = sort(unique(round.([zoom_lo, zoom_hi]; digits = 2)))
            # Mark the exact source region on the full log-scale panel so the
            # reader can immediately connect it to the magnified inset.
            plot!(stacked,
                  Shape([39.5, 50.5, 50.5, 39.5],
                        [zoom_ylims[1], zoom_ylims[1],
                         zoom_ylims[2], zoom_ylims[2]]);
                  subplot = 2, fillalpha = 0, linecolor = :red,
                  linewidth = 2, label = "")
            plot!(stacked, zoom.T, zoom.pep_best;
                  inset = (2, bbox(0.04, 0.04, 0.82, 0.42, :top, :right)),
                  subplot = 5, color = zoom_color, marker = :circle,
                  markersize = 5.6, markerstrokewidth = 0.25,
                  linewidth = 2, label = "", legend = false,
                  xlims = (38, 60), ylims = zoom_ylims,
                  xticks = ([40, 50], ["", ""]),
                  yticks = (zoom_yticks, latex_tick.(zoom_yticks)),
                  yscale = :identity, framestyle = :box,
                  foreground_color_border = :red,
                  tickfontsize = 8, grid = true,
                  background_color_inside = :white)
            plot!(stacked, zoom.T, zoom.pep_analytic;
                  subplot = 5, color = zoom_color, linestyle = :dash,
                  marker = :diamond, markersize = 4.0,
                  markerstrokewidth = 0.25, linewidth = 2, label = "")
            at50 = zoom[zoom.T .== 50, :]
            best50, analytic50 = first(at50.pep_best), first(at50.pep_analytic)
            improvement = 100 * (1 - best50 / analytic50)
            plot!(stacked, [52, 52], [best50, analytic50];
                  subplot = 5, color = :black, linewidth = 1,
                  marker = :none, label = "")
            for value in (best50, analytic50)
                plot!(stacked, [51.6, 52.4], [value, value];
                      subplot = 5, color = :black, linewidth = 1,
                      marker = :none, label = "")
            end
            annotate!(stacked, 56, (best50 + analytic50) / 2,
                      text(@sprintf("%.0f%% lower", improvement), 8, :black);
                      subplot = 5)
        end
    end
    save_pdf_png(stacked, outdir, "best_PEP_and_alpha_vs_T_split")
end

function write_summary(outdir, run, raw, shape, counts, complete, clean, sampled,
                       best, kappa, partial, design_run, design_rows, design_matches,
                       analytic_mode, analytic_run, analytic_rows,
                       analytic_matches, include_uncertified,
                       mark_uncertified, colormap,
                       normalized, sweep_ymax, sweep_xmax, plot_horizons,
                       split_horizon_panels,
                       plot_style)
    selected_best = best[best.kappa .== kappa, :]
    open(joinpath(outdir, "plot_summary.txt"), "w") do io
        println(io, "run = ", relpath(run, joinpath(@__DIR__, "..")))
        println(io, "raw_rows = ", nrow(raw))
        println(io, "expected_rows = ", shape.expected_rows)
        println(io, "expected_alpha_per_group = ", shape.nalpha)
        println(io, "observed_groups = ", nrow(counts))
        println(io, "complete_groups = ", length(complete))
        println(io, "valid_rows_in_complete_groups = ", nrow(clean))
        println(io, "selected_kappa = ", kappa)
        println(io, "matplotlib_colormap = ", colormap)
        println(io, "normalized_worst_cases = ", normalized)
        println(io, "alpha_sweep_normalized = false")
        println(io, "sweep_ymax = ", something(sweep_ymax, "auto"))
        println(io, "sweep_xmax = ", something(sweep_xmax, "auto"))
        println(io, "T_plot_horizons = ", something(plot_horizons, "all"))
        println(io, "split_horizon_panels = ", split_horizon_panels)
        println(io, "titles_enabled = ", plot_style.show_titles)
        println(io, "title_font_size = ", plot_style.title_font_size)
        println(io, "axis_font_size = ", plot_style.axis_font_size)
        println(io, "legend_font_size = ", plot_style.legend_font_size)
        println(io, "tick_font_size = ", plot_style.tick_font_size)
        println(io, "star_size = ", plot_style.star_size)
        println(io, "line_width = ", plot_style.line_width)
        normalized && println(io,
            "worst_case_normalization = W(alpha,T) / W(0,T), " *
            "W(0,T) = rho^2 / (2(1-L_A)^2), rho = R_u; " *
            "the alpha-sweep figure is excluded")
        println(io, "partial = ", partial)
        println(io, "design_run = ", design_run === nothing ? "none" : relpath(design_run, joinpath(@__DIR__, "..")))
        println(io, "valid_FOM_design_rows = ", design_rows)
        println(io, "matched_FOM_design_groups = ", design_matches)
        println(io, "design_only_groups = ", count(isnan.(best.pep_best_sampled) .&
                                                    (best.best_source .== "FOM design")))
        println(io, "analytical_mode = ", analytic_mode)
        println(io, "analytical_from_PEP_run = ", analytic_run === nothing ? "none" : relpath(analytic_run, joinpath(@__DIR__, "..")))
        println(io, "valid_analytical_from_PEP_rows = ", analytic_rows)
        println(io, "matched_analytical_from_PEP_groups = ", analytic_matches)
        println(io, "include_uncertified_analytic = ", include_uncertified)
        println(io, "mark_uncertified_analytic = ", mark_uncertified)
        println(io, "selected_best_point_rows = ", nrow(selected_best))
        println(io, "uncertified_analytic_points = ", count(.!selected_best.analytic_certified .& isfinite.(selected_best.pep_analytic)))
        include_uncertified && println(io, "UNCERTIFIED rows are finite solver estimates, not certified PEP bounds; see the diagnostics CSV beside the analytical Arrow file.")
        println(io, "best-point source counts = ", combine(groupby(selected_best, :best_source), nrow))
        if analytic_mode == :recovered
            println(io, "analytical curves use recovered A/B at alpha_best, the paper corollary, and a fresh class-PEP solve.")
            println(io, "nearest-sweep analytical diagnostics are retained in sampled_best.csv but are not plotted.")
        elseif analytic_mode == :theorem1
            println(io, "analytical curves use the sharp uniform Theorem-1 step size and a fresh class-PEP solve.")
            println(io, "Theorem 1 gives the same alpha as the former formula; its performance here is evaluated at the exact alpha rather than at the nearest sweep sample.")
        else
            println(io, "analytical curves use the former norm-bound formula; performance is the nearest sampled sweep value.")
        end
    end
    CSV.write(joinpath(outdir, "group_coverage.csv"), sort(counts[counts.kappa .== kappa, :], GROUP_KEY))
    CSV.write(joinpath(outdir, "sampled_best.csv"), sort(sampled[sampled.kappa .== kappa, :], GROUP_KEY))
    CSV.write(joinpath(outdir, "best_points.csv"), sort(selected_best, GROUP_KEY))
end

function main()
    run, opts = parse_args(ARGS)
    colormap = String(opts["colormap"])
    normalized = get(opts, "normalize-worst-cases", false)
    # The alpha-sweep figure always shows raw regret, so its default clipping
    # is independent of normalization used by the other figures.
    sweep_ymax = parse_sweep_ymax(opts["sweep-ymax"])
    sweep_xmax = haskey(opts, "sweep-xmax") ? parse_positive_plot_value(opts, "sweep-xmax") : nothing
    plot_style = parse_plot_style(opts)
    # Validate before loading potentially large result tables.
    matplotlib_palette(colormap, 2)
    run_parent = isdir(run) ? run : dirname(run)
    outdir = abspath(get(opts, "out", joinpath(run_parent, "paper_figures")))
    mkpath(outdir)

    raw = load_run(run)
    required = [:T, :norm_A, :norm_B, :kappa, :alpha, :pep,
                :mu_f, :L_f, :R_x, :R_u]
    missingcols = filter(c -> !hasproperty(raw, c), required)
    isempty(missingcols) || error("sweep data lack columns $missingcols")

    shape = expected_shape(run, raw)
    complete_raw, counts, complete = complete_sweep_rows(raw, shape.nalpha)
    partial = shape.expected_rows !== nothing && nrow(raw) != shape.expected_rows
    partial && !get(opts, "allow-partial", false) &&
        error("run is partial ($(nrow(raw))/$(shape.expected_rows) rows); pass --allow-partial to plot only complete alpha groups")
    clean = clean_rows(complete_raw, [:alpha, :pep])
    sampled = sampled_best_table(clean)
    nrow(sampled) > 0 || error("no complete, valid alpha groups available")

    sibling_design = joinpath(run_parent, "merged_design.arrow")
    design_run = if get(opts, "sweep-only", false)
        nothing
    elseif haskey(opts, "design")
        abspath(opts["design"])
    elseif isfile(sibling_design)
        @info "using sibling merged design file" design = sibling_design
        abspath(sibling_design)
    else
        nothing
    end
    design = design_run === nothing ? nothing :
             fom_design_table(load_run(design_run);
                              include_parameters = get(opts, "include-design-only", false))
    best, design_matches = select_best_points(sampled, design;
        allow_partial_design = get(opts, "allow-partial-design", false),
        include_design_only = get(opts, "include-design-only", false))
    design_rows = design === nothing ? 0 : nrow(design)

    analytic_mode = parse_analytic_mode(opts["analytic-mode"])
    include_uncertified = get(opts, "include-uncertified-analytic", false)
    mark_uncertified = include_uncertified || get(opts, "mark-uncertified-analytic", false)
    include_uncertified && analytic_mode != :theorem1 &&
        error("--include-uncertified-analytic currently requires --analytic-mode theorem1")
    mark_uncertified && !include_uncertified &&
        error("--mark-uncertified-analytic requires --include-uncertified-analytic")
    analytic_run = nothing
    analytic_rows = 0
    if analytic_mode == :former
        best, analytic_matches = attach_former_analytical(best)
    else
        sibling_candidates = if analytic_mode == :recovered
            [joinpath(run_parent, "analytic_from_pep", "merged.arrow"),
             joinpath(run_parent, "merged_analytic_recovered.arrow"),
             joinpath(run_parent, "merged_analytic.arrow")]
        else
            [joinpath(run_parent, "analytic_theorem1", "merged.arrow"),
             joinpath(run_parent, "merged_analytic_theorem1.arrow")]
        end
        sibling_analytic = findfirst(isfile, sibling_candidates)
        analytic_run = if get(opts, "sweep-only", false)
            nothing
        elseif haskey(opts, "analytic")
            abspath(opts["analytic"])
        elseif sibling_analytic !== nothing
            path = sibling_candidates[sibling_analytic]
            @info "using analytical PEP results" analytic_mode analytical = path
            abspath(path)
        else
            nothing
        end
        analytical = analytic_run === nothing ? nothing :
                     pep_analytic_table(load_run(analytic_run);
                                        include_uncertified)
        if include_uncertified
            best, _ = add_analytic_only_sampled(best, sampled, analytical)
        end
        source = analytic_mode == :recovered ? "recovered A/B + PEP" :
                 "Theorem 1 uniform bound + PEP"
        best, analytic_matches = attach_pep_analytical(
            best, analytical, source)
        analytic_rows = analytical === nothing ? 0 : nrow(analytical)
    end

    # Keep the entire alpha-sweep figure, including its optimum stars, in raw
    # regret units even when all summary/best-point figures are normalized.
    sweep_best = best
    if normalized
        sampled = normalize_worst_case_table(
            sampled, [:pep_best_sampled, :pep_analytic_sampled])
        best = normalize_worst_case_table(
            best, [:pep_best_sampled, :pep_analytic_sampled,
                   :pep_best, :pep_analytic])
    end

    target_kappa = haskey(opts, "kappa") ? parse(Float64, opts["kappa"]) : maximum(best.kappa)
    kappa = nearest_value(best.kappa, target_kappa; label = "kappa", tolerance = 1e-10)
    kappa === nothing && error("no kappa available")
    horizons = csv_ints(opts["horizons"])
    plot_horizons = lowercase(strip(opts["T-plot-horizons"])) == "all" ?
        nothing : sort(unique(csv_ints(opts["T-plot-horizons"])))
    if plot_horizons !== nothing &&
       (isempty(plot_horizons) || any(T -> T <= 0, plot_horizons))
        error("--T-plot-horizons needs positive comma-separated integers or 'all'")
    end
    curve_T = parse(Int, opts["curve-T"])
    norm_A = csv_floats(opts["norm-a"])
    split_horizon_panels = opts["split-horizon-panels"]

    write_summary(outdir, run, raw, shape, counts, complete, clean, sampled, best,
                  kappa, partial, design_run, design_rows, design_matches,
                  analytic_mode, analytic_run, analytic_rows, analytic_matches,
                  include_uncertified, mark_uncertified, colormap, normalized,
                  sweep_ymax, sweep_xmax, plot_horizons, split_horizon_panels,
                  plot_style)
    paper_figure2(best, horizons, kappa, outdir;
                  partial, colormap, normalized, mark_uncertified, plot_style)
    analytic_matches > 0 &&
        relative_improvement_figure(best, horizons, kappa, outdir;
                                    colormap, normalized, mark_uncertified,
                                    plot_style)
    paper_figure3(clean, sweep_best, curve_T, norm_A, kappa, outdir;
                  partial, colormap, sweep_ymax, sweep_xmax, plot_style)
    horizon_figures(best, norm_A, kappa, outdir;
                    partial, colormap, normalized, mark_uncertified,
                    plot_horizons, split_horizon_panels,
                    plot_style)

    @info "paper-style figures written" outdir raw_rows = nrow(raw) complete_groups = length(complete) selected_kappa = kappa colormap normalized alpha_sweep_normalized = false sweep_ymax plot_style partial design_run design_matches analytic_mode analytic_run analytic_matches
end

main()
