# Evaluate an exact analytical step size with a fresh worst-case PEP
# certificate. Two methods share the same filtering, parallel execution,
# crash-safe storage, and output schema:
#
#   recovered: recover a worst-case realization at alpha_best and
#       evaluate the plant-specific corollary on those matrices;
#   theorem1 (default): evaluate the sharp uniform Theorem-1 threshold directly from
#       the norm bounds, then solve PEP at that exact step size.
#
# In recovered mode, for every finite FOM winner this script:
#   1. solves the PEP again at alpha_best and retains its primal Gram matrix;
#   2. recovers explicit matrices A and B from all forward/adjoint oracle pairs;
#   3. evaluates the paper's corollary with G=(I-A)\B and
#      A'PA-P=-I to obtain alpha_analytic;
#   4. solves the original class PEP at alpha_analytic.
#
# Results are crash-safe and resumable. Each completed group is written to an
# Arrow shard and each recovered realization is saved as JLD2.
#
# Usage from Julia/:
#   julia --project=. --threads=auto scripts/evaluate_analytic_from_pep.jl \
#       results/full_YYYYMMDD_HHMMSS/merged_design.arrow --method recovered
#
#   julia --project=. --threads=auto scripts/evaluate_analytic_from_pep.jl \
#       results/full_YYYYMMDD_HHMMSS/merged_design.arrow --method theorem1
#
# Useful laptop subsets:
#   ... --horizons 10,15,20 --kappa 0 --norm-a 0.2,0.4,0.6,0.8
#   ... --paper-subset                # union needed by default paper plots
#   ... --paper-subset --paper-horizons 5
#   ... --smoke                       # first unfinished group only

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

include(joinpath(@__DIR__, "..", "src", "PEPControl.jl"))
include(joinpath(@__DIR__, "..", "src", "PEPWorstCase.jl"))
using .PEPControl, .PEPWorstCase
using Arrow, DataFrames, Dates, JLD2, LinearAlgebra, Printf

const PD = PEPControl.PEPDesign
const GROUP_KEY = (:T, :norm_A, :kappa)

function parse_args(args)
    isempty(args) && error(
        "usage: evaluate_analytic_from_pep.jl DESIGN_ARROW [options]")
    design = abspath(args[1])
    opts = Dict{String,Any}(
        "batch-rows" => "1",
        "factor-rtol" => "1e-8",
        "factor-atol" => "1e-10",
        "cert-tol" => "1e-7",
        # Recovery needs the primal Gram matrix. Mosek can stop with
        # SLOW_PROGRESS around a few 1e-6 on these facially reduced PEPs;
        # retain the exact gap in the output and accept up to 1e-5.
        "gap-tol" => "1e-5",
        "facial-reduction" => "explicit",
        "method" => "theorem1",
    )
    i = 2
    while i <= length(args)
        arg = args[i]
        if arg == "--smoke"
            opts["limit"] = "1"
            i += 1
        elseif arg == "--paper-subset"
            opts["paper-subset"] = true
            i += 1
        elseif arg == "--retain-uncertified-analytic" || arg == "--only-unusable-design"
            opts[arg[3:end]] = true
            i += 1
        elseif startswith(arg, "--")
            i < length(args) || error("$arg expects a value")
            opts[arg[3:end]] = args[i + 1]
            i += 2
        else
            error("unexpected argument $arg")
        end
    end
    return design, opts
end

csv_floats(value) = parse.(Float64,
    filter(!isempty, strip.(split(value, ','))))
csv_ints(value) = parse.(Int, filter(!isempty, strip.(split(value, ','))))

function design_file(path)
    isfile(path) && return path
    isdir(path) || error("design input does not exist: $path")
    merged = joinpath(path, "merged.arrow")
    isfile(merged) || error("design directory contains no merged.arrow: $path")
    return merged
end

function fom_winners(path; only_unusable = false)
    frame = DataFrame(Arrow.Table(design_file(path)); copycols = true)
    required = [:T, :norm_A, :kappa, :mu_f, :L_f, :norm_B, :R_x, :R_u,
                :optimizer, :alpha_best, :pep_best]
    missingcols = filter(column -> !hasproperty(frame, column), required)
    isempty(missingcols) || error("design data lack columns $missingcols")
    mask = lowercase.(String.(frame.optimizer)) .== "fom"
    usable = isfinite.(frame.alpha_best) .& (frame.alpha_best .> 0) .&
             isfinite.(frame.pep_best)
    if only_unusable
        mask .&= .!usable .& (frame.norm_A .< 1)
    else
        mask .&= usable
    end
    finite_fom = frame[mask, :]
    nrow(finite_fom) > 0 || error("design data contain no finite FOM result")

    winner_indices = Int[]
    for group in groupby(finite_fom, collect(GROUP_KEY); sort = true)
        parent_rows = parentindices(group)[1]
        push!(winner_indices, parent_rows[argmin(group.pep_best)])
    end
    return sort(finite_fom[winner_indices, :], collect(GROUP_KEY))
end

function nearest_values(available, requested, label)
    values = sort(unique(Float64.(available)))
    selected = Float64[]
    for target in requested
        value = values[argmin(abs.(values .- target))]
        value == target || @warn "snapped requested $label to design grid" requested = target selected = value
        push!(selected, value)
    end
    return unique(selected)
end

function filter_winners(winners, opts)
    selected = winners
    if get(opts, "paper-subset", false)
        target_kappa = if haskey(opts, "kappa")
            only(nearest_values(selected.kappa,
                                [parse(Float64, opts["kappa"])], "kappa"))
        else
            maximum(selected.kappa)
        end
        at_kappa = selected[selected.kappa .== target_kappa, :]
        paper_horizons = haskey(opts, "paper-horizons") ?
                         csv_ints(opts["paper-horizons"]) : [10, 15, 20]
        requested_norm_A = haskey(opts, "paper-norm-a") ?
                           csv_floats(opts["paper-norm-a"]) :
                           [0.2, 0.4, 0.6, 0.8]
        paper_norm_A = nearest_values(
            at_kappa.norm_A, requested_norm_A, "norm_A")
        mask = (selected.kappa .== target_kappa) .&
               (in.(selected.T, Ref(paper_horizons)) .|
                in.(selected.norm_A, Ref(paper_norm_A)))
        selected = selected[mask, :]
    end
    if haskey(opts, "horizons")
        horizons = csv_ints(opts["horizons"])
        selected = selected[in.(selected.T, Ref(horizons)), :]
    end
    if haskey(opts, "kappa") && nrow(selected) > 0
        requested = parse(Float64, opts["kappa"])
        value = nearest_values(selected.kappa, [requested], "kappa")
        selected = selected[in.(selected.kappa, Ref(value)), :]
    end
    if haskey(opts, "norm-a") && nrow(selected) > 0
        requested = csv_floats(opts["norm-a"])
        values = nearest_values(selected.norm_A, requested, "norm_A")
        selected = selected[in.(selected.norm_A, Ref(values)), :]
    end
    if haskey(opts, "limit")
        limit = parse(Int, opts["limit"])
        limit >= 1 || error("--limit must be positive")
        selected = first(selected, min(limit, nrow(selected)))
    end
    nrow(selected) > 0 || error("no FOM design groups match the requested filters")
    return selected
end

const ResultRow = @NamedTuple{
    T::Int, n_steps::Int, kappa::Float64, mu_f::Float64, L_f::Float64,
    norm_A::Float64, norm_B::Float64, R_x::Float64, R_u::Float64,
    alpha_best::Float64, pep_best_design::Float64, pep_best_resolved::Float64,
    alpha_analytic::Float64, pep_analytic::Float64,
    gram_rank::Int, gram_threshold::Float64, gram_min_eigenvalue::Float64,
    gram_relative_error::Float64,
    recovered_norm_A::Float64, recovered_norm_B::Float64,
    operator_A_fit_residual::Float64, operator_B_fit_residual::Float64,
    spectral_radius_A::Float64, norm_G::Float64, norm_P::Float64,
    L_G::Float64, L_tilde::Float64, zeta::Float64, omega::Float64,
    lyapunov_residual::Float64, lyapunov_iterations::Int,
    best_relative_gap::Float64, analytic_relative_gap::Float64,
    best_scalar_violation::Float64, analytic_scalar_violation::Float64,
    best_cone_violation::Float64, analytic_cone_violation::Float64,
    best_solver_status::String, analytic_solver_status::String,
    realization_file::String, status::String, solve_time::Float64
}

function failed_result(row, message, elapsed)
    nan = NaN
    return ResultRow((
        Int(row.T), Int(row.T) + 1, Float64(row.kappa), Float64(row.mu_f),
        Float64(row.L_f), Float64(row.norm_A), Float64(row.norm_B),
        Float64(row.R_x), Float64(row.R_u), Float64(row.alpha_best),
        Float64(row.pep_best), nan, nan, nan,
        0, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
        nan, nan, nan, nan, nan, 0,
        nan, nan, nan, nan, nan, nan,
        "", "", "", "ERROR: " * first(message, 500), elapsed))
end

function artifact_name(row)
    a = replace(@sprintf("%.12g", Float64(row.norm_A)), "+" => "p")
    k = replace(@sprintf("%.12g", Float64(row.kappa)), "+" => "p")
    return @sprintf("T%03d_A%s_kappa%s.jld2", Int(row.T), a, k)
end

function evaluate_theorem1_group(row, instance, compiled;
                                 cert_tol, gap_tol, facial_reduction,
                                 retain_uncertified,
                                 started)
    nan = NaN
    alpha = theorem1_uniform_step_size(
        instance.norm_A, instance.norm_B, instance.L_f)
    solution = PD.solve_pep(
        compiled, fill(alpha, nsteps(instance)); warn = false,
        cert_tol, gap_tol, facial_reduction)
    if !solution.diagnostics.certified && !retain_uncertified
        error("PEP at the Theorem-1 step size is not certified: " *
              "status=$(solution.status), gap=$(solution.diagnostics.relative_gap)")
    end
    isfinite(solution.obj) && solution.obj > 0 ||
        error("PEP at the Theorem-1 step size has no finite positive estimate")
    diagnostics = solution.diagnostics
    return ResultRow((
        T = instance.T,
        n_steps = nsteps(instance),
        kappa = Float64(row.kappa),
        mu_f = instance.mu_f,
        L_f = instance.L_f,
        norm_A = instance.norm_A,
        norm_B = instance.norm_B,
        R_x = instance.R_x,
        R_u = instance.R_u,
        alpha_best = Float64(row.alpha_best),
        pep_best_design = Float64(row.pep_best),
        pep_best_resolved = nan,
        alpha_analytic = alpha,
        pep_analytic = solution.obj,
        gram_rank = 0,
        gram_threshold = nan,
        gram_min_eigenvalue = nan,
        gram_relative_error = nan,
        recovered_norm_A = nan,
        recovered_norm_B = nan,
        operator_A_fit_residual = nan,
        operator_B_fit_residual = nan,
        spectral_radius_A = nan,
        norm_G = nan,
        norm_P = nan,
        L_G = nan,
        L_tilde = nan,
        zeta = nan,
        omega = nan,
        lyapunov_residual = nan,
        lyapunov_iterations = 0,
        best_relative_gap = nan,
        analytic_relative_gap = diagnostics.relative_gap,
        best_scalar_violation = nan,
        analytic_scalar_violation = diagnostics.max_scalar_violation,
        best_cone_violation = nan,
        analytic_cone_violation = diagnostics.max_cone_violation,
        best_solver_status = "NOT_RUN",
        analytic_solver_status = string(solution.status),
        realization_file = "",
        status = solution.diagnostics.certified ? "OK" : "UNCERTIFIED",
        solve_time = time() - started,
    ))
end

function evaluate_group(row, outdir, artifact_lock;
                        factor_rtol, factor_atol, cert_tol, gap_tol,
                        facial_reduction, method, retain_uncertified)
    started = time()
    try
        instance = ControlInstance(
            T = Int(row.T), norm_A = Float64(row.norm_A),
            norm_B = Float64(row.norm_B), R_x = Float64(row.R_x),
            R_u = Float64(row.R_u), mu_f = Float64(row.mu_f),
            L_f = Float64(row.L_f))
        model, _ = build_pep_control(instance)
        compiled = PD.compile(model)
        if method == :theorem1
            return evaluate_theorem1_group(
                row, instance, compiled; cert_tol, gap_tol,
                facial_reduction, retain_uncertified, started)
        end

        alpha_best = Float64(row.alpha_best)
        eta_best = fill(alpha_best, nsteps(instance))
        best_solution = PD.solve_pep(
            compiled, eta_best; warn = false, cert_tol, gap_tol,
            facial_reduction)
        best_solution.diagnostics.certified || error(
            "PEP at alpha_best is not certified: status=$(best_solution.status), " *
            "gap=$(best_solution.diagnostics.relative_gap)")

        gram = factor_gram(best_solution.G; rtol = factor_rtol,
                           atol = factor_atol)
        gram.rank > 0 || error("best-alpha PEP returned a rank-zero Gram matrix")
        operators = extract_control_operators(model, eta_best, gram.factor)
        analytical = corollary_step_size(operators.A, operators.B, instance.L_f)
        isfinite(analytical.alpha) && analytical.alpha > 0 ||
            error("analytical step size is not finite and positive")

        eta_analytic = fill(analytical.alpha, nsteps(instance))
        analytic_solution = PD.solve_pep(
            compiled, eta_analytic; warn = false, cert_tol, gap_tol,
            facial_reduction)
        analytic_solution.diagnostics.certified || error(
            "PEP at alpha_analytic is not certified: status=$(analytic_solution.status), " *
            "gap=$(analytic_solution.diagnostics.relative_gap)")

        relative_artifact = joinpath("realizations", artifact_name(row))
        absolute_artifact = joinpath(outdir, relative_artifact)
        lock(artifact_lock) do
            JLD2.jldsave(absolute_artifact;
                T = instance.T, kappa = Float64(row.kappa),
                norm_A_bound = instance.norm_A, norm_B_bound = instance.norm_B,
                alpha_best, pep_best = best_solution.obj,
                alpha_analytic = analytical.alpha,
                pep_analytic = analytic_solution.obj,
                A = operators.A, B = operators.B,
                G = analytical.G, P = analytical.P,
                operator_A_fit_residual = operators.residual_A,
                operator_B_fit_residual = operators.residual_B,
                gram_eigenvalues = gram.eigenvalues,
                gram_threshold = gram.threshold,
                gram_relative_error = gram.relative_error,
                best_diagnostics = best_solution.diagnostics,
                analytic_diagnostics = analytic_solution.diagnostics)
        end

        bd = best_solution.diagnostics
        ad = analytic_solution.diagnostics
        return ResultRow((
            instance.T, nsteps(instance), Float64(row.kappa), instance.mu_f,
            instance.L_f, instance.norm_A, instance.norm_B, instance.R_x,
            instance.R_u, alpha_best, Float64(row.pep_best), best_solution.obj,
            analytical.alpha, analytic_solution.obj,
            gram.rank, gram.threshold, gram.min_eigenvalue,
            gram.relative_error, opnorm(operators.A), opnorm(operators.B),
            operators.residual_A, operators.residual_B,
            analytical.spectral_radius_A, analytical.norm_G,
            opnorm(analytical.P), analytical.L_G, analytical.L_tilde,
            analytical.zeta, analytical.omega, analytical.lyapunov_residual,
            analytical.lyapunov_iterations,
            bd.relative_gap, ad.relative_gap,
            bd.max_scalar_violation, ad.max_scalar_violation,
            bd.max_cone_violation, ad.max_cone_violation,
            string(best_solution.status), string(analytic_solution.status),
            relative_artifact, "OK", time() - started))
    catch err
        message = sprint(showerror, err, catch_backtrace())
        @error "analytical-from-PEP group failed" T = row.T norm_A = row.norm_A kappa = row.kappa exception = (err, catch_backtrace())
        return failed_result(row, message, time() - started)
    end
end

function main()
    design, opts = parse_args(ARGS)
    input_file = design_file(design)
    input_parent = dirname(input_file)
    method = Symbol(lowercase(strip(opts["method"])))
    method in (:recovered, :theorem1) ||
        error("--method must be recovered or theorem1")
    retain_uncertified = get(opts, "retain-uncertified-analytic", false)
    only_unusable = get(opts, "only-unusable-design", false)
    (retain_uncertified || only_unusable) && method != :theorem1 &&
        error("uncertified and unusable-design options require --method theorem1")
    default_output = method == :recovered ? "analytic_from_pep" :
                     "analytic_theorem1"
    outdir = abspath(get(opts, "out", joinpath(input_parent, default_output)))
    method == :recovered && mkpath(joinpath(outdir, "realizations"))
    mkpath(outdir)

    factor_rtol = parse(Float64, opts["factor-rtol"])
    factor_atol = parse(Float64, opts["factor-atol"])
    cert_tol = parse(Float64, opts["cert-tol"])
    gap_tol = parse(Float64, opts["gap-tol"])
    facial_reduction = Symbol(lowercase(opts["facial-reduction"]))
    facial_reduction in (:none, :explicit) ||
        error("--facial-reduction must be none or explicit")
    batch_rows = parse(Int, opts["batch-rows"])
    batch_rows >= 1 || error("--batch-rows must be positive")

    winners = filter_winners(fom_winners(input_file; only_unusable), opts)
    array = array_context()
    task_dir = array_task_dir(outdir, array.task_id, array.task_count)
    mkpath(task_dir)
    assigned = array_partition(collect(eachrow(winners)),
                               array.task_id, array.task_count)
    stored = load_results(task_dir)
    successful = nrow(stored) == 0 ? stored :
        stored[in.(stored.status, Ref(("OK", "UNCERTIFIED"))), :]
    done = Set{Tuple}((Int(row.T), Float64(row.norm_A), Float64(row.kappa))
                      for row in eachrow(successful))
    todo = [row for row in assigned
            if (Int(row.T), Float64(row.norm_A), Float64(row.kappa)) ∉ done]

    config = (
        design = input_file, method, selected_groups = nrow(winners),
        factor_rtol, factor_atol, cert_tol, gap_tol, facial_reduction,
        retain_uncertified, only_unusable,
        interpretation = method == :recovered ?
            "recover A,B at alpha_best; compute corollary alpha; " *
            "evaluate alpha in original class PEP" :
            "compute the exact Theorem-1 uniform alpha from norm bounds; " *
            "evaluate alpha in original class PEP")
    # Refresh the manifest on resumed runs so it describes the current
    # selection (for example, a smoke check followed by --paper-subset).
    array.task_id == 0 && write_config(outdir, merge(config, (
        array_task_count = array.task_count,)))
    writer = ResultWriter(ResultRow, task_dir; batch = batch_rows)
    artifact_lock = ReentrantLock()
    BLAS.set_num_threads(1)

    @info "evaluate_analytic_from_pep" method design = input_file outdir task_dir task_id = array.task_id task_count = array.task_count selected = nrow(winners) assigned = length(assigned) resumed = length(done) remaining = length(todo) threads = Threads.nthreads()
    finished = Threads.Atomic{Int}(0)
    failures = Threads.Atomic{Int}(0)
    started = time()
    pforeach(todo) do row
        result = evaluate_group(row, outdir, artifact_lock;
            factor_rtol, factor_atol, cert_tol, gap_tol, facial_reduction,
            method, retain_uncertified)
        record!(writer, result)
        if startswith(result.status, "ERROR")
            Threads.atomic_add!(failures, 1)
        end
        count = Threads.atomic_add!(finished, 1) + 1
        @info @sprintf("[%5.1f%%] T=%d ||A||=%.3f kappa=%.3g status=%s (%.1f s)",
                      100count / max(length(todo), 1), result.T, result.norm_A,
                      result.kappa, result.status, time() - started)
    end

    new_rows = close_writer!(writer)
    failures[] == 0 || error(
        "$(failures[]) analytical group(s) failed; successful shards were kept " *
        "and the failed groups will be retried when this task is resubmitted")
    open(joinpath(task_dir, "_SUCCESS"), "w") do io
        println(io, "completed=", Dates.now())
    end
    merged = array.task_count == 1 ? merge_parts(outdir) : nothing
    @info "evaluate_analytic_from_pep: done" method task_id = array.task_id new_rows merged walltime = time() - started
end

main()
